! -*- f90 -*-

!連続変換する場合は@NCMPを使って成分数を明示する必要がある。
!平成16年2月4日(水)-rオプションで(不完全ながら)動径分布関数も出力でき
!るようになった。
!
!outputtype
! 1: NGPH
! 2: WGPH
! 3: MDVW
! 4: RDF
! 5: Pair Interaction (@PAIR)
! 6: NEAR  ... 各分子から近い順に10個まで出力する。
! 7: 最近接4分子のなすO-O-O角度を分布の形式で出力する。(場合の数は分子数の6倍)
! 8: 水素結合で連続する4分子のなす二面角の分布。(そのような組みあわせは非常にたくさんあるが…)
! 9: 動径ヒストグラム(規格化しない動径分布関数)
!10: tetrahedrality。分布を得たい場合はDistribTool/dt0.plで変換すればいい。
!11: 第n隣接分子までの距離。
!12: q6local(5A/3A)
!13: pair + ohlen(@PAOH)
subroutine traj2(stream, outputtype, rarg)
  use common_module
  use rigid_module
  use flex_module
  use grid_module
  use mol_module
  use box_module
  use site_module
  use monitor_module
  use distrib_module
  use standard_interaction_module
  use rdf_module
  use cutoff_module
  use rigid_setup_module
  use flex_setup_module
  use nvde_module
  use spce_module
  use tip4p_module
  use interaction2_module
  use error_module
  use neighbors_module
  use tetrahedrality_module
  use random_module
  use q6local2

  implicit none
  integer, intent(IN) :: stream
  integer, intent(IN) :: outputtype
  real(kind=8), intent(IN) :: rarg

  type(sInteraction)  :: ww(MAXCOMBI)
  type(sAddress)      :: addr(MAXCOMPO)
  type(sRigid)        :: rigid( MAXCOMPO )
  type(sFlex)         :: flex( MAXCOMPO )
  type(sMol)          :: mol( MAXCOMPO )
  type(sSite)         :: site
  type(sBox)          :: box
  type(sGrid)         :: grid
  type(sRDF)          :: rdf(MAXCOMBI)
  type(sCutOff)       :: co
  type(sStdInt)       :: si(MAXCOMPO)
  type(sNeighbors)    :: nei(MAXMOL)   ! for tetrahedrality
  type(sNeighbors)    :: nei3(MAXMOL)  ! for dihedral angle
  type(sHistogram), pointer :: hist !three-body-angle
  type(pRandom), pointer :: random
  type(sQ6System), pointer :: q6s
  type(vector3)       :: cell
#ifdef MEMUSAGE
  type(sMonitor_integer)      :: mon(12)
#endif
  !type(sAdjacency),allocatable :: hbond(:,:)
  !
  !temporary neighbor list
  !
  type(sNeighbors) :: newnei(17)
  integer          :: newneiidx(17)
  integer          :: nmodified
  !
  !三体角の余弦
  !
  real(kind=8)     :: costh
  !
  !柔軟分子の質点の質量
  !
  real(kind=8) :: molmass(MAXSITE)
  type(vector3) :: com( MAXMOL, MAXCOMPO )
  character(len=5) :: tag
  character(len=8) :: lastid
  integer          :: ncompo, compo, incompo, ncombi
  integer          :: mdvwsite
  integer          :: i,j,k,l,ii,jj,kk,ll
  integer          :: nmol
  integer          :: combi( MAXCOMBI ), combj( MAXCOMBI )
  real(kind=8)     :: cutoff
  integer          :: nloop
  real(kind=8), allocatable :: ep1(:), ohlen(:)
  integer          :: pair
  real(kind=8)     :: dx,dy,dz,ox,oy,oz
  real(kind=8)     :: dx1,dy1,dz1,dx2,dy2,dz2, dd1, dd2
  logical          :: isWater
  integer          :: err
  incompo = 0
  nloop = 0

  call new(box)
  call new(random,5678)

#ifdef MEMUSAGE
  call new( mon(1), MON_MAXIMUM, MAXGRIDX,&
       MON_HALT, "MAXGRIDX                  ", 0 )
  call new( mon(2), MON_MAXIMUM, MAXGRIDY,&
       MON_HALT, "MAXGRIDY                  ", 0 )
  call new( mon(3), MON_MAXIMUM, MAXGRIDZ,&
       MON_HALT, "MAXGRIDZ                  ", 0 )
  call new( mon(4), MON_MAXIMUM, MAXINGRID,&
       MON_HALT, "MAXINGRID                 ", 0 )
  call new( mon(5), MON_MAXIMUM, MAXNEIBORCELL,&
       MON_HALT, "MAXNEIBORCELL             ", 0 )
  call new( mon(6), MON_MAXIMUM, MAXMOL,&
       MON_HALT, "MAXMOL                    ", 0 )
  call new( mon(7), MON_MAXIMUM, MAXsite_total,&
       MON_HALT, "MAXSITE                   ", 0 )
  call new( mon(8), MON_MAXIMUM, MAXPAIR,&
       MON_HALT, "MAXPAIR                   ", 0 )
  call new( mon(9), MON_MAXIMUM, MAXCELL,&
       MON_HALT, "MAXCELL                   ", 0 )
  call new( mon(10), MON_MAXIMUM, MAXPARTNER,&
       MON_HALT, "MAXPARTNER                ", 0 )
#endif

  do
     ncompo  = 0
     call new(site, MAXsite_total)
     do
        read(STDIN,'(a5)',END=999) tag
        if(tag.eq."@ID08")then
           read(STDIN,*) lastid
           cycle
        endif
        if(tag.eq."@NCMP")then
           read(STDIN,*) incompo
           if(incompo.gt.MAXCOMPO)then
              call die( error_too_many_components )
           endif
           cycle
        endif
        if(tag.eq.'@NX4A' .or. tag.eq.'@WTG3' .or. tag.eq.'@AR3A' .or. tag .eq. '@APC1' )then
           if(lastid.eq."")then
              call die( error_anonymous_molecule )
           endif
           ncompo = ncompo + 1
           if(ncompo.gt.MAXCOMPO)then
              call die( error_too_many_components )
           endif
           mol(ncompo)%id = lastid
           err = Rigid_Setup( rigid(ncompo), lastid, .not. FIXED, mol(ncompo), si(ncompo) )
           if ( err .eq. error_id_not_found ) then
              err = Flex_Setup( lastid, .not. FIXED, mol(ncompo), si(ncompo), molmass )
              if ( err .ne. error_none ) then
                 write( STDERR, * ) "Unknown molecule type: ",lastid
                 call die( err, "traj2 1" )
              endif
           endif
           if ( .not. mol(ncompo)%isRigid ) then
              call Flex_ReadAR3A(flex(ncompo),mol(ncompo),site,.true.,STDIN)
              call Flex_SetMass(flex(ncompo),mol(ncompo),molmass)
           endif
           if ( mol(ncompo)%isRigid ) then
              call Rigid_ReadNX4A(rigid(ncompo),mol(ncompo),site, STDIN )
           endif
           !
           !指定された成分数を読みこんだら入力をうちきる。
           !
           if(ncompo.eq.incompo)then
              exit
           endif
           lastid=""
           cycle
        endif
        call Box_Loader(box,STDIN,tag)
        call CutOff_Loader(co,STDIN,tag)
     enddo
999  continue

     write(STDERR,*) "Number of components : ", ncompo
     write(STDERR,*) "Output type          : ", outputtype
     write(STDERR,*) "Cutoff length        : ", co%out

     if ( ncompo == 0 ) goto 9999

     nloop = nloop + 1

     if ( outputtype .eq. 1 .or. outputtype .eq. 2 .or. &
          outputtype .eq. 4 .or. outputtype .eq. 5 .or. &
          outputtype .eq. 6 .or. outputtype .eq. 7 .or. &
          outputtype .eq. 8 .or. outputtype .eq. 9 .or. &
          outputtype .eq.10 .or. outputtype == 11 .or.  &
          outputtype == 12  .or. outputtype == 13       &
          ) then
        cutoff = 4d0
        if ( outputtype .eq. 4 ) then
           cutoff = 16d0
        endif
        if ( outputtype .eq. 5 .or. outputtype== 9 .or. outputtype == 13 ) then
           cutoff = co%out
        endif
        if ( outputtype .eq. 6 .or. outputtype .eq. 7 .or. outputtype .eq. 10 .or. &
             outputtype == 11 .or. outputtype == 12 ) then
           cutoff = 6d0
        endif
        if ( active( box ) ) then
           if ( cutoff .eq. 0d0 ) then
              write(STDERR,*) "Cutoff length must be specified."
              stop
           endif
           if ( box%size%vec(1)*0.5d0 < cutoff ) cutoff = box%size%vec(1)*0.5d0
           if ( box%size%vec(1)*0.5d0 < cutoff ) cutoff = box%size%vec(1)*0.5d0
           if ( box%size%vec(1)*0.5d0 < cutoff ) cutoff = box%size%vec(1)*0.5d0
        endif
        !
        !同種分子間の相互作用(分子数をnとするとn(n+1)/2対)
        !
        do compo=1,ncompo
           !
           !第2パラメータを省略した場合は自己相互作用とみなす。
           !
           !DECではoptionalの扱いがうまくいかないらしいので明示的に-1を与える。
           !
           call interaction_initialize(ww(compo),mol(compo)%nmol,-1)
           combi( compo ) = compo
           combj( compo ) = compo
        enddo
        !
        !異種分子間の相互作用(分子数をn,mとするとnm対)
        !
        ncombi = ncompo
        do i=1,ncompo
           do j=i+1,ncompo
              ncombi = ncombi + 1
              combi( ncombi ) = i
              combj( ncombi ) = j
              call interaction_initialize(ww(ncombi),mol(i)%nmol,mol(j)%nmol)
           enddo
        enddo

        if ( outputtype .eq. 4 ) then
           !rdfの初期化
           if ( .not. associated(rdf(1)%hist) ) then
              do k=1,ncombi
                 i = combi( k )
                 j = combj( k )
                 call rdf_initialize( rdf(k), mol(i)%nsite, mol(i)%name, mol(j)%nsite, mol(j)%name, nint(cutoff/0.05d0), 0.05d0 )
              enddo
           endif
        endif

        call site_resetforce(site)
        do compo=1,ncompo
           if ( mol(compo)%isRigid ) then
              !
              !剛体の位置と向きから、サイトの位置を決める。
              !
              call rigid_setsiteposition(rigid(compo),mol(compo),site)
           else
              !
              !柔軟分子のサイト位置から、重心の位置を決める。
              !
              call flex_setcom( flex( compo ), mol( compo ), site )
           endif
        enddo
        !
        !重心座標だけを直列化
        !(重心を別配列で計算するのであれば、TIP4Pに第5のサイトを与える必要はないのでは？)
        !
        do compo=1, ncompo
           nmol = mol(compo)%nmol
           if ( mol(compo)%isRigid ) then
              com( 1:nmol, compo ) = rigid( compo )%mol(1:nmol)%com
           else
              do i=1,nmol
                 k = mol( compo )%offset + i*mol( compo )%nsite
                 com( i, compo )%vec(1) = site%x(k)
                 com( i, compo )%vec(2) = site%y(k)
                 com( i, compo )%vec(3) = site%z(k)
              enddo
           endif
        enddo
        do i=1,ncompo
           call Address_Initialize(addr(i),mol(i)%nmol)
        enddo
        if ( cutoff .ne. 0d0 ) then
           call Grid_SetAdaptive(grid, cutoff)
           call Grid_Update( grid, box, 0d0 )
        endif
        if ( cutoff .ne. 0d0 .and. useful( grid ) ) then
           do i=1,ncombi
              call Interaction_Reset(ww(i))
           enddo
           !
           !分子をセルに帰属させる。
           !
           do i=1,ncompo
              call Address_Assign(addr(i),box,grid,mol(i)%nmol,com(1,i))
           enddo
           !
           !同種分子間相互作用対リストの作成
           !
           do i=1,ncompo
              call Grid_NeighborList(grid,grid%homo,ww(i),addr(i),addr(i))
           enddo
           !
           !異種分子間相互作用対リストの作成
           !
           do i=ncompo+1,ncombi
              call Grid_NeighborList(grid,grid%hetero,ww(i),&
                   addr(combi(i)),addr(combj(i)))
           enddo
           !do j=1,ncombi
           !   call interaction_alltoall(ww(j))
           !enddo
           do k=1,ncombi
              i=combi(k)
              j=combj(k)
              call Interaction_Compress2(ww(k),box%mode.ne.BOX_NONE,&
                   box ,cutoff**2, com(1,i),com(1,j),.false.)
           enddo
        else ! useful(grid)
           write(STDERR,*) "GRID IS NOT USEFUL."
           do k=1,ncombi
              i=combi(k)
              j=combj(k)
              call interaction_alltoall(ww(k))
              !
              !カットオフ外の対を除いて、相互作用対リストを圧縮する。
              !相互作用対リストを圧縮して、あとの力計算処理を高速化する
              !(ベクトル機では有効)
              !
              call interaction_compress( ww(k), com(1,i), com(1,j), .false., cutoff, box )
           enddo
        endif
        if ( outputtype .eq. 1 .or. outputtype .eq. 2 ) then
           do i=1,ncompo
              !
              !水-水の場合のみ、水素結合グラフを出力できる。他の場合は未定。
              !
              isWater = .false.
              if ( mol(i)%id(1:8) .eq. "TIP4P   " ) isWater = .true.
              if ( mol(i)%id(1:8) .eq. "ST2_____" ) isWater = .true.
              if ( mol(i)%id(1:8) .eq. NvdE_ID08  ) isWater = .true.
              if ( mol(i)%id(1:8) .eq. SPCE_ID08  ) isWater = .true.
              if ( isWater ) then
                 call tip4p_SaveGraph(ww(i),mol(i),site,stream,outputtype)
              endif
           enddo
        elseif ( outputtype .eq. 4 ) then
           do k=1,ncombi
              i = combi( k )
              j = combj( k )
              call RadialDistrib(ww(k),mol(i),mol(j),site,rdf(k))
           enddo
        elseif ( outputtype .eq. 5 ) then
           write(STDOUT,'("@NCMP")')
           write(STDOUT,*) ncompo
           do k=1,ncombi
              i = combi( k )
              j = combj( k )
              call CutOff_Dimmer2( co, ww(k), com(1,i), com(1,j) )
              allocate( ep1(ww(k)%npair) )
              allocate( ohlen(ww(k)%npair) )
              call pairpotential(ww(k),mol(i),mol(j),site,ep1,.true.,si(i),si(j) )
              !成分番号(0から始まる)
              write(STDOUT,'("@CMP2")')
              write(STDOUT,*) i-1,j-1
              !分子対の番号と相互作用(J/mol)
              write(STDOUT,'("@PAIR")')
              write(STDOUT,*) mol(i)%nmol, mol(j)%nmol
              do pair=1,ww(k)%npair
                 ii = ww(k)%pair_i(pair)
                 jj = ww(k)%pair_j(pair)
                 dx = com(ii,i)%vec(1) - com(jj,j)%vec(1)
                 dy = com(ii,i)%vec(2) - com(jj,j)%vec(2)
                 dz = com(ii,i)%vec(3) - com(jj,j)%vec(3)
                 ox = ww(k)%ox(pair)
                 oy = ww(k)%oy(pair)
                 oz = ww(k)%oz(pair)
                 dx = dx-ox
                 dy = dy-oy
                 dz = dz-oz
                 write(STDOUT,'(2(i5,1x),2(e24.17,1x))') ii-1,jj-1,ep1(pair) * I2J, sqrt(dx**2 + dy**2 + dz**2)
                 !write(STDOUT,*) pair,ep1(pair) * I2J
              enddo
              write(STDOUT,*) -1,-1,0.0d0,0d0
              deallocate( ep1 )
           enddo
        elseif ( outputtype .eq. 13 ) then
           write(STDOUT,'("@NCMP")')
           write(STDOUT,*) ncompo
           do k=1,ncombi
              i = combi( k )
              j = combj( k )
              call CutOff_Dimmer2( co, ww(k), com(1,i), com(1,j) )
              allocate( ep1(ww(k)%npair) )
              allocate( ohlen(ww(k)%npair) )
              call pairpotential(ww(k),mol(i),mol(j),site,ep1,.true.,si(i),si(j) )
              call stdohlength(ww(k),mol(i),mol(j),si(i),si(j),site,ohlen )
              !成分番号(0から始まる)
              write(STDOUT,'("@CMP2")')
              write(STDOUT,*) i-1,j-1
              !分子対の番号と相互作用(J/mol)
              write(STDOUT,'("@PAOH")')
              write(STDOUT,*) mol(i)%nmol, mol(j)%nmol
              do pair=1,ww(k)%npair
                 ii = ww(k)%pair_i(pair)
                 jj = ww(k)%pair_j(pair)
                 dx = com(ii,i)%vec(1) - com(jj,j)%vec(1)
                 dy = com(ii,i)%vec(2) - com(jj,j)%vec(2)
                 dz = com(ii,i)%vec(3) - com(jj,j)%vec(3)
                 ox = ww(k)%ox(pair)
                 oy = ww(k)%oy(pair)
                 oz = ww(k)%oz(pair)
                 dx = dx-ox
                 dy = dy-oy
                 dz = dz-oz
                 write(STDOUT,'(2(i5,1x),3(e24.17,1x))') ii-1,jj-1,ep1(pair) * I2J, sqrt(dx**2 + dy**2 + dz**2), ohlen(pair)
                 !write(STDOUT,*) pair,ep1(pair) * I2J
              enddo
              write(STDOUT,*) -1,-1,0.0d0,0d0,0d0
              deallocate( ep1 )
              deallocate( ohlen )
           enddo
        else if ( outputtype .eq. 6 ) then
           !
           !ある分子の近傍10分子を記録しておく。
           !1成分のみ。
           !traj2 -p | near.plと同じ結果を得た。(平成17年2月3日(木))
           !
           call neighbors_all2( mol(1)%nmol, nei, ww(1), com(1,1), 10, ORDER_BY_DISTANCE )
           call writetag( stream, "@NEAR" )
           write(stream,*) mol(1)%nmol, 10
           do i=1, mol(1)%nmol
              write(stream, '(10i5)') ( nei(i)%near(j)-1, j=1,10 )
           enddo
        else if ( outputtype .eq. 7 ) then
           !
           !ある分子の近傍4分子の成す角度の分布を求める。
           !1成分のみ。
           !
           call new( hist, 21, 0.1d0, -1d0 )
           call neighbors_all2( mol(1)%nmol, nei, ww(1), com(1,1), 4, ORDER_BY_DISTANCE )
           do i=1, mol(1)%nmol
              do j=1, nei(i)%nnear
                 do k=j+1, nei(i)%nnear
                    costh = nei(i)%dx(j)*nei(i)%dx(k) + nei(i)%dy(j)*nei(i)%dy(k) + nei(i)%dz(j)*nei(i)%dz(k)
                    costh = costh / sqrt ( nei(i)%dd(j) * nei(i)%dd(k) )
                    call histogram_accumulate( hist, costh, 1d0 )
                 enddo
              enddo
           enddo
           call histogram_save( hist, "@HST2", stream )
           call histogram_done( hist )
        else if ( outputtype .eq. 8 ) then
           !
           !ある分子の近傍4分子の成す角度の分布を求める。
           !1成分のみ。
           !
           call histogram_initialize( hist, 21, 0.1d0, -1d0 )
           !allocate( hbond(mol(1)%nmol, mol(1)%nmol) )
           !hbond(:,:)%dd = 0d0
           call neighbors_all2( mol(1)%nmol, nei3, ww(1), com(1,1), 0, ORDER_NONE, 3d0 )
           do i=1, mol(1)%nmol
              do jj=1, nei3(i)%nnear
                 j = nei3(i)%near(jj)
                 do kk=1, nei3(j)%nnear
                    k = nei3(j)%near(kk)
                    if ( k .ne. i ) then
                       dx1 = nei3(i)%dy(jj) * nei3(j)%dz(kk) - nei3(i)%dz(jj) * nei3(j)%dy(kk)
                       dy1 = nei3(i)%dz(jj) * nei3(j)%dx(kk) - nei3(i)%dx(jj) * nei3(j)%dz(kk)
                       dz1 = nei3(i)%dx(jj) * nei3(j)%dy(kk) - nei3(i)%dy(jj) * nei3(j)%dx(kk)
                       dd1 = dx1**2 + dy1**2 + dz1**2
                       do ll=1, nei3(k)%nnear
                          l = nei3(k)%near(ll)
                          if ( l .ne. j .and. l < i ) then
                             dx2 = nei3(j)%dy(kk) * nei3(k)%dz(ll) - nei3(j)%dz(kk) * nei3(k)%dy(ll)
                             dy2 = nei3(j)%dz(kk) * nei3(k)%dx(ll) - nei3(j)%dx(kk) * nei3(k)%dz(ll)
                             dz2 = nei3(j)%dx(kk) * nei3(k)%dy(ll) - nei3(j)%dy(kk) * nei3(k)%dx(ll)
                             dd2 = dx2**2 + dy2**2 + dz2**2
                             costh = dx1*dx2 + dy1*dy2 + dz1*dz2
                             costh = costh / sqrt ( dd1 * dd2 )
                             !write(STDERR,*) i,j,k,l,costh
                             call histogram_accumulate( hist, costh, 1d0 )
                          endif
                       enddo
                    endif
                 enddo
              enddo
           enddo
           call histogram_save( hist, "@HST2", stream )
           call histogram_done( hist )
        elseif ( outputtype .eq. 9 ) then
           call new( hist, 61, 0.1d0, 0d0 )
           do k=1,ncombi
              i = combi( k )
              j = combj( k )
              do pair=1,ww(k)%npair
                 ii = ww(k)%pair_i(pair)
                 jj = ww(k)%pair_j(pair)
                 dx = com(ii,i)%vec(1) - com(jj,j)%vec(1)
                 dy = com(ii,i)%vec(2) - com(jj,j)%vec(2)
                 dz = com(ii,i)%vec(3) - com(jj,j)%vec(3)
                 ox = ww(k)%ox(pair)
                 oy = ww(k)%oy(pair)
                 oz = ww(k)%oz(pair)
                 dx = dx-ox
                 dy = dy-oy
                 dz = dz-oz
                 call histogram_accumulate( hist, sqrt(dx**2 + dy**2 + dz**2), 1d0 )
              enddo
           enddo
           call histogram_save( hist, "@HST2", stream )
           call histogram_done( hist )
        else if ( outputtype .eq. 10 ) then
           !
           !四面体度。
           !
           call writetag( stream, "@TETR" )
           write( stream,* ) mol(1)%nmol
           call neighbors_all2( mol(1)%nmol, nei, ww(1), com(1,1), 4, ORDER_BY_DISTANCE )
           do i=1, mol(1)%nmol
              write( stream,* ) tetrahedrality( nei(i)%dx, nei(i)%dy, nei(i)%dz, nei(i)%dd )
           enddo
        else if ( outputtype .eq. 11 ) then
           !
           !ある分子の第n近傍分子までの距離
           !
           j = rarg
           call neighbors_all2( mol(1)%nmol, nei, ww(1), com(1,1), j, ORDER_BY_DISTANCE )
           call writetag( stream, "@NEAX" )
           write(stream,*) mol(1)%nmol, j
           do i=1, mol(1)%nmol
              write(stream, fmt0 ) sqrt(nei(i)%dd(j))
           enddo
        else if ( outputtype .eq. 12 ) then
           !
           !Q6Local
           !
           call new( q6s )
           call sQ6System_prepareall( q6s, mol(1)%nmol, ww(1), com(1,1) )
           call writetag( stream, "@Q6P6" );
           write(6,*) mol(1)%nmol, 5d0, 3d0
           do i=1, mol(1)%nmol
              write(6,*) q6s%q6(i)
           enddo
           call done( q6s )
        endif
     endif
     if ( outputtype .eq. 3 ) then
        call site_resetforce(site)
        write(STDERR,*) "type 3"
        do compo=1,ncompo
           if ( mol(compo)%isRigid ) then
              do i=1,mol(compo)%nmol
                 call box_renormalize( box, rigid( compo )%mol(i)%com, cell )
              enddo
              !
              !剛体の位置と向きから、サイトの位置を決める。
              !
              call rigid_setsiteposition(rigid(compo),mol(compo),site)
           else
              !
              !柔軟分子のサイト位置から、重心の位置を決める。
              !
              call flex_setcom( flex( compo ), mol( compo ), site )
           endif
        enddo
        !call save(box,stream, BOX_BOX3)
        !write(stream,'("@MDVW")')
        call save(box,stream, BOX_MDVW)
        mdvwsite = 0
        do compo=1,ncompo
           mdvwsite = mdvwsite + Mol_CountSiteMDVW( mol(compo) )
        enddo
        write(stream,*) mdvwsite
        do compo=1,ncompo
           if ( mol(compo)%isRigid ) then
              call save( rigid(compo), mol(compo), stream, RIGID_MDVW )
           else
              call save( flex(compo), mol(compo), site, stream, FLEX_MDVW)
           endif
        enddo
     endif
     !
     !Reset site list
     !
     call new(site, MAXSITE_TOTAL)
  enddo
9999 continue
  if ( outputtype .eq. 4 ) then
     do k=1,ncombi
        i = combi( k )
        j = combj( k )
        if ( i .eq. j ) then
           call rdf_normalize( rdf(k), nloop, box, mol(i) )
        else
           call rdf_normalize( rdf(k), nloop, box, mol(i), mol(j) )
        endif
        !call rdf_show( rdf(k) )
     enddo
  endif
end subroutine traj2

! -*- f90 -*-
#define DEBUG

!ToDo
!Bookkeeping&Grid to accelerate calculations
!Q6 implementation
!Verification run
!
#define OP

program mctest
  use common_module
  use site_module
  use rigid_module
  use box_module
  use error_module
  use cutoff_module
  use standard_interaction_module
  use vector_module
  use book_module
  use grid_module
  use tetrahedrality_module

  use oplsmeoh_module
  use tip4p_module
  use tip5p_module
  use pairinfo_module
  use quat_module
  use interaction_module
  use mol_module
  use interaction2_module
  use monitor_module
#ifdef OP
  use neighbors_module
  use stack_module
#endif
  implicit none 
#ifdef OP
  type(sNeighbors) :: nei(MAXMOL)  ! for tetrahedrality
  type(sNeighbors) :: nei3(MAXMOL) ! for Q6
  type(sNeighbors) :: nei6(MAXMOL) ! for Q6
#endif
  type(sBox)  :: box
  type(sCutoff) :: co
  type(vector3) :: com( MAXMOL, MAXCOMPO ), comt(1)
  integer :: combi( MAXCOMBI ), combj( MAXCOMBI )
  character(len=5) :: tag
  character(len=8) :: lastid

  type(sInteraction)  :: wwt
  type(sInteraction)  :: ww( MAXCOMBI )
  type(sSite) :: site
  type(sMol) :: mol( MAXCOMPO ),molt
  type(sRigid) :: rigid( MAXCOMPO ),rigidt
  type(sStdInt) :: si(MAXCOMPO),sit
  type(sPairInfo)    :: pairinfo(MAXCOMBI),pairtest(MAXCOMBI)
  type(sBook) :: book
  type(sGrid) :: grid
  type(sAddress) :: addr(MAXCOMPO)
  integer          :: i,j,k,loop,nloop
  real(kind=8)     :: epsum,oldsum
  real(kind=8), allocatable :: random(:)
  real(kind=8), allocatable :: delta(:)
  !
  !@NCMPで明示的に指定される成分数
  !
  integer :: incompo
  !
  !成分の数(とりあえず上限は2)
  !
  integer      :: ncompo
  !
  !組み合わせの数
  !
  integer      :: ncombi
  integer :: dof
  !
  !分子の総数
  !
  integer      :: nmol
  !
  !エネルギー差分
  !
  real(kind=8) :: deltae
  !
  !微小回転の大きさ
  !
  real(kind=8) :: deltar
  !
  !微小回転前の四元子
  !
  real(kind=8) :: olda,oldb,oldc,oldd
  !
  !微小回転後の四元子
  !
  real(kind=8) :: qa,qb,qc,qd
  !
  !各軸回りの微小回転角度
  !
  real(kind=8) :: rotx,roty,rotz
  !
  !内側ループ変数
  !
  integer      :: trial
  !
  !動かす分子の番号
  !
  integer      :: target
  !
  !微小並進量
  !
  real(kind=8) :: deltat, deltax, deltay, deltaz
  !
  !default radius to cut off interaction
  !
  real(kind=8) :: rc
  !
  !pair energy after trial
  !
  real(kind=8),dimension(:),allocatable :: ep !,vr
  !
  !true if accepted
  !
  logical      :: accept
  !
  !beta in 1/K
  !
  real(kind=8) :: beta
  !
  !metropolis's ratio
  !
  real(kind=8) :: ratio
  !
  !number of trials in a loop ( usually set to twice of number of
  !molecules )
  !
  integer      :: ntrial
  !
  !Tetrahedrality and their average
  !
  real(kind=8) :: tetra(MAXMOL), sumtetra,newtetra
#ifdef OP
  type(indexedreal8stack) :: stack
  !
  !temporary neighbor list
  !
  type(sNeighbors) :: newnei(17)
  integer          :: newneiidx(17)
  integer          :: nmodified
  type(sNeighbors) :: newnei3(99)
  integer          :: newnei3idx(99)
  integer          :: nei3modified
  type(sNeighbors) :: newnei6(99)
  integer          :: newnei6idx(99)
  integer          :: nei6modified
  type(sAdjacency),allocatable :: hbond(:,:)
  real(kind=8)     :: q6( MAXMOL ), newq6, sumq6
  type(real8stack) :: relx, rely, relz
  type(indexedreal8stack) :: q6stack, r6a
#endif
  integer          :: ii,jj,kk
  external qx
  real(kind=8)     :: qx
  !for debug
  real(kind=8) :: px,py,pz,pp

  !
  !interval between saving snapshots
  !
  integer          :: snpo

  integer      :: seed
  integer      :: compo

  !write(6,*) box%mode
  !call new(nose)
  call new(site, MAXsite_total)
  !call new(an)
  !call new(av)
  co%mode = 0
  box%mode = 0
  write(STDERR,*) co%mode, box%mode
  call new(book)
  call new(grid)
  call new(box)
  call new(co)

  beta = 1d0/400d0
  dof=0d0
  ncompo=0
  incompo=0
  seed = 12345
  lastid="        "
  snpo = 0
  !read config
  do
     read(STDIN,'(a5)',END=999) tag
     write(STDERR,*) tag
     if(tag.eq."@SNPO")then
        read(STDIN,*) snpo
        cycle
     endif
     if(tag.eq."@MCLP")then
        read(STDIN,*) nloop
        cycle
     endif
     if(tag.eq."@ID08")then
        read(STDIN,*) lastid
        cycle
     endif
     if(tag.eq."@NCMP")then
        read(STDIN,*) incompo
        if(incompo.gt.MAXCOMPO)then
           call die( error_too_many_components, "mctest 1" )
        endif
        cycle
     endif
     if(tag.eq.'@NX4A')then
        if(lastid.eq."        ")then
           write(STDERR,*) "Error: Anonymous data"
           stop
        endif
        ncompo = ncompo + 1
        if(ncompo.gt.MAXCOMPO)then
           write(STDERR,*) "Number of componets exceeded the&
                & limit."
           stop
        endif
        mol(ncompo)%id=lastid
        if(lastid.eq."TIP4P   ")then
           call new(mol(ncompo),5,6, RIGID_MODE, .not. FIXED )
           call Rigid_TIP4P_Constructor(rigid(ncompo))
           call TIP4P_SetInteraction(si(ncompo))
        else if(lastid.eq."TIP5P   ")then
           call new(mol(ncompo),6,6, RIGID_MODE, .not. FIXED )
           call Rigid_TIP5P_Constructor(rigid(ncompo))
           call TIP5P_SetInteraction(si(ncompo))
        else if (lastid.eq."OPLSMEOH")then
           call new(mol(ncompo),4,6, RIGID_MODE, .not. FIXED )
           call Rigid_OPLSMeOH_Constructor(rigid(ncompo))
           call OPLSMeOH_SetInteraction(si(ncompo))
        else
           write(STDERR,*) "Unknown rigid body: ",lastid
           stop
        endif
        call Rigid_ReadNX4A(rigid(ncompo),mol(ncompo),site&
             & ,STDIN)
        dof = dof + Mol_DoF(mol(ncompo))
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
     call Cutoff_Loader(co,STDIN,tag)
     call Grid_Loader(grid,STDIN,tag)
     call Book_Loader(book,STDIN,tag)
  enddo
999 continue
#ifdef OP
  allocate( hbond(mol(1)%nmol, mol(1)%nmol) )
#endif
  do i=1,mol(1)%nmol
     call qnormalize( &
          rigid(1)%mol(i)%quat(1)%vec(1),&
          rigid(1)%mol(i)%quat(1)%vec(2),&
          rigid(1)%mol(i)%quat(1)%vec(3),&
          rigid(1)%mol(i)%quat(1)%vec(4) )
  enddo
  if( grid%mode .ne. GRID_NONE)then
     if(.not.active(book))then
        call Book_Initialize(book, BOOK_FIX, 0d0, 1)
     endif
  endif
  !
  !箱の大きさが指定されていてCutOffが指定されていない場合
  !
  if(box%mode.ne.BOX_NONE.and.co%mode.eq.CUTOFF_NONE)then
     co%mode=CUTOFF_POLY
     rc=dmin1(box%size%vec(1),box%size%vec(2),box%size%vec(3))*0.5d0
     call cutoff_initialize(co,rc-2d0,rc)
  endif
  !
  !いちおう変な設定をチェックしておく。
  !
  if( active(book) .and. grid%mode .eq. GRID_NONE)then
     write(STDERR,*) "ERROR: GRID IS NECESSARY FOR BOOKING."
     call die( 0, "mctest 2" )
  endif
  if( active(book) .and. co%mode .eq. CUTOFF_NONE)then
     write(STDERR,*) "ERROR: CUTOFF IS NECESSARY FOR BOOKING."
     call die( 0, "mctest 3" )
  endif
  !
  !同種分子間の相互作用(分子数をnとするとn(n+1)/2対)
  !
  do j=1,ncompo
     !
     !第2パラメータを省略した場合は自己相互作用とみなす。
     !
     !DECではoptionalの扱いがうまくいかないらしいので明示的に-1を与える。
     !
     call interaction_initialize(ww(j),mol(j)%nmol,-1)
     combi( j ) = j
     combj( j ) = j
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
  write(*,*) "j",ncombi,ncompo
  !
  !系が小さい場合を想定;全分子同士が相互作用する(グリッド分割しない)
  !
  if( active(book) )then
     if( book%margin .le. 0d0 )then
        book%mode = BOOK_FIX
     endif
     if(book%mode .eq. BOOK_AUTOINTERVAL)then
        call die( 0, "Main 27" )
        !
        !まだ考えていない。帳簿法はほとんど使わないから無視しても構わ
        !ないのだが・・・!!!
        !
        !com0(1:water%nmol)=rwater%mol(1:water%nmol)%com
        !call Book_ResetInterval(bk,co%Out+bk%margin,water%nmol,com0)
     endif
     !
     !初期エネルギー計算のために、とりあえず挿入(あとで除去しろ)
     !
     do j=1,ncombi
        call interaction_alltoall(ww(j))
     enddo
  else
     !
     !系が小さい場合を想定;全分子同士が相互作用する(グリッド分割しない)
     !
     do j=1,ncombi
        call interaction_alltoall(ww(j))
     enddo
  endif
  if( grid%mode.ne.GRID_NONE )then
     do i=1,ncompo
        call Address_Initialize( addr(i),mol(i)%nmol )
     enddo
  endif
  !
  !make interaction table
  !
  do j=1, ncompo
     com( 1:mol(j)%nmol, j ) = rigid( j )%mol(1:mol(j)%nmol)%com
  enddo
  !
  !set all site position
  !
  do j=1,ncompo
     call rigid_setsiteposition(rigid(j),mol(j),site)
  enddo
  pairinfo(1)%ep(:,:) = 0d0
  !pairinfo(1)%vr(:,:) = 0d0
  epsum = 0d0
  do k=1,ncombi
     i = combi( k )
     j = combj( k )
     if ( co%mode.ne.CUTOFF_NONE) then
        call Interaction_Compress(ww(k),com(1,i),com(1,j),&
             .false. ,CutOff_OuterLimit(co),box)
     else
        call Interaction_Compress(ww(k),com(1,i),com(1,j),&
             .false. ,-1d0,box)
     endif
     !
     !Smooth cutoffによる減衰係数を事前計算
     !
     !call CutOff_Dimmer(co, ww(k), site, mol(i), mol(j) )
     call CutOff_Dimmer2( co, ww(k), com(1,i), com(1,j) )
     !
     !LJ+Coulombの標準相互作用の場合
     !
     call pairenergy(ww(k),mol(i),mol(j),site,pairinfo(k), .false. ,si(i),si(j))
     epsum = epsum + pairinfo(k)%epsum
  enddo
  !
  !Count total number of molecules
  !
  nmol = 0
  do j=1,ncompo
     nmol = nmol + mol(j)%nmol
  enddo
  !
  !トータルエネルギーが./Mixtureの結果と微妙に違う場合がある。
  !何か初期化を忘れているかもしれぬ。注意を払うべし。
  !
  write(STDOUT,*) epsum*I2J/(dble(nmol)*1000d0)
  !
  !Trial専用に設けた分子セットの初期化を行う。
  !TIP4Pを決めうちする。
  call new( molt,5,6, RIGID_MODE, .not. FIXED )
  call rigid_tip4p_constructor( rigidt )
  call tip4p_setinteraction( sit )
  call mol_allocate(molt, site, 1, RIGID_MODE )
  call interaction_initialize( wwt, mol(1)%nmol, molt%nmol )
  call interaction_alltoall( wwt )
  allocate( ep(wwt%npair0) )
  !allocate( vr(wwt%npair0) )
  allocate( delta(mol(1)%nmol) )
  delta(:) = 1d-1
  !
  !主loop
  !
  ntrial = 2*mol(1)%nmol
  allocate(random(ntrial*8))
  call mtrngi( seed )
  do loop=1,nloop
     !
     ![0,1)乱数
     !
     call mtrndv(random,ntrial*8)
     !
     !とりあえず、1MC stepに一回、グリッドを再設計する。
     !グリッドのサイズにはすこし余裕をみておく。
     !
     if ( grid%mode .ne. GRID_NONE ) then
        !
        !grid and addr are local in this block
        !
        !
        !Serialize the center-of-molecules.
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
        !
        !Decide the size of grid
        !
        call Grid_Update( grid, box, book%outer )
        !
        !Assign molecules to the cell
        !
        do i=1,ncompo
           call Address_Assign( addr(i), box, grid, mol(i)%nmol, com(1,i) )
        enddo
        !
        !Flush pairlist
        !
        do k=1,ncombi
           call Interaction_Reset(ww(k))
           i=combi(k)
           j=combj(k)
           !
           !分子間相互作用対リストの作成
           !
           if ( ww(k)%isomol ) then
              call Grid_NeighborList(grid,grid%homo,ww(k),addr(i),addr(j))
           else
              call Grid_NeighborList(grid,grid%hetero,ww(k),addr(i),addr(j))
           endif
           call Interaction_Compress2(ww(k), active( box ),&
                box ,co%out**2, com(1,i),com(1,j),.false.)
        enddo

        if( book%mode.eq.BOOK_AUTOINTERVAL )then
           !
           !帳簿のマージン半径内に分子がどのくらい侵入してきたかを
           !チェックして、次に帳簿を更新すべきタイミングを決定する。
           !
           !call Book_AssumeInterval( book, ww%nmol_i, com0, Book_BestMargin(book,ww%nmol_i,com0) )
        endif
        book%countdown=book%interval
     endif
#ifdef OP
     !
     !ある分子の近傍16分子を記録しておく。
     !
     do k=1,ncombi
        i=combi(k)
        j=combj(k)
        !modified from neighbors_all 2012-05-28
        call neighbors_all2( mol(i)%nmol, nei, ww(k), com(1,i), 16 )
        hbond(:,:)%dd = 0d0
        call neighbors_all2( mol(i)%nmol, nei3, ww(k), com(1,i), 0, 3d0, hbond )
        call neighbors_all2( mol(i)%nmol, nei6, ww(k), com(1,i), 0, 6d0 )
     enddo
     !
     !全分子のtetrahedralityを計算しておく
     !
     do compo=1, ncompo
        sumtetra = 0d0
        do i=1, mol(compo)%nmol
           !write(6,*) i,nei(i)%nnear,nei(i)%near(3),nei(i)%dx(3),nei(i)%dy(3),nei(i)%dz(3)
           newtetra = tetrahedrality( nei(i)%dx, nei(i)%dy, nei(i)%dz, nei(i)%dd )
           !px(1:4) = nei(i)%dx(1:4)
           !py(1:4) = nei(i)%dy(1:4)
           !pz(1:4) = nei(i)%dz(1:4)
           !pp(1:4) = 1d0/sqrt(nei(i)%dd(1:4))
           !px(1:4) = px(1:4) * pp(1:4)
           !py(1:4) = py(1:4) * pp(1:4)
           !pz(1:4) = pz(1:4) * pp(1:4)
           !write(6,*) "Q6TEST", qx( 6, 4, px, py, pz )
#ifdef DEBUG
           if ( 0.000001d0< dabs( ( newtetra - tetra(i))/newtetra ) ) write(6,*) "TETRA", i, tetra(i), newtetra
#endif
           tetra(i) = newtetra
           sumtetra = sumtetra + newtetra
           !write(6,*) i, sumtetra, newtetra
        enddo
        write(6,*) "tetra:", sumtetra / mol(compo)%nmol
     enddo
     !
     !全分子のq6tを計算しておく
     !
     do compo=1, 1
        sumq6 = 0d0
        !
        !for all moecules
        !
        do i=1, mol(compo)%nmol
           !
           !for all neighbor molecules of i within 6A ( including i itself )
           !
           call new( relx, 200 )
           call new( rely, 200 )
           call new( relz, 200 )
           nei6(i)%near( nei6(i)%nnear+1 ) = i
           !write(6,*) "Q6::", nei6(i)%nnear
           call new( r6a, nei6(i)%nnear+1 )
           do j=1, nei6(i)%nnear
              if ( nei6(i)%dd(j) < 6d0**2 ) then
                 call push( r6a, nei6(i)%near(j), 0d0 )
              endif
           enddo
           call push( r6a, i, 0d0 )
           !write(6,*) "<6A pair", r6a%depth
           do j=1, r6a%depth
              jj = r6a%index(j)
              do k=j+1, r6a%depth
                 kk = r6a%index(k)
                 !
                 !if they are bound
                 !
                 if ( hbond( jj, kk )%dd .ne. 0d0 .and. hbond( jj, kk )%dd .lt. 3d0**2 ) then
                    !write(6,*) i,j,k,jj,kk, hbond(jj,kk)%dd,relx%depth
                    !
                    !record relative position vectors
                    !
                    pp = 1d0 / sqrt( hbond( jj, kk )%dd )
                    px = hbond( jj, kk )%dx * pp
                    py = hbond( jj, kk )%dy * pp
                    pz = hbond( jj, kk )%dz * pp
                    call push( relx, px )
                    call push( relx,-px )
                    call push( rely, py )
                    call push( rely,-py )
                    call push( relz, pz )
                    call push( relz,-pz )
                 endif
              enddo
           enddo
           newq6 = qx( 6, relx%depth, relx%value, rely%value, relz%value )
#ifdef DEBUG
           !write(6,*) "Q6:", relx%depth, newq6
           if ( 0.000001d0< dabs( ( newq6 - q6(i))/newq6 ) ) write(6,*) "Q6", i, q6(i), newq6
#endif
           q6(i) = newq6
           sumq6 = sumq6 + q6(i)
        enddo
        write(6,*) "q6:", sumq6 / mol(compo)%nmol
     enddo
#endif

     

     do trial=0,ntrial-1
        !trial
        !
        !分子を選ぶ
        !
        target = random( trial*8 + 1 )*mol(1)%nmol + 1
        deltar = delta(target)
        deltat = delta(target)
        !write(6,*) "target",target
        rotx   = deltar * ( random( trial*8 + 2 ) - 0.5d0 )
        roty   = deltar * ( random( trial*8 + 3 ) - 0.5d0 )
        rotz   = deltar * ( random( trial*8 + 4 ) - 0.5d0 )
        call qinfrotator( rotx,roty,rotz, qa,qb,qc,qd )
        olda   = rigid(1)%mol( target )%quat(1)%vec(1)
        oldb   = rigid(1)%mol( target )%quat(1)%vec(2)
        oldc   = rigid(1)%mol( target )%quat(1)%vec(3)
        oldd   = rigid(1)%mol( target )%quat(1)%vec(4)
        call qadd( qa,qb,qc,qd, olda,oldb,oldc,oldd )
        !write(STDOUT,*) olda,oldb,oldc,oldd
        !write(STDOUT,*) qa,qb,qc,qd
        rigidt%mol( 1 )%quat( 1 )%vec(1) = qa
        rigidt%mol( 1 )%quat( 1 )%vec(2) = qb
        rigidt%mol( 1 )%quat( 1 )%vec(3) = qc
        rigidt%mol( 1 )%quat( 1 )%vec(4) = qd
        deltax = deltat * ( random( trial*8 + 5 ) - 0.5d0 )
        deltay = deltat * ( random( trial*8 + 6 ) - 0.5d0 )
        deltaz = deltat * ( random( trial*8 + 7 ) - 0.5d0 )
        comt(1)%vec(1) = deltax
        comt(1)%vec(2) = deltay
        comt(1)%vec(3) = deltaz
        comt(1)%vec(:) = comt(1)%vec(:) + rigid(1)%mol( target )%com%vec(:)
        rigidt%mol( 1 )%com = comt(1)
        !
        !com()のほうも新しい座標にしておく
        !(Tetrahedralityの計算に必要)
        !
        com( target, 1 ) = comt(1)
        !write(STDOUT,*) comt
        call rigid_setsiteposition(rigidt,molt,site)
        !
        !まず、rc+delta以内の隣接分子の表を作る。Compressが使えるはず。
        !本当はこの段階でtargetをpairlistから除外しておいた方がよい。
        !
        !
     !グリッド分割する場合は、グリッドの大きさを、相互作用距離+最大移
     !動距離よりも大きくとっておく。こうすれば、分子が移動した結果隣の
     !グリッドに移動してしまっても、移動前の隣接分子が突然見えなくなっ
     !てしまうことはない。
        !
        !グリッドを使わない場合も、半径(cutoff距離
        ! +移動距離)内の分子の一覧を作る。
        !
        if ( co%mode.ne.CUTOFF_NONE) then
           call Interaction_Compress(wwt,com(1,1),comt,&
                .false. ,CutOff_OuterLimit(co)+deltat,box)
        else
           call Interaction_Compress(wwt,com(1,1),comt,&
                .false. ,-1d0,box)
        endif
        !
        !Smooth cutoffによる減衰係数を事前計算
        !
        !call CutOff_Dimmer(co, ww(k), site, mol(i), mol(j) )
        call CutOff_Dimmer2( co, wwt, com(1,1), comt )
        !
        !LJ+Coulombの標準相互作用の場合。
        !
        call pairpotential(wwt,mol(1),molt,site,ep,.false., si(1), sit )
        !
        !Umbrella potentialのような外場はここで計算する。
        !
        !---------------------------------------------------------------------
#ifdef OP
        call neighbors_displacement( nei, target, deltax, deltay, deltaz, nmodified, newnei, newneiidx )
        call neighbors_displacement( nei3, target, deltax, deltay, deltaz, nei3modified, newnei3, newnei3idx )
        call neighbors_displacement( nei6, target, deltax, deltay, deltaz, nei6modified, newnei6, newnei6idx )
        !
        !hbondの内容は上書きする(rejectした場合は復旧必要)
        !
        do i=1, newnei3( 1 )%nnear
           j = newnei3( 1 )%near(i)
           ! target-->j
           hbond( target, j )%dx = newnei3(1)%dx(i)
           hbond( target, j )%dy = newnei3(1)%dy(i)
           hbond( target, j )%dz = newnei3(1)%dz(i)
           hbond( target, j )%dd = newnei3(1)%dd(i)
           ! j-->target
           hbond( j, target )%dx =-newnei3(1)%dx(i)
           hbond( j, target )%dy =-newnei3(1)%dy(i)
           hbond( j, target )%dz =-newnei3(1)%dz(i)
           hbond( j, target )%dd = newnei3(1)%dd(i)
        enddo
        call new( stack, nmodified )
        do i=1, nmodified
           call push( stack, newneiidx( i ), tetrahedrality( newnei( i )%dx, newnei( i )%dy, newnei( i )%dz, newnei( i )%dd ) )
        enddo

#ifdef DEBUG
        !write(6,*) "NMODIFIED", nModified
#endif

        !
        !q6 recalculation
        !差分計算にしないと、まともに再計算していたのでは遅すぎる。
        !
        call new( q6stack, mol(1)%nmol )
        do ii=1, nei6modified
           i = newnei6idx( ii )
           !
           !for all neighbor molecules of i within 6A ( including i itself )
           !
           newnei6(ii)%near( newnei6(ii)%nnear+1 ) = i
           call new( relx, 200 )
           call new( rely, 200 )
           call new( relz, 200 )
           call new( r6a, newnei6(ii)%nnear+1 )
           do j=1, newnei6(ii)%nnear
              if ( newnei6(ii)%dd(j) < 6d0**2 ) then
                 call push( r6a, newnei6(ii)%near(j), 0d0 )
              endif
           enddo
           call push( r6a, i, 0d0 )
           do j=1, r6a%depth
              jj = r6a%index(j)
              do k=j+1, r6a%depth
                 kk = r6a%index(k)
                 !
                 !if they are bound
                 !
                 if ( hbond( jj, kk )%dd .ne. 0d0 .and. hbond( jj, kk )%dd .lt. 3d0**2 ) then
                    !
                    !record relative position vectors
                    !
                    pp = 1d0 / sqrt( hbond( jj, kk )%dd )
                    px = hbond( jj, kk )%dx * pp
                    py = hbond( jj, kk )%dy * pp
                    pz = hbond( jj, kk )%dz * pp
                    call push( relx, px )
                    call push( relx,-px )
                    call push( rely, py )
                    call push( rely,-py )
                    call push( relz, pz )
                    call push( relz,-pz )
                 endif
              enddo
           enddo
           call push( q6stack, i, qx( 6, relx%depth, relx%value, rely%value, relz%value ) )
        enddo
#endif
        !
        !---------------------------------------------------------------------
        !
        deltae = 0d0
        do k=1,wwt%npair
           i = wwt%pair_i(k)
           if ( i /= target ) then
              deltae = deltae + ep(k) - pairinfo(1)%ep(i,target)
           endif
        enddo
        write(STDOUT,*) deltat, deltar, deltae*I2J*J2K
        !
        !新しい分子配置をうけいれるかどうかを判定する。
        !
        accept = .false.
        if ( deltae < 0d0 ) then
           accept = .true.
        else
           ratio = dexp(-( deltae * I2J * J2K )*beta)
           if ( random( trial*8 + 8 ) < ratio ) then
              accept = .true.
              !write(STDERR,*) "Ratio",ratio,deltae * I2J * J2K,"dE"
           endif
        endif
        !
        !うけいれた場合
        !
        if ( accept ) then
           write(6,*) "accept"
           !
           !座標を更新する。
           !
           rigid(1)%mol( target ) = rigidt%mol(1)
           !
           !com(), Site Positionも1分子分だけ更新する。
           !
           com( target, 1 ) = rigidt%mol(1)%com
           i=mol(1)%offset + (target-1)* mol(1)%nsite
           j=molt%offset
           do k=1,mol(1)%nsite
              site%x(i+k)=site%x(j+k)
              site%y(i+k)=site%y(j+k)
              site%z(i+k)=site%z(j+k)
           enddo
           !
           !相互作用表を更新する。
           !
           epsum = epsum + deltae
           pairinfo%epsum = epsum
           do k=1,wwt%npair
              i = wwt%pair_i(k)
              if ( i /= target ) then
                 pairinfo(1)%ep(i,target) = ep(k)
                 pairinfo(1)%ep(target,i) = ep(k)
              endif
           enddo
           !
           !移動量を更新する。
           !
           delta( target ) = delta( target ) * 1.05
#ifdef OP
           !
           !四面体度を更新する
           !
           do i=1, stack%depth
              j = stack%index(i)
              sumtetra = sumtetra - tetra(j) + stack%value(i)
              tetra(j) = stack%value(i)
           enddo
           !
           !q6を更新する
           !
           do i=1, q6stack%depth
              j = q6stack%index(i)
              sumq6 = sumq6 - q6(j) + q6stack%value(i)
              q6(j) = q6stack%value(i)
           enddo
           !
           !隣接表を更新する。
           !
           do i=1, nmodified
              j = newneiidx(i)
              nei(j) = newnei(i)
           enddo
           do i=1, nei3modified
              j = newnei3idx(i)
              nei3(j) = newnei3(i)
           enddo
           do i=1, nei6modified
              j = newnei6idx(i)
              nei6(j) = newnei6(i)
           enddo
#endif /*OP*/
        else
           write(6,*) "reject"
           !
           !com()の内容をもとにもどす。
           !
           com( target, 1 )%vec(:) = rigid(1)%mol( target )%com%vec(:)
           delta( target ) = delta( target ) * 0.95
#ifdef OP
           !
           !hbondの復旧
           !
           do i=1, nei3( target )%nnear
              j = nei3( target )%near(i)
              ! target-->j
              hbond( target, j )%dx = nei3( target )%dx(i)
              hbond( target, j )%dy = nei3( target )%dy(i)
              hbond( target, j )%dz = nei3( target )%dz(i)
              hbond( target, j )%dd = nei3( target )%dd(i)
              ! j-->target
              hbond( j, target )%dx =-nei3( target )%dx(i)
              hbond( j, target )%dy =-nei3( target )%dy(i)
              hbond( j, target )%dz =-nei3( target )%dz(i)
              hbond( j, target )%dd = nei3( target )%dd(i)
           enddo
#endif
        endif
     enddo
     !
     !calculate HB network
     !

     !
     !Save configuration
     !@SNPOで指定すべき。
     !
     if( 0 < snpo )then
        if(mod(loop,snpo).eq.0)then
           call Box_Save( box, STDOUT, BOX_BOX3 )
           call Rigid_Save( rigid(1), mol(1), STDOUT, RIGID_NX4A )
        endif
     endif
     !
     !確認のため、全エネルギーを計算しなおして照合する。
     !
     write(6,*) "epsum1",epsum
#ifdef DEBUG
     oldsum = epsum
     !write(6,*) "epsum1",epsum
     !
     !make interaction table
     !
     do j=1, ncompo
        !
        !ちゃんと動いていれば、そもそもcomの内容はrigidのcomの内容と一致するはず。
        !
        do i=1, mol(j)%nmol
           do k=1,3
              if ( com( i, j )%vec(k) /= rigid( j )%mol( i )%com%vec(k) ) then
                 write(STDERR,*) "COM", j,i,k
              endif
           enddo
        enddo
     enddo
     !
     !set all site position
     !
     do j=1,ncompo
        call rigid_setsiteposition(rigid(j),mol(j),site)
     enddo
     pairtest(1)%ep(:,:) = 0d0
     !pairtest(1)%vr(:,:) = 0d0
     epsum = 0d0
     do k=1,ncombi
        i = combi( k )
        j = combj( k )
        if ( co%mode.ne.CUTOFF_NONE) then
           call Interaction_Compress(ww(k),com(1,i),com(1,j),&
                .false. ,CutOff_OuterLimit(co),box)
        else
           call Interaction_Compress(ww(k),com(1,i),com(1,j),&
                .false. ,-1d0,box)
        endif
        !
        !Smooth cutoffによる減衰係数を事前計算
        !
        !call CutOff_Dimmer(co, ww(k), site, mol(i), mol(j) )
        call CutOff_Dimmer2( co, ww(k), com(1,i), com(1,j) )
        !
        !LJ+Coulombの標準相互作用の場合
        !
        call pairenergy(ww(k),mol(i),mol(j),site ,pairtest(k), .false., si(i),si(j))
        epsum = epsum + pairtest(k)%epsum
     enddo
     do i=1,mol(1)%nmol
        do j=1,mol(1)%nmol
           if( pairtest(1)%ep(i,j) /= 0d0 &
                .or. pairinfo(1)%ep(i,j) /= 0d0 ) then
              if ( dabs(&
                   (pairinfo(1)%ep(i,j)-pairtest(1)%ep(i,j))/&
                   (pairinfo(1)%ep(i,j)+pairtest(1)%ep(i,j)))&
                   > 0.0000000001d0 ) then
                 write(STDERR,*) i,j,pairinfo(1)%ep(i,j),pairtest(1)&
                      & %ep(i,j)
              endif
           endif
        enddo
     enddo
     if ( dabs((epsum - oldsum)/epsum) > 0.00000001) write(6,*)&
          & "epsum",epsum,oldsum
#endif
  enddo
  deallocate( random )
  deallocate( ep )
  !deallocate( vr )
  deallocate( delta )
  !最終配置の出力、乱数の内部状態の出力
end program mctest

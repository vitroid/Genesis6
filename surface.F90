! -*- f90 -*-
#undef YAPLOT
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
program main
!subroutine traj2(stream, outputtype)
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
  use error_module
  use rigid_setup_module
  use flex_setup_module
  use interaction_module
 
  implicit none
  !integer, intent(IN) :: stream
  !integer, intent(IN) :: outputtype
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
  !
  !柔軟分子の質点の質量
  !
  real(kind=8) :: molmass(MAXSITE)
  type(vector3) :: com( MAXMOL, MAXCOMPO )
  character(len=5) :: tag
  character(len=8) :: lastid
  integer          :: ncompo, compo, incompo, ncombi
  integer          :: mdvwsite
  integer          :: i,j,k,ii,jj
  integer          :: nmol
  integer          :: combi( MAXCOMBI ), combj( MAXCOMBI )
  real(kind=8)     :: cutoff
  integer          :: nloop
  real(kind=8), allocatable :: ep1(:)
  real(kind=8)     :: vrmat(3,3)
  integer          :: pair
  integer          :: stream
  integer          :: outputtype
  integer          :: slice(256+200)
  integer          :: ibin, binmin,binmax, min,max,sum,mid
  real(kind=8)     :: bin
  integer          :: minmol, surfmol(MAXCOMPO),surfsite(MAXCOMPO),minsite
  real(kind=8)     :: surfx(MAXCOMPO), surfy(MAXCOMPO), surfz(MAXCOMPO)
  real(kind=8)     :: dx,dy,dz,minsx,minsy,minsz,gx,gy,gz,cx,cy,cz,mindist
  real(kind=8)     :: ex,ey,ez, lastd, d
  integer          :: ix,iy,m,s,check
  logical          :: converged
  integer          :: nearmol( 500, MAXCOMPO ), nnearmol( MAXCOMPO )
  real(kind=8)     :: nearmolx( 500, MAXCOMPO )
  real(kind=8)     :: nearmoly( 500, MAXCOMPO )
  real(kind=8)     :: nearmolz( 500, MAXCOMPO )
  real(kind=8)     :: dmin, dmindz,deltaz, dmi(MAXCOMPO)
  integer          :: is,js,dmincompo,dminmol,lastcompo
  integer          :: err

  stream = STDOUT
  outputtype = 6
  incompo = 0
  nloop = 0

  call new(box)
  do
     ncompo  = 0
     call new(site, MAXsite_total)
     do
        read(STDIN,'(a5)',END=999) tag
        !write(STDERR,*) tag
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
        if(tag.eq.'@NX4A' .or. tag.eq.'@WTG3' .or. tag.eq.'@AR3A')then
           if(lastid.eq."")then
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
           err = rigid_setup( rigid(ncompo), lastid, FIXED, mol(ncompo), si(ncompo) )
           if ( err .ne. error_id_not_found ) then
              err = Flex_Setup( lastid, .not. FIXED, mol(ncompo), si(ncompo), molmass )
              if ( err .ne. error_none ) then
                 call die( err, "surface 1" )
              endif
           endif
           if ( .not. mol(ncompo)%isRigid ) then
              call Flex_ReadAR3A(flex(ncompo),mol(ncompo),site,.true.,STDIN)
              call Flex_SetMass(flex(ncompo),mol(ncompo),molmass)
           endif
           if ( mol(ncompo)%isRigid ) then
              call Rigid_ReadNX4A(rigid(ncompo),mol(ncompo),site ,STDIN)
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

     cutoff = 5d0
     if ( active( box ) ) then
        if ( cutoff .eq. 0d0 ) then
           write(STDERR,*) "Cutoff length must be specified."
           stop
        endif
        if ( box%size%vec(1)*0.5d0 < cutoff ) cutoff = box%size%vec(1)*0.5d0
        if ( box%size%vec(1)*0.5d0 < cutoff ) cutoff = box%size%vec(1)*0.5d0
        if ( box%size%vec(1)*0.5d0 < cutoff ) cutoff = box%size%vec(1)*0.5d0
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
           do i=1,nmol
              dx = com( i, compo)%vec(1)
              dy = com( i, compo)%vec(2)
              dz = com( i, compo)%vec(3)+24d0
              call box_renormalize0(box, dx,dy,dz,cx,cy,cz)
#ifdef YAPLOT
              write(11,'(a1,i4)')          '@', compo+2
              write(11,'(a1,1(f12.5,1x))') 'r', 1.0
              write(11,'(a1,3(f12.5,1x))') 'o', dx,dy,dz
#endif
           enddo
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
     call Grid_SetAdaptive(grid, cutoff)
     call Grid_Update( grid, box, 0d0 )
     !
     !分子をセルに帰属させる。
     !
     do i=1,ncompo
        call Address_Assign(addr(i),box,grid,mol(i)%nmol,com(1,i))
     enddo
     !
     !ここから表面探索専用
     !まず、てきとうなxy断面を決め、格子分割する。
     !そこから、一番近い水と二酸化炭素を捜す。
     !二酸化炭素の方が遠い場合は、xy断面をその二酸化炭素の場所に移し、再度
     !一番近い水を捜して、その中点を探索開始点とする。
     !この最初の過程だけは、格子分割は使わないで大域的探索を行う。
     !
     !   成分1に注目する。z方向の密度分布を平滑化し、中央値をさがす。
     slice(:) = 0d0
     do i=1,mol(1)%nmol
        bin = com(i,1)%vec(3) / box%size%vec(3)
        bin = bin - dnint( bin )
        ibin = bin * 256 + 128 + 100
        slice(ibin) = slice(ibin) + 1
     enddo
     do i=1,100
        slice(256+100+i) = slice(i)
        slice(101-i) = slice(256+101-i)
     enddo
     max = 0
     min = 999999
     do ibin=1,256
        sum = 0
        do i=-20,20
           sum = sum + slice(ibin+i+100)
        enddo
        !write(6,*) ibin, slice(ibin+100), dble(sum) / 41.0
        if ( max < sum ) then
           max = sum
           binmax = ibin
        endif
        if ( sum < min ) then
           min = sum
           binmin = ibin
        endif
     enddo
     !write(6,*) "minmax:", binmin,binmax
     mid = ( max + min ) / 2
     ibin = binmin
     do 
        sum = 0
        do i=-20,20
           sum = sum + slice(ibin+i+100)
        enddo
        if ( mid < sum ) then
           goto 99
        endif
        ibin = ibin + 1
        if ( ibin .eq. 256 ) ibin = 1
     enddo
99   continue
     !write(6,*) "binmid:", ibin
     !これで探索基準面が決まった。
     write(6,'("@GRD2")')
     write(6,*) 32,32
     do ix=1,32
        gx = dble(ix) * box%size%vec(1) / 32d0
        do iy=1,32
           gy = dble(iy) * box%size%vec(2) / 32d0
           gz = dble(ibin) * box%size%vec(3) / 256d0
           !find nearest voronoi division
           !find nearest molecule
           dmin = 99999
           do compo=1,ncompo
              do i=1,mol(compo)%nmol
                 !dx = rigid(compo)%mol(i)%com%vec(1) - gx
                 !dy = rigid(compo)%mol(i)%com%vec(2) - gy
                 !dz = rigid(compo)%mol(i)%com%vec(3) - gz
                 dx = com( i, compo)%vec(1) - gx
                 dy = com( i, compo)%vec(2) - gy
                 dz = com( i, compo)%vec(3) - gz
                 call box_renormalize0(box, dx,dy,dz,cx,cy,cz)
                 !do is=1, mol(compo)%nsite
                    ex = dx !+ rigid(compo)%mol(i)%intrax(is)
                    ey = dy !+ rigid(compo)%mol(i)%intray(is)
                    ez = dz !+ rigid(compo)%mol(i)%intraz(is)
                    if ( ex**2 + ey**2 + ez**2 < dmin ) then
                       dmin = ex**2 + ey**2 + ez**2
                       dmincompo = compo
                       dminmol = i
                    endif
                 !enddo
              enddo
           enddo
           if ( dmincompo .eq. 1 ) then
              deltaz = +0.1d0
           else
              deltaz = -0.1d0
           endif
           lastcompo = dmincompo
           !write(6,*) dmincompo, dminmol, dmin, gz
98         continue
           gz = gz + deltaz
           dmin = 99999
           dmi(:) = 9999
           do compo=1,ncompo
              do i=1,mol(compo)%nmol
                 !dx = rigid(compo)%mol(i)%com%vec(1) - gx
                 !dy = rigid(compo)%mol(i)%com%vec(2) - gy
                 !dz = rigid(compo)%mol(i)%com%vec(3) - gz
                 dx = com( i, compo)%vec(1) - gx
                 dy = com( i, compo)%vec(2) - gy
                 dz = com( i, compo)%vec(3) - gz
                 call box_renormalize0(box, dx,dy,dz,cx,cy,cz)
                 !do is=1, mol(compo)%nsite
                    ex = dx !+ rigid(compo)%mol(i)%intrax(is)
                    ey = dy !+ rigid(compo)%mol(i)%intray(is)
                    ez = dz !+ rigid(compo)%mol(i)%intraz(is)
                    if ( ex**2 + ey**2 + ez**2 < dmin ) then
                       dmin = ex**2 + ey**2 + ez**2
                       dmincompo = compo
                       dminmol = i
                    endif
                    if ( ex**2 + ey**2 + ez**2 < dmi(compo) ) then
                       dmi(compo) = ex**2 + ey**2 + ez**2
                    endif
                 !enddo
              enddo
           enddo
           !write(6,*) "@",gz, dmincompo, dminmol, dmi(1),dmi(2), mol(1)%nmol,mol(2)%nmol
           if ( dmincompo .eq. lastcompo ) goto 98
           write(6,'(2(i4,1x),3(f8.4))') ix,iy,gx,gy,gz
           dx = gx
           dy = gy
           dz = gz+24d0
           call box_renormalize0(box, dx,dy,dz,cx,cy,cz)
#ifdef YAPLOT
           write(11,'(a1,i4)')          '@', 2
           write(11,'(a1,1(f12.5,1x))') 'r', 0.3
           write(11,'(a1,3(f12.5,1x))') 'o', dx,dy,dz
#endif
        enddo
     enddo
#ifdef YAPLOT
     write(11,*)
#endif
  enddo
9999 continue
end program main
!end subroutine traj2

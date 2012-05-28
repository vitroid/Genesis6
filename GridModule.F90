! -*- f90 -*-
!#define DEBUG
module grid_module
  use common_module
  use vector_module
  use error_module
  implicit none
#define mAddr(g,x,y,z) (((z)*g%ndivy+(y))*g%ndivx+(x)+1)
!共通部分
  
  !分子群ごとの属性
  type sAddress
     sequence
     !maximum number of residents in a cell
     integer :: nmax
     integer,dimension(MAXmol) :: cell
     !those two will not be used in Water2
     integer,dimension(MAXCELL) :: nresident
     integer,dimension(MAXINGRID,MAXCELL) :: resident
  end type sAddress
  
  type sNeighborCell
     sequence
     !vectorを2倍すると原点に戻る場合は.true.
     !logical,dimension(MAXNEIBORCELL) :: doublezero
     integer :: n,n0
     !以下の配列は、aclを生成する時にのみ必要。
     integer,dimension(MAXNEIBORCELL) :: x,y,z
     !vectorを2倍すると原点に戻るようなベクトルは別個に表にしておく。
     integer,dimension(MAXNEIBORCELL) :: x0,y0,z0

     !adjacent cell list for each cell
     integer,dimension(MAXCELL*MAXNEIBORCELL) :: celli,cellj,cellp&
          & ,cell0i,cell0j,cell0p
  end type sNeighborCell

  type sGrid
     sequence
     !実際のセルの大きさ
     type(vector3) :: cell
     !要求されたセルの大きさ(Adaptiveの場合)
     type(vector3) :: requested
     ! 分割数固定モードか、距離を固定して分割数を動的に決めるか
     integer :: mode
     integer :: ndivx,ndivy,ndivz
     ! 隣接セルの対表。同種だと13通り、異種だと26通り。
     type(sNeighborCell) :: homo,hetero
  end type sGrid
  integer,parameter :: GRID_NONE=0,GRID_FIX=1,GRID_ADAPTIVE=2
 ! 別名の定義
  interface new
     module procedure Grid_Constructor
  end interface

  interface load
     module procedure grid_loader
  end interface

  interface save
     module procedure grid_save
  end interface

  interface load_binary
     module procedure grid_binaryloader
  end interface

  interface save_binary
     module procedure grid_savebinary
  end interface

  interface active
     module procedure grid_isactive
  end interface

  interface useful
     module procedure grid_isuseful
  end interface

contains
  !とりあえず意味のある値を入れておく。あとで変更してもよい。
  subroutine Grid_Constructor(g)
    type(sGrid),intent(INOUT) :: g
    g%mode=GRID_NONE
    g%ndivx=1
    g%ndivy=1
    g%ndivz=1
    g%cell%vec(:) = 0d0
    g%requested%vec(:) = 0d0
    return
  end subroutine Grid_Constructor
  
  !
  !boxの変化に応じて、分割の仕方を動的に変更する。
  !
  subroutine Grid_Update( g, box, range )
    use box_module
    type(sGrid),  intent(inout) :: g
    type(sBox),   intent(in)    :: box
    real(kind=8), intent(IN)    :: range
    !
    !range is given just for backward compatibility.
    !
    !
    !local
    !
    logical      :: fUpdated
    real(kind=8) :: r

    r = range
    if(g%mode.eq.GRID_ADAPTIVE)then
       fUpdated=Grid_ResetDivision(g,box)
       r = 0d0
    else
       fUpdated=Grid_ResetCellSize(g,box)
    endif
#ifdef DEBUG
    write(STDERR,*) "GRID:", g%ndivx, g%ndivy, g%ndivz
#endif
    if(fUpdated)then
       call NeighborCell_ListHomo(g%homo,g,r)
       call NeighborCell_ListHetero(g%hetero,g,r)
       call NeighborCell_MakeAcl(g%homo,g)
       call NeighborCell_MakeAcl(g%hetero,g)
    endif
  end subroutine Grid_Update

  function Grid_IsActive( g )
    type(sGrid),intent(IN) :: g
    logical :: Grid_IsActive
    Grid_IsActive = g%mode.ne.GRID_NONE
  end function Grid_IsActive

  function Grid_IsUseful( g )
    type(sGrid),intent(IN) :: g
    logical :: Grid_IsUseful
    Grid_IsUseful = ( 3 < g%ndivx .or. 3 < g%ndivy .or. 3 < g%ndivz )
  end function Grid_IsUseful

  function Grid_ReadVOXN(g,file)
    integer,intent(IN) ::file
    type(sGrid),intent(OUT) :: g
    logical :: Grid_ReadVOXN
    read(file,*) g%ndivx,g%ndivy,g%ndivz
    Grid_ReadVOXN=(g%ndivx > 0)
    if(g%ndivx > 0)then
       g%mode = GRID_FIX
    else
       g%mode = GRID_NONE
    endif
    if(g%ndivx > MAXGRIDX .or. g%ndivy > MAXGRIDY .or. g&
         & %ndivz > MAXGRIDZ)then
       write(STDERR,*) "ERROR: TOO FINE GRID ",g%ndivx,g%ndivy,g%ndivz
       call die( 0, "Grid 1" )
    endif
#ifdef VERBOSE
    write(STDERR,*) "@VOXN",g%ndivx,g%ndivy,g%ndivz
#endif
    return
  end function Grid_ReadVOXN
  
  function Grid_ReadBinaryVOXN(g,file)
    integer,intent(IN) ::file
    type(sGrid),intent(OUT) :: g
    logical :: Grid_ReadBinaryVOXN
    read(file) g%ndivx,g%ndivy,g%ndivz
    Grid_ReadBinaryVOXN=(g%ndivx > 0)
    if(g%ndivx > 0)then
       g%mode = GRID_FIX
    else
       g%mode = GRID_NONE
    endif
    if(g%ndivx > MAXGRIDX .or. g%ndivy > MAXGRIDY .or. g&
         & %ndivz > MAXGRIDZ)then
       write(STDERR,*) "ERROR: TOO FINE GRID ",g%ndivx,g%ndivy,g%ndivz
       call die( 0, "Grid 2" )
    endif
#ifdef VERBOSE
    write(STDERR,*) "@VOXN",g%ndivx,g%ndivy,g%ndivz
#endif
    return
  end function Grid_ReadBinaryVOXN
  
  subroutine Grid_SetAdaptive(g, size)
    type(sGrid),intent(OUT) :: g
    real(kind=8), intent(IN) :: size
    g%mode          = GRID_ADAPTIVE
    g%requested%vec(1:3) = size
  end subroutine Grid_SetAdaptive
    
!VOXA形式の場合は、隣接セルが26個になるようにセル分割数のほうを調節す
!る。当面、分割数はxyzすべて同じに限るので、立方体以外のシステムでは効
!率が悪化する。(隣接セル表キャッシュの都合上)
  function Grid_ReadBinaryVOXA(g,file)
    integer,intent(IN) ::file
    type(sGrid),intent(OUT) :: g
    logical :: Grid_ReadBinaryVOXA
    read(file) g%requested%vec(1:3)
    Grid_ReadBinaryVOXA=(g%requested%vec(1) > 0)
    if(g%requested%vec(1) > 0)then
       g%mode = GRID_ADAPTIVE
    else
       g%mode = GRID_NONE
    endif
#ifdef VERBOSE
    write(STDERR,*) "@VOXA",g%requested%vec(1:3)
#endif
    return
  end function Grid_ReadBinaryVOXA
  
  function Grid_ReadVOXA(g,file)
    integer,intent(IN) ::file
    type(sGrid),intent(OUT) :: g
    logical :: Grid_ReadVOXA
    read(file,*) g%requested%vec(1:3)
    Grid_ReadVOXA=(g%requested%vec(1) > 0)
    if(g%requested%vec(1) > 0)then
       g%mode = GRID_ADAPTIVE
    else
       g%mode = GRID_NONE
    endif
#ifdef VERBOSE
    write(STDERR,*) "@VOXA",g%requested%vec(1:3)
#endif
    return
  end function Grid_ReadVOXA
  
  subroutine Grid_WriteBinaryVOXN(g,file)
    integer,intent(IN) ::file
    type(sGrid),intent(IN) :: g
    write(file) "@VOXN"
    write(file) g%ndivx,g%ndivy,g%ndivz
    return
  end subroutine Grid_WriteBinaryVOXN
  
  subroutine Grid_WriteBinaryVOXA(g,file)
    integer,intent(IN) ::file
    type(sGrid),intent(IN) :: g
    write(file) "@VOXA"
    write(file) g%requested%vec(1:3)
    return
  end subroutine Grid_WriteBinaryVOXA
  
  subroutine Grid_WriteVOXN(g,file)
    integer,intent(IN) ::file
    type(sGrid),intent(IN) :: g
    write(file,1) 
    write(file,*) g%ndivx,g%ndivy,g%ndivz
1   format("@VOXN")
    return
  end subroutine Grid_WriteVOXN
  
  subroutine Grid_WriteVOXA(g,file)
    integer,intent(IN) ::file
    type(sGrid),intent(IN) :: g
    write(file,1)
    write(file,*) g%requested%vec(1:3)
1   format( "@VOXA")
    return
  end subroutine Grid_WriteVOXA
  
  subroutine Grid_SetCellSize(g,b,x,y,z)
    use box_module
    real(kind=8),intent(in) :: x,y,z
    type(sGrid),intent(OUT) :: g
    type(sBox),intent(in) :: b
    logical :: result
    g%requested%vec(1)=x
    g%requested%vec(2)=y
    g%requested%vec(3)=z
    result=Grid_ResetDivision(g,b)
#ifdef VERBOSE
    write(STDERR,*) "@VOXA",g%cell%vec(1:3)
#endif
  end subroutine Grid_SetCellSize
  
!あるセルに隣接するセルの表を生成する。
!Homo(geneous)な系では、セルAに対してセルBから及ぼす力を計算すれば、そ
!の逆は計算しなくてよいので、セル対応表も半分だけ準備すればよい。
  subroutine NeighborCell_ListHomo(nc,g,maxradius)
    type(sGrid),intent(IN) :: g
    type(sNeighborCell),intent(OUT) :: nc
    real(kind=8),intent(in) :: maxradius
    integer :: jx,jy,jz
    integer :: pjx,pjy,pjz
    integer :: mjx,mjy,mjz
    real(kind=8)  :: maxradius2
    integer :: nmaxx,nmaxy,nmaxz
    integer :: nminx,nminy,nminz
    real(kind=8)  :: ddx,ddy,ddz
    real(kind=8)  :: dix,diy,diz
  !integer :: kk
  !格子の対称性を使って計算量を減らそうと思ったが、意外に難しい。とい
  !うのは、周期境界をまたぐことで、セル対a:bとb:aが同一になるケース(a
  !→bベクトルの2倍が0ベクトルと相同な場合。)があり、その場合には分子
  !対i,jを作る際に重複を避けるような工夫が必要になってくる。カットオフ
  !がセルサイズよりも十分短い場合にはあまりこういうケースは考えなくて
  !よいので、ひとまずややこしいながらも実装しておく。2種類作って速度を
  !比較するか。

  !この配列はその場で確保した方がいい。
    integer,dimension(0:MAXGRIDX-1,0:MAXGRIDY-1,0:MAXGRIDZ-1) :: mark
    !kk=0
    mark(0:g%ndivx-1,0:g%ndivy-1,0:g%ndivz-1)=0
    if(maxradius == 0d0)then
       !最近接セル26個を無条件に追加。
       nmaxx = 1
       nmaxy = 1
       nmaxz = 1
       nminx = -1
       nminy = -1
       nminz = -1
       if ( nmaxx > (g%ndivx-1)/2 ) nmaxx=(g%ndivx-1)/2
       if ( nmaxy > (g%ndivy-1)/2 ) nmaxy=(g%ndivy-1)/2
       if ( nmaxz > (g%ndivz-1)/2 ) nmaxz=(g%ndivz-1)/2
       if ( nminx < -g%ndivx/2 )    nminx=-g%ndivx/2
       if ( nminy < -g%ndivy/2 )    nminy=-g%ndivy/2
       if ( nminz < -g%ndivz/2 )    nminz=-g%ndivz/2
       nc%n=0
       nc%n0=0
       do jz=nminz,nmaxz
          diz=iabs(jz)-1
          if(diz < 0)diz=0
          pjz = jz
          if(pjz < 0)pjz=pjz+g%ndivz
          mjz = g%ndivz - pjz
          if(pjz == 0)mjz=0
          do jy=nminy,nmaxy
             diy=iabs(jy)-1
             if(diy < 0)diy=0
             pjy = jy
             if(pjy < 0)pjy=pjy+g%ndivy
             mjy = g%ndivy - pjy
             if(pjy == 0)mjy=0           
             do jx=nminx,nmaxx
                dix=iabs(jx)-1
                if(dix < 0)dix=0
                pjx = jx
                if(pjx < 0)pjx=pjx+g%ndivx
                mjx = g%ndivx - pjx
                if(pjx == 0)mjx=0
                if(mark(pjx,pjy,pjz) == 0)then
                   mark(mjx,mjy,mjz)=mark(mjx,mjy,mjz)+1
                   mark(pjx,pjy,pjz)=mark(pjx,pjy,pjz)+1
                   if(mark(pjx,pjy,pjz).ne.1)then
                      nc%n0 = nc%n0+1
#ifdef VERBOSE
                      if ( nc%n0 > MAXNEIBORCELL ) then
                         write(STDERR,*) "ERROR: TOO MANY NEIGHBOR&
                              & CELLS"&
                              & ,nc%n0
                         call die( 0, "Grid 3" )
                      endif
#endif
                      nc%x0(nc%n0) = pjx
                      nc%y0(nc%n0) = pjy
                      nc%z0(nc%n0) = pjz
                   else
                      nc%n = nc%n+1
#ifdef VERBOSE
                      if ( nc%n > MAXNEIBORCELL ) then
                         write(STDERR,*) "ERROR: TOO MANY NEIGHBOR&
                              & CELLS"&
                              & ,nc%n
                         call die( 0, "Grid 4" )
                      endif
#endif
                      nc%x(nc%n) = pjx
                      nc%y(nc%n) = pjy
                      nc%z(nc%n) = pjz
                   endif
                endif
             enddo
          enddo
       enddo
    else
       maxradius2=maxradius*maxradius
       nmaxx = int(maxradius/g%cell%vec(1)+1d0)
       nmaxy = int(maxradius/g%cell%vec(2)+1d0)
       nmaxz = int(maxradius/g%cell%vec(3)+1d0)
       nminx = -nmaxx
       nminy = -nmaxy
       nminz = -nmaxz
       !write(6,*) nminx,nmaxx
       !write(6,*) nminy,nmaxy
       !write(6,*) nminz,nmaxz
       if ( nmaxx > (g%ndivx-1)/2 ) nmaxx=(g%ndivx-1)/2
       if ( nmaxy > (g%ndivy-1)/2 ) nmaxy=(g%ndivy-1)/2
       if ( nmaxz > (g%ndivz-1)/2 ) nmaxz=(g%ndivz-1)/2
       if ( nminx < -g%ndivx/2 )    nminx=-g%ndivx/2
       if ( nminy < -g%ndivy/2 )    nminy=-g%ndivy/2
       if ( nminz < -g%ndivz/2 )    nminz=-g%ndivz/2
       !write(6,*) nminx,nmaxx
       !write(6,*) nminy,nmaxy
       !write(6,*) nminz,nmaxz
       nc%n=0
       nc%n0=0
       !mark(0,0,0)=.true.
       do jz=nminz,nmaxz
          diz=iabs(jz)-1
          if(diz < 0)diz=0
          ddz=(diz*g%cell%vec(3))**2
          pjz = jz
          if(pjz < 0)pjz=pjz+g%ndivz
          mjz = g%ndivz - pjz
          if(pjz == 0)mjz=0
          if(ddz < maxradius2)then
             do jy=nminy,nmaxy
                diy=iabs(jy)-1
                if(diy < 0)diy=0
                ddy=(diy*g%cell%vec(2))**2
                pjy = jy
                if(pjy < 0)pjy=pjy+g%ndivy
                mjy = g%ndivy - pjy
                if(pjy == 0)mjy=0           
                if(ddy+ddz < maxradius2)then
                   do jx=nminx,nmaxx
                      dix=iabs(jx)-1
                      if(dix < 0)dix=0
                      ddx=(dix*g%cell%vec(1))**2
                      pjx = jx
                      if(pjx < 0)pjx=pjx+g%ndivx
                      mjx = g%ndivx - pjx
                      if(pjx == 0)mjx=0
                      if(ddx+ddy+ddz < maxradius2)then
                         if(mark(pjx,pjy,pjz) == 0)then
                            mark(mjx,mjy,mjz)=mark(mjx,mjy,mjz)+1
                            mark(pjx,pjy,pjz)=mark(pjx,pjy,pjz)+1
                            !kk=kk+1
                            if(mark(pjx,pjy,pjz).ne.1)then
                               nc%n0 = nc%n0+1
#ifdef VERBOSE
                               if(nc%n0 > MAXNEIBORCELL)then
                                  write(STDERR,*) "ERROR: TOO MANY&
                                       & NEIGHBOR CELLS",nc%n0
                                  call die( 0, "Grid 5" )
                               endif
#endif
                               nc%x0(nc%n0) = pjx
                               nc%y0(nc%n0) = pjy
                               nc%z0(nc%n0) = pjz
                            else
                               nc%n = nc%n+1
#ifdef VERBOSE
                               if(nc%n > MAXNEIBORCELL)then
                                  write(STDERR,*) "ERROR: TOO MANY&
                                       & NEIGHBOR CELLS",nc%n
                                  call die( 0, "Grid 6" )
                               endif
#endif
                               nc%x(nc%n) = pjx
                               nc%y(nc%n) = pjy
                               nc%z(nc%n) = pjz
                            endif
                         endif
                      endif
                   enddo
                endif
             enddo
          endif
       enddo
    endif
#ifdef VERBOSE
    write(STDERR,*) "Number of adjacent cells:",nc%n,"+",nc%n0
#endif
  end subroutine NeighborCell_ListHomo

  subroutine Address_Initialize(a,n)
    !number of molecules
    integer,intent(in) :: n
    type(sAddress),intent(INOUT) :: a
  end subroutine Address_Initialize
  
  function Grid_ResetCellSize(g,b)
    use box_module
    type(sGrid),intent(INOUT) :: g
    type(sBox),intent(IN)     :: b
    type(vector3)             :: old
    logical                   :: Grid_ResetCellSize
    old=g%cell
    g%cell%vec(1) = b%size%vec(1)/g%ndivx
    g%cell%vec(2) = b%size%vec(2)/g%ndivy
    g%cell%vec(3) = b%size%vec(3)/g%ndivz
    Grid_ResetCellSize=((old%vec(1) /= g%cell%vec(1)) .or. (old%vec(2) /= g%cell%vec(2)) .or. (old%vec(3) /= g%cell%vec(3)))
    return
  end function Grid_ResetCellSize

  function Grid_ResetDivision(g,b)
    use box_module
    type(sGrid),intent(INOUT) :: g
    type(sBox),intent(IN) :: b
    integer :: oldx,oldy,oldz
    logical :: Grid_ResetDivision
    oldx=g%ndivx
    oldy=g%ndivy
    oldz=g%ndivz
  !分割数は切り捨てで計算。
    g%ndivx = b%size%vec(1)/g%requested%vec(1)
    g%ndivy = b%size%vec(2)/g%requested%vec(2)
    g%ndivz = b%size%vec(3)/g%requested%vec(3)
    Grid_ResetDivision = Grid_ResetCellSize(g,b)
  end function Grid_ResetDivision

  subroutine Address_Assign(a,b,g,n,com)
    use box_module
    use interaction_module
    type(sBox),intent(in) :: b
    type(sGrid),intent(IN) :: g
    type(sAddress),intent(INOUT) :: a
    integer,intent(in) :: n
    type(vector3),dimension(*),intent(IN) :: com
    integer :: jx,jy,jz,c
    integer :: i,nneib
    type(vector3)         :: tmp,ci,delta
  !これをちゃんと0にしとかないと、あとで問題が生じる。
    a%resident(:,:)=0
    a%nresident(:)=0
    a%nmax=0
    ci = g%cell
    call invert(ci)
#ifdef VERBOSE
    write(STDERR,*) g%cell%vec(1),g%cell%vec(2),g%cell%vec(3)
#endif
    do i=1,n
       delta = com(i)
       call renormalize(b,delta,tmp)
       jx=dnint(delta%vec(1)*ci%vec(1))
       if(jx < 0)jx = jx+g%ndivx
       jy=dnint(delta%vec(2)*ci%vec(2))
       if(jy < 0)jy = jy+g%ndivy
       jz=dnint(delta%vec(3)*ci%vec(3))
       if(jz < 0)jz = jz+g%ndivz
       a%cell(i)=mAddr(g,jx,jy,jz)
#ifdef VERBOSE
       !write(STDERR,*) i,jx,jy,jz,mAddr(g,jx,jy,jz)
#endif
    enddo
    !平成１２年４月１０日(月)ループを分割した。前半はvectorize可能。
    do i=1,n
       c =a%cell(i)
       nneib = a%nresident(c)
       nneib=nneib+1
       a%nresident(c)=nneib
       if(nneib > a%nmax)a%nmax=nneib
       a%resident(nneib,c)=i
#ifdef QUITEVERBOSE
       write(STDERR,*) i,c,nneib
#endif
    enddo
#ifdef VERBOSE
    write(STDERR,*) 'Maximum number of molecules in a cell=',a%nmax
#endif
  end subroutine Address_Assign

#ifdef MPI
  subroutine Grid_NeighborList(NPROCS,MYRANK,g,nc,in,ai,aj)
#else
  subroutine Grid_NeighborList(g,nc,in,ai,aj)
#endif
    use interaction_module
    type(sGrid),intent(IN) :: g
    type(sAddress),intent(IN) :: ai,aj
    type(sNeighborCell),intent(IN) :: nc
    type(sInteraction),intent(OUT) :: in
    integer :: cellpair
    integer :: i,j,k,p3i,ii,jj,cell1,cell2,p1,n1
    integer :: ng
    !number of molecules in adjacent cells
    integer,dimension(:),allocatable :: nmac_i
    !storage point of molecules
    integer,dimension(:),allocatable :: lvpos
    integer,dimension(MAXINGRID*MAXNEIBORCELL,MAXCELL) ::&
         & neighbormollist
    integer,dimension(MAXINGRID*MAXNEIBORCELL,MAXCELL) ::&
         & neighbormollist0
#ifdef MPI
    integer,intent(in) :: NPROCS,MYRANK
    integer :: from,to
#else
    integer,parameter :: from=0
#endif
    ng=g%ndivx*g%ndivy*g%ndivz
#ifdef MPI
    !aclも必要な範囲しか作らない。
    from=(ng*MYRANK)/NPROCS
    to  =(ng*(MYRANK+1))/NPROCS
    ng  = to-from
#endif
#ifdef VERBOSE
    write(STDERR,*) "Number of cells",ng
    write(STDERR,*) "Number of neighbor cells",nc%n
#endif
  !基本にたちかえって設計しなおす。
  !まず、各セルごとの、隣接分子の個数nmac_iを数える。
  !storeposは、隣接分子リストを並列生成するのに必要。
    !safely deallocated inside this routine.
    allocate(nmac_i(ng))
    nmac_i(:) = 0
    do cellpair=1,nc%n
       do i=1,ng
        !jの算出がもしスピードを低下させるなら、逆ベクトルのaclの表を
        !作って別計算にすることも可能。
          j=nc%cellj((i-1)*nc%n+cellpair)
          nmac_i(i)=nmac_i(i)+aj%nresident(j)
        !write(STDERR,*) cellpair,i,j,nmac_i(i),aj%nresident(j)
       enddo
    enddo
    !neighbormollist0(:,:)=0
    neighbormollist0(1:aj%nmax*nc%n,1:ng)=0
    do jj=1,aj%nmax
!OCL VECTOR,NOVREC
       do k=1,nc%n*ng
          cell1=nc%celli(k)
          cell2=nc%cellj(k)
          p3i=jj+(nc%cellp(k)-1)*aj%nmax
          neighbormollist0(p3i,cell1)=aj%resident(jj,cell2)
       enddo
    enddo
    do cell1=1,ng
       i=0
       do j=1,aj%nmax*nc%n
          if(neighbormollist0(j,cell1).ne.0)then
             i=i+1
             neighbormollist(i,cell1)=neighbormollist0(j,cell1)
          endif
       enddo
    enddo
  !縮小した隣接分子リストをもとに、list vectorを生成する。
    !safely deallocated inside this routine.
    allocate(lvpos(ng))
  !list vectorは分子i側の順序で作成する。
    k=0
    do i=1,ng
       lvpos(i)=k
       k=k+nmac_i(i)*ai%nresident(i+from)
    enddo
#ifdef VERBOSE
    write(STDERR,*) "Vector length(1)",k,"(",MAXPAIR,")"
    if(k > MAXPAIR)then
       write(STDERR,*) "Too many pairs"
       call die( 0, "Grid 7" )
    endif
#endif
  !平成１２年５月２２日(月)一旦セルごとの隣接分子表を作成する。並列機
  !の場合はこのリストを分割する。
  !ここまでの処理は並列機では重複して実行されることになる。あまり時間
  !がかかるようだと、この方式のメリットは減る。
  !倍とっておけば十分だろう。
    in%npair0=k
  !平成１２年５月２２日(月)今回の改造でこの部分がシンプルになったのは
  !収穫。内側の2重ループはlist vectorを準備しておいて1重化できそう。
    do cell1=1,ng
       p1=lvpos(cell1)
       n1=nmac_i(cell1)
       do ii=1,ai%nresident(cell1+from)
          i=ai%resident(ii,cell1+from)
          j=(ii-1)*n1+p1
!平成１２年５月２２日(月)なんでこんなところが遅いんだ？ネストが深すぎ
!る？ここは確かに遅い。2重ループを1重にできないか？
!slow 1745,455
!OCL VECTOR,NOVREC
          do jj=1,nmac_i(cell1)
             in%pair_i0(j+jj)=i
             in%pair_j0(j+jj)=neighbormollist(jj,cell1)
          enddo
       enddo
    enddo
    !keep neighbormollist()
    deallocate(lvpos)
    !n0の系列(2倍すると0ベクトルになるような相対セルベクトルの場合)は別
  !処理するしかない。これもよい方法があるだろうか?
  !一般にはこのようなセル対は少ないはずだが、分子ラベルを比較して登録
  !するかどうか判断する必要があるので、より時間がかかりうる。

  !平成１２年４月２４日(月)やっとこの部分が最も時間のかかる領域になっ
  !た(他の部分はベクトル化の方針が立った)

  !平成１２年５月２２日(月)この部分を並列機用にするには、各CPUにできる
  !だけ均一に処理を分散する必要がある。
  !nmac_i,neighbormollistを再利用する。
    nmac_i(:) = 0
    do cellpair=1,nc%n0
       do i=1,ng
          !jの算出がもしスピードを低下させるなら、逆ベクトルのaclの表を
          !作って別計算にすることも可能。
          j=nc%cell0j((i-1)*nc%n0+cellpair)
          nmac_i(i)=nmac_i(i)+aj%nresident(j)
#ifdef QUITEVERBOSE
          write(STDERR,*) "N-zero", i,nmac_i(i)
#endif
       enddo
    enddo
    neighbormollist0(1:aj%nmax*nc%n0,1:ng)=0
  !neighbormollist0(:,:)=0
    do jj=1,aj%nmax
!OCL VECTOR,NOVREC
       do k=1,nc%n0*ng
          cell1=nc%cell0i(k)
          cell2=nc%cell0j(k)
          p3i=jj+(nc%cell0p(k)-1)*aj%nmax
          neighbormollist0(p3i,cell1)=aj%resident(jj,cell2)
       enddo
    enddo
    do cell1=1,ng
       i=0
       do j=1,aj%nmax*nc%n0
          if(neighbormollist0(j,cell1).ne.0)then
             i=i+1
             neighbormollist(i,cell1)=neighbormollist0(j,cell1)
          endif
       enddo
    enddo
    k=in%npair0
    do cell1=1,ng
       do ii=1,ai%nresident(cell1+from)
          i=ai%resident(ii,cell1+from)
          !slow 568
#ifdef QUITEVERBOSE
          write(STDERR,*) ng, cell1, nmac_i(cell1)
#endif
          do jj=1,nmac_i(cell1)
             j=neighbormollist(jj,cell1)
             if(j > i)then
                k=k+1
                in%pair_i0(k)=i
                in%pair_j0(k)=j
             endif
          enddo
       enddo
    enddo
#ifdef VERBOSE
    write(STDERR,*) "Vector length(2)",k
    if(k > MAXPAIR)then
       write(STDERR,*) "Too many pairs"
       call die( 0, "Grid 8" )
    endif
#endif
    
    deallocate(nmac_i)
    in%npair0=k
  !以上でlist vectorがまずできた。力保管場所表を生成する前に、対距離を
  !測定して、遠い要素を除去する。この操作はvectorizeできるはず。

  !加算があとじゃないとうまくvector化できないのかもしれない。
  !平成１２年４月２５日(火)収集拡散方式において、newkの初期値を1にする
  !とsegmentation faultするのは翻訳系のバグのような気がする。
  !平成１２年４月２５日(火)原因発見。newkの最終値はpairの個数+1になっ
  !ているので1引かなければいけない。
  !newk=0
  !収集拡散の計算は、ある単一の配列に対してしか行えないような気配があ
  !る。そこで、それにあわせてコードする。
    return
  end subroutine Grid_NeighborList

#ifdef MPI
  subroutine NeighborCell_ListHetero(NPROCS,MYRANK,nc,g,maxradius)
#else
  subroutine NeighborCell_ListHetero(nc,g,maxradius)
#endif
    type(sGrid),intent(IN) :: g
    type(sNeighborCell),intent(OUT) :: nc
    real(kind=8),intent(in) :: maxradius
    integer :: jx,jy,jz
    integer :: pjx,pjy,pjz
    integer :: mjx,mjy,mjz
    real(kind=8)  :: maxradius2
    integer :: nmaxx,nmaxy,nmaxz
    integer :: nminx,nminy,nminz
    real(kind=8)  :: ddx,ddy,ddz
    real(kind=8)  :: dix,diy,diz
#ifdef MPI
    integer :: from,to
    integer,intent(IN) :: NPROCS,MYRANK
#endif
    if(maxradius == 0d0)then
       !最近接セル26個を無条件に追加。
       nmaxx = 1
       nmaxy = 1
       nmaxz = 1
       nminx = -1
       nminy = -1
       nminz = -1
       if ( nmaxx > (g%ndivx-1)/2 ) nmaxx=(g%ndivx-1)/2
       if ( nmaxy > (g%ndivy-1)/2 ) nmaxy=(g%ndivy-1)/2
       if ( nmaxz > (g%ndivz-1)/2 ) nmaxz=(g%ndivz-1)/2
       if ( nminx < -g%ndivx/2 )    nminx=-g%ndivx/2
       if ( nminy < -g%ndivy/2 )    nminy=-g%ndivy/2
       if ( nminz < -g%ndivz/2 )    nminz=-g%ndivz/2
       nc%n=0
       !always zero in case of hetero-interaction
       nc%n0=0
       do jz=nminz,nmaxz
          diz=iabs(jz)-1
          if(diz < 0)diz=0
          pjz = jz
          if(pjz < 0)pjz=pjz+g%ndivz
          mjz = g%ndivz - pjz
          if(pjz == 0)mjz=0
          do jy=nminy,nmaxy
             diy=iabs(jy)-1
             if(diy < 0)diy=0
             pjy = jy
             if(pjy < 0)pjy=pjy+g%ndivy
             mjy = g%ndivy - pjy
             if(pjy == 0)mjy=0           
             do jx=nminx,nmaxx
                dix=iabs(jx)-1
                if(dix < 0)dix=0
                pjx = jx
                if(pjx < 0)pjx=pjx+g%ndivx
                mjx = g%ndivx - pjx
                if(pjx == 0)mjx=0
                nc%n = nc%n+1
#ifdef VERBOSE
                if ( nc%n > MAXNEIBORCELL ) then
                   write(STDERR,*) "ERROR: TOO MANY NEIGHBOR&
                        & CELLS"&
                        & ,nc%n
                   call die( 0, "Grid 9" )
                endif
#endif
                nc%x(nc%n) = pjx
                nc%y(nc%n) = pjy
                nc%z(nc%n) = pjz
             enddo
          enddo
       enddo
    else
       maxradius2=maxradius*maxradius
       nmaxx = int(maxradius/g%cell%vec(1)+1d0)
       nmaxy = int(maxradius/g%cell%vec(2)+1d0)
       nmaxz = int(maxradius/g%cell%vec(3)+1d0)
       nminx = -nmaxx
       nminy = -nmaxy
       nminz = -nmaxz
       !write(6,*) nminx,nmaxx
       !write(6,*) nminy,nmaxy
       !write(6,*) nminz,nmaxz
       if(nmaxx > (g%ndivx-1)/2)nmaxx=(g%ndivx-1)/2
       if(nmaxy > (g%ndivy-1)/2)nmaxy=(g%ndivy-1)/2
       if(nmaxz > (g%ndivz-1)/2)nmaxz=(g%ndivz-1)/2
       if(nminx < -g%ndivx/2)nminx=-g%ndivx/2
       if(nminy < -g%ndivy/2)nminy=-g%ndivy/2
       if(nminz < -g%ndivz/2)nminz=-g%ndivz/2
       !write(6,*) nminx,nmaxx
       !write(6,*) nminy,nmaxy
       !write(6,*) nminz,nmaxz
       nc%n=0
       !always zero
       nc%n0=0
       do jz=nminz,nmaxz
          diz=iabs(jz)-1
          if(diz < 0)diz=0
          ddz=(diz*g%cell%vec(3))**2
          pjz = jz
          if(pjz < 0)pjz=pjz+g%ndivz
          mjz = g%ndivz - pjz
          if(pjz == 0)mjz=0
          if(ddz < maxradius2)then
             do jy=nminy,nmaxy
                diy=iabs(jy)-1
                if(diy < 0)diy=0
                ddy=(diy*g%cell%vec(2))**2
                pjy = jy
                if(pjy < 0)pjy=pjy+g%ndivy
                mjy = g%ndivy - pjy
                if(pjy == 0)mjy=0           
                if(ddy+ddz < maxradius2)then
                   do jx=nminx,nmaxx
                      dix=iabs(jx)-1
                      if(dix < 0)dix=0
                      ddx=(dix*g%cell%vec(1))**2
                      pjx = jx
                      if(pjx < 0)pjx=pjx+g%ndivx
                      mjx = g%ndivx - pjx
                      if(pjx == 0)mjx=0
                      if(ddx+ddy+ddz < maxradius2)then
                         nc%n = nc%n+1
#ifdef VERBOSE
                         if(nc%n > MAXNEIBORCELL)then
                            write(STDERR,*) "ERROR: TOO MANY NEIGHBOR&
                                 & CELLS",nc%n
                            call die( 0, "Grid 10" )
                         endif
#endif
                         nc%x(nc%n) = pjx
                         nc%y(nc%n) = pjy
                         nc%z(nc%n) = pjz
                      endif
                   enddo
                endif
             enddo
          endif
       enddo
    endif
#ifdef VERBOSE
    write(STDERR,*) "Number of adjacent cells:",nc%n,"+",nc%n0
#endif
    return
  end subroutine NeighborCell_ListHetero

#ifdef MPI
  subroutine NeighborCell_MakeAcl(NPROCS,MYRANK,nc,g)
#else
  subroutine NeighborCell_MakeAcl(nc,g)
#endif
    type(sGrid),intent(IN) :: g
    type(sNeighborCell),intent(INOUT) :: nc
#ifdef MPI
    integer :: from,to
    integer,intent(IN) :: NPROCS,MYRANK
#endif
    integer :: ng
    real(kind=8) :: dx,dy,dz
    integer :: i,j,k
    integer :: jx,jy,jz,ix,iy,iz
    integer,dimension(:),allocatable :: lvx,lvy,lvz
    ng = g%ndivx*g%ndivy*g%ndivz
    !safely deallocated inside this routine.
    allocate(lvx(ng))
    allocate(lvy(ng))
    allocate(lvz(ng))
    i=0
    do iz=0,g%ndivz-1
       do iy=0,g%ndivy-1
          do ix=0,g%ndivx-1
             !i = mAddr(g,ix,iy,iz)
             i=i+1
             lvx(i)=ix
             lvy(i)=iy
             lvz(i)=iz
          enddo
       enddo
    enddo
#ifdef MPI
    from=(ng*MYRANK)/NPROCS
    to  =(ng*(MYRANK+1))/NPROCS
    ng=to-from
    write(STDERR,*) "MPI",MYRANK,"/",NPROCS,":",from,to,"(",ng,")"
    lvx(1:ng) = lvx(from+1:to)
    lvy(1:ng) = lvy(from+1:to)
    lvz(1:ng) = lvz(from+1:to)
#endif
    do j=1,nc%n
       jx = nc%x(j)
       jy = nc%y(j)
       jz = nc%z(j)
       do i=1,ng
          ix = lvx(i)
          iy = lvy(i)
          iz = lvz(i)
          dz = mod(iz+jz,g%ndivz)
          dy = mod(iy+jy,g%ndivy)
          dx = mod(ix+jx,g%ndivx)
          k = mAddr(g,dx,dy,dz)
          nc%celli((i-1)*nc%n+j) = i
          nc%cellj((i-1)*nc%n+j) = k
          nc%cellp((i-1)*nc%n+j) = j
       enddo
    enddo
    do j=1,nc%n0
       jx = nc%x0(j)
       jy = nc%y0(j)
       jz = nc%z0(j)
       do i=1,ng
          ix = lvx(i)
          iy = lvy(i)
          iz = lvz(i)
          dz = mod(iz+jz,g%ndivz)
          dy = mod(iy+jy,g%ndivy)
          dx = mod(ix+jx,g%ndivx)
          k = mAddr(g,dx,dy,dz)
          nc%cell0i((i-1)*nc%n0+j) = i
          nc%cell0j((i-1)*nc%n0+j) = k
          nc%cell0p((i-1)*nc%n0+j) = j
       enddo
    enddo
    deallocate(lvx)
    deallocate(lvy)
    deallocate(lvz)
  end subroutine NeighborCell_MakeAcl

  subroutine Grid_Loader(g,file,tag)
    integer,intent(IN) ::file
    type(sGrid),intent(OUT) :: g
    character(len=5),intent(IN) :: tag
    logical :: result
    if(tag == "@VOXN")then
       result = Grid_ReadVOXN(g,file)
    endif
    if(tag == "@VOXA")then
       result = Grid_ReadVOXA(g,file)
    endif
  end subroutine Grid_Loader
  
  subroutine Grid_BinaryLoader(g,file,tag)
    integer,intent(IN) ::file
    type(sGrid),intent(OUT) :: g
    character(len=5),intent(IN) :: tag
    logical :: result
    if(tag == "@VOXN")then
       result = Grid_ReadBinaryVOXN(g,file)
    endif
    if(tag == "@VOXA")then
       result = Grid_ReadBinaryVOXA(g,file)
    endif
  end subroutine Grid_BinaryLoader
  
  subroutine Grid_Save(g,file)
    integer,intent(IN) ::file
    type(sGrid),intent(IN) :: g
    if(g%mode  ==  GRID_FIX)then
       call Grid_WriteVOXN(g,file)
    else
       if(g%mode  == GRID_ADAPTIVE)then
          call Grid_WriteVOXA(g,file)
       endif
    endif
  end subroutine Grid_Save
  
  subroutine Grid_SaveBinary(g,file)
    integer,intent(IN) ::file
    type(sGrid),intent(IN) :: g
    if(g%mode  ==  GRID_FIX)then
       call Grid_WriteBinaryVOXN(g,file)
    else
       if(g%mode  ==  GRID_ADAPTIVE)then
          call Grid_WriteBinaryVOXA(g,file)
       endif
    endif
  end subroutine Grid_SaveBinary
end module grid_module

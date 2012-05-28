! -*- f90 -*-
module distrib_module
  implicit none
  type sHistogram
     integer :: nbin
     real(kind=8) :: binwidth, minvalue
     real(kind=8), pointer :: bins(:)
  end type sHistogram


  type sHistogramPtr
     type( sHistogram ), pointer :: p
  end type sHistogramPtr

  !
  !依存関係をへらすため、sRDFだけここで定義する。
  !
  type sRDF
     type(sHistogramPtr), pointer :: hist(:,:) ! サイト種ごとに別個に
     integer, pointer :: itype(:), jtype(:) !  サイトの種類。
     integer :: nitype, njtype
     character(len=8), pointer :: iname(:), jname(:) !種類ごとの名前
     real(kind=8), pointer :: imult(:), jmult(:) !種類ごとの重複度
  end type sRDF

  interface new
     module procedure histogram_initialize
  end interface

  interface done
     module procedure histogram_done
  end interface

contains


  subroutine histogram_initialize( hist, nbin, binwidth, minvalue )
    use error_module
    type(sHistogram), pointer       :: hist
    integer                         :: nbin
    real(kind=8)                    :: binwidth, minvalue
    if ( associated( hist ) ) then
       call die( 0, "sHistogram is already allocated." )
    endif
    allocate( hist )
    hist%nbin = nbin
    hist%minvalue = minvalue
    hist%binwidth = binwidth
    allocate( hist%bins( nbin ) )
    hist%bins(:) = 0d0
  end subroutine histogram_initialize

  subroutine histogram_dup( hist, from )
    type(sHistogram), pointer     :: hist
    type(sHistogram), intent(IN)  :: from
    call new( hist, from%nbin, from%binwidth, from%minvalue )
    hist%bins(:) = from%bins(:)
  end subroutine histogram_dup
    

  subroutine histogram_done( hist )
    type(sHistogram), pointer :: hist
    deallocate( hist%bins )
    deallocate( hist )
  end subroutine histogram_done
  !
  !重心間動径分布関数の積算
  !
  subroutine comrdf_accumulate(hist,nsite,com,box)
    use interaction_module
    use box_module
    use vector_module
    use grid_module
    type(sHistogram), intent(inout) :: hist
    type(sBox),       intent(in)    :: box
    integer,          intent(in)    :: nsite
    type(vector3),    intent(in)    :: com(nsite)
    !
    !local
    !
    type(sInteraction) :: ww
    real(kind=8)       :: maxvalue
    real(kind=8)       :: dx,dy,dz,r
    integer            :: i,j,k,bin
    real(kind=8)       :: size
    type(sGrid)        :: grid
    type(sAddress)     :: addr

    maxvalue = hist%minvalue + hist%binwidth * hist%nbin
    !pair listの生成
    call interaction_initialize( ww, nsite, -1 )
    !
    !グリッド分割する場合
    !
    size = maxvalue
    call Grid_SetAdaptive(grid, size) 
    call Grid_Update( grid, box, 0d0 )
    ! 試しに分割して、3x3x3分割以下なら分割する意味はない。
    if( useful( grid ) )then
       !
       !各分子が所属するセルの表を格納するaddrを準備
       !
       call Address_Initialize(addr,nsite)
       !
       !相互作用対リストを初期化する。
       !
       call Interaction_Reset( ww )
       !
       !分子をセルに割りあてる。
       !
       call Address_Assign(addr,box,grid,nsite,com)
       !
       !隣接するセルの表から、相互作用対リストを生成する
       !
       !
       !同種分子間相互作用対リストの作成
       !
       call Grid_NeighborList(grid,grid%homo,ww,addr,addr)
       !
       !カットオフ外の対を除いて、相互作用対リストを圧縮する。
       !
       call Interaction_Compress2(ww,box%mode.ne.BOX_NONE,&
               box ,maxvalue**2, com,com,.false.)
    else
       !
       !全ての分子はほかの全ての分子と相互作用しうる。
       !
       call interaction_alltoall(ww)
       !
       !カットオフ外の対を除いて、相互作用対リストを圧縮する。
       !相互作用対リストを圧縮して、あとの力計算処理を高速化する
       !(ベクトル機では有効)
       !
       call interaction_compress( ww, com, com, .false., maxvalue,&
            & box )
    endif
    !積算
    do k=1,ww%npair
       i = ww%pair_i(k)
       j = ww%pair_j(k)
       dx = com(i)%vec(1) - com(j)%vec(1) - ww%ox(k)
       dy = com(i)%vec(2) - com(j)%vec(2) - ww%oy(k)
       dz = com(i)%vec(3) - com(j)%vec(3) - ww%oz(k)
       r  = sqrt( dx**2 + dy**2 + dz**2 )
       call histogram_accumulate( hist, r, 1d0 )
    enddo
  end subroutine comrdf_accumulate

  subroutine histogram_accumulate( hist, value, multi )
    type(sHistogram), intent(inout) :: hist
    real(kind=8)    , intent(in)    :: value
    real(kind=8)    , intent(in)    :: multi
    integer  :: bin
    bin = nint( ( value - hist%minvalue ) / hist%binwidth ) + 1
    if ( 1 <= bin .and. bin <= hist%nbin ) then
       hist%bins(bin) = hist%bins(bin) + multi
    endif
  end subroutine histogram_accumulate

  subroutine histogram_save( hist, tag, file )
    use common_module
    type(sHistogram), intent(in) :: hist
    character(len=5), intent(IN) :: tag
    integer, intent(IN) :: file

    integer      :: i
    call writetag( file, tag )
    write(file,*) hist%nbin, hist%binwidth, hist%minvalue
    do i=1, hist%nbin
       write(file,*) hist%bins(i)
    enddo
  end subroutine histogram_save

  subroutine histogram_load( hist, file )
    use common_module
    type(sHistogram), pointer :: hist
    integer, intent(IN) :: file

    integer :: nbin
    real(kind=8) :: binwidth, minvalue

    integer      :: i
    read(file,*) nbin, binwidth, minvalue
    call new( hist, nbin, binwidth, minvalue )
    do i=1, nbin
       read(file,*) hist%bins(i)
    enddo
  end subroutine histogram_load

  subroutine histogram_show( hist, file )
    type(sHistogram), intent(in) :: hist
    integer, intent(IN) :: file

    integer      :: i
    do i=1, hist%nbin
       write(file,*) ( i - 1 ) * hist%binwidth + hist%minvalue, hist%bins(i)
    enddo
  end subroutine histogram_show

  subroutine comrdf_normalize( hist, density, naccum )
    use physconst_module
    type(sHistogram), intent(inout) :: hist
    real(kind=8),       intent(in) :: density ! 1AA**3あたりの分子数
    integer,            intent(in) :: naccum
    !
    !local
    !
    integer :: i
    real(kind=8) :: radius,x,volume,navg
    do i=1, hist%nbin
       radius = (i-1) * hist%binwidth + hist%minvalue
       !
       !半径r厚さ2xのうす皮の体積は
       ! 4/3 pi ( (r+x)**3 - (r-x)**3 ) = 4/3 pi 2 ( 3 r**2 x + x**3 )
       !
       x      = hist%binwidth * 0.5d0
       volume = 4d0/3d0 * pi * 2d0 * ( 3d0 * radius**2 * x + x**3 )
       !
       !薄皮に含まれる粒子数の期待値
       !
       navg   = density * volume
       hist%bins(i) = hist%bins(i) / ( naccum * navg )
       !test
       write(6,*) radius, hist%bins(i)
    enddo
  end subroutine comrdf_normalize

!
!tgtがからむ三体角を計算し、histogramに足しこむ。
!
subroutine partial_accumulate( hist, set, nei, nei2, amount, verbose )
  use neighbors_module
  !use distrib_module
  use error_module
  use set_module
  !
  !頻度分布
  !
  type(sHistogram), pointer :: hist
  !
  !変更する必要のある分子の集合
  !
  type(sSet), intent(IN) :: set
  !
  !最近接4分子、5Aまでの隣接分子のリスト
  !
  type(sNeighbors), intent(IN) :: nei(*), nei2(*)
  !
  !頻度分布への加算量
  !
  real(kind=8), intent(IN) :: amount
  logical, intent(IN) :: verbose
  integer :: i,j,k,ii
  real(kind=8) :: costh
  do ii=1, set%size
     i = set%set(ii)
     do j=1, nei(i)%nnear
        do k=j+1, nei(i)%nnear
           !
           ! O-O-O角の余弦を計算する。
           !
           costh = nei(i)%dx(j)*nei(i)%dx(k) + nei(i)%dy(j)*nei(i)%dy(k) + nei(i)%dz(j)*nei(i)%dz(k)
           costh = costh / sqrt ( nei(i)%dd(j) * nei(i)%dd(k) )
           call histogram_accumulate( hist, costh, amount )
           if ( verbose ) write(STDERR,*) "<", nei(i)%near(j), 0, nei(i)%near(k), costh
        enddo
     enddo
  enddo
end subroutine partial_accumulate

subroutine partial_accumulate_indexed( hist, set, nei, nei2, amount, indexsize, index, verbose )
  use neighbors_module
  !use distrib_module
  use error_module
  use set_module
  type(sHistogram), pointer :: hist
  type(sSet), intent(IN) :: set
  type(sNeighbors), intent(IN) :: nei(*), nei2(*)
  real(kind=8), intent(IN) :: amount
  integer, intent(IN) :: indexsize
  integer, intent(IN) :: index(*)
  logical, intent(IN) :: verbose

  integer :: i,j,k,ii
  real(kind=8) :: costh
  do i=1, indexsize
     ii = index(i)
     if ( exist( set, ii ) ) then
        do j=1, nei(i)%nnear
           do k=j+1, nei(i)%nnear
              !
              ! O-O-O角の余弦を計算する。
              !
              costh = nei(i)%dx(j)*nei(i)%dx(k) + nei(i)%dy(j)*nei(i)%dy(k) + nei(i)%dz(j)*nei(i)%dz(k)
              costh = costh / sqrt ( nei(i)%dd(j) * nei(i)%dd(k) )
              call histogram_accumulate( hist, costh, amount )
              if ( verbose ) write(STDERR,*) ">", nei(i)%near(j), 0, nei(i)%near(k), costh
           enddo
        enddo
     endif
  enddo
end subroutine partial_accumulate_indexed

function kullback_div( hist1, hist2 )
  real(kind=8) :: kullback_div
  type(sHistogram), pointer :: hist1, hist2
  real(kind=8) :: div, sum1, sum2, p, q
  integer :: i
  sum1 = 0d0
  sum2 = 0d0
  do i=1, hist1%nbin
     sum1 = sum1 + hist1%bins(i)
     sum2 = sum2 + hist2%bins(i)
  enddo
  div = 0d0
  do i=1, hist1%nbin
     p = hist1%bins(i) / sum1
     q = hist2%bins(i) / sum2
     if ( q == 0 ) q = p
     if ( p /= 0 ) then
        div = div + p * log ( p / q )
        !write(6,*) p, log ( p / q )
     endif
  enddo
  kullback_div = div
end function kullback_div

end module distrib_module

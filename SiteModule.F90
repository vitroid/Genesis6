! -*- f90 -*-
!
!相互作用点(サイト)を一括して扱うモジュール。
!
!剛体分子は、分子全体を一つの剛体として扱い、運動は剛体の回転と並進と
!して記述する。一方、分子間相互作用は、剛体内座標に固定された相互作用
!点(サイト)同士の相対座標で計算される。Genesis6では、系に含まれる全分
!子の全サイトを、単一の配列として扱うことで、最も計算時間がかかる力の
!計算の処理をできるだけ効率化できるように設計している。
!
!柔軟分子の場合は、相互作用点一点一点が質量を持つものと考える。このた
!め、剛体分子に含まれるサイトと柔軟分子のサイトでは、属性が異なる(柔軟
!分子の質点は、一つ一つがいわば分子としてふるまうので、質量や速度、加
!速度といった属性を保持しなければならない。剛体分子の場合は剛体全体の
!回転速度や加速度といった属性はあるがサイトにはそのような属性は不要)。
!柔軟分子特有の属性はFlexModule.F90で付与している。剛体分子の属性は
!RigidModule.F90で扱う。
!
!相互作用をカットオフする場合、一般的には分子の重心からの距離を基準に
!して、相互作用をスクリーニングする。この場合、重心に対しても力が働く
!ことになる。(φ(r)=s(r)φ0(r)をrで微分して力を求める)そこで、TIP4Pな
!どの剛体分子を扱う場合には、重心サイトを仮想的なサイトとして追加し、
!重心に働く力を他のサイト同様に算出し、剛体運動に還元させている。
!
module site_module
  use common_module
  implicit none

  type sSite
     sequence
     !全サイトの個数
     integer :: nsite
     !後半は剛体分子のsiteが入る。
     integer :: nflexsite
     real(kind=8),pointer :: x(:), y(:), z(:)
     real(kind=8),pointer :: fx(:), fy(:), fz(:)
  end type sSite
  interface new
     module procedure Site_Constructor
  end interface
contains
  subroutine Site_constructor(si, n)
    type(sSite)        :: si
    integer,intent(IN) :: n
    !avoid gfortran error
    !if ( associated( si%x ) ) then
    if ( allocated( si%x ) ) then
       deallocate( si%x )
       deallocate( si%y )
       deallocate( si%z )
       deallocate( si%fx )
       deallocate( si%fy )
       deallocate( si%fz )
    endif
    !write(STDERR,*) "SI ALLOC",n
    allocate(si%x(n))
    allocate(si%y(n))
    allocate(si%z(n))
    allocate(si%fx(n))
    allocate(si%fy(n))
    allocate(si%fz(n))
    call Site_initialize(si)
  end subroutine Site_constructor
  
  subroutine Site_initialize(si)
    type(sSite) :: si
    si%nsite=0
    si%nflexsite=0
  end subroutine Site_initialize
  
  subroutine Site_resetforce(si)
    type(sSite) :: si
    !write(STDERR,*) "SI",si%nsite
    si%fx(1:si%Nsite)=0d0
    si%fy(1:si%Nsite)=0d0
    si%fz(1:si%Nsite)=0d0
  end subroutine Site_resetforce
end module site_module

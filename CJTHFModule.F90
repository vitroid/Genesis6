! -*- f90 -*-

!
!THFの構造と相互作用。Chandrasekhar/Jorgensenのパラメータによる5サイトモデル
!(so called CJ Flat model)
!ID is CJTHF___
!JCP77,5073(1982), JACS 103, 3976(1981)
!
module cjthf_module
  use physconst_module
  implicit none
  real(kind=8),private,parameter :: QC1 = 0.25d0
  real(kind=8),private,parameter :: QC2 = 0d0
  real(kind=8),private,parameter :: QO  =-0.500d0
  !real(kind=8),private,parameter :: SIGC = (7290000.0 / 1825.0)**(1d0/6d0)
  !real(kind=8),private,parameter :: SIGO = (500000.0 /  600.0)**(1d0/6d0)
  !real(kind=8),private,parameter :: EPSC = 1825d0 / SIGC**6 / 4d0 * CA
  !real(kind=8),private,parameter :: EPSO = 600d0 / SIGO**6 / 4d0 * CA
  real(kind=8),private,parameter :: SIGC = 3.98331002404014142825d0
  real(kind=8),private,parameter :: SIGO = 3.06763104238899231290d0
  real(kind=8),private,parameter :: EPSC = 0.47789215157775979392d0
  real(kind=8),private,parameter :: EPSO = 0.75312001838671926989d0
  real(kind=8),private,parameter :: MASSC = 12.01d0
  real(kind=8),private,parameter :: MASSO = 15.9949d0
  real(kind=8),private,parameter :: MASSH = 1.00787d0
  real(kind=8),private,parameter :: COC   = 111.0
  real(kind=8),private,parameter :: OCC   = 109.4
  real(kind=8),private,parameter :: CO    = 1.411
  real(kind=8),private,parameter :: CC    = 1.529
  integer, parameter             :: UATHFSITE = 6
  ! United-Atom -- 最後のサイトは重心、名前は空白にしておく。
  character(len=8),parameter     :: uathfName(UATHFSITE) = (/"O ", "Ca", "Cb", "Cb", "Ca", "  "/)
  real(kind=8),private,parameter :: mass(UATHFSITE) = (/MASSO, MASSC+2*MASSH, MASSC+2*MASSH, MASSC+2*MASSH, MASSC+2*MASSH, 0d0/)

contains
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine cjthf_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, UATHFSITE, LJ_COULOMB)
    !charges [Q]
    si%param(1,1) = QO
    si%param(1,2) = QC1
    si%param(1,3) = QC2
    si%param(1,4) = QC2
    si%param(1,5) = QC1
    si%param(1,6) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = EPSO
    si%param(2,2) = EPSC
    si%param(2,3) = EPSC
    si%param(2,4) = EPSC
    si%param(2,5) = EPSC
    si%param(2,6) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = SIGO
    si%param(3,2) = SIGC
    si%param(3,3) = SIGC
    si%param(3,4) = SIGC
    si%param(3,5) = SIGC
    si%param(3,6) = 0d0
  end subroutine cjthf_setinteraction

  subroutine Rigid_CJTHF_Constructor(r)
    use physconst_module
    use rigid_module
    type(sRigid),intent(INOUT) :: r
    type(vector3) :: first, second, third
    real(kind=8) :: comx,comy,comz
    integer :: i

    r%molx(1)  = 0d0
    r%moly(1)  = 0d0
    r%molz(1)  = 0d0

    r%molx(2)  = CO * dsin( COC*0.5d0*deg2rad )
    r%moly(2)  = CO * dcos( COC*0.5d0*deg2rad )
    r%molz(2)  = 0d0

    r%molx(5) =-CO * dsin( COC*0.5d0*deg2rad )
    r%moly(5) = CO * dcos( COC*0.5d0*deg2rad )
    r%molz(5) = 0d0

    r%molx(3)  = r%molx(2) + CC * dcos( ( 270d0 - OCC - COC*0.5d0 ) * deg2rad )
    r%moly(3)  = r%moly(2) + CC * dsin( ( 270d0 - OCC - COC*0.5d0 ) * deg2rad )
    r%molz(3)  = 0d0

    r%molx(4)  = r%molx(5) - CC * dcos( ( 270d0 - OCC - COC*0.5d0 ) * deg2rad )
    r%moly(4)  = r%moly(5) + CC * dsin( ( 270d0 - OCC - COC*0.5d0 ) * deg2rad )
    r%molz(4)  = 0d0

    !write(6,*) r%molx(3)-r%molx(4)
    !result: 1.52906 平成16年2月4日(水)

    comx = 0d0
    comy = 0d0
    comz = 0d0
    r%mass = 0d0
    do i=1,UATHFSITE
       comx = comx + r%molx(i) * mass(i)
       comy = comy + r%moly(i) * mass(i)
       comz = comz + r%molz(i) * mass(i)
       r%mass = r%mass + mass(i)
    enddo
    r%massi= 1d0/r%mass
    comx = comz * r%massi
    comy = comy * r%massi
    comz = comz * r%massi
    do i=1,UATHFSITE
       r%molx(i) = r%molx(i) - comx
       r%moly(i) = r%moly(i) - comy
       r%molz(i) = r%molz(i) - comz
    enddo

    r%molx(6)  = 0d0
    r%moly(6)  = 0d0
    r%molz(6)  = 0d0
    !comx = 0d0
    !comy = 0d0
    !comz = 0d0
    !do i=1,UATHFSITE
    !   comx = comx + r%molx(i) * mass(i)
    !   comy = comy + r%moly(i) * mass(i)
    !   comz = comz + r%molz(i) * mass(i)
    !enddo
    !write(STDERR,*) "CHECK: ", comx,comy,comz

    r%Ixx = 0d0
    r%Iyy = 0d0
    r%Izz = 0d0
    do i=1,UATHFSITE
       r%Ixx = r%Ixx + mass(i)*(r%moly(i)**2 + r%molz(i)**2)
       r%Iyy = r%Iyy + mass(i)*(r%molz(i)**2 + r%molx(i)**2)
       r%Izz = r%Izz + mass(i)*(r%molx(i)**2 + r%moly(i)**2)
    enddo
    r%ixxi = 1d0/r%ixx
    r%iyyi = 1d0/r%iyy
    r%izzi = 1d0/r%izz

    !write(6,*) sqrt( ( r%molx(11) - r%molx(1) )**2 + ( r%moly(11) - r%moly(1) )**2 + ( r%molz(11) - r%molz(1) )**2 )
    !write(6,*) sqrt( ( r%molx(11) - r%molx(8) )**2 + ( r%moly(11) - r%moly(8) )**2 + ( r%molz(11) - r%molz(8) )**2 )
    !write(6,*) sqrt( ( r%molx(5) - r%molx(8) )**2 + ( r%moly(5) - r%moly(8) )**2 + ( r%molz(5) - r%molz(8) )**2 )
    !write(6,*) sqrt( ( r%molx(5) - r%molx(2) )**2 + ( r%moly(5) - r%moly(2) )**2 + ( r%molz(5) - r%molz(2) )**2 )
    !write(6,*) sqrt( ( r%molx(1) - r%molx(2) )**2 + ( r%moly(1) - r%moly(2) )**2 + ( r%molz(1) - r%molz(2) )**2 )
    !do i=1,14
    !   write(6,10) "t", r%molx(i), r%moly(i), r%molz(i), i
    !enddo
    !10format(a1,1x,3(f12.5),i5)
    !stop
  end subroutine Rigid_CJTHF_Constructor

end module cjthf_module

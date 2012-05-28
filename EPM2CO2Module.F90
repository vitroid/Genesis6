! -*- f90 -*-
module epm2co2_module
  use common_module
  use physconst_module
  use interaction_module
  implicit none
  ! J.Phys.Chem.B105, 12092 (2001)
  real(kind=8), parameter :: Q_EPM2CO2_C =+0.6512d0
  real(kind=8), parameter :: Q_EPM2CO2_O =-0.3256d0
  real(kind=8), parameter :: EPS_EPM2CO2_C =(28.129/J2K * 1d-3)
  real(kind=8), parameter :: SIG_EPM2CO2_C =2.757d0
  real(kind=8), parameter :: EPS_EPM2CO2_O =(80.507/J2K * 1d-3)
  real(kind=8), parameter :: SIG_EPM2CO2_O =3.033d0
  !重心を含めて4点
  integer, parameter :: EPM2CO2SITE = 4
  !
  character(len=8),parameter :: co2Name(EPM2CO2SITE) = (/"C", "O", "O", " "/)

contains
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine epm2co2_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, EPM2CO2SITE, LJ_COULOMB)
    !charges [Q]
    si%param(1,1) = Q_EPM2CO2_C
    si%param(1,2) = Q_EPM2CO2_O
    si%param(1,3) = Q_EPM2CO2_O
    si%param(1,4) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = EPS_EPM2CO2_C
    si%param(2,2) = EPS_EPM2CO2_O
    si%param(2,3) = EPS_EPM2CO2_O
    si%param(2,4) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = SIG_EPM2CO2_C
    si%param(3,2) = SIG_EPM2CO2_O
    si%param(3,3) = SIG_EPM2CO2_O
    si%param(3,4) = 0d0
  end subroutine epm2co2_setinteraction


  subroutine Rigid_EPM2CO2_Constructor(r)
    use rigid_module
    type(sRigid),intent(INOUT) :: r
    real(kind=8),dimension(4) :: mass
    integer :: i
    mass(1)=12d0
    mass(2)=16d0
    mass(3)=16d0
    mass(4)=0d0
    r%molx(1)=0d0
    r%molx(2)=0d0
    r%molx(3)=0d0
    r%molx(4)=0d0
    r%moly(1)=0d0
    r%moly(2)=0d0
    r%moly(3)=0d0
    r%moly(4)=0d0
    r%molz(1)=0d0
    r%molz(2)=-1.149d0
    r%molz(3)=+1.149d0
    r%molz(4)=0d0
    r%Ixx = 0d0
    r%Iyy = 0d0
    r%Izz = 0d0
    r%mass = 0d0
    do i=1,4
       r%Ixx = r%Ixx + mass(i)*(r%moly(i)**2 + r%molz(i)**2)
       r%Iyy = r%Iyy + mass(i)*(r%molz(i)**2 + r%molx(i)**2)
       r%Izz = r%Izz + mass(i)*(r%molx(i)**2 + r%moly(i)**2)
       r%mass = r%mass + mass(i)
    enddo
    r%ixxi = 1d0/r%ixx
    r%iyyi = 1d0/r%iyy
    r%izzi = 0d0
    r%massi= 1d0/r%mass
    return
  end subroutine Rigid_EPM2CO2_Constructor

end module epm2co2_module

! -*- f90 -*-
module dumbbell_module
  use common_module
  use physconst_module
  use interaction_module
  implicit none
  ! J.Phys.Chem.B105, 12092 (2001)
  real(kind=8), parameter :: EPS_DUMBBELL = 0d0
  real(kind=8), parameter :: SIG_DUMBBELL = 1d0
  !重心を含めて4点
  integer, parameter :: DUMBBELLSITE = 3
  character(len=8), parameter :: DUMBBELL_ID08="DUMBBELL"
  !
  character(len=8),parameter :: DumbbellName(DUMBBELLSITE) = (/"A", "B", " "/)

contains
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine dumbbell_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, DUMBBELLSITE, LJ_COULOMB )
    !charges [Q]
    si%param(1,1) = 0d0
    si%param(1,2) = 0d0
    si%param(1,3) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = 0d0
    si%param(2,2) = 0d0
    si%param(2,3) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = 1d0
    si%param(3,2) = 1d0
    si%param(3,3) = 0d0
  end subroutine dumbbell_setinteraction


  subroutine Rigid_DUMBBELL_Constructor(r)
    use rigid_module
    type(sRigid),intent(INOUT) :: r
    real(kind=8),dimension(4) :: mass
    integer :: i
    mass(1)=1d0
    mass(2)=1d0
    mass(3)=0d0
    r%molx(1)=0d0
    r%molx(2)=0d0
    r%molx(3)=0d0
    r%moly(1)=0d0
    r%moly(2)=0d0
    r%moly(3)=0d0
    r%molz(1)=-0.5d0
    r%molz(2)=+0.5d0
    r%molz(3)=0d0
    r%Ixx = 0d0
    r%Iyy = 0d0
    r%Izz = 0d0
    r%mass = 0d0
    do i=1,3
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
  end subroutine Rigid_DUMBBELL_Constructor

end module dumbbell_module

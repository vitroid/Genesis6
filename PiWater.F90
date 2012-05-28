! -*- f90 -*-
module piwater_module
  use common_module
  use physconst_module
  use interaction_module
  use tip4p_constant_module
  implicit none
  ! Use TIP4P parameters but the number of sites is 3 and the angle is pi.
  !重心を含めて4点
  integer, parameter :: PIWATERSITE = 4
  !use WaterName in spce_module
  !character(len=8),parameter :: co2Name(PIWATERSITE) = (/"C", "O", "O", " "/)
  character(len=8), parameter :: piwater_id08 = "PI_WATER"

contains
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine piwater_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, PIWATERSITE, LJ_COULOMB)
    !charges [Q]
    si%param(1,1) = WM
    si%param(1,2) = WH
    si%param(1,3) = WH
    si%param(1,4) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = EPSWO
    si%param(2,2) = 0d0
    si%param(2,3) = 0d0
    si%param(2,4) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = SIGWO
    si%param(3,2) = 0d0
    si%param(3,3) = 0d0
    si%param(3,4) = 0d0
  end subroutine piwater_setinteraction


  subroutine Rigid_PIWATER_Constructor(r)
    use rigid_module
    type(sRigid),intent(INOUT) :: r
    real(kind=8),dimension(4) :: mass
    integer :: i
    mass(1)=16d0
    mass(2)=1d0
    mass(3)=1d0
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
    r%molz(2)=-0.9572d0
    r%molz(3)=+0.9572d0
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
  end subroutine Rigid_PIWATER_Constructor

end module piwater_module

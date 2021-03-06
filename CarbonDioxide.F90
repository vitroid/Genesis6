! -*- f90 -*-
module co2_module
  use common_module
  use physconst_module
  use interaction_module
  implicit none
  ! J.Phys.Chem.B105, 12093 (2001)
  real(kind=8), parameter :: Q_CO2_C =-0.3256d0
  real(kind=8), parameter :: Q_CO2_O =+0.6512d0
  real(kind=8), parameter :: EPS_CO2_C =(28.129/J2K * 1d-3)
  real(kind=8), parameter :: SIG_CO2_C =2.757d0
  real(kind=8), parameter :: EPS_CO2_O =(80.507/J2K * 1d-3)
  real(kind=8), parameter :: SIG_CO2_O =3.033d0
  !$B=E?4$r4^$a$F(B4$BE@(B
  integer, parameter :: CO2SITE = 4

contains
  !
  !StdInteraction$B$r;H$&>l9g$N=i4|2=(B
  !
  subroutine co2_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, CO2SITE, LJ_COULOMB)
    !charges [Q]
    si%param(1,1) = Q_CO2_C
    si%param(1,2) = Q_CO2_O
    si%param(1,3) = Q_CO2_O
    si%param(1,4) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = EPS_CO2_C
    si%param(2,2) = EPS_CO2_O
    si%param(2,3) = EPS_CO2_O
    si%param(2,4) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = SIG_CO2_C
    si%param(3,2) = SIG_CO2_O
    si%param(3,3) = SIG_CO2_O
    si%param(3,4) = 0d0
  end subroutine co2_setinteraction

end module co2_tip4p_module

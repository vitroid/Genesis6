! -*- f90 -*-
module argon_module
  use physconst_module
  use interaction_module
  implicit none
  !Argon's LJ parameters (Theory of Simple Liquid P.53)
  !sigma = 3.405A, eps = 119.8K = 0.997730kJ/mol
  real(kind=8), parameter :: EPSAR =(119.8/J2K * 1d-3)
  real(kind=8), parameter :: SIGAR =3.405
  !United atom model of methane by Jorgensen JCP 89, 3742 (1988)
  !id should be LJME____.
  real(kind=8), parameter :: EPSME =(0.2940 * CA)
  real(kind=8), parameter :: SIGME =3.730
  !Ideal Gas
  real(kind=8), parameter :: EPSIG =0d0
  real(kind=8), parameter :: SIGIG =1d0
  !ReduceAr
  real(kind=8), parameter :: EPSRA =1d0
  real(kind=8), parameter :: SIGRA =1d0
  !重心を含めて2点
  integer, parameter :: ARGONSITE = 2
  !
  character(len=8),parameter :: argonName(ARGONSITE) = (/"Ar", "  "/)
  character(len=8),parameter :: uaMethaneName(ARGONSITE) = (/"Me", "  "/)
  character(len=8),parameter :: IdealGasName(ARGONSITE) = (/"X", " "/)
  character(len=8) :: LJAR_ID08     = 'LJAR    '
  character(len=8) :: LJME_ID08     = 'LJME____'
  character(len=8) :: IDEALGAS_ID08 = 'IDEALGAS'
  character(len=8) :: REDUCEAR_ID08 = 'REDUCEAR'

contains

  subroutine Argon_GetMass( mass )
    implicit none
    real(kind=8) :: mass(ARGONSITE)
    mass(1) = 18d0
  end subroutine Argon_GetMass

  subroutine UAMethane_GetMass( mass )
    implicit none
    real(kind=8) :: mass(ARGONSITE)
    mass(1) = 16d0
  end subroutine UAMethane_GetMass

  subroutine IdealGas_GetMass( mass )
    implicit none
    real(kind=8) :: mass(ARGONSITE)
    mass(1) = 1d0
  end subroutine IdealGas_GetMass
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine ljatom_setinteraction(si,eps,sig)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    real(kind=8),intent(in)     :: eps,sig
    call si_allocate(si, ARGONSITE, LJ_COULOMB)
    !charges [Q]
    si%param(1,1) = 0d0
    si%param(1,2) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = eps
    si%param(2,2) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = sig
    si%param(3,2) = 0d0
  end subroutine ljatom_setinteraction

  subroutine argon_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call ljatom_setinteraction(si,EPSAR,SIGAR)
  end subroutine argon_setinteraction

  subroutine uamethane_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call ljatom_setinteraction(si,EPSME,SIGME)
  end subroutine uamethane_setinteraction

  subroutine idealgas_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call ljatom_setinteraction(si,EPSIG,SIGIG)
  end subroutine idealgas_setinteraction

  subroutine reducear_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call ljatom_setinteraction(si,EPSRA,SIGRA)
  end subroutine reducear_setinteraction

end module argon_module

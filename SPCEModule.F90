! -*- f90 -*-
module spce_module
  use common_module
  use physconst_module
  use interaction_module
  implicit none
  real(kind=8), parameter :: theta=109.5d0*PI/180d0
  real(kind=8), parameter :: l1=1d0
  real(kind=8), parameter :: hmass=1d0,omass=16d0
  ! J.Phys.Chem.B105, 12093 (2001)
  real(kind=8), parameter :: Q_SPCE_O =-0.8476d0
  real(kind=8), parameter :: Q_SPCE_H =+0.4238d0
  real(kind=8), parameter :: EPS_SPCE_O =(78.208d0/J2K * 1d-3)
  real(kind=8), parameter :: SIG_SPCE_O =3.166d0
  !重心を含めて4点
  integer, parameter :: SPCESITE = 4
  integer, parameter :: DUMBWATERSITE = 4
  character(len=8),parameter     :: waterName(SPCESITE) = (/"O", "H", "H", " "/)
  character(len=8), parameter :: spce_id08 = "SPC_E   "
  character(len=8), parameter :: DUMBWATER_ID08 = "DUMBWATR"
  
contains
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine spce_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, SPCESITE, LJ_COULOMB)
    !charges [Q]
    si%param(1,1) = Q_SPCE_O
    si%param(1,2) = Q_SPCE_H
    si%param(1,3) = Q_SPCE_H
    si%param(1,4) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = EPS_SPCE_O
    si%param(2,2) = 0d0
    si%param(2,3) = 0d0
    si%param(2,4) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = SIG_SPCE_O
    si%param(3,2) = 0d0
    si%param(3,3) = 0d0
    si%param(3,4) = 0d0
  end subroutine spce_setinteraction

  subroutine Rigid_SPCE_Constructor(r)
    use physconst_module
    use rigid_module
    type(sRigid),intent(INOUT) :: r
    real(kind=8) :: zshift
    integer :: i
    ! 4 sites; OHH and center of mass
    r%molx(1)=0d0
    r%moly(1)=0d0
    r%molz(1)=0d0
    r%molx(2)=0d0
    r%moly(2)=l1*dsin(0.5d0*theta)
    r%molz(2)=l1*dcos(0.5d0*theta)
    r%molx(3)=0d0
    r%moly(3)=-l1*dsin(0.5d0*theta)
    r%molz(3)=l1*dcos(0.5d0*theta)
    zshift=r%molz(2)*2d0*hmass/(2d0*hmass+omass)
    do i=1,3
       r%molz(i)=r%molz(i)-zshift
    enddo
    r%molx(4)=0d0
    r%moly(4)=0d0
    r%molz(4)=0d0
    r%Ixx = omass*(r%moly(1)**2+r%molz(1)**2)+hmass*(r%moly(2)**2+r&
         & %molz(2)**2)*2d0
    r%Iyy = omass*(r%molz(1)**2+r%molx(1)**2)+hmass*(r%molz(2)**2+r&
         & %molx(2)**2)*2d0
    r%Izz = omass*(r%molx(1)**2+r%moly(1)**2)+hmass*(r%molx(2)**2+r&
         & %moly(2)**2)*2d0
    r%ixxi = 1d0/r%ixx
    r%iyyi = 1d0/r%iyy
    r%izzi = 1d0/r%izz
    r%mass = 2d0*hmass + omass
    r%massi= 1d0/r%mass
  end subroutine Rigid_SPCE_Constructor

  subroutine SPCE_GetMass( mass )
    real(kind=8) :: mass(*)
    mass(1) = omass
    mass(2) = hmass
    mass(3) = hmass
    mass(4) = 0d0
  end subroutine SPCE_GetMass

  subroutine Flex_Rattle_Register_SPCE(flex, mol, rattle)
    use flex_module
    use rattle_module
    use mol_module
    implicit none
    type(sFlex), intent(INOUT) :: flex
    type(sMol),  intent(IN)    :: mol
    type(sRattle), intent(OUT)  :: rattle
    integer :: nsite,i,m,base
    nsite = mol%nmol*mol%nsite 
    call rattle_constructor( rattle, nsite, nsite*3 )
    !
    !サイトを全部登録して、内部番号をうけとる。
    !
    do i=1, nsite
       flex%xinfo(i)%rattlesite = rattle_register_site(rattle, i+mol%offset, flex%mass(i) )
    enddo
    !
    !拘束を登録する。
    !
    do m=1, mol%nmol
       base = mol%offset + (m-1) * mol%nsite
       call rattle_register_pair( rattle, base+1, base+2, l1 )
       call rattle_register_pair( rattle, base+1, base+3, l1 )
       call rattle_register_pair( rattle, base+2, base+3, l1 * 2d0 * sin( theta / 2d0 ) )
    enddo
  end subroutine Flex_Rattle_Register_SPCE

  !
  !StdInteractionを使う場合の初期化
  !
  subroutine dumbwater_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, DUMBWATERSITE, LJ_COULOMB)
    !charges [Q]
    si%param(1,1) = 0d0
    si%param(1,2) = 0d0
    si%param(1,3) = 0d0
    si%param(1,4) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = 0d0
    si%param(2,2) = 0d0
    si%param(2,3) = 0d0
    si%param(2,4) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = 0d0
    si%param(3,2) = 0d0
    si%param(3,3) = 0d0
    si%param(3,4) = 0d0
  end subroutine dumbwater_setinteraction


end module spce_module

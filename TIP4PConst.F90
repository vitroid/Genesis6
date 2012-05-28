module tip4p_constant_module
  use physconst_module
  implicit none
  real(kind=8), parameter :: AA = (600000d0*CA * 100d0)
  real(kind=8), parameter :: BB = (-610d0*CA * 100d0)
  !     corresponding eps and sig are: 0.1550kcal/mol & 3.1535AA
  real(kind=8), parameter :: WM = (-1.04d0)
  real(kind=8), parameter :: WH = (+0.52d0)
  real(kind=8), parameter :: WWMMright = (WM*WM*EE**2*NA/(4d0*PI*EPS&
       & *1d-10*1000d0) * 100d0)
  real(kind=8), parameter :: WWMM = (WM*WM*COEFF * 100d0)
  real(kind=8), parameter :: WWMH = (WM*WH*COEFF * 100d0)
  real(kind=8), parameter :: WWHH = (WH*WH*COEFF * 100d0)
  !TIP4P Water's LJ parameter (decomposed)
  real(kind=8), parameter :: EPSWO_org =(0.1550d0 * CA)
  real(kind=8), parameter :: SIGWO_org =3.1535d0
  real(kind=8), parameter :: EPSWO =(0.15504166666666666666666d0 * CA)
  real(kind=8), parameter :: SIGWO =3.15357794197649532464d0
  !重心を含めて5点
  integer, parameter :: TIP4PSITE = 5
  !
  character(len=8),parameter :: tip4pName(TIP4PSITE)=(/"O", " ", "H", "H", " "/)
end module tip4p_constant_module

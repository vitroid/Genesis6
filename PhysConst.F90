! -*- f90 -*-
module physconst_module
  !use common_module
  real(kind=8),parameter  :: TOL=1d-12,TOL2=TOL*2d0
  real(kind=8) Gear5_0,Gear5_1,Gear5_3,Gear5_4,Gear5_5
  real(kind=8) Gear4_0,Gear4_2,Gear4_3,Gear4_4
  parameter(Gear5_0= (3d0/16d0))
  parameter(Gear5_1= (251d0/360d0))
  parameter(Gear5_3= (11d0/18d0))
  parameter(Gear5_4= (1d0/6d0))
  parameter(Gear5_5= (1d0/60d0))
  parameter(Gear4_0= (251d0/720d0))
  parameter(Gear4_2= (11d0/12d0))
  parameter(Gear4_3= (1d0/3d0))
  parameter(Gear4_4= (1d0/24d0))
  
  real(kind=8), parameter :: PIright =  3.14159265358979324d0
  real(kind=8), parameter :: PItanaka = 3.1415926535897932d0
  real(kind=8), parameter :: PI = PItanaka
  real(kind=8), parameter :: deg2rad = PI/180d0
  real(kind=8), parameter :: EE = 1.6021892d-19
  real(kind=8), parameter :: NA = 6.0225d23
  !real(kind=8), parameter :: NA = 1d0

  real(kind=8), parameter :: kb = 1.380066d-23
  !real(kind=8), parameter :: kb = 1d0
  real(kind=8), parameter :: rgas=na*kb
  real(kind=8), parameter :: CAright = 4.186d0
  !American Calorie
  real(kind=8), parameter :: CAtanaka = 4.184d0
  real(kind=8), parameter :: CA = CAtanaka
  real(kind=8), parameter :: EPSright = 0.885418782d-11
  real(kind=8), parameter :: EPStanaka = (EE*EE*NA*1d7/(4d0*PI*CA *332.17752d0))
  real(kind=8), parameter :: EPS = EPStanaka
  ! angstrom -> meter
  real(kind=8), parameter :: ANG2M = 1d-10
  !real(kind=8), parameter :: ANG2M = 1d0

  ! m -> km, J -> kJ, etc
  real(kind=8), parameter :: X2KX  = 1d-3
  !real(kind=8), parameter :: X2KX  = 1d0
  real(kind=8), parameter :: KX2X  = 1d0 / X2KX
  real(kind=8), parameter :: COEFFright = (EE*EE*NA/(4d0*PI*EPS*ANG2M * 1000d0))
  real(kind=8), parameter :: COEFFtanaka = (CA*332.17752d0)
  real(kind=8), parameter :: COEFF = COEFFtanaka

  !Internal <--> Joule/mol
  real(kind=8), parameter :: I2J=10d0
  !real(kind=8), parameter :: I2J=1d0
  real(kind=8), parameter :: J2I=1.0d0/I2J
  !Pascal <--> Atm
  real(kind=8), parameter :: a2p=101326D0
  !real(kind=8), parameter :: a2p=1d0
  real(kind=8), parameter :: p2a=(1D0/a2p)
  !Joule/mol <--> Kelvin ( = 1/rgas )
  real(kind=8), parameter :: J2K=0.120273d0
  
  real(kind=8), parameter :: ROT0 =(- 561.677d0  *100d0)
  real(kind=8), parameter :: ROT1 =( 1289.25d0   *100d0)
  real(kind=8), parameter :: ROT2 =(- 955.296d0  *100d0)
  real(kind=8), parameter :: ROT3 =(  291.073d0  *100d0)
  real(kind=8), parameter :: ROT4 =(-  31.6509d0 *100d0)

end module physconst_module

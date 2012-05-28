! -*- f90 -*-

!
!THFの構造と相互作用。Platheのパラメータによる13サイトモデル
!(so called Faller model)
!ID is FLATFTHF
!
module flatfthf_module
  implicit none
  real(kind=8),private,parameter :: QC1 = 0.228d0
  real(kind=8),private,parameter :: QC2 = 0.061d0
  real(kind=8),private,parameter :: QO  =-0.577d0
  real(kind=8),private,parameter :: QH  = 0d0
  real(kind=8),private,parameter :: SIGC = 3.06d0
  real(kind=8),private,parameter :: SIGO = 1.93d0
  real(kind=8),private,parameter :: SIGH = 2.43d0
  real(kind=8),private,parameter :: EPSC = 0.290d0
  real(kind=8),private,parameter :: EPSO = 0.509d0
  real(kind=8),private,parameter :: EPSH = 0.200d0
  real(kind=8),private,parameter :: MASSC = 12.01d0
  real(kind=8),private,parameter :: MASSO = 15.9949d0
  real(kind=8),private,parameter :: MASSH = 1.00787d0
  real(kind=8),private,parameter :: COC   = 111.2 + 1.59 !+ 2.76
  real(kind=8),private,parameter :: OCC   = 106.1 + 1.59 !+ 2.76
  real(kind=8),private,parameter :: CCC   = 101.4 !+ 2.76
  real(kind=8),private,parameter :: OCH   = 109.15
  real(kind=8),private,parameter :: CO    = 1.428
  real(kind=8),private,parameter :: CH    = 1.115
  real(kind=8),private,parameter :: CC    = 1.536
  integer, parameter             :: THFSITE = 14
  character(len=8),parameter     :: thfName(THFSITE) = (/"O", "C", "H", "H", "C", "H", "H", "C", "H", "H", "C", "H", "H", " "/)
  real(kind=8),private,parameter :: mass(THFSITE) = (/MASSO, MASSC, MASSH, MASSH, MASSC, MASSH, MASSH, MASSC, MASSH, MASSH, MASSC, MASSH, MASSH, 0d0/)

contains
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine flatfthf_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, THFSITE, LJ_COULOMB)
    !charges [Q]
    si%param(1,1) = QO
    si%param(1,2) = QC1
    si%param(1,3) = QH
    si%param(1,4) = QH
    si%param(1,5) = QC2
    si%param(1,6) = QH
    si%param(1,7) = QH
    si%param(1,8) = QC2
    si%param(1,9) = QH
    si%param(1,10) = QH
    si%param(1,11) = QC1
    si%param(1,12) = QH
    si%param(1,13) = QH
    si%param(1,14) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = EPSO
    si%param(2,2) = EPSC
    si%param(2,3) = EPSH
    si%param(2,4) = EPSH
    si%param(2,5) = EPSC
    si%param(2,6) = EPSH
    si%param(2,7) = EPSH
    si%param(2,8) = EPSC
    si%param(2,9) = EPSH
    si%param(2,10) = EPSH
    si%param(2,11) = EPSC
    si%param(2,12) = EPSH
    si%param(2,13) = EPSH
    si%param(2,14) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = SIGO
    si%param(3,2) = SIGC
    si%param(3,3) = SIGH
    si%param(3,4) = SIGH
    si%param(3,5) = SIGC
    si%param(3,6) = SIGH
    si%param(3,7) = SIGH
    si%param(3,8) = SIGC
    si%param(3,9) = SIGH
    si%param(3,10) = SIGH
    si%param(3,11) = SIGC
    si%param(3,12) = SIGH
    si%param(3,13) = SIGH
    si%param(3,14) = 0d0
  end subroutine flatfthf_setinteraction

  subroutine Rigid_FLATFTHF_Constructor(r)
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

    r%molx(11) =-CO * dsin( COC*0.5d0*deg2rad )
    r%moly(11) = CO * dcos( COC*0.5d0*deg2rad )
    r%molz(11) = 0d0

    r%molx(5)  = r%molx(2) + CC * dcos( ( 270d0 - OCC - COC*0.5d0 ) * deg2rad )
    r%moly(5)  = r%moly(2) + CC * dsin( ( 270d0 - OCC - COC*0.5d0 ) * deg2rad )
    r%molz(5)  = 0d0

    r%molx(8)  = r%molx(11) - CC * dcos( ( 270d0 - OCC - COC*0.5d0 ) * deg2rad )
    r%moly(8)  = r%moly(11) + CC * dsin( ( 270d0 - OCC - COC*0.5d0 ) * deg2rad )
    r%molz(8)  = 0d0

    first%vec(1) = r%molx(1) - r%molx(2)
    first%vec(2) = r%moly(1) - r%moly(2)
    first%vec(3) = r%molz(1) - r%molz(2)
    second%vec(1) = r%molx(5) - r%molx(2)
    second%vec(2) = r%moly(5) - r%moly(2)
    second%vec(3) = r%molz(5) - r%molz(2)
    call thirdvector1( first, second, OCH*deg2rad, OCH*deg2rad, third )
    call vector_scale1( third, CH )
    r%molx(3) = third%vec(1) + r%molx(2)
    r%moly(3) = third%vec(2) + r%moly(2)
    r%molz(3) = third%vec(3) + r%molz(2)
    r%molx(4) = third%vec(1) + r%molx(2)
    r%moly(4) = third%vec(2) + r%moly(2)
    r%molz(4) =-third%vec(3) + r%molz(2)

    first%vec(1) = r%molx(2) - r%molx(5)
    first%vec(2) = r%moly(2) - r%moly(5)
    first%vec(3) = r%molz(2) - r%molz(5)
    second%vec(1) = r%molx(8) - r%molx(5)
    second%vec(2) = r%moly(8) - r%moly(5)
    second%vec(3) = r%molz(8) - r%molz(5)
    call thirdvector1( first, second, OCH*deg2rad, OCH*deg2rad, third )
    call vector_scale1( third, CH )
    r%molx(6) = third%vec(1) + r%molx(5)
    r%moly(6) = third%vec(2) + r%moly(5)
    r%molz(6) = third%vec(3) + r%molz(5)
    r%molx(7) = third%vec(1) + r%molx(5)
    r%moly(7) = third%vec(2) + r%moly(5)
    r%molz(7) =-third%vec(3) + r%molz(5)

    first%vec(1) = r%molx(5) - r%molx(8)
    first%vec(2) = r%moly(5) - r%moly(8)
    first%vec(3) = r%molz(5) - r%molz(8)
    second%vec(1) = r%molx(11) - r%molx(8)
    second%vec(2) = r%moly(11) - r%moly(8)
    second%vec(3) = r%molz(11) - r%molz(8)
    call thirdvector1( first, second, OCH*deg2rad, OCH*deg2rad, third )
    call vector_scale1( third, CH )
    r%molx(9) = third%vec(1) + r%molx(8)
    r%moly(9) = third%vec(2) + r%moly(8)
    r%molz(9) = third%vec(3) + r%molz(8)
    r%molx(10) = third%vec(1) + r%molx(8)
    r%moly(10) = third%vec(2) + r%moly(8)
    r%molz(10) =-third%vec(3) + r%molz(8)

    first%vec(1) = r%molx(8) - r%molx(11)
    first%vec(2) = r%moly(8) - r%moly(11)
    first%vec(3) = r%molz(8) - r%molz(11)
    second%vec(1) = r%molx(1) - r%molx(11)
    second%vec(2) = r%moly(1) - r%moly(11)
    second%vec(3) = r%molz(1) - r%molz(11)
    call thirdvector1( first, second, OCH*deg2rad, OCH*deg2rad, third )
    call vector_scale1( third, CH )
    r%molx(12) = third%vec(1) + r%molx(11)
    r%moly(12) = third%vec(2) + r%moly(11)
    r%molz(12) = third%vec(3) + r%molz(11)
    r%molx(13) = third%vec(1) + r%molx(11)
    r%moly(13) = third%vec(2) + r%moly(11)
    r%molz(13) =-third%vec(3) + r%molz(11)

    comx = 0d0
    comy = 0d0
    comz = 0d0
    r%mass = 0d0
    do i=1,THFSITE
       comx = comx + r%molx(i) * mass(i)
       comy = comy + r%moly(i) * mass(i)
       comz = comz + r%molz(i) * mass(i)
       r%mass = r%mass + mass(i)
    enddo
    r%massi= 1d0/r%mass
    comx = comz * r%massi
    comy = comy * r%massi
    comz = comz * r%massi
    do i=1,THFSITE
       r%molx(i) = r%molx(i) - comx
       r%moly(i) = r%moly(i) - comy
       r%molz(i) = r%molz(i) - comz
    enddo
    r%molx(14) = 0d0
    r%moly(14) = 0d0
    r%molz(14) = 0d0

    !comx = 0d0
    !comy = 0d0
    !comz = 0d0
    !do i=1,THFSITE
    !   comx = comx + r%molx(i) * mass(i)
    !   comy = comy + r%moly(i) * mass(i)
    !   comz = comz + r%molz(i) * mass(i)
    !enddo
    !write(STDERR,*) "CHECK: ", comx,comy,comz

    r%Ixx = 0d0
    r%Iyy = 0d0
    r%Izz = 0d0
    do i=1,THFSITE
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
  end subroutine Rigid_FLATFTHF_Constructor

end module flatfthf_module
!
!http://alphix6.mpip-mainz.mpg.de/~mplathe/forcefields/thf.html
!
!Forcefield for tetrahydrofurane
!Keys: tetrahydrofurane, THF 
!Short description
!explicit atom force-field for THF 
!Use/purpose
!Thermodynamic properties of THF at ambient conditions 
!Source
!R. Faller, H. Schmitz, O. Biermann, and F. Muller-Plathe ``Automatic Parameterization of Force Fields for Liquids by Simplex Optimization'', J. Comp. Chem. 20, 1009-1017 (1999). 
!Parameters
!for use with the YASP MD program
!
!Atoms:
!         mass        epsilon              sigma       q
!         [amu]       [kJ/mol]              [nm]      [e]
!C        12.01        0.290               0.306       0.228/0.061 (*)
!O        15.9949      0.509               0.193      -0.577
!H         1.00787     0.200               0.243       0.0
!
!(*) the lower value is used for the carbons not connected directly to
!the oxygen
!Note: The sigmas of O and H have been inadvertently
!interchanged in Table 3 of the source.
!
!
!Bond lengths (constrained using SHAKE)
!            length[nm]
!CO            0.1428
!CH            0.1115
!CC            0.1536
!
!harmonic bond angles: (all force constants: k = 450kJ/mol)
!              phi0         
!COC           111.2
!OCC           106.1
!CCC           101.4
!OCH           109.0/109.3
!C3C2H         111.0/113.2
!C2C3H         110.4(2x)/112.8(2x)
!C3C4H         110.4(2x)/113.7(2x)
!HC2H          108.2
!HC3H          108.1
!If there are more than one angle per type they are distributed to the
!different hydrogens for sake of consistency.
!
!Last modified: Fri Nov 10 08:35:42 MET 2000 

!平面剛体分子にするために
!内角の和が540度に足りない。
!111.2 + 106.1*2 + 101.4*2 = 526.2 = 540 - 2.76*5
!
!そこで内角をすこし大きくする。

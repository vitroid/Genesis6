! -*- f90 -*-

!
!Methaneの構造と相互作用。SSPKのパラメータによる5サイトモデル
!Yaxiong Sun, David Spellmeyer, David A. Pearlman, Peter Kollman;
!J. Am. Chem. Soc.; 1992; 114(17); 6798-6801.
!ID is SSPKMET_
!
module sspkmethane_module
  use physconst_module
  implicit none
  !論文には2通りの電荷の値が掲載されている。ここではESPの方を使用。
  !このパラメータはおかしい。分子半径が過小評価される。もしかしたらσ
  !が半分になっているかもしれない。何にせよ、Referenceをちゃんと掲載し
  !ていないのでチェックのしようがない。平成15年11月17日(月)
  real(kind=8),private,parameter :: QC  =-0.464d0
  real(kind=8),private,parameter :: QH  = 0.116d0
  real(kind=8),private,parameter :: SIGC = 1.9082d0
  real(kind=8),private,parameter :: SIGH = 1.4872d0
  real(kind=8),private,parameter :: EPSC = 0.1094 * CA
  real(kind=8),private,parameter :: EPSH = 0.0157 * CA
  real(kind=8),private,parameter :: MASSC = 12.01d0
  real(kind=8),private,parameter :: MASSH = 1.00787d0
  ! http://michele.usc.edu/105a/bonding/structure.html
  ! これで本当に良いかどうかは不明。SSPK論文には長さが明記されていない。
  real(kind=8),private,parameter :: CH    = 1.10d0
  integer, parameter             :: SSPKMETHANESITE = 6
  character(len=8),parameter     :: methaneName(SSPKMETHANESITE) = (/"C", "H", "H", "H", "H", " "/)

contains
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine sspkmethane_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, SSPKMETHANESITE, LJ_COULOMB)
    !charges [Q]
    si%param(1,1) = QC
    si%param(1,2) = QH
    si%param(1,3) = QH
    si%param(1,4) = QH
    si%param(1,5) = QH
    si%param(1,6) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = EPSC
    si%param(2,2) = EPSH
    si%param(2,3) = EPSH
    si%param(2,4) = EPSH
    si%param(2,5) = EPSH
    si%param(2,6) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = SIGC
    si%param(3,2) = SIGH
    si%param(3,3) = SIGH
    si%param(3,4) = SIGH
    si%param(3,5) = SIGH
    si%param(3,6) = 0d0
  end subroutine sspkmethane_setinteraction

  subroutine Rigid_SSPKMETHANE_Constructor(r)
    use physconst_module
    use rigid_module
    type(sRigid),intent(INOUT) :: r
    real(kind=8),dimension(SSPKMETHANESITE) :: mass
    real(kind=8) :: delta
    real(kind=8) :: comx,comy,comz
    integer :: i
    mass(1)=MASSC
    mass(2)=MASSH
    mass(3)=MASSH
    mass(4)=MASSH
    mass(5)=MASSH
    mass(6)=0d0

    delta = CH / sqrt(3d0)

    r%molx(1)  = 0d0
    r%moly(1)  = 0d0
    r%molz(1)  = 0d0

    r%molx(2)  = -delta
    r%moly(2)  = -delta
    r%molz(2)  = -delta

    r%molx(3)  = -delta
    r%moly(3)  = +delta
    r%molz(3)  = +delta

    r%molx(4)  = +delta
    r%moly(4)  = -delta
    r%molz(4)  = +delta

    r%molx(5)  = +delta
    r%moly(5)  = +delta
    r%molz(5)  = -delta

    r%molx(6) = 0d0
    r%moly(6) = 0d0
    r%molz(6) = 0d0

    !確認のため重心位置を計算する。
    comx = 0d0
    comy = 0d0
    comz = 0d0
    r%mass = 0d0
    do i=1,SSPKMETHANESITE
       comx = comx + r%molx(i) * mass(i)
       comy = comy + r%moly(i) * mass(i)
       comz = comz + r%molz(i) * mass(i)
       r%mass = r%mass + mass(i)
    enddo
    r%massi= 1d0/r%mass
    comx = comz * r%massi
    comy = comy * r%massi
    comz = comz * r%massi
    write(STDERR,*) comx,comy,comz
    do i=1,SSPKMETHANESITE
       r%molx(i) = r%molx(i) - comx
       r%moly(i) = r%moly(i) - comy
       r%molz(i) = r%molz(i) - comz
    enddo

    r%Ixx = 0d0
    r%Iyy = 0d0
    r%Izz = 0d0
    do i=1,4
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
  end subroutine Rigid_SSPKMETHANE_Constructor

end module sspkmethane_module

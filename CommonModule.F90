! -*- f90 -*-
module common_module
! -*- f90 -*-
!ここで指定している定数については、実際の配列の使用量との差を簡単にモ
!ニタできるようにしておくとよい。配列は固定的にとらないと、なかなかベ
!クトル機で速度が出ないので、いたしかたない。

!------------------- Variables for grid tesselation ---------------------------
!Number of grids in one direction. 16 for 4096 particle system
!ACTEST/GRID_PR3を使う場合は、これを大きくしてもさほどメモリを消費しな
!い。
!integer,parameter :: MAXGRID=7-1    !for 4096/rc7AA
!integer,parameter :: MAXGRID=4-1    !for 512/rc6AA
!integer,parameter :: MAXGRID=13-1  !for 4096 traj2ngph
integer,parameter :: MAXGRIDX=16  !for si512 1.8A
integer,parameter :: MAXGRIDY=16  !for si512 1.8A
integer,parameter :: MAXGRIDZ=16  !for si512 1.8A
!integer,parameter :: MAXGRID=1-1           
!integer,parameter :: MAXGRID=4-1           
!Number of molecules in a grid cell. 60 for 4096/(8x8x8)
!integer,parameter :: MAXINGRID=30!60
!integer,parameter :: MAXINGRID=60!512rc8
integer,parameter :: MAXINGRID=500
!Number of neighboring cells. 300 for rc=9.0
!For mixture
integer,parameter :: MAXNEIBORCELL=27
!integer,parameter :: MAXNEIBORCELL=300

!------------------- Variables for system definition --------------------------
!Maximum number of molecules in a system.
!integer,parameter :: MAXmol=4097
!integer,parameter :: MAXmol=64
integer,parameter :: MAXmol=2049
!integer,parameter :: MAXmol=513
!integer,parameter :: MAXmol=769
!integer,parameter :: MAXmol=217
!Maximum number of sites in a molecule. 5 for TIP4P water. (including CoM)
!6 for TIP5P
!14 for THF
integer,parameter :: MAXsite=14
!Maximum number of interaction partner. 300 for rc=9.0
!セル分割を粗くする場合、この数字をもっと大きくとる必要がある。
!integer,parameter :: MAXpartner=217!300
!integer,parameter :: MAXpartner=50!300
integer,parameter :: MAXpartner=800

!------------------- Variables for constraint ---------------------------------
!Maximum number of constraint pair in a molecule. 5 for EtOH
integer,parameter :: MAXconst=5
!Maximum number of partner sites constrained with a site.
integer,parameter :: MAXconst_partner=3

!------------------- Variables for Ewald's sum --------------------------------
!Maximum number of vectors
integer,parameter :: MAXK=150

!------------------- Automatically calculated constants.  ---------------------
integer,parameter :: MAXsite_total=MAXsite*MAXmol
integer,parameter :: MAXconst_pair=MAXmol*MAXconst
integer,parameter :: MAXconst_site=MAXmol*MAXsite
!
!DEC F90ではここで宣言するとInternal Errorをひきおこすため、
!InteractionModule.F90で再宣言している。将来は動的に割当てるべき。
!
integer,parameter :: Maxpair = MAXmol*MAXpartner
!integer,parameter :: MAXpair=60000!for 512
!integer,parameter :: MAXpair=60000!for 512rc7
!integer,parameter :: MAXpair=200000!for 768 rc=9
!integer,parameter :: MAXpair=10000!for 96 rc=7
!integer,parameter :: MAXpair=1000000!for 4096 rc=7
!integer,parameter :: MAXpair=300000
!integer,parameter :: MAXpair=30000
!integer,parameter :: MAXpair=140000!for 4096 traj2ngph
  integer,parameter :: MAXCELL=MAXGRIDX * MAXGRIDY * MAXGRIDZ + 1
  integer :: stdin,stdout,trajin,trajout,vstrout
  parameter(stdin=5,stdout=6,trajin=12,trajout=13,vstrout=11)
  !
  !成分数と、相互作用の組合せの数の最大値
  !
  integer, parameter :: MAXCOMPO=3, MAXCOMBI=MAXCOMPO*(MAXCOMPO+1)/2

  character(len=100) ::      fmt0 = '(99(e25.17))'
  character(len=100) ::      fmt1 = '(i11,1x,99(e25.17))'
  character(len=100) ::      fmt2 = '(2(i11,1x),99(e25.17))'

contains
  subroutine writetag( file, tag )
    implicit none
    integer, intent(IN)          :: file
    character(len=5), intent(IN) :: tag
    write( file, '(a5)' ) tag
  end subroutine writetag

end module common_module

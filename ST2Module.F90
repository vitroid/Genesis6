! -*- f90 -*-
!Stillinger and Rahman ST2 See archive #00028.
!
!
module st2_module
  use common_module
  use physconst_module
  use interaction_module
  implicit none
  real(kind=8), parameter :: theta=109.47d0*PI/180d0
  real(kind=8), parameter :: l1=1d0
  real(kind=8), parameter :: l2=0.8d0
  real(kind=8), parameter :: hmass=1d0,omass=16d0
  ! J.Phys.Chem.B105, 12093 (2001)
  real(kind=8), parameter :: Q_ST2_M =-0.2357
  real(kind=8), parameter :: Q_ST2_H =+0.2357
  real(kind=8), parameter :: EPS_ST2_O =0.52605e-21*NA*1d-3
  !from Stillinger and Rahman's original Paper JCP 60, 1545 (1974)
  !real(kind=8), parameter :: EPS_ST2_O =7.5750d-2 * CA
  real(kind=8), parameter :: SIG_ST2_O =3.100d0
  ! JCP 60,1545(1974)
  real(kind=8), parameter :: RL = 2.0160d0
  real(kind=8), parameter :: RU = 3.1287d0
  !重心を含めて6点
  integer, parameter :: ST2SITE = 6
  character(len=8),parameter     :: st2Name(ST2SITE) = (/"O", "H", "H", "M", "M", " "/)
  
contains
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine st2_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, ST2SITE, LJ_COULOMB)
    !charges [Q]
    si%param(1,1) = 0d0
    si%param(1,2) = Q_ST2_H
    si%param(1,3) = Q_ST2_H
    si%param(1,4) = Q_ST2_M
    si%param(1,5) = Q_ST2_M
    si%param(1,6) = 0d0
    !LJ-eps [kJ/mol]
    !write(STDERR,*) EPS_ST2_O
    si%param(2,1) = EPS_ST2_O
    si%param(2,2) = 0d0
    si%param(2,3) = 0d0
    si%param(2,4) = 0d0
    si%param(2,5) = 0d0
    si%param(2,6) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = SIG_ST2_O
    si%param(3,2) = 0d0
    si%param(3,3) = 0d0
    si%param(3,4) = 0d0
    si%param(3,5) = 0d0
    si%param(3,6) = 0d0
  end subroutine st2_setinteraction

  subroutine Rigid_ST2_Constructor(r)
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
    r%molx(4)=l2*dsin(0.5d0*theta)
    r%moly(4)=0d0
    r%molz(4)=-l2*dcos(0.5d0*theta)
    r%molx(5)=-l2*dsin(0.5d0*theta)
    r%moly(5)=0d0
    r%molz(5)=-l2*dcos(0.5d0*theta)
    zshift=r%molz(2)*2d0*hmass/(2d0*hmass+omass)
    do i=1,ST2SITE-1
       r%molz(i)=r%molz(i)-zshift
    enddo
    r%molx(6)=0d0
    r%moly(6)=0d0
    r%molz(6)=0d0
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
  end subroutine Rigid_ST2_Constructor

  subroutine ST2_GetMass( mass )
    real(kind=8) :: mass(*)
    mass(1) = omass
    mass(2) = hmass
    mass(3) = hmass
    mass(4) = 0d0
    mass(5) = 0d0
    mass(6) = 0d0
  end subroutine ST2_GetMass

  subroutine force_st2(iv,mi,mj,site,ep1,vrmat)
    use mol_module
    use site_module
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    real(kind=8),intent(inout) :: ep1(*)
    real(kind=8),intent(OUT)   :: vrmat(3,3)

    real(kind=8), parameter :: QMM = ( Q_ST2_M * Q_ST2_M * COEFF * 100d0)
    real(kind=8), parameter :: QHM = ( Q_ST2_H * Q_ST2_M * COEFF * 100d0)
    real(kind=8), parameter :: QHH = ( Q_ST2_H * Q_ST2_H * COEFF * 100d0)
    real(kind=8), parameter :: AA = ( 4d0*EPS_ST2_O*SIG_ST2_O**12 * 100d0)
    real(kind=8), parameter :: BB = (-4d0*EPS_ST2_O*SIG_ST2_O** 6 * 100d0)

    real(kind=8) :: fxi1,fyi1,fzi1
    real(kind=8) :: fxi2,fyi2,fzi2
    real(kind=8) :: fxi3,fyi3,fzi3
    real(kind=8) :: fxi4,fyi4,fzi4
    real(kind=8) :: fxi5,fyi5,fzi5
    real(kind=8) :: fxi6,fyi6,fzi6
    real(kind=8) :: fxj1,fyj1,fzj1
    real(kind=8) :: fxj2,fyj2,fzj2
    real(kind=8) :: fxj3,fyj3,fzj3
    real(kind=8) :: fxj4,fyj4,fzj4
    real(kind=8) :: fxj5,fyj5,fzj5
    real(kind=8) :: fxj6,fyj6,fzj6
    real(kind=8) :: dx,dy,dz,dx00,dy00,dz00
    real(kind=8) :: dhx,dhy,dhz
    real(kind=8) :: d6,d12,e0,f0,radiusi
    real(kind=8) :: xi,yi,zi,xj,yj,zj
    real(kind=8) forcex,forcey,forcez
    real(kind=8) ox,oy,oz
    real(kind=8) cx,cy,cz
    real(kind=8) ep
    real(kind=8) dd,hh,dh
    real(kind=8) :: vrx,vry,vrz
    integer :: i,j,k
    integer :: ii,jj
    real(kind=8) :: r12, s12, ds12 !, s12d
    real(kind=8) :: dhc, hhc, eplj
    real(kind=8) :: dhcx,dhcy,dhcz
    real(kind=8) :: c,d
#ifdef VERBOSE
    integer :: count
#endif
    !integer :: kk
#ifdef VPOPTIMIZE
    !0番目の要素は、loop内のifをへらすためのトリック
    real(kind=8),dimension(:),allocatable :: fxi,fyi,fzi,fxj,fyj,fzj
    !real(kind=8) :: r_r0,r_r1,r_r0s,r_r1s,dr
    allocate(fxi(0:ST2SITE*(iv%npair+1)))
    allocate(fyi(0:ST2SITE*(iv%npair+1)))
    allocate(fzi(0:ST2SITE*(iv%npair+1)))
    allocate(fxj(0:ST2SITE*(iv%npair+1)))
    allocate(fyj(0:ST2SITE*(iv%npair+1)))
    allocate(fzj(0:ST2SITE*(iv%npair+1)))
    fxi(:)=0d0
    fyi(:)=0d0
    fzi(:)=0d0
    fxj(:)=0d0
    fyj(:)=0d0
    fzj(:)=0d0
#endif /*VPOPTIMIZE*/
#ifdef VERBOSE
    count=0
#endif
    ep1(1:iv%npair) = 0d0
    vrmat(:,:) = 0d0
    !kk=0
  !VPP用のベクトル化制御行
  !OCL VECTOR,REPEAT(1000000)
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+ST2SITE)
       dy00 = site%y(ii+ST2SITE)
       dz00 = site%z(ii+ST2SITE)
       dx00 = dx00 - site%x(jj+ST2SITE)
       dy00 = dy00 - site%y(jj+ST2SITE)
       dz00 = dz00 - site%z(jj+ST2SITE)
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       dx00 = dx00-ox
       dy00 = dy00-oy
       dz00 = dz00-oz
       hh   = iv%eratio( k )
       dh   = iv%fratio( k )
       !
       !ST2/BNC特有の、近距離カットオフ関数
       !
       d12  = (dx00**2 + dy00**2 + dz00**2)
       if ( d12 < ru**2 ) then
          if ( d12 < rl**2 ) then
             s12 = 0d0
             ds12 = 0d0
          else
             r12  = sqrt(d12)
             s12  = (r12 - rl)**2*(3d0*ru - rl - 2d0*r12)/(ru-rl)**3
             ds12 = 6d0*(r12-rl)*(ru-r12)/(ru-rl)**3
             !
             !エネルギーが収束しない原因は、ds12がスムーズな関数でないせいかもしれないので、まずは
             !遅いがスムーズな関数に置きかえる。しかしやはりエネルギーの保存が悪い。ということは、式の導出に問題があると考えるべき。
             !
             !c = 1d0/(2d0*PI)
             !d = -rl/(ru-rl)
             !s12 = -c*sin( (r12 - rl)/(ru-rl) * 2.0* PI ) + c*r12*2d0*pi/(ru-rl) + d
             !ds12 = -c*2.0*PI/(ru-rl)*cos( (r12 - rl)/(ru-rl) * 2.0* PI ) + c*2d0*pi/(ru-rl)
             !微分が正しいことを確認平成16年2月12日(木)
             !r12 = r12 - 0.00001d0
             !s12d  = (r12 - rl)**2*(3d0*ru - rl - 2d0*r12)/(ru-rl)**3
             !write(6,*) (s12 - s12d)/0.00001d0, ds12
             ds12 = -ds12 / r12
          endif
       else
          s12 = 1d0
          ds12 = 0d0
       endif
       !s12 = 1d0
       !ds12 = 0d0
       dhc = dh * s12 + hh * ds12
       hhc = hh * s12
       !write(11,'(5(e17.10,1x))') sqrt(d12), hhc,dhc*sqrt(d12),hh,dh*sqrt(d12)
       dhcx  = dhc*dx00
       dhcy  = dhc*dy00
       dhcz  = dhc*dz00
       dhx  = dh*dx00
       dhy  = dh*dy00
       dhz  = dh*dz00
#ifdef VERBOSE
          count=count+1
#endif
          cx = -ox
          cy = -oy
          cz = -oz
          fxi1 = 0d0
          fyi1 = 0d0
          fzi1 = 0d0
          fxi2 = 0d0
          fyi2 = 0d0
          fzi2 = 0d0
          fxi3 = 0d0
          fyi3 = 0d0
          fzi3 = 0d0
          fxi4 = 0d0
          fyi4 = 0d0
          fzi4 = 0d0
          fxi5 = 0d0
          fyi5 = 0d0
          fzi5 = 0d0
          fxi6 = 0d0
          fyi6 = 0d0
          fzi6 = 0d0
          fxj1 = 0d0
          fyj1 = 0d0
          fzj1 = 0d0
          fxj2 = 0d0
          fyj2 = 0d0
          fzj2 = 0d0
          fxj3 = 0d0
          fyj3 = 0d0
          fzj3 = 0d0
          fxj4 = 0d0
          fyj4 = 0d0
          fzj4 = 0d0
          fxj5 = 0d0
          fyj5 = 0d0
          fzj5 = 0d0
          fxj6 = 0d0
          fyj6 = 0d0
          fzj6 = 0d0
  
          ep=0d0
          !サイトの順序をTIP4Pから変更し、OHHMMとする。
          !HHMMの相互作用が16通り生じる。
          xi = site%x(ii+2)
          yi = site%y(ii+2)
          zi = site%z(ii+2)
          xj = site%x(jj+2)
          yj = site%y(jj+2)
          zj = site%z(jj+2)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHH*radiusi
          f0 = QHH*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi2 = fxi2 + forcex
          fyi2 = fyi2 + forcey
          fzi2 = fzi2 + forcez
          fxj2 = fxj2 + forcex
          fyj2 = fyj2 + forcey
          fzj2 = fzj2 + forcez
          xi = site%x(ii+2)
          yi = site%y(ii+2)
          zi = site%z(ii+2)
          xj = site%x(jj+3)
          yj = site%y(jj+3)
          zj = site%z(jj+3)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHH*radiusi
          f0 = QHH*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi2 = fxi2 + forcex
          fyi2 = fyi2 + forcey
          fzi2 = fzi2 + forcez
          fxj3 = fxj3 + forcex
          fyj3 = fyj3 + forcey
          fzj3 = fzj3 + forcez
          xi = site%x(ii+3)
          yi = site%y(ii+3)
          zi = site%z(ii+3)
          xj = site%x(jj+2)
          yj = site%y(jj+2)
          zj = site%z(jj+2)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHH*radiusi
          f0 = QHH*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi3 = fxi3 + forcex
          fyi3 = fyi3 + forcey
          fzi3 = fzi3 + forcez
          fxj2 = fxj2 + forcex
          fyj2 = fyj2 + forcey
          fzj2 = fzj2 + forcez
          xi = site%x(ii+3)
          yi = site%y(ii+3)
          zi = site%z(ii+3)
          xj = site%x(jj+3)
          yj = site%y(jj+3)
          zj = site%z(jj+3)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHH*radiusi
          f0 = QHH*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi3 = fxi3 + forcex
          fyi3 = fyi3 + forcey
          fzi3 = fzi3 + forcez
          fxj3 = fxj3 + forcex
          fyj3 = fyj3 + forcey
          fzj3 = fzj3 + forcez
  
  
          xi = site%x(ii+4)
          yi = site%y(ii+4)
          zi = site%z(ii+4)
          xj = site%x(jj+4)
          yj = site%y(jj+4)
          zj = site%z(jj+4)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QMM*radiusi
          f0 = QMM*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi4 = fxi4 + forcex
          fyi4 = fyi4 + forcey
          fzi4 = fzi4 + forcez
          fxj4 = fxj4 + forcex
          fyj4 = fyj4 + forcey
          fzj4 = fzj4 + forcez
          xi = site%x(ii+4)
          yi = site%y(ii+4)
          zi = site%z(ii+4)
          xj = site%x(jj+5)
          yj = site%y(jj+5)
          zj = site%z(jj+5)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QMM*radiusi
          f0 = QMM*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi4 = fxi4 + forcex
          fyi4 = fyi4 + forcey
          fzi4 = fzi4 + forcez
          fxj5 = fxj5 + forcex
          fyj5 = fyj5 + forcey
          fzj5 = fzj5 + forcez
          xi = site%x(ii+5)
          yi = site%y(ii+5)
          zi = site%z(ii+5)
          xj = site%x(jj+4)
          yj = site%y(jj+4)
          zj = site%z(jj+4)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QMM*radiusi
          f0 = QMM*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi5 = fxi5 + forcex
          fyi5 = fyi5 + forcey
          fzi5 = fzi5 + forcez
          fxj4 = fxj4 + forcex
          fyj4 = fyj4 + forcey
          fzj4 = fzj4 + forcez
          xi = site%x(ii+5)
          yi = site%y(ii+5)
          zi = site%z(ii+5)
          xj = site%x(jj+5)
          yj = site%y(jj+5)
          zj = site%z(jj+5)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QMM*radiusi
          f0 = QMM*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi5 = fxi5 + forcex
          fyi5 = fyi5 + forcey
          fzi5 = fzi5 + forcez
          fxj5 = fxj5 + forcex
          fyj5 = fyj5 + forcey
          fzj5 = fzj5 + forcez
  
  
  
          xi = site%x(ii+2)
          yi = site%y(ii+2)
          zi = site%z(ii+2)
          xj = site%x(jj+4)
          yj = site%y(jj+4)
          zj = site%z(jj+4)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          f0 = QHM*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi2 = fxi2 + forcex
          fyi2 = fyi2 + forcey
          fzi2 = fzi2 + forcez
          fxj4 = fxj4 + forcex
          fyj4 = fyj4 + forcey
          fzj4 = fzj4 + forcez
          xi = site%x(ii+2)
          yi = site%y(ii+2)
          zi = site%z(ii+2)
          xj = site%x(jj+5)
          yj = site%y(jj+5)
          zj = site%z(jj+5)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          f0 = QHM*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi2 = fxi2 + forcex
          fyi2 = fyi2 + forcey
          fzi2 = fzi2 + forcez
          fxj5 = fxj5 + forcex
          fyj5 = fyj5 + forcey
          fzj5 = fzj5 + forcez
          xi = site%x(ii+3)
          yi = site%y(ii+3)
          zi = site%z(ii+3)
          xj = site%x(jj+4)
          yj = site%y(jj+4)
          zj = site%z(jj+4)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          f0 = QHM*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi3 = fxi3 + forcex
          fyi3 = fyi3 + forcey
          fzi3 = fzi3 + forcez
          fxj4 = fxj4 + forcex
          fyj4 = fyj4 + forcey
          fzj4 = fzj4 + forcez
          xi = site%x(ii+3)
          yi = site%y(ii+3)
          zi = site%z(ii+3)
          xj = site%x(jj+5)
          yj = site%y(jj+5)
          zj = site%z(jj+5)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          f0 = QHM*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi3 = fxi3 + forcex
          fyi3 = fyi3 + forcey
          fzi3 = fzi3 + forcez
          fxj5 = fxj5 + forcex
          fyj5 = fyj5 + forcey
          fzj5 = fzj5 + forcez
          xi = site%x(ii+4)
          yi = site%y(ii+4)
          zi = site%z(ii+4)
          xj = site%x(jj+2)
          yj = site%y(jj+2)
          zj = site%z(jj+2)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          f0 = QHM*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi4 = fxi4 + forcex
          fyi4 = fyi4 + forcey
          fzi4 = fzi4 + forcez
          fxj2 = fxj2 + forcex
          fyj2 = fyj2 + forcey
          fzj2 = fzj2 + forcez
          xi = site%x(ii+4)
          yi = site%y(ii+4)
          zi = site%z(ii+4)
          xj = site%x(jj+3)
          yj = site%y(jj+3)
          zj = site%z(jj+3)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          f0 = QHM*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi4 = fxi4 + forcex
          fyi4 = fyi4 + forcey
          fzi4 = fzi4 + forcez
          fxj3 = fxj3 + forcex
          fyj3 = fyj3 + forcey
          fzj3 = fzj3 + forcez
          xi = site%x(ii+5)
          yi = site%y(ii+5)
          zi = site%z(ii+5)
          xj = site%x(jj+2)
          yj = site%y(jj+2)
          zj = site%z(jj+2)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          f0 = QHM*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi5 = fxi5 + forcex
          fyi5 = fyi5 + forcey
          fzi5 = fzi5 + forcez
          fxj2 = fxj2 + forcex
          fyj2 = fyj2 + forcey
          fzj2 = fzj2 + forcez
          xi = site%x(ii+5)
          yi = site%y(ii+5)
          zi = site%z(ii+5)
          xj = site%x(jj+3)
          yj = site%y(jj+3)
          zj = site%z(jj+3)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          f0 = QHM*radiusi*dd
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi5 = fxi5 + forcex
          fyi5 = fyi5 + forcey
          fzi5 = fzi5 + forcez
          fxj3 = fxj3 + forcex
          fyj3 = fyj3 + forcey
          fzj3 = fzj3 + forcez

          xi = site%x(ii+1)
          yi = site%y(ii+1)
          zi = site%z(ii+1)
          xj = site%x(jj+1)
          yj = site%y(jj+1)
          zj = site%z(jj+1)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          d6 = dd*dd*dd
          d12= d6*d6
          e0 = AA*d12 + BB*d6
          f0 = dd*(AA*12d0*d12 + BB*6d0*d6)
          eplj = e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi1 = fxi1 + forcex
          fyi1 = fyi1 + forcey
          fzi1 = fzi1 + forcez
          fxj1 = fxj1 + forcex
          fyj1 = fyj1 + forcey
          fzj1 = fzj1 + forcez
  
          fxi1=fxi1*hh
          fyi1=fyi1*hh
          fzi1=fzi1*hh
          fxi2=fxi2*hhc
          fyi2=fyi2*hhc
          fzi2=fzi2*hhc
          fxi3=fxi3*hhc
          fyi3=fyi3*hhc
          fzi3=fzi3*hhc
          fxi4=fxi4*hhc
          fyi4=fyi4*hhc
          fzi4=fzi4*hhc
          fxi5=fxi5*hhc
          fyi5=fyi5*hhc
          fzi5=fzi5*hhc
          fxi6=fxi6*hh
          fyi6=fyi6*hh
          fzi6=fzi6*hh
  
          fxj1=fxj1*hh
          fyj1=fyj1*hh
          fzj1=fzj1*hh
          fxj2=fxj2*hhc
          fyj2=fyj2*hhc
          fzj2=fzj2*hhc
          fxj3=fxj3*hhc
          fyj3=fyj3*hhc
          fzj3=fzj3*hhc
          fxj4=fxj4*hhc
          fyj4=fyj4*hhc
          fzj4=fzj4*hhc
          fxj5=fxj5*hhc
          fyj5=fyj5*hhc
          fzj5=fzj5*hhc
          fxj6=fxj6*hh
          fyj6=fyj6*hh
          fzj6=fzj6*hh
          dhx=dhx*eplj
          dhy=dhy*eplj
          dhz=dhz*eplj
          dhcx=dhcx*ep
          dhcy=dhcy*ep
          dhcz=dhcz*ep
          fxi6=fxi6+dhx+dhcx
          fyi6=fyi6+dhy+dhcy
          fzi6=fzi6+dhz+dhcz
          fxj6=fxj6+dhx+dhcx
          fyj6=fyj6+dhy+dhcy
          fzj6=fzj6+dhz+dhcz
          ep1(k)=ep*hhc + eplj*hh 
          vrx = (fxi1+fxi2+fxi3+fxi4+fxi5+fxi6)
          vry = (fyi1+fyi2+fyi3+fyi4+fyi5+fyi6)
          vrz = (fzi1+fzi2+fzi3+fzi4+fzi5+fzi6)
          vrmat(1,1) = vrmat(1,1) + dx00 * vrx
          vrmat(1,2) = vrmat(1,2) + dx00 * vry
          vrmat(1,3) = vrmat(1,3) + dx00 * vrz
          vrmat(2,1) = vrmat(2,1) + dy00 * vrx
          vrmat(2,2) = vrmat(2,2) + dy00 * vry
          vrmat(2,3) = vrmat(2,3) + dy00 * vrz
          vrmat(3,1) = vrmat(3,1) + dz00 * vrx
          vrmat(3,2) = vrmat(3,2) + dz00 * vry
          vrmat(3,3) = vrmat(3,3) + dz00 * vrz
          !vrsum=vrsum+vrx+vry+vrz
          !write(6,*) i,j,k,ep
          !     非vector計算機の場合には、この場で力の集計をした方がてっとり
          !     早い。
#ifdef VPOPTIMIZE
          !     ......計算結果をlist vectorに格納する。
          fxi(k*ST2SITE+0)=fxi1
          fxi(k*ST2SITE+1)=fxi2
          fxi(k*ST2SITE+2)=fxi3
          fxi(k*ST2SITE+3)=fxi4
          fxi(k*ST2SITE+4)=fxi5
          fxi(k*ST2SITE+5)=fxi6
          fyi(k*ST2SITE+0)=fyi1
          fyi(k*ST2SITE+1)=fyi2
          fyi(k*ST2SITE+2)=fyi3
          fyi(k*ST2SITE+3)=fyi4
          fyi(k*ST2SITE+4)=fyi5
          fyi(k*ST2SITE+5)=fyi6
          fzi(k*ST2SITE+0)=fzi1
          fzi(k*ST2SITE+1)=fzi2
          fzi(k*ST2SITE+2)=fzi3
          fzi(k*ST2SITE+3)=fzi4
          fzi(k*ST2SITE+4)=fzi5
          fzi(k*ST2SITE+5)=fzi6
  
          fxj(k*ST2SITE+0)=-fxj1
          fxj(k*ST2SITE+1)=-fxj2
          fxj(k*ST2SITE+2)=-fxj3
          fxj(k*ST2SITE+3)=-fxj4
          fxj(k*ST2SITE+4)=-fxj5
          fxj(k*ST2SITE+5)=-fxj6
          fyj(k*ST2SITE+0)=-fyj1
          fyj(k*ST2SITE+1)=-fyj2
          fyj(k*ST2SITE+2)=-fyj3
          fyj(k*ST2SITE+3)=-fyj4
          fyj(k*ST2SITE+4)=-fyj5
          fyj(k*ST2SITE+5)=-fyj6
          fzj(k*ST2SITE+0)=-fzj1
          fzj(k*ST2SITE+1)=-fzj2
          fzj(k*ST2SITE+2)=-fzj3
          fzj(k*ST2SITE+3)=-fzj4
          fzj(k*ST2SITE+4)=-fzj5
          fzj(k*ST2SITE+5)=-fzj6
#else 
          site%fx(ii+1)=site%fx(ii+1)+fxi1
          site%fx(ii+2)=site%fx(ii+2)+fxi2
          site%fx(ii+3)=site%fx(ii+3)+fxi3
          site%fx(ii+4)=site%fx(ii+4)+fxi4
          site%fx(ii+5)=site%fx(ii+5)+fxi5
          site%fx(ii+6)=site%fx(ii+6)+fxi6
          site%fx(jj+1)=site%fx(jj+1)-fxj1
          site%fx(jj+2)=site%fx(jj+2)-fxj2
          site%fx(jj+3)=site%fx(jj+3)-fxj3
          site%fx(jj+4)=site%fx(jj+4)-fxj4
          site%fx(jj+5)=site%fx(jj+5)-fxj5
          site%fx(jj+6)=site%fx(jj+6)-fxj6
  
          site%fy(ii+1)=site%fy(ii+1)+fyi1
          site%fy(ii+2)=site%fy(ii+2)+fyi2
          site%fy(ii+3)=site%fy(ii+3)+fyi3
          site%fy(ii+4)=site%fy(ii+4)+fyi4
          site%fy(ii+5)=site%fy(ii+5)+fyi5
          site%fy(ii+6)=site%fy(ii+6)+fyi6
          site%fy(jj+1)=site%fy(jj+1)-fyj1
          site%fy(jj+2)=site%fy(jj+2)-fyj2
          site%fy(jj+3)=site%fy(jj+3)-fyj3
          site%fy(jj+4)=site%fy(jj+4)-fyj4
          site%fy(jj+5)=site%fy(jj+5)-fyj5
          site%fy(jj+6)=site%fy(jj+6)-fyj6
  
          site%fz(ii+1)=site%fz(ii+1)+fzi1
          site%fz(ii+2)=site%fz(ii+2)+fzi2
          site%fz(ii+3)=site%fz(ii+3)+fzi3
          site%fz(ii+4)=site%fz(ii+4)+fzi4
          site%fz(ii+5)=site%fz(ii+5)+fzi5
          site%fz(ii+6)=site%fz(ii+6)+fzi6
          site%fz(jj+1)=site%fz(jj+1)-fzj1
          site%fz(jj+2)=site%fz(jj+2)-fzj2
          site%fz(jj+3)=site%fz(jj+3)-fzj3
          site%fz(jj+4)=site%fz(jj+4)-fzj4
          site%fz(jj+5)=site%fz(jj+5)-fzj5
          site%fz(jj+6)=site%fz(jj+6)-fzj6
          !f90  -O6 2.724u 0.020s 0:03.05 89.8% 0+11k 6+13io 0pf+0w
          !kf90 -O6 2.751u 0.014s 0:03.12 88.4% 0+11k 7+11io 0pf+0w
          !この部分を除去すると
          !kf90 -O6 1.427u 0.021s 0:01.77 81.3% 0+11k 7+11io 0pf+0w
#endif /*VPOPTIMIZE*/
       !endif
    enddo
#ifdef VPOPTIMIZE
    !     各サイトに加わる力を集計する。Vectorizeのためのテクニック
    !iv%partner_iの値が0の要素は集計しないようにしているので、最初に0ク
    !リアが必要。
    !lとiのループの順序を逆にするとさらに遅くなる。
    !xyzの処理を分割すると若干速くなる(しかし可読性がおちる)
    !この部分の処理が1/3以上の時間を消費。
    !相互作用が全てのさいとに力をおよぼすとはかぎらない。
    !内側のloopは手で展開しておいてもいいのではないか？
  
    !平成１２年５月２０日(土)この部分遅い。その上、並列化しても、計算量
    !が減らない(0要素が増えるだけ)ので、効率を落としている。
    !この方法自体を見直すのが一つ。
    !相互作用表をできるだけ均質にする方法が一つ。
    !たぶん、均質性を仮定している限り、将来もここの速度が問題になる(溶液
    !などを扱う場合)
    do l=1,iv%maxpartner_i
       do ii=1,mi%nsite*mi%nmol
          !i=(ii-1)/5+1
          !s=mod(ii-1,5)+1
          i=mi%lvmol(ii)
          s=mi%lvsite(ii)
          k=iv%partner_i(i,l)
          !if(k.ne.0)then
          site%fx(ii+mi%offset)=site%fx(ii+mi%offset)+fxi(k*ST2SITE+s-1)
          site%fy(ii+mi%offset)=site%fy(ii+mi%offset)+fyi(k*ST2SITE+s-1)
          site%fz(ii+mi%offset)=site%fz(ii+mi%offset)+fzi(k*ST2SITE+s-1)
          !endif
       enddo
    enddo
    do l=1,iv%maxpartner_j
       do jj=1,mj%nsite*mj%nmol
          !j=(jj-1)/5+1
          !s=mod(jj-1,5)+1
          j=mj%lvmol(jj)
          s=mj%lvsite(jj)
          k=iv%partner_j(j,l)
          !平成１２年４月２４日(月)ifがないほうが2割速い。
          !0要素を用いるトリックは結構使えるかも。
          !if(k.ne.0)then
          site%fx(jj+mj%offset)=site%fx(jj+mj%offset)+fxj(k*ST2SITE+s-1)
          site%fy(jj+mj%offset)=site%fy(jj+mj%offset)+fyj(k*ST2SITE+s-1)
          site%fz(jj+mj%offset)=site%fz(jj+mj%offset)+fzj(k*ST2SITE+s-1)
          !endif
       enddo
    enddo
    deallocate(fxi)
    deallocate(fyi)
    deallocate(fzi)
    deallocate(fxj)
    deallocate(fyj)
    deallocate(fzj)
#endif /*VPOPTIMIZE*/
    !write(6,*) kk
#ifdef VERBOSE
    write(STDERR,*) "INTERACTION PAIRS: ",count
#endif
  end subroutine force_st2


  subroutine pairpotential_st2(iv,mi,mj,site,ep1 )
    use mol_module
    use site_module
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    real(kind=8),intent(inout) :: ep1(*)

    real(kind=8), parameter :: QMM = ( Q_ST2_M * Q_ST2_M * COEFF * 100d0)
    real(kind=8), parameter :: QHM = ( Q_ST2_H * Q_ST2_M * COEFF * 100d0)
    real(kind=8), parameter :: QHH = ( Q_ST2_H * Q_ST2_H * COEFF * 100d0)
    real(kind=8), parameter :: AA = ( 4d0*EPS_ST2_O*SIG_ST2_O**12 * 100d0)
    real(kind=8), parameter :: BB = (-4d0*EPS_ST2_O*SIG_ST2_O** 6 * 100d0)

    real(kind=8) :: dx,dy,dz,dx00,dy00,dz00
    real(kind=8) :: dhx,dhy,dhz
    real(kind=8) :: d6,d12,e0,radiusi
    real(kind=8) :: xi,yi,zi,xj,yj,zj
    real(kind=8) ox,oy,oz
    real(kind=8) cx,cy,cz
    real(kind=8) ep
    real(kind=8) dd,hh,dh
    integer :: i,j,k
    integer :: ii,jj
    real(kind=8) :: r12, s12
    real(kind=8) :: hhc, eplj
    real(kind=8) :: c,d
#ifdef VERBOSE
    integer :: count
#endif
    !integer :: kk
#ifdef VPOPTIMIZE
    !0番目の要素は、loop内のifをへらすためのトリック
    real(kind=8),dimension(:),allocatable :: fxi,fyi,fzi,fxj,fyj,fzj
    !real(kind=8) :: r_r0,r_r1,r_r0s,r_r1s,dr
#endif /*VPOPTIMIZE*/
#ifdef VERBOSE
    count=0
#endif
    ep1(1:iv%npair) = 0d0
    !kk=0
  !VPP用のベクトル化制御行
  !OCL VECTOR,REPEAT(1000000)
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+ST2SITE)
       dy00 = site%y(ii+ST2SITE)
       dz00 = site%z(ii+ST2SITE)
       dx00 = dx00 - site%x(jj+ST2SITE)
       dy00 = dy00 - site%y(jj+ST2SITE)
       dz00 = dz00 - site%z(jj+ST2SITE)
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       dx00 = dx00-ox
       dy00 = dy00-oy
       dz00 = dz00-oz
       hh   = iv%eratio( k )
       !
       !ST2/BNC特有の、近距離カットオフ関数
       !
       d12  = (dx00**2 + dy00**2 + dz00**2)
       if ( d12 < ru**2 ) then
          if ( d12 < rl**2 ) then
             s12 = 0d0
          else
             r12  = sqrt(d12)
             s12  = (r12 - rl)**2*(3d0*ru - rl - 2d0*r12)/(ru-rl)**3
          endif
       else
          s12 = 1d0
       endif
       hhc = hh * s12
#ifdef VERBOSE
          count=count+1
#endif
          cx = -ox
          cy = -oy
          cz = -oz
          ep=0d0
          !サイトの順序をTIP4Pから変更し、OHHMMとする。
          !HHMMの相互作用が16通り生じる。
          xi = site%x(ii+2)
          yi = site%y(ii+2)
          zi = site%z(ii+2)
          xj = site%x(jj+2)
          yj = site%y(jj+2)
          zj = site%z(jj+2)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHH*radiusi
          ep = ep + e0
          xi = site%x(ii+2)
          yi = site%y(ii+2)
          zi = site%z(ii+2)
          xj = site%x(jj+3)
          yj = site%y(jj+3)
          zj = site%z(jj+3)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHH*radiusi
          ep = ep + e0
          xi = site%x(ii+3)
          yi = site%y(ii+3)
          zi = site%z(ii+3)
          xj = site%x(jj+2)
          yj = site%y(jj+2)
          zj = site%z(jj+2)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHH*radiusi
          ep = ep + e0
          xi = site%x(ii+3)
          yi = site%y(ii+3)
          zi = site%z(ii+3)
          xj = site%x(jj+3)
          yj = site%y(jj+3)
          zj = site%z(jj+3)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHH*radiusi
          ep = ep + e0
  
          xi = site%x(ii+4)
          yi = site%y(ii+4)
          zi = site%z(ii+4)
          xj = site%x(jj+4)
          yj = site%y(jj+4)
          zj = site%z(jj+4)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QMM*radiusi
          ep = ep + e0
          xi = site%x(ii+4)
          yi = site%y(ii+4)
          zi = site%z(ii+4)
          xj = site%x(jj+5)
          yj = site%y(jj+5)
          zj = site%z(jj+5)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QMM*radiusi
          ep = ep + e0
          xi = site%x(ii+5)
          yi = site%y(ii+5)
          zi = site%z(ii+5)
          xj = site%x(jj+4)
          yj = site%y(jj+4)
          zj = site%z(jj+4)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QMM*radiusi
          ep = ep + e0
          xi = site%x(ii+5)
          yi = site%y(ii+5)
          zi = site%z(ii+5)
          xj = site%x(jj+5)
          yj = site%y(jj+5)
          zj = site%z(jj+5)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QMM*radiusi
          ep = ep + e0
  
          xi = site%x(ii+2)
          yi = site%y(ii+2)
          zi = site%z(ii+2)
          xj = site%x(jj+4)
          yj = site%y(jj+4)
          zj = site%z(jj+4)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          ep = ep + e0
          xi = site%x(ii+2)
          yi = site%y(ii+2)
          zi = site%z(ii+2)
          xj = site%x(jj+5)
          yj = site%y(jj+5)
          zj = site%z(jj+5)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          ep = ep + e0
          xi = site%x(ii+3)
          yi = site%y(ii+3)
          zi = site%z(ii+3)
          xj = site%x(jj+4)
          yj = site%y(jj+4)
          zj = site%z(jj+4)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          ep = ep + e0
          xi = site%x(ii+3)
          yi = site%y(ii+3)
          zi = site%z(ii+3)
          xj = site%x(jj+5)
          yj = site%y(jj+5)
          zj = site%z(jj+5)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          ep = ep + e0
          xi = site%x(ii+4)
          yi = site%y(ii+4)
          zi = site%z(ii+4)
          xj = site%x(jj+2)
          yj = site%y(jj+2)
          zj = site%z(jj+2)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          ep = ep + e0
          xi = site%x(ii+4)
          yi = site%y(ii+4)
          zi = site%z(ii+4)
          xj = site%x(jj+3)
          yj = site%y(jj+3)
          zj = site%z(jj+3)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          ep = ep + e0
          xi = site%x(ii+5)
          yi = site%y(ii+5)
          zi = site%z(ii+5)
          xj = site%x(jj+2)
          yj = site%y(jj+2)
          zj = site%z(jj+2)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          ep = ep + e0
          xi = site%x(ii+5)
          yi = site%y(ii+5)
          zi = site%z(ii+5)
          xj = site%x(jj+3)
          yj = site%y(jj+3)
          zj = site%z(jj+3)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = QHM*radiusi
          ep = ep + e0

          xi = site%x(ii+1)
          yi = site%y(ii+1)
          zi = site%z(ii+1)
          xj = site%x(jj+1)
          yj = site%y(jj+1)
          zj = site%z(jj+1)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          d6 = dd*dd*dd
          d12= d6*d6
          e0 = AA*d12 + BB*d6
          eplj = e0
          ep1(k)=ep*hhc + eplj*hh 
    enddo
  end subroutine pairpotential_st2

end module st2_module

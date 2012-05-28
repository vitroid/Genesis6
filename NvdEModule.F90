! -*- f90 -*-


module nvde_module
  use common_module
  use physconst_module
  use interaction_module
  implicit none
  !Nada Model
  real(kind=8), parameter :: nvde_oh   =   0.980d0
  real(kind=8), parameter :: nvde_om   =   0.230d0
  real(kind=8), parameter :: nvde_ol   =   0.8892d0
  real(kind=8), parameter :: nvde_hoh  = 108.00d0 * PI / 180d0
  real(kind=8), parameter :: nvde_lol  = 111.0d0  * PI / 180d0
  real(kind=8), parameter :: nvde_epso =  85.9766d0 * rgas * 1d-3
  real(kind=8), parameter :: nvde_epsh =  13.8817d0 * rgas * 1d-3
  real(kind=8), parameter :: nvde_sigo =   3.115d0
  real(kind=8), parameter :: nvde_sigh =   0.673d0
  real(kind=8), parameter :: nvde_qh   =   0.477d0
  real(kind=8), parameter :: nvde_qm   =  -0.866d0
  real(kind=8), parameter :: nvde_ql   =  -0.044d0
  real(kind=8), parameter :: nvde_hmass=1d0, nvde_omass=16d0
  !real(kind=8), parameter, private :: epsoh = sqrt( nvde_epso * nvde_epsh )
  real(kind=8), parameter, private :: epsoh = 34.54708914250229563036d0 * rgas * 1d-3
  real(kind=8), parameter, private :: sigoh = ( nvde_sigo + nvde_sigh ) * 0.5d0
  
  real(kind=8), parameter, private :: QQMM = nvde_qm * nvde_qm * COEFF * 100d0
  real(kind=8), parameter, private :: QQMH = nvde_qm * nvde_qh * COEFF * 100d0
  real(kind=8), parameter, private :: QQML = nvde_qm * nvde_ql * COEFF * 100d0
  real(kind=8), parameter, private :: QQHH = nvde_qh * nvde_qh * COEFF * 100d0
  real(kind=8), parameter, private :: QQHL = nvde_qh * nvde_ql * COEFF * 100d0
  real(kind=8), parameter, private :: QQLL = nvde_ql * nvde_ql * COEFF * 100d0
  
  real(kind=8), parameter, private :: AAOO = ( 4d0 * nvde_epso * nvde_sigo**12 * 100d0)
  real(kind=8), parameter, private :: BBOO = (-4d0 * nvde_epso * nvde_sigo** 6 * 100d0)
  real(kind=8), parameter, private :: AAHH = ( 4d0 * nvde_epsh * nvde_sigh**12 * 100d0)
  real(kind=8), parameter, private :: BBHH = (-4d0 * nvde_epsh * nvde_sigh** 6 * 100d0)
  real(kind=8), parameter, private :: AAOH = ( 4d0 * epsoh * sigoh**12 * 100d0)
  real(kind=8), parameter, private :: BBOH = (-4d0 * epsoh * sigoh** 6 * 100d0)
  !real(kind=8), parameter :: CC5P = 0.2410d0
  !comparison test with TIP4P
  !real(kind=8), parameter :: AA5P = ( 600000d0*CA * 100d0)
  !real(kind=8), parameter :: BB5P = (-610d0*CA * 100d0)
  !real(kind=8), parameter :: CC5P = 0.52d0
  !real(kind=8), parameter :: MM5P = ( CC5P*CC5P*COEFF * 100d0)
  !real(kind=8), parameter :: HM5P = (-CC5P*CC5P*COEFF * 100d0)
  !real(kind=8), parameter :: HH5P = ( CC5P*CC5P*COEFF * 100d0)
  !重心を含めて7点
  integer, parameter :: NVDESITE = 7
  character(len=8),parameter :: nvdeName(NVDESITE)=(/"O", "H", "H", "M", "L", "L", " "/)
  real(kind=8),private,parameter :: mass(NVDESITE) = (/nvde_omass, nvde_hmass, nvde_hmass, 0d0, 0d0, 0d0, 0d0/)
  character(len=8), parameter :: nvde_id08 = "NVDE____"

contains
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine nvde_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si

    call si_allocate(si, NVDESITE, LJ_COULOMB)
    !site order: OHHMLL&CoM
    !charges [Q]
    si%param(1,1) = 0d0
    si%param(1,2) = nvde_qh
    si%param(1,3) = nvde_qh
    si%param(1,4) = nvde_qm
    si%param(1,5) = nvde_ql
    si%param(1,6) = nvde_ql
    si%param(1,7) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = nvde_epso
    si%param(2,2) = nvde_epsh
    si%param(2,3) = nvde_epsh
    si%param(2,4) = 0d0
    si%param(2,5) = 0d0
    si%param(2,6) = 0d0
    si%param(2,7) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = nvde_sigo
    si%param(3,2) = nvde_sigh
    si%param(3,3) = nvde_sigh
    si%param(3,4) = 0d0
    si%param(3,5) = 0d0
    si%param(3,6) = 0d0
    si%param(3,7) = 0d0
  end subroutine nvde_setinteraction

  subroutine Rigid_NvdE_Constructor(r)
    use rigid_module
    type(sRigid),intent(INOUT) :: r
    real(kind=8) :: comx,comy,comz
    integer :: i
    ! 7 sites; OHHMLL and center of mass
    r%molx(1) = 0d0
    r%moly(1) = 0d0
    r%molz(1) = 0d0
    r%molx(2) = 0d0
    r%moly(2) = nvde_oh * dsin( 0.5d0 * nvde_hoh )
    r%molz(2) = nvde_oh * dcos( 0.5d0 * nvde_hoh )
    r%molx(3) = 0d0
    r%moly(3) =-nvde_oh * dsin( 0.5d0 * nvde_hoh )
    r%molz(3) = nvde_oh * dcos( 0.5d0 * nvde_hoh )
    r%molx(4) = 0d0
    r%moly(4) = 0d0
    r%molz(4) = nvde_om 
    r%molx(5) = nvde_ol * dsin( 0.5d0 * nvde_lol )
    r%moly(5) = 0d0
    r%molz(5) =-nvde_ol * dcos( 0.5d0 * nvde_lol )
    r%molx(6) =-nvde_ol * dsin( 0.5d0 * nvde_lol )
    r%moly(6) = 0d0
    r%molz(6) =-nvde_ol * dcos( 0.5d0 * nvde_lol )
    comx = 0d0
    comy = 0d0
    comz = 0d0
    r%mass = 0d0
    do i=1,NVDESITE-1
       comx = comx + r%molx(i) * mass(i)
       comy = comy + r%moly(i) * mass(i)
       comz = comz + r%molz(i) * mass(i)
       r%mass = r%mass + mass(i)
    enddo
    r%massi= 1d0/r%mass
    comx = comz * r%massi
    comy = comy * r%massi
    comz = comz * r%massi
    do i=1,NVDESITE - 1
       r%molx(i) = r%molx(i) - comx
       r%moly(i) = r%moly(i) - comy
       r%molz(i) = r%molz(i) - comz
    enddo
    r%molx(NVDESITE) = 0d0
    r%moly(NVDESITE) = 0d0
    r%molz(NVDESITE) = 0d0

    r%Ixx = 0d0
    r%Iyy = 0d0
    r%Izz = 0d0
    do i=1,NVDESITE
       r%Ixx = r%Ixx + mass(i)*(r%moly(i)**2 + r%molz(i)**2)
       r%Iyy = r%Iyy + mass(i)*(r%molz(i)**2 + r%molx(i)**2)
       r%Izz = r%Izz + mass(i)*(r%molx(i)**2 + r%moly(i)**2)
    enddo
    r%ixxi = 1d0/r%ixx
    r%iyyi = 1d0/r%iyy
    r%izzi = 1d0/r%izz
  end subroutine Rigid_NvdE_Constructor


  subroutine force_nvde(iv,mi,mj,site,ep1,vrmat)
    use site_module
    use mol_module
    use cutoff_module
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    real(kind=8),intent(out) :: vrmat(3,3)
    real(kind=8), intent(INOUT) :: ep1(*)
    real(kind=8) :: fxi1,fyi1,fzi1
    real(kind=8) :: fxi2,fyi2,fzi2
    real(kind=8) :: fxi3,fyi3,fzi3
    real(kind=8) :: fxi4,fyi4,fzi4
    real(kind=8) :: fxi5,fyi5,fzi5
    real(kind=8) :: fxi6,fyi6,fzi6
    real(kind=8) :: fxi7,fyi7,fzi7
    real(kind=8) :: fxj1,fyj1,fzj1
    real(kind=8) :: fxj2,fyj2,fzj2
    real(kind=8) :: fxj3,fyj3,fzj3
    real(kind=8) :: fxj4,fyj4,fzj4
    real(kind=8) :: fxj5,fyj5,fzj5
    real(kind=8) :: fxj6,fyj6,fzj6
    real(kind=8) :: fxj7,fyj7,fzj7
#ifdef SHOWFORCE
    real(kind=8) :: fxs,fys,fzs,fpair
#endif
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
    !integer :: kk
#ifdef VERBOSE
    integer :: count
#endif
    !integer :: kk
#ifdef VPOPTIMIZE
    integer :: l,s
    !0番目の要素は、loop内のifをへらすためのトリック
    real(kind=8),dimension(:),allocatable :: fxi,fyi,fzi,fxj,fyj,fzj
    real(kind=8) :: r_r0,r_r1,r_r0s,r_r1s,dr
    allocate(fxi(0:NVDESITE*(iv%npair+1)))
    allocate(fyi(0:NVDESITE*(iv%npair+1)))
    allocate(fzi(0:NVDESITE*(iv%npair+1)))
    allocate(fxj(0:NVDESITE*(iv%npair+1)))
    allocate(fyj(0:NVDESITE*(iv%npair+1)))
    allocate(fzj(0:NVDESITE*(iv%npair+1)))
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
       !平成16年5月19日(水)write文を入れるとAddress Errorを回避できる。なぜ？
       i = iv%pair_i(k)
       j = iv%pair_j(k)
#ifdef VPOPTIMIZE
#else
       if ( i.eq.1 .and. j.eq.0 ) then
          !絶対実行されない文だが効果がなぜかある。
          write(STDERR,*) i,j,k,iv%Npair,mi%offset,mj%offset
       endif
#endif
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+NVDESITE)
       dy00 = site%y(ii+NVDESITE)
       dz00 = site%z(ii+NVDESITE)
       dx00 = dx00 - site%x(jj+NVDESITE)
       dy00 = dy00 - site%y(jj+NVDESITE)
       dz00 = dz00 - site%z(jj+NVDESITE)
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       dx00 = dx00-ox
       dy00 = dy00-oy
       dz00 = dz00-oz
       hh   = iv%eratio( k )
       dh   = iv%fratio( k )
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
       fxi7 = 0d0
       fyi7 = 0d0
       fzi7 = 0d0

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
       fxj7 = 0d0
       fyj7 = 0d0
       fzj7 = 0d0
       
       ep=0d0
       !MM
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
       e0 = QQMM*radiusi
       f0 = QQMM*radiusi*dd
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

       !HH
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
       e0 = QQHH*radiusi
       f0 = QQHH*radiusi*dd
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = e0 + AAHH*d12 + BBHH*d6
       f0 = f0 + dd*(AAHH*12d0*d12 + BBHH*6d0*d6)
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
       e0 = QQHH*radiusi
       f0 = QQHH*radiusi*dd
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = e0 + AAHH*d12 + BBHH*d6
       f0 = f0 + dd*(AAHH*12d0*d12 + BBHH*6d0*d6)
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
       e0 = QQHH*radiusi
       f0 = QQHH*radiusi*dd
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = e0 + AAHH*d12 + BBHH*d6
       f0 = f0 + dd*(AAHH*12d0*d12 + BBHH*6d0*d6)
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
       e0 = QQHH*radiusi
       f0 = QQHH*radiusi*dd
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = e0 + AAHH*d12 + BBHH*d6
       f0 = f0 + dd*(AAHH*12d0*d12 + BBHH*6d0*d6)
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

       !LL
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
       e0 = QQLL*radiusi
       f0 = QQLL*radiusi*dd
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

       xi = site%x(ii+5)
       yi = site%y(ii+5)
       zi = site%z(ii+5)
       xj = site%x(jj+6)
       yj = site%y(jj+6)
       zj = site%z(jj+6)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQLL*radiusi
       f0 = QQLL*radiusi*dd
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi5 = fxi5 + forcex
       fyi5 = fyi5 + forcey
       fzi5 = fzi5 + forcez
       fxj6 = fxj6 + forcex
       fyj6 = fyj6 + forcey
       fzj6 = fzj6 + forcez

       xi = site%x(ii+6)
       yi = site%y(ii+6)
       zi = site%z(ii+6)
       xj = site%x(jj+5)
       yj = site%y(jj+5)
       zj = site%z(jj+5)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQLL*radiusi
       f0 = QQLL*radiusi*dd
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi6 = fxi6 + forcex
       fyi6 = fyi6 + forcey
       fzi6 = fzi6 + forcez
       fxj5 = fxj5 + forcex
       fyj5 = fyj5 + forcey
       fzj5 = fzj5 + forcez

       xi = site%x(ii+6)
       yi = site%y(ii+6)
       zi = site%z(ii+6)
       xj = site%x(jj+6)
       yj = site%y(jj+6)
       zj = site%z(jj+6)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQLL*radiusi
       f0 = QQLL*radiusi*dd
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi6 = fxi6 + forcex
       fyi6 = fyi6 + forcey
       fzi6 = fzi6 + forcez
       fxj6 = fxj6 + forcex
       fyj6 = fyj6 + forcey
       fzj6 = fzj6 + forcez

       !MH
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
       e0 = QQMH*radiusi
       f0 = QQMH*radiusi*dd
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
       e0 = QQMH*radiusi
       f0 = QQMH*radiusi*dd
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
       e0 = QQMH*radiusi
       f0 = QQMH*radiusi*dd
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
       e0 = QQMH*radiusi
       f0 = QQMH*radiusi*dd
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

       !ML
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
       e0 = QQML*radiusi
       f0 = QQML*radiusi*dd
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

       xi = site%x(ii+6)
       yi = site%y(ii+6)
       zi = site%z(ii+6)
       xj = site%x(jj+4)
       yj = site%y(jj+4)
       zj = site%z(jj+4)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQML*radiusi
       f0 = QQML*radiusi*dd
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi6 = fxi6 + forcex
       fyi6 = fyi6 + forcey
       fzi6 = fzi6 + forcez
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
       e0 = QQML*radiusi
       f0 = QQML*radiusi*dd
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

       xi = site%x(ii+4)
       yi = site%y(ii+4)
       zi = site%z(ii+4)
       xj = site%x(jj+6)
       yj = site%y(jj+6)
       zj = site%z(jj+6)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQML*radiusi
       f0 = QQML*radiusi*dd
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi4 = fxi4 + forcex
       fyi4 = fyi4 + forcey
       fzi4 = fzi4 + forcez
       fxj6 = fxj6 + forcex
       fyj6 = fyj6 + forcey
       fzj6 = fzj6 + forcez

       !HL
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
       e0 = QQHL*radiusi
       f0 = QQHL*radiusi*dd
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

       xi = site%x(ii+6)
       yi = site%y(ii+6)
       zi = site%z(ii+6)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQHL*radiusi
       f0 = QQHL*radiusi*dd
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi6 = fxi6 + forcex
       fyi6 = fyi6 + forcey
       fzi6 = fzi6 + forcez
       fxj2 = fxj2 + forcex
       fyj2 = fyj2 + forcey
       fzj2 = fzj2 + forcez

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
       e0 = QQHL*radiusi
       f0 = QQHL*radiusi*dd
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

       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+6)
       yj = site%y(jj+6)
       zj = site%z(jj+6)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQHL*radiusi
       f0 = QQHL*radiusi*dd
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi2 = fxi2 + forcex
       fyi2 = fyi2 + forcey
       fzi2 = fzi2 + forcez
       fxj6 = fxj6 + forcex
       fyj6 = fyj6 + forcey
       fzj6 = fzj6 + forcez

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
       e0 = QQHL*radiusi
       f0 = QQHL*radiusi*dd
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

       xi = site%x(ii+6)
       yi = site%y(ii+6)
       zi = site%z(ii+6)
       xj = site%x(jj+3)
       yj = site%y(jj+3)
       zj = site%z(jj+3)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQHL*radiusi
       f0 = QQHL*radiusi*dd
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi6 = fxi6 + forcex
       fyi6 = fyi6 + forcey
       fzi6 = fzi6 + forcez
       fxj3 = fxj3 + forcex
       fyj3 = fyj3 + forcey
       fzj3 = fzj3 + forcez

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
       e0 = QQHL*radiusi
       f0 = QQHL*radiusi*dd
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

       xi = site%x(ii+3)
       yi = site%y(ii+3)
       zi = site%z(ii+3)
       xj = site%x(jj+6)
       yj = site%y(jj+6)
       zj = site%z(jj+6)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQHL*radiusi
       f0 = QQHL*radiusi*dd
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi3 = fxi3 + forcex
       fyi3 = fyi3 + forcey
       fzi3 = fzi3 + forcez
       fxj6 = fxj6 + forcex
       fyj6 = fyj6 + forcey
       fzj6 = fzj6 + forcez

       !OO
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
       e0 = AAOO*d12 + BBOO*d6
       f0 = dd*(AAOO*12d0*d12 + BBOO*6d0*d6)
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi1 = fxi1 + forcex
       fyi1 = fyi1 + forcey
       fzi1 = fzi1 + forcez
       fxj1 = fxj1 + forcex
       fyj1 = fyj1 + forcey
       fzj1 = fzj1 + forcez

       !OH
       xi = site%x(ii+1)
       yi = site%y(ii+1)
       zi = site%z(ii+1)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = AAOH*d12 + BBOH*d6
       f0 = dd*(AAOH*12d0*d12 + BBOH*6d0*d6)
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi1 = fxi1 + forcex
       fyi1 = fyi1 + forcey
       fzi1 = fzi1 + forcez
       fxj2 = fxj2 + forcex
       fyj2 = fyj2 + forcey
       fzj2 = fzj2 + forcez

       xi = site%x(ii+1)
       yi = site%y(ii+1)
       zi = site%z(ii+1)
       xj = site%x(jj+3)
       yj = site%y(jj+3)
       zj = site%z(jj+3)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = AAOH*d12 + BBOH*d6
       f0 = dd*(AAOH*12d0*d12 + BBOH*6d0*d6)
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi1 = fxi1 + forcex
       fyi1 = fyi1 + forcey
       fzi1 = fzi1 + forcez
       fxj3 = fxj3 + forcex
       fyj3 = fyj3 + forcey
       fzj3 = fzj3 + forcez

       xi = site%x(ii+3)
       yi = site%y(ii+3)
       zi = site%z(ii+3)
       xj = site%x(jj+1)
       yj = site%y(jj+1)
       zj = site%z(jj+1)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = AAOH*d12 + BBOH*d6
       f0 = dd*(AAOH*12d0*d12 + BBOH*6d0*d6)
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi3 = fxi3 + forcex
       fyi3 = fyi3 + forcey
       fzi3 = fzi3 + forcez
       fxj1 = fxj1 + forcex
       fyj1 = fyj1 + forcey
       fzj1 = fzj1 + forcez

       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+1)
       yj = site%y(jj+1)
       zj = site%z(jj+1)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = AAOH*d12 + BBOH*d6
       f0 = dd*(AAOH*12d0*d12 + BBOH*6d0*d6)
       ep = ep + e0
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi2 = fxi2 + forcex
       fyi2 = fyi2 + forcey
       fzi2 = fzi2 + forcez
       fxj1 = fxj1 + forcex
       fyj1 = fyj1 + forcey
       fzj1 = fzj1 + forcez

       !if ( k .eq. 1 ) then
       !   write( STDOUT, * ) fxi1,fxi2,fxi3,fxi4,fxi5,fxi6,fxi7
       !   write( STDOUT, * ) fyi1,fyi2,fyi3,fyi4,fyi5,fyi6,fyi7
       !   write( STDOUT, * ) fzi1,fzi2,fzi3,fzi4,fzi5,fzi6,fzi7
       !   write( STDOUT, * ) fxj1,fxj2,fxj3,fxj4,fxj5,fxj6,fxj7
       !   write( STDOUT, * ) fyj1,fyj2,fyj3,fyj4,fyj5,fyj6,fyj7
       !   write( STDOUT, * ) fzj1,fzj2,fzj3,fzj4,fzj5,fzj6,fzj7
       !endif
       fxi1=fxi1*hh
       fyi1=fyi1*hh
       fzi1=fzi1*hh
       fxi2=fxi2*hh
       fyi2=fyi2*hh
       fzi2=fzi2*hh
       fxi3=fxi3*hh
       fyi3=fyi3*hh
       fzi3=fzi3*hh
       fxi4=fxi4*hh
       fyi4=fyi4*hh
       fzi4=fzi4*hh
       fxi5=fxi5*hh
       fyi5=fyi5*hh
       fzi5=fzi5*hh
       fxi6=fxi6*hh
       fyi6=fyi6*hh
       fzi6=fzi6*hh
       fxi7=fxi7*hh
       fyi7=fyi7*hh
       fzi7=fzi7*hh
       
       fxj1=fxj1*hh
       fyj1=fyj1*hh
       fzj1=fzj1*hh
       fxj2=fxj2*hh
       fyj2=fyj2*hh
       fzj2=fzj2*hh
       fxj3=fxj3*hh
       fyj3=fyj3*hh
       fzj3=fzj3*hh
       fxj4=fxj4*hh
       fyj4=fyj4*hh
       fzj4=fzj4*hh
       fxj5=fxj5*hh
       fyj5=fyj5*hh
       fzj5=fzj5*hh
       fxj6=fxj6*hh
       fyj6=fyj6*hh
       fzj6=fzj6*hh
       fxj7=fxj7*hh
       fyj7=fyj7*hh
       fzj7=fzj7*hh
       dhx=dhx*ep
       dhy=dhy*ep
       dhz=dhz*ep
       fxi7=fxi7+dhx
       fyi7=fyi7+dhy
       fzi7=fzi7+dhz
       fxj7=fxj7+dhx
       fyj7=fyj7+dhy
       fzj7=fzj7+dhz
       ep1(k) = ep*hh
       !vrx = dx00*(fxi1+fxi2+fxi3+fxi4+fxi5)
       !vry = dy00*(fyi1+fyi2+fyi3+fyi4+fyi5)
       !vrz = dz00*(fzi1+fzi2+fzi3+fzi4+fzi5)
       !vrsum=vrsum+vrx+vry+vrz
       vrx = (fxi1+fxi2+fxi3+fxi4+fxi5+fxi6+fxi7)
       vry = (fyi1+fyi2+fyi3+fyi4+fyi5+fyi6+fyi7)
       vrz = (fzi1+fzi2+fzi3+fzi4+fzi5+fzi6+fzi7)
       vrmat(1,1) = vrmat(1,1) + dx00 * vrx
       vrmat(1,2) = vrmat(1,2) + dx00 * vry
       vrmat(1,3) = vrmat(1,3) + dx00 * vrz
       vrmat(2,1) = vrmat(2,1) + dy00 * vrx
       vrmat(2,2) = vrmat(2,2) + dy00 * vry
       vrmat(2,3) = vrmat(2,3) + dy00 * vrz
       vrmat(3,1) = vrmat(3,1) + dz00 * vrx
       vrmat(3,2) = vrmat(3,2) + dz00 * vry
       vrmat(3,3) = vrmat(3,3) + dz00 * vrz
       !if( i.eq.1 .and. j.lt.10 )then
       !   write(STDOUT,*) i,j,vrx,dx00,(site%x(ii+kk),kk=1,5)
       !endif
       !write(6,*) i,j,k,ep
       !
       ! 非vector計算機の場合には、この場で力の集計をした方がてっとり
       !     早い。
#ifdef VPOPTIMIZE
       !     ......計算結果をlist vectorに格納する。
       fxi(k*NVDESITE+0)=fxi1
       fxi(k*NVDESITE+1)=fxi2
       fxi(k*NVDESITE+2)=fxi3
       fxi(k*NVDESITE+3)=fxi4
       fxi(k*NVDESITE+4)=fxi5
       fxi(k*NVDESITE+5)=fxi6
       fxi(k*NVDESITE+6)=fxi7

       fyi(k*NVDESITE+0)=fyi1
       fyi(k*NVDESITE+1)=fyi2
       fyi(k*NVDESITE+2)=fyi3
       fyi(k*NVDESITE+3)=fyi4
       fyi(k*NVDESITE+4)=fyi5
       fyi(k*NVDESITE+5)=fyi6
       fyi(k*NVDESITE+6)=fyi7

       fzi(k*NVDESITE+0)=fzi1
       fzi(k*NVDESITE+1)=fzi2
       fzi(k*NVDESITE+2)=fzi3
       fzi(k*NVDESITE+3)=fzi4
       fzi(k*NVDESITE+4)=fzi5
       fzi(k*NVDESITE+5)=fzi6
       fzi(k*NVDESITE+6)=fzi7
       
       fxj(k*NVDESITE+0)=-fxj1
       fxj(k*NVDESITE+1)=-fxj2
       fxj(k*NVDESITE+2)=-fxj3
       fxj(k*NVDESITE+3)=-fxj4
       fxj(k*NVDESITE+4)=-fxj5
       fxj(k*NVDESITE+5)=-fxj6
       fxj(k*NVDESITE+6)=-fxj7

       fyj(k*NVDESITE+0)=-fyj1
       fyj(k*NVDESITE+1)=-fyj2
       fyj(k*NVDESITE+2)=-fyj3
       fyj(k*NVDESITE+3)=-fyj4
       fyj(k*NVDESITE+4)=-fyj5
       fyj(k*NVDESITE+5)=-fyj6
       fyj(k*NVDESITE+6)=-fyj7

       fzj(k*NVDESITE+0)=-fzj1
       fzj(k*NVDESITE+1)=-fzj2
       fzj(k*NVDESITE+2)=-fzj3
       fzj(k*NVDESITE+3)=-fzj4
       fzj(k*NVDESITE+4)=-fzj5
       fzj(k*NVDESITE+5)=-fzj6
       fzj(k*NVDESITE+6)=-fzj7
#else 
       site%fx(ii+1)=site%fx(ii+1)+fxi1
       site%fx(ii+2)=site%fx(ii+2)+fxi2
       site%fx(ii+3)=site%fx(ii+3)+fxi3
       site%fx(ii+4)=site%fx(ii+4)+fxi4
       site%fx(ii+5)=site%fx(ii+5)+fxi5
       site%fx(ii+6)=site%fx(ii+6)+fxi6
       site%fx(ii+7)=site%fx(ii+7)+fxi7
       site%fx(jj+1)=site%fx(jj+1)-fxj1
       site%fx(jj+2)=site%fx(jj+2)-fxj2
       site%fx(jj+3)=site%fx(jj+3)-fxj3
       site%fx(jj+4)=site%fx(jj+4)-fxj4
       site%fx(jj+5)=site%fx(jj+5)-fxj5
       site%fx(jj+6)=site%fx(jj+6)-fxj6
       site%fx(jj+7)=site%fx(jj+7)-fxj7
       
       site%fy(ii+1)=site%fy(ii+1)+fyi1
       site%fy(ii+2)=site%fy(ii+2)+fyi2
       site%fy(ii+3)=site%fy(ii+3)+fyi3
       site%fy(ii+4)=site%fy(ii+4)+fyi4
       site%fy(ii+5)=site%fy(ii+5)+fyi5
       site%fy(ii+6)=site%fy(ii+6)+fyi6
       site%fy(ii+7)=site%fy(ii+7)+fyi7
       site%fy(jj+1)=site%fy(jj+1)-fyj1
       site%fy(jj+2)=site%fy(jj+2)-fyj2
       site%fy(jj+3)=site%fy(jj+3)-fyj3
       site%fy(jj+4)=site%fy(jj+4)-fyj4
       site%fy(jj+5)=site%fy(jj+5)-fyj5
       site%fy(jj+6)=site%fy(jj+6)-fyj6
       site%fy(jj+7)=site%fy(jj+7)-fyj7

       site%fz(ii+1)=site%fz(ii+1)+fzi1
       site%fz(ii+2)=site%fz(ii+2)+fzi2
       site%fz(ii+3)=site%fz(ii+3)+fzi3
       site%fz(ii+4)=site%fz(ii+4)+fzi4
       site%fz(ii+5)=site%fz(ii+5)+fzi5
       site%fz(ii+6)=site%fz(ii+6)+fzi6
       site%fz(ii+7)=site%fz(ii+7)+fzi7
       site%fz(jj+1)=site%fz(jj+1)-fzj1
       site%fz(jj+2)=site%fz(jj+2)-fzj2
       site%fz(jj+3)=site%fz(jj+3)-fzj3
       site%fz(jj+4)=site%fz(jj+4)-fzj4
       site%fz(jj+5)=site%fz(jj+5)-fzj5
       site%fz(jj+6)=site%fz(jj+6)-fzj6
       site%fz(jj+7)=site%fz(jj+7)-fzj7
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
          site%fx(ii+mi%offset)=site%fx(ii+mi%offset)+fxi(k*NVDESITE&
               & +s-1)
          site%fy(ii+mi%offset)=site%fy(ii+mi%offset)+fyi(k*NVDESITE&
               & +s-1)
          site%fz(ii+mi%offset)=site%fz(ii+mi%offset)+fzi(k*NVDESITE&
               & +s-1)
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
          site%fx(jj+mj%offset)=site%fx(jj+mj%offset)+fxj(k*NVDESITE&
               & +s-1)
          site%fy(jj+mj%offset)=site%fy(jj+mj%offset)+fyj(k*NVDESITE&
               & +s-1)
          site%fz(jj+mj%offset)=site%fz(jj+mj%offset)+fzj(k*NVDESITE&
               & +s-1)
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
  end subroutine force_nvde


  subroutine potential_nvde(iv,mi,mj,site,ep1)
    use site_module
    use mol_module
    use cutoff_module
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    real(kind=8), intent(INOUT) :: ep1(*)
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
    !integer :: kk
#ifdef VERBOSE
    integer :: count
#endif
    !integer :: kk
#ifdef VERBOSE
    count=0
#endif
    ep1(1:iv%npair) = 0d0
    !kk=0
    !VPP用のベクトル化制御行
    !OCL VECTOR,REPEAT(1000000)
    do k=1,iv%Npair
       !平成16年5月19日(水)write文を入れるとAddress Errorを回避できる。なぜ？
       i = iv%pair_i(k)
       j = iv%pair_j(k)
#ifdef VPOPTIMIZE
#else
       if ( i.eq.1 .and. j.eq.0 ) then
          !絶対実行されない文だが効果がなぜかある。
          write(STDERR,*) i,j,k,iv%Npair,mi%offset,mj%offset
       endif
#endif
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+NVDESITE)
       dy00 = site%y(ii+NVDESITE)
       dz00 = site%z(ii+NVDESITE)
       dx00 = dx00 - site%x(jj+NVDESITE)
       dy00 = dy00 - site%y(jj+NVDESITE)
       dz00 = dz00 - site%z(jj+NVDESITE)
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       dx00 = dx00-ox
       dy00 = dy00-oy
       dz00 = dz00-oz
       hh   = iv%eratio( k )
#ifdef VERBOSE
       count=count+1
#endif
       cx = -ox
       cy = -oy
       cz = -oz
       ep=0d0
       !MM
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
       e0 = QQMM*radiusi
       ep = ep + e0

       !HH
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
       e0 = QQHH*radiusi
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = e0 + AAHH*d12 + BBHH*d6
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
       e0 = QQHH*radiusi
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = e0 + AAHH*d12 + BBHH*d6
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
       e0 = QQHH*radiusi
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = e0 + AAHH*d12 + BBHH*d6
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
       e0 = QQHH*radiusi
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = e0 + AAHH*d12 + BBHH*d6
       ep = ep + e0

       !LL
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
       e0 = QQLL*radiusi
       ep = ep + e0

       xi = site%x(ii+5)
       yi = site%y(ii+5)
       zi = site%z(ii+5)
       xj = site%x(jj+6)
       yj = site%y(jj+6)
       zj = site%z(jj+6)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQLL*radiusi
       ep = ep + e0

       xi = site%x(ii+6)
       yi = site%y(ii+6)
       zi = site%z(ii+6)
       xj = site%x(jj+5)
       yj = site%y(jj+5)
       zj = site%z(jj+5)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQLL*radiusi
       ep = ep + e0

       xi = site%x(ii+6)
       yi = site%y(ii+6)
       zi = site%z(ii+6)
       xj = site%x(jj+6)
       yj = site%y(jj+6)
       zj = site%z(jj+6)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQLL*radiusi
       ep = ep + e0

       !MH
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
       e0 = QQMH*radiusi
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
       e0 = QQMH*radiusi
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
       e0 = QQMH*radiusi
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
       e0 = QQMH*radiusi
       ep = ep + e0

       !ML
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
       e0 = QQML*radiusi
       ep = ep + e0

       xi = site%x(ii+6)
       yi = site%y(ii+6)
       zi = site%z(ii+6)
       xj = site%x(jj+4)
       yj = site%y(jj+4)
       zj = site%z(jj+4)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQML*radiusi
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
       e0 = QQML*radiusi
       ep = ep + e0

       xi = site%x(ii+4)
       yi = site%y(ii+4)
       zi = site%z(ii+4)
       xj = site%x(jj+6)
       yj = site%y(jj+6)
       zj = site%z(jj+6)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQML*radiusi
       ep = ep + e0

       !HL
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
       e0 = QQHL*radiusi
       ep = ep + e0

       xi = site%x(ii+6)
       yi = site%y(ii+6)
       zi = site%z(ii+6)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQHL*radiusi
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
       e0 = QQHL*radiusi
       ep = ep + e0

       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+6)
       yj = site%y(jj+6)
       zj = site%z(jj+6)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQHL*radiusi
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
       e0 = QQHL*radiusi
       ep = ep + e0

       xi = site%x(ii+6)
       yi = site%y(ii+6)
       zi = site%z(ii+6)
       xj = site%x(jj+3)
       yj = site%y(jj+3)
       zj = site%z(jj+3)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQHL*radiusi
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
       e0 = QQHL*radiusi
       ep = ep + e0

       xi = site%x(ii+3)
       yi = site%y(ii+3)
       zi = site%z(ii+3)
       xj = site%x(jj+6)
       yj = site%y(jj+6)
       zj = site%z(jj+6)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQHL*radiusi
       ep = ep + e0

       !OO
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
       e0 = AAOO*d12 + BBOO*d6
       ep = ep + e0

       !OH
       xi = site%x(ii+1)
       yi = site%y(ii+1)
       zi = site%z(ii+1)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = AAOH*d12 + BBOH*d6
       ep = ep + e0

       xi = site%x(ii+1)
       yi = site%y(ii+1)
       zi = site%z(ii+1)
       xj = site%x(jj+3)
       yj = site%y(jj+3)
       zj = site%z(jj+3)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = AAOH*d12 + BBOH*d6
       ep = ep + e0

       xi = site%x(ii+3)
       yi = site%y(ii+3)
       zi = site%z(ii+3)
       xj = site%x(jj+1)
       yj = site%y(jj+1)
       zj = site%z(jj+1)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = AAOH*d12 + BBOH*d6
       ep = ep + e0

       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+1)
       yj = site%y(jj+1)
       zj = site%z(jj+1)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = AAOH*d12 + BBOH*d6
       ep = ep + e0

       ep1(k) = ep*hh
    enddo
  end subroutine potential_nvde


end module nvde_module

! -*- f90 -*-

!1次元一時配列の確保と開放。これは若干効果がある。

!力の計算のsubroutine化。kf90を使えばちゃんとinline展開してくれるので
!速度の低下は防げる。
!平成１２年５月９日(火)ワークステーションのkapコンパイラを使わないと
!inline展開できないのだが、kapを使うと結果が大きくことなってくるようだ。
!当面保留する。
!vppでは当然ながら速度は低下しない。
!#define TEST12
#undef TEST12

!EwaldのコードはTEST12の方にしか挿入していない(手抜き)
#ifdef EWALD
#define TEST12
#endif

!2次元配列を利用した力の集計。
!小さいsystemでは効果があるかもしれない。もし効果があるならifで
!分岐させるべき。vppでは未check。
#undef TEST11

!2次元配列の要素の入れ替え。入れかえてもwsではほとんど速度は変化しない。
#undef TEST10
#ifdef TEST10
#define fxi(x,y) fxiq(y,x)
#define fyi(x,y) fyiq(y,x)
#define fzi(x,y) fziq(y,x)
#define fxj(x,y) fxjq(y,x)
#define fyj(x,y) fyjq(y,x)
#define fzj(x,y) fzjq(y,x)
#endif

module tip4p_module
  use physconst_module
  use interaction_module
  implicit none
contains
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine tip4p_setinteraction(si)
    use tip4p_constant_module
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, TIP4PSITE, LJ_COULOMB)
    !site order: OMHH&CoM
    !charges [Q]
    si%param(1,1) = 0d0
    si%param(1,2) = WM
    si%param(1,3) = WH
    si%param(1,4) = WH
    si%param(1,5) = 0d0
    !LJ-eps [kJ/mol]
    !write(STDERR,*) EPSWO
    si%param(2,1) = EPSWO
    si%param(2,2) = 0d0
    si%param(2,3) = 0d0
    si%param(2,4) = 0d0
    si%param(2,5) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = SIGWO
    si%param(3,2) = 0d0
    si%param(3,3) = 0d0
    si%param(3,4) = 0d0
    si%param(3,5) = 0d0
  end subroutine tip4p_setinteraction

  subroutine Rigid_TIP4P_Constructor(r)
    use rigid_module
    type(sRigid),intent(INOUT) :: r
    r%Ixx=(1.75618411509879202d+00)
    r%Iyy=(6.10236519209664374d-01)
    r%Izz=(1.14594759588912742d+00)
    r%ixxi = 1d0/r%ixx
    r%iyyi = 1d0/r%iyy
    r%izzi = 1d0/r%izz
    r%mass = 18d0
    r%massi= 1d0/r%mass
    r%molx(1)=0d0
    r%molx(2)=0d0
    r%molx(3)=0d0
    r%molx(4)=0d0
    r%molx(5)=0d0
    r%moly(1)=0.00000000000000000d+00
    r%moly(2)= 0.00000000000000000d+00
    r%moly(3)= 7.56950327263661182d-01
    r%moly(4)=-7.56950327263661182d-01
    r%molz(5)=0d0
    r%molz(1)=-6.50980307353661025d-02
    r%molz(2)= 8.49019692646338919d-02 
    r%molz(3)=5.20784245882928820d-01
    r%molz(4)=5.20784245882928820d-01
    r%molz(5)=0d0
    return
  end subroutine Rigid_TIP4P_Constructor

#ifdef EWALD
  subroutine interaction_force_tip4p(ewald,iv,mi,mj,site,ep1,vrmat)
#else
  subroutine interaction_force_tip4p(iv,mi,mj,site,ep1,vrmat)
#endif
    use site_module
    use mol_module
    use cutoff_module
    use tip4p_constant_module
#ifdef EWALD
    use ewald_module
    type(sEwald),intent(in) :: ewald
#endif
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
    real(kind=8) :: fxj1,fyj1,fzj1
    real(kind=8) :: fxj2,fyj2,fzj2
    real(kind=8) :: fxj3,fyj3,fzj3
    real(kind=8) :: fxj4,fyj4,fzj4
    real(kind=8) :: fxj5,fyj5,fzj5
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
#ifdef TEST10
    real(kind=8),dimension(0:MAXpair,TIP4PSITE) :: fxiq,fyiq,fziq&
         & ,fxjq,fyjq,fzjq
#else
    real(kind=8),dimension(:),allocatable :: fxi,fyi,fzi,fxj,fyj,fzj
#ifdef TEST11
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm1,fym1,fzm1
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm2,fym2,fzm2
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm3,fym3,fzm3
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm4,fym4,fzm4
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm5,fym5,fzm5
#endif
#endif /*VPOPTIMIZE*/
    real(kind=8) :: r_r0,r_r1,r_r0s,r_r1s,dr
    allocate(fxi(0:TIP4PSITE*(iv%npair+1)))
    allocate(fyi(0:TIP4PSITE*(iv%npair+1)))
    allocate(fzi(0:TIP4PSITE*(iv%npair+1)))
    allocate(fxj(0:TIP4PSITE*(iv%npair+1)))
    allocate(fyj(0:TIP4PSITE*(iv%npair+1)))
    allocate(fzj(0:TIP4PSITE*(iv%npair+1)))
    fxi(:)=0d0
    fyi(:)=0d0
    fzi(:)=0d0
    fxj(:)=0d0
    fyj(:)=0d0
    fzj(:)=0d0
#ifdef TEST11
#endif /*TEST11*/
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
       !write(STDERR,*) i,j,k,iv%Npair,mi%offset,mj%offset
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+TIP4PSITE)
       dy00 = site%y(ii+TIP4PSITE)
       dz00 = site%z(ii+TIP4PSITE)
       dx00 = dx00 - site%x(jj+TIP4PSITE)
       dy00 = dy00 - site%y(jj+TIP4PSITE)
       dz00 = dz00 - site%z(jj+TIP4PSITE)
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
       
       ep=0d0
#ifdef TEST12
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+2)-site%x(jj+2)+cx,site&
            & %y(ii+2)-site%y(jj+2)+cy,site%z(ii+2)-site%z(jj+2)+cz&
            & ,WWMM,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+2)-site%x(jj+2)+cx,site%y(ii+2)&
            & -site%y(jj+2)+cy,site%z(ii+2)-site%z(jj+2)+cz,WWMM,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi2 = fxi2 + forcex
       fyi2 = fyi2 + forcey
       fzi2 = fzi2 + forcez
       fxj2 = fxj2 + forcex
       fyj2 = fyj2 + forcey
       fzj2 = fzj2 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+2)-site%x(jj+3)+cx,site&
            & %y(ii+2)-site%y(jj+3)+cy,site%z(ii+2)-site%z(jj+3)+cz&
            & ,WWMH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+2)-site%x(jj+3)+cx,site%y(ii+2)&
            & -site%y(jj+3)+cy,site%z(ii+2)-site%z(jj+3)+cz,WWMH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi2 = fxi2 + forcex
       fyi2 = fyi2 + forcey
       fzi2 = fzi2 + forcez
       fxj3 = fxj3 + forcex
       fyj3 = fyj3 + forcey
       fzj3 = fzj3 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+2)-site%x(jj+4)+cx,site&
            & %y(ii+2)-site%y(jj+4)+cy,site%z(ii+2)-site%z(jj+4)+cz&
            & ,WWMH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+2)-site%x(jj+4)+cx,site%y(ii+2)&
            & -site%y(jj+4)+cy,site%z(ii+2)-site%z(jj+4)+cz,WWMH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi2 = fxi2 + forcex
       fyi2 = fyi2 + forcey
       fzi2 = fzi2 + forcez
       fxj4 = fxj4 + forcex
       fyj4 = fyj4 + forcey
       fzj4 = fzj4 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+3)-site%x(jj+2)+cx,site&
            & %y(ii+3)-site%y(jj+2)+cy,site%z(ii+3)-site%z(jj+2)+cz&
            & ,WWMH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+3)-site%x(jj+2)+cx,site%y(ii+3)&
            & -site%y(jj+2)+cy,site%z(ii+3)-site%z(jj+2)+cz,WWMH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi3 = fxi3 + forcex
       fyi3 = fyi3 + forcey
       fzi3 = fzi3 + forcez
       fxj2 = fxj2 + forcex
       fyj2 = fyj2 + forcey
       fzj2 = fzj2 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+3)-site%x(jj+3)+cx,site&
            & %y(ii+3)-site%y(jj+3)+cy,site%z(ii+3)-site%z(jj+3)+cz&
            & ,WWHH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+3)-site%x(jj+3)+cx,site%y(ii+3)&
            & -site%y(jj+3)+cy,site%z(ii+3)-site%z(jj+3)+cz,WWHH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi3 = fxi3 + forcex
       fyi3 = fyi3 + forcey
       fzi3 = fzi3 + forcez
       fxj3 = fxj3 + forcex
       fyj3 = fyj3 + forcey
       fzj3 = fzj3 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+3)-site%x(jj+4)+cx,site&
            & %y(ii+3)-site%y(jj+4)+cy,site%z(ii+3)-site%z(jj+4)+cz&
            & ,WWHH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+3)-site%x(jj+4)+cx,site%y(ii+3)&
            & -site%y(jj+4)+cy,site%z(ii+3)-site%z(jj+4)+cz,WWHH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi3 = fxi3 + forcex
       fyi3 = fyi3 + forcey
       fzi3 = fzi3 + forcez
       fxj4 = fxj4 + forcex
       fyj4 = fyj4 + forcey
       fzj4 = fzj4 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+4)-site%x(jj+2)+cx,site&
            & %y(ii+4)-site%y(jj+2)+cy,site%z(ii+4)-site%z(jj+2)+cz&
            & ,WWMH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+4)-site%x(jj+2)+cx,site%y(ii+4)&
            & -site%y(jj+2)+cy,site%z(ii+4)-site%z(jj+2)+cz,WWMH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi4 = fxi4 + forcex
       fyi4 = fyi4 + forcey
       fzi4 = fzi4 + forcez
       fxj2 = fxj2 + forcex
       fyj2 = fyj2 + forcey
       fzj2 = fzj2 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+4)-site%x(jj+3)+cx,site&
            & %y(ii+4)-site%y(jj+3)+cy,site%z(ii+4)-site%z(jj+3)+cz&
            & ,WWHH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+4)-site%x(jj+3)+cx,site%y(ii+4)&
            & -site%y(jj+3)+cy,site%z(ii+4)-site%z(jj+3)+cz,WWHH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi4 = fxi4 + forcex
       fyi4 = fyi4 + forcey
       fzi4 = fzi4 + forcez
       fxj3 = fxj3 + forcex
       fyj3 = fyj3 + forcey
       fzj3 = fzj3 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+4)-site%x(jj+4)+cx,site&
            & %y(ii+4)-site%y(jj+4)+cy,site%z(ii+4)-site%z(jj+4)+cz&
            & ,WWHH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+4)-site%x(jj+4)+cx,site%y(ii+4)&
            & -site%y(jj+4)+cy,site%z(ii+4)-site%z(jj+4)+cz,WWHH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi4 = fxi4 + forcex
       fyi4 = fyi4 + forcey
       fzi4 = fzi4 + forcez
       fxj4 = fxj4 + forcex
       fyj4 = fyj4 + forcey
       fzj4 = fzj4 + forcez
       call LJ_Force(site%x(ii+1)-site%x(jj+1)+cx,site%y(ii+1)-site&
            & %y(jj+1)+cy,site%z(ii+1)-site%z(jj+1)+cz,AA,BB,e0&
            & ,forcex,forcey,forcez)
       ep = ep + e0
       fxi1 = fxi1 + forcex
       fyi1 = fyi1 + forcey
       fzi1 = fzi1 + forcez
       fxj1 = fxj1 + forcex
       fyj1 = fyj1 + forcey
       fzj1 = fzj1 + forcez
#else /*TEST12*/
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
       e0 = WWMM*radiusi
       f0 = WWMM*radiusi*dd
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
       e0 = WWMH*radiusi
       f0 = WWMH*radiusi*dd
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
       e0 = WWMH*radiusi
       f0 = WWMH*radiusi*dd
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
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = WWMH*radiusi
       f0 = WWMH*radiusi*dd
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
       e0 = WWMH*radiusi
       f0 = WWMH*radiusi*dd
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
       e0 = WWHH*radiusi
       f0 = WWHH*radiusi*dd
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
       e0 = WWHH*radiusi
       f0 = WWHH*radiusi*dd
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
       xj = site%x(jj+3)
       yj = site%y(jj+3)
       zj = site%z(jj+3)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = WWHH*radiusi
       f0 = WWHH*radiusi*dd
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
       e0 = WWHH*radiusi
       f0 = WWHH*radiusi*dd
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
#endif /*TEST12*/

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
       dhx=dhx*ep
       dhy=dhy*ep
       dhz=dhz*ep
       fxi5=fxi5+dhx
       fyi5=fyi5+dhy
       fzi5=fzi5+dhz
       fxj5=fxj5+dhx
       fyj5=fyj5+dhy
       fzj5=fzj5+dhz
       ep1(k) = ep*hh
       !vrx = dx00*(fxi1+fxi2+fxi3+fxi4+fxi5)
       !vry = dy00*(fyi1+fyi2+fyi3+fyi4+fyi5)
       !vrz = dz00*(fzi1+fzi2+fzi3+fzi4+fzi5)
       !vrsum=vrsum+vrx+vry+vrz
       vrx = (fxi1+fxi2+fxi3+fxi4+fxi5)
       vry = (fyi1+fyi2+fyi3+fyi4+fyi5)
       vrz = (fzi1+fzi2+fzi3+fzi4+fzi5)
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
#ifdef SHOWFORCE
       !見苦しいが、対力はここで表示させざるをえない。
       !TIP4P以外は未実装
       !対力は反発の時正、引力のとき負。
       fxs=fxi1+fxi2+fxi3+fxi4+fxi5
       fys=fyi1+fyi2+fyi3+fyi4+fyi5
       fzs=fzi1+fzi2+fzi3+fzi4+fzi5
       fpair=dsqrt(fxs**2+fys**2+fzs**2)
       if(dx00*fxs+dy00*fys+dz00*fzs.lt.0) fpair=-fpair
       !出力は原則として配列を0から数える。
       !write(STDOUT,*) i-1,j-1,dsqrt(dx00**2+dy00**2+dz00**2),fpair
#endif
       !write(6,*) i,j,k,ep
       !
       ! 非vector計算機の場合には、この場で力の集計をした方がてっとり
       !     早い。
#ifdef VPOPTIMIZE
       !     ......計算結果をlist vectorに格納する。
       fxi(k*TIP4PSITE+0)=fxi1
       fxi(k*TIP4PSITE+1)=fxi2
       fxi(k*TIP4PSITE+2)=fxi3
       fxi(k*TIP4PSITE+3)=fxi4
       fxi(k*TIP4PSITE+4)=fxi5
       fyi(k*TIP4PSITE+0)=fyi1
       fyi(k*TIP4PSITE+1)=fyi2
       fyi(k*TIP4PSITE+2)=fyi3
       fyi(k*TIP4PSITE+3)=fyi4
       fyi(k*TIP4PSITE+4)=fyi5
       fzi(k*TIP4PSITE+0)=fzi1
       fzi(k*TIP4PSITE+1)=fzi2
       fzi(k*TIP4PSITE+2)=fzi3
       fzi(k*TIP4PSITE+3)=fzi4
       fzi(k*TIP4PSITE+4)=fzi5
       
       fxj(k*TIP4PSITE+0)=-fxj1
       fxj(k*TIP4PSITE+1)=-fxj2
       fxj(k*TIP4PSITE+2)=-fxj3
       fxj(k*TIP4PSITE+3)=-fxj4
       fxj(k*TIP4PSITE+4)=-fxj5
       fyj(k*TIP4PSITE+0)=-fyj1
       fyj(k*TIP4PSITE+1)=-fyj2
       fyj(k*TIP4PSITE+2)=-fyj3
       fyj(k*TIP4PSITE+3)=-fyj4
       fyj(k*TIP4PSITE+4)=-fyj5
       fzj(k*TIP4PSITE+0)=-fzj1
       fzj(k*TIP4PSITE+1)=-fzj2
       fzj(k*TIP4PSITE+2)=-fzj3
       fzj(k*TIP4PSITE+3)=-fzj4
       fzj(k*TIP4PSITE+4)=-fzj5
#else 
       site%fx(ii+1)=site%fx(ii+1)+fxi1
       site%fx(ii+2)=site%fx(ii+2)+fxi2
       site%fx(ii+3)=site%fx(ii+3)+fxi3
       site%fx(ii+4)=site%fx(ii+4)+fxi4
       site%fx(ii+5)=site%fx(ii+5)+fxi5
       site%fx(jj+1)=site%fx(jj+1)-fxj1
       site%fx(jj+2)=site%fx(jj+2)-fxj2
       site%fx(jj+3)=site%fx(jj+3)-fxj3
       site%fx(jj+4)=site%fx(jj+4)-fxj4
       site%fx(jj+5)=site%fx(jj+5)-fxj5
       
       site%fy(ii+1)=site%fy(ii+1)+fyi1
       site%fy(ii+2)=site%fy(ii+2)+fyi2
       site%fy(ii+3)=site%fy(ii+3)+fyi3
       site%fy(ii+4)=site%fy(ii+4)+fyi4
       site%fy(ii+5)=site%fy(ii+5)+fyi5
       site%fy(jj+1)=site%fy(jj+1)-fyj1
       site%fy(jj+2)=site%fy(jj+2)-fyj2
       site%fy(jj+3)=site%fy(jj+3)-fyj3
       site%fy(jj+4)=site%fy(jj+4)-fyj4
       site%fy(jj+5)=site%fy(jj+5)-fyj5

       site%fz(ii+1)=site%fz(ii+1)+fzi1
       site%fz(ii+2)=site%fz(ii+2)+fzi2
       site%fz(ii+3)=site%fz(ii+3)+fzi3
       site%fz(ii+4)=site%fz(ii+4)+fzi4
       site%fz(ii+5)=site%fz(ii+5)+fzi5
       site%fz(jj+1)=site%fz(jj+1)-fzj1
       site%fz(jj+2)=site%fz(jj+2)-fzj2
       site%fz(jj+3)=site%fz(jj+3)-fzj3
       site%fz(jj+4)=site%fz(jj+4)-fzj4
       site%fz(jj+5)=site%fz(jj+5)-fzj5
       !f90  -O6 2.724u 0.020s 0:03.05 89.8% 0+11k 6+13io 0pf+0w
       !kf90 -O6 2.751u 0.014s 0:03.12 88.4% 0+11k 7+11io 0pf+0w
       !この部分を除去すると
       !kf90 -O6 1.427u 0.021s 0:01.77 81.3% 0+11k 7+11io 0pf+0w
#endif /*VPOPTIMIZE*/
       !endif
    enddo
#ifdef VPOPTIMIZE
#ifdef TEST11
    fxm1(:,:)=0d0
    fym1(:,:)=0d0
    fzm1(:,:)=0d0
    fxm2(:,:)=0d0
    fym2(:,:)=0d0
    fzm2(:,:)=0d0
    fxm3(:,:)=0d0
    fym3(:,:)=0d0
    fzm3(:,:)=0d0
    fxm4(:,:)=0d0
    fym4(:,:)=0d0
    fzm4(:,:)=0d0
    fxm5(:,:)=0d0
    fym5(:,:)=0d0
    fzm5(:,:)=0d0
!OCL NOVREC
    do k=1,iv%npair
       i=iv%pair_i(k)
       j=iv%pair_j(k)
       fxm1(j,i)=fxj(k*5)
       fym1(j,i)=fyj(k*TIP4PSITE)
       fzm1(j,i)=fzj(k*TIP4PSITE)
       
       fxm2(j,i)=fxj(k*TIP4PSITE+1)
       fym2(j,i)=fyj(k*TIP4PSITE+1)
       fzm2(j,i)=fzj(k*TIP4PSITE+1)
       
       fxm3(j,i)=fxj(k*TIP4PSITE+2)
       fym3(j,i)=fyj(k*TIP4PSITE+2)
       fzm3(j,i)=fzj(k*TIP4PSITE+2)
       
       fxm4(j,i)=fxj(k*TIP4PSITE+3)
       fym4(j,i)=fyj(k*TIP4PSITE+3)
       fzm4(j,i)=fzj(k*TIP4PSITE+3)
       
       fxm5(j,i)=fxj(k*TIP4PSITE+4)
       fym5(j,i)=fyj(k*TIP4PSITE+4)
       fzm5(j,i)=fzj(k*TIP4PSITE+4)
    enddo
  !異種グループ間相互作用の場合は、一旦iに加わる力を集計する。
    if(.not.iv%isomol)then
       do i=1,mi%nmol
          !OCL NOVREC
          do j=1,mj%nmol
             site%fx((j-1)*TIP4PSITE+mj%offset+1)=site%fx((j-1)&
                  & *TIP4PSITE+mj%offset+1)+fxm1(j,i)
             site%fx((j-1)*TIP4PSITE+mj%offset+2)=site%fx((j-1)&
                  & *TIP4PSITE+mj%offset+2)+fxm2(j,i)
             site%fx((j-1)*TIP4PSITE+mj%offset+3)=site%fx((j-1)&
                  & *TIP4PSITE+mj%offset+3)+fxm3(j,i)
             site%fx((j-1)*TIP4PSITE+mj%offset+4)=site%fx((j-1)&
                  & *TIP4PSITE+mj%offset+4)+fxm4(j,i)
             site%fx((j-1)*TIP4PSITE+mj%offset+5)=site%fx((j-1)&
                  & *TIP4PSITE+mj%offset+5)+fxm5(j,i)
             site%fy((j-1)*TIP4PSITE+mj%offset+1)=site%fy((j-1)&
                  & *TIP4PSITE+mj%offset+1)+fym1(j,i)
             site%fy((j-1)*TIP4PSITE+mj%offset+2)=site%fy((j-1)&
                  & *TIP4PSITE+mj%offset+2)+fym2(j,i)
             site%fy((j-1)*TIP4PSITE+mj%offset+3)=site%fy((j-1)&
                  & *TIP4PSITE+mj%offset+3)+fym3(j,i)
             site%fy((j-1)*TIP4PSITE+mj%offset+4)=site%fy((j-1)&
                  & *TIP4PSITE+mj%offset+4)+fym4(j,i)
             site%fy((j-1)*TIP4PSITE+mj%offset+5)=site%fy((j-1)&
                  & *TIP4PSITE+mj%offset+5)+fym5(j,i)
             site%fz((j-1)*TIP4PSITE+mj%offset+1)=site%fz((j-1)&
                  & *TIP4PSITE+mj%offset+1)+fzm1(j,i)
             site%fz((j-1)*TIP4PSITE+mj%offset+2)=site%fz((j-1)&
                  & *TIP4PSITE+mj%offset+2)+fzm2(j,i)
             site%fz((j-1)*TIP4PSITE+mj%offset+3)=site%fz((j-1)&
                  & *TIP4PSITE+mj%offset+3)+fzm3(j,i)
             site%fz((j-1)*TIP4PSITE+mj%offset+4)=site%fz((j-1)&
                  & *TIP4PSITE+mj%offset+4)+fzm4(j,i)
             site%fz((j-1)*TIP4PSITE+mj%offset+5)=site%fz((j-1)&
                  & *TIP4PSITE+mj%offset+5)+fzm5(j,i)
          enddo
       enddo
       fxm1(:,:)=0d0
       fym1(:,:)=0d0
       fzm1(:,:)=0d0
       fxm2(:,:)=0d0
       fym2(:,:)=0d0
       fzm2(:,:)=0d0
       fxm3(:,:)=0d0
       fym3(:,:)=0d0
       fzm3(:,:)=0d0
       fxm4(:,:)=0d0
       fym4(:,:)=0d0
       fzm4(:,:)=0d0
       fxm5(:,:)=0d0
       fym5(:,:)=0d0
       fzm5(:,:)=0d0
    endif
  !OCL NOVREC
    do k=1,iv%npair
       i=iv%pair_i(k)
       j=iv%pair_j(k)
       fxm1(i,j)=fxi(k*TIP4PSITE)
       fym1(i,j)=fyi(k*TIP4PSITE)
       fzm1(i,j)=fzi(k*TIP4PSITE)
       
       fxm2(i,j)=fxi(k*TIP4PSITE+1)
       fym2(i,j)=fyi(k*TIP4PSITE+1)
       fzm2(i,j)=fzi(k*TIP4PSITE+1)
       
       fxm3(i,j)=fxi(k*TIP4PSITE+2)
       fym3(i,j)=fyi(k*TIP4PSITE+2)
       fzm3(i,j)=fzi(k*TIP4PSITE+2)
       
       fxm4(i,j)=fxi(k*TIP4PSITE+3)
       fym4(i,j)=fyi(k*TIP4PSITE+3)
       fzm4(i,j)=fzi(k*TIP4PSITE+3)
       
       fxm5(i,j)=fxi(k*TIP4PSITE+4)
       fym5(i,j)=fyi(k*TIP4PSITE+4)
       fzm5(i,j)=fzi(k*TIP4PSITE+4)
    enddo
    do j=1,mj%nmol
       !OCL NOVREC
       do i=1,mi%nmol
          site%fx((i-1)*TIP4PSITE+mi%offset+1)=site%fx((i-1)&
               & *TIP4PSITE+mi%offset+1)+fxm1(i,j)
          site%fx((i-1)*TIP4PSITE+mi%offset+2)=site%fx((i-1)&
               & *TIP4PSITE+mi%offset+2)+fxm2(i,j)
          site%fx((i-1)*TIP4PSITE+mi%offset+3)=site%fx((i-1)&
               & *TIP4PSITE+mi%offset+3)+fxm3(i,j)
          site%fx((i-1)*TIP4PSITE+mi%offset+4)=site%fx((i-1)&
               & *TIP4PSITE+mi%offset+4)+fxm4(i,j)
          site%fx((i-1)*TIP4PSITE+mi%offset+5)=site%fx((i-1)&
               & *TIP4PSITE+mi%offset+5)+fxm5(i,j)
          site%fy((i-1)*TIP4PSITE+mi%offset+1)=site%fy((i-1)&
               & *TIP4PSITE+mi%offset+1)+fym1(i,j)
          site%fy((i-1)*TIP4PSITE+mi%offset+2)=site%fy((i-1)&
               & *TIP4PSITE+mi%offset+2)+fym2(i,j)
          site%fy((i-1)*TIP4PSITE+mi%offset+3)=site%fy((i-1)&
               & *TIP4PSITE+mi%offset+3)+fym3(i,j)
          site%fy((i-1)*TIP4PSITE+mi%offset+4)=site%fy((i-1)&
               & *TIP4PSITE+mi%offset+4)+fym4(i,j)
          site%fy((i-1)*TIP4PSITE+mi%offset+5)=site%fy((i-1)&
               & *TIP4PSITE+mi%offset+5)+fym5(i,j)
          site%fz((i-1)*TIP4PSITE+mi%offset+1)=site%fz((i-1)&
               & *TIP4PSITE+mi%offset+1)+fzm1(i,j)
          site%fz((i-1)*TIP4PSITE+mi%offset+2)=site%fz((i-1)&
               & *TIP4PSITE+mi%offset+2)+fzm2(i,j)
          site%fz((i-1)*TIP4PSITE+mi%offset+3)=site%fz((i-1)&
               & *TIP4PSITE+mi%offset+3)+fzm3(i,j)
          site%fz((i-1)*TIP4PSITE+mi%offset+4)=site%fz((i-1)&
               & *TIP4PSITE+mi%offset+4)+fzm4(i,j)
          site%fz((i-1)*TIP4PSITE+mi%offset+5)=site%fz((i-1)&
               & *TIP4PSITE+mi%offset+5)+fzm5(i,j)
       enddo
    enddo
    deallocate(fxi)
    deallocate(fyi)
    deallocate(fzi)
    deallocate(fxj)
    deallocate(fyj)
    deallocate(fzj)
#else /*TEST11*/
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
          site%fx(ii+mi%offset)=site%fx(ii+mi%offset)+fxi(k*TIP4PSITE&
               & +s-1)
          site%fy(ii+mi%offset)=site%fy(ii+mi%offset)+fyi(k*TIP4PSITE&
               & +s-1)
          site%fz(ii+mi%offset)=site%fz(ii+mi%offset)+fzi(k*TIP4PSITE&
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
          site%fx(jj+mj%offset)=site%fx(jj+mj%offset)+fxj(k*TIP4PSITE&
               & +s-1)
          site%fy(jj+mj%offset)=site%fy(jj+mj%offset)+fyj(k*TIP4PSITE&
               & +s-1)
          site%fz(jj+mj%offset)=site%fz(jj+mj%offset)+fzj(k*TIP4PSITE&
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
#endif /*TEST11*/
#endif /*VPOPTIMIZE*/
    !write(6,*) kk
#ifdef VERBOSE
    write(STDERR,*) "INTERACTION PAIRS: ",count
#endif
  end subroutine interaction_force_tip4p

  subroutine interaction_yaphb_tip4p(iv,mi,mj,site,r,file)
    use site_module
    use mol_module
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    real(kind=8),intent(in) :: r
    integer,intent(in) :: file
    real(kind=8) :: dx,dy,dz,dx00,dy00,dz00
    real(kind=8) :: rr2
    !real(kind=8) :: dhx,dhy,dhz
    real(kind=8) :: xi,yi,zi,xj,yj,zj
    real(kind=8) ox,oy,oz
    real(kind=8) cx,cy,cz
    !real(kind=8) ep
    !real(kind=8) dd !,hh,dh
    integer :: i,j,k
    integer :: ii,jj
    !integer :: kk
    !kk=0
    !VPP用のベクトル化制御行
!OCL VECTOR,REPEAT(1000000)
    rr2=r*r
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+5)
       dy00 = site%y(ii+5)
       dz00 = site%z(ii+5)
       dx00 = dx00 - site%x(jj+5)
       dy00 = dy00 - site%y(jj+5)
       dz00 = dz00 - site%z(jj+5)
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       dx00 = dx00-ox
       dy00 = dy00-oy
       dz00 = dz00-oz
       cx = -ox
       cy = -oy
       cz = -oz
       
       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+3)
       yj = site%y(jj+3)
       zj = site%z(jj+3)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       if(dx**2 + dy**2 + dz**2.lt.rr2)then
          write(file,3) xi+cx,yi+cy,zi+cz,xj,yj,zj
       endif
       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+4)
       yj = site%y(jj+4)
       zj = site%z(jj+4)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       if(dx**2 + dy**2 + dz**2.lt.rr2)then
          write(file,3) xi+cx,yi+cy,zi+cz,xj,yj,zj
       endif
       xi = site%x(ii+3)
       yi = site%y(ii+3)
       zi = site%z(ii+3)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       if(dx**2 + dy**2 + dz**2.lt.rr2)then
          write(file,3) xi+cx,yi+cy,zi+cz,xj,yj,zj
       endif
       xi = site%x(ii+4)
       yi = site%y(ii+4)
       zi = site%z(ii+4)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       if(dx**2 + dy**2 + dz**2.lt.rr2)then
          write(file,3) xi+cx,yi+cy,zi+cz,xj,yj,zj
       endif
    enddo
3   format("l",6(1x,f7.2))
  end subroutine interaction_yaphb_tip4p

  subroutine interaction_ngph_tip4p(iv,mi,mj,site,r,file)
    use mol_module
    use site_module
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    real(kind=8),intent(in) :: r
    integer,intent(in) :: file
    real(kind=8) :: dx,dy,dz,dx00,dy00,dz00
    real(kind=8) :: rr2
    !real(kind=8) :: dhx,dhy,dhz
    real(kind=8) :: xi,yi,zi,xj,yj,zj
    real(kind=8) ox,oy,oz
    real(kind=8) cx,cy,cz
    !real(kind=8) ep
    !real(kind=8) dd !,hh,dh
    integer :: i,j,k
    integer :: ii,jj
    !integer :: kk
    !kk=0
!VPP用のベクトル化制御行
!OCL VECTOR,REPEAT(1000000)
    rr2=r*r
    write(file,2)
2   format("@NGPH")
    write(file,*) mi%nmol
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+5)
       dy00 = site%y(ii+5)
       dz00 = site%z(ii+5)
       dx00 = dx00 - site%x(jj+5)
       dy00 = dy00 - site%y(jj+5)
       dz00 = dz00 - site%z(jj+5)
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       dx00 = dx00-ox
       dy00 = dy00-oy
       dz00 = dz00-oz
       cx = -ox
       cy = -oy
       cz = -oz
       
       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+3)
       yj = site%y(jj+3)
       zj = site%z(jj+3)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       if(dx**2 + dy**2 + dz**2.lt.rr2)then
          write(file,3) j-1,i-1
       endif
       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+4)
       yj = site%y(jj+4)
       zj = site%z(jj+4)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       if(dx**2 + dy**2 + dz**2.lt.rr2)then
          write(file,3) j-1,i-1
       endif
       xi = site%x(ii+3)
       yi = site%y(ii+3)
       zi = site%z(ii+3)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       if(dx**2 + dy**2 + dz**2.lt.rr2)then
          write(file,3) i-1,j-1
       endif
       xi = site%x(ii+4)
       yi = site%y(ii+4)
       zi = site%z(ii+4)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       if(dx**2 + dy**2 + dz**2.lt.rr2)then
          write(file,3) i-1,j-1
       endif
    enddo
    write(file,3) -1,-1
3   format(2(i4,1x))
  end subroutine interaction_ngph_tip4p

!平成13年8月3日(金)WGPH形式はNGPH形式に、腕番号を付加したもの
  subroutine interaction_wgph_tip4p(iv,mi,mj,site,r,file)
    use site_module
    use mol_module
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    real(kind=8),intent(in) :: r
    integer,intent(in) :: file
    real(kind=8) :: dx,dy,dz,dx00,dy00,dz00
    real(kind=8) :: rr2
    !real(kind=8) :: dhx,dhy,dhz
    real(kind=8) :: xi,yi,zi,xj,yj,zj
    real(kind=8) ox,oy,oz
    real(kind=8) cx,cy,cz
    !real(kind=8) ep
    !real(kind=8) dd !,hh,dh
    integer :: i,j,k
    integer :: ii,jj
    !integer :: kk
    !kk=0
!VPP用のベクトル化制御行
!OCL VECTOR,REPEAT(1000000)
    rr2=r*r
    write(file,2)
2   format("@WGPH")
    write(file,*) mi%nmol
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+5)
       dy00 = site%y(ii+5)
       dz00 = site%z(ii+5)
       dx00 = dx00 - site%x(jj+5)
       dy00 = dy00 - site%y(jj+5)
       dz00 = dz00 - site%z(jj+5)
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       dx00 = dx00-ox
       dy00 = dy00-oy
       dz00 = dz00-oz
       cx = -ox
       cy = -oy
       cz = -oz
       
       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+3)
       yj = site%y(jj+3)
       zj = site%z(jj+3)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       if(dx**2 + dy**2 + dz**2.lt.rr2)then
          write(file,3) j-1,i-1,0
       endif
       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+4)
       yj = site%y(jj+4)
       zj = site%z(jj+4)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       if(dx**2 + dy**2 + dz**2.lt.rr2)then
          write(file,3) j-1,i-1,1
       endif
       xi = site%x(ii+3)
       yi = site%y(ii+3)
       zi = site%z(ii+3)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       if(dx**2 + dy**2 + dz**2.lt.rr2)then
          write(file,3) i-1,j-1,0
       endif
       xi = site%x(ii+4)
       yi = site%y(ii+4)
       zi = site%z(ii+4)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       if(dx**2 + dy**2 + dz**2.lt.rr2)then
          write(file,3) i-1,j-1,1
       endif
    enddo
    write(file,3) -1,-1
3   format(3(i4,1x))
  end subroutine interaction_wgph_tip4p

  !
  !MonteCarlo用のエネルギー計算算程。
  !
  subroutine energy_tip4p(iv,mi,mj,site,epsum,vrmat)
    use tip4p_constant_module
    use site_module
    use mol_module
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    real(kind=8),intent(out) :: epsum,vrmat(3,3)
    real(kind=8) :: fxi1,fyi1,fzi1
    real(kind=8) :: fxi2,fyi2,fzi2
    real(kind=8) :: fxi3,fyi3,fzi3
    real(kind=8) :: fxi4,fyi4,fzi4
    real(kind=8) :: fxi5,fyi5,fzi5
    real(kind=8) :: fxj1,fyj1,fzj1
    real(kind=8) :: fxj2,fyj2,fzj2
    real(kind=8) :: fxj3,fyj3,fzj3
    real(kind=8) :: fxj4,fyj4,fzj4
    real(kind=8) :: fxj5,fyj5,fzj5
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
    !0番目の要素は、loop内のifをへらすためのトリック
#ifdef TEST10
    real(kind=8),dimension(0:MAXpair,TIP4PSITE) :: fxiq,fyiq,fziq&
         & ,fxjq,fyjq,fzjq
#else
    real(kind=8),dimension(:),allocatable :: fxi,fyi,fzi,fxj,fyj,fzj
#ifdef TEST11
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm1,fym1,fzm1
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm2,fym2,fzm2
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm3,fym3,fzm3
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm4,fym4,fzm4
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm5,fym5,fzm5
#endif
#endif /*VPOPTIMIZE*/
    real(kind=8) :: r_r0,r_r1,r_r0s,r_r1s,dr
    allocate(fxi(0:TIP4PSITE*(iv%npair+1)))
    allocate(fyi(0:TIP4PSITE*(iv%npair+1)))
    allocate(fzi(0:TIP4PSITE*(iv%npair+1)))
    allocate(fxj(0:TIP4PSITE*(iv%npair+1)))
    allocate(fyj(0:TIP4PSITE*(iv%npair+1)))
    allocate(fzj(0:TIP4PSITE*(iv%npair+1)))
    fxi(:)=0d0
    fyi(:)=0d0
    fzi(:)=0d0
    fxj(:)=0d0
    fyj(:)=0d0
    fzj(:)=0d0
#ifdef TEST11
#endif /*TEST11*/
#endif /*VPOPTIMIZE*/
#ifdef VERBOSE
    count=0
#endif
    epsum=0d0
    !kk=0
    !VPP用のベクトル化制御行
    !OCL VECTOR,REPEAT(1000000)
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       !write(STDERR,*) i,j,k,iv%Npair,mi%offset,mj%offset
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+TIP4PSITE)
       dy00 = site%y(ii+TIP4PSITE)
       dz00 = site%z(ii+TIP4PSITE)
       dx00 = dx00 - site%x(jj+TIP4PSITE)
       dy00 = dy00 - site%y(jj+TIP4PSITE)
       dz00 = dz00 - site%z(jj+TIP4PSITE)
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
       
       ep=0d0
#ifdef TEST12
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+2)-site%x(jj+2)+cx,site&
            & %y(ii+2)-site%y(jj+2)+cy,site%z(ii+2)-site%z(jj+2)+cz&
            & ,WWMM,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+2)-site%x(jj+2)+cx,site%y(ii+2)&
            & -site%y(jj+2)+cy,site%z(ii+2)-site%z(jj+2)+cz,WWMM,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi2 = fxi2 + forcex
       fyi2 = fyi2 + forcey
       fzi2 = fzi2 + forcez
       fxj2 = fxj2 + forcex
       fyj2 = fyj2 + forcey
       fzj2 = fzj2 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+2)-site%x(jj+3)+cx,site&
            & %y(ii+2)-site%y(jj+3)+cy,site%z(ii+2)-site%z(jj+3)+cz&
            & ,WWMH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+2)-site%x(jj+3)+cx,site%y(ii+2)&
            & -site%y(jj+3)+cy,site%z(ii+2)-site%z(jj+3)+cz,WWMH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi2 = fxi2 + forcex
       fyi2 = fyi2 + forcey
       fzi2 = fzi2 + forcez
       fxj3 = fxj3 + forcex
       fyj3 = fyj3 + forcey
       fzj3 = fzj3 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+2)-site%x(jj+4)+cx,site&
            & %y(ii+2)-site%y(jj+4)+cy,site%z(ii+2)-site%z(jj+4)+cz&
            & ,WWMH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+2)-site%x(jj+4)+cx,site%y(ii+2)&
            & -site%y(jj+4)+cy,site%z(ii+2)-site%z(jj+4)+cz,WWMH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi2 = fxi2 + forcex
       fyi2 = fyi2 + forcey
       fzi2 = fzi2 + forcez
       fxj4 = fxj4 + forcex
       fyj4 = fyj4 + forcey
       fzj4 = fzj4 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+3)-site%x(jj+2)+cx,site&
            & %y(ii+3)-site%y(jj+2)+cy,site%z(ii+3)-site%z(jj+2)+cz&
            & ,WWMH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+3)-site%x(jj+2)+cx,site%y(ii+3)&
            & -site%y(jj+2)+cy,site%z(ii+3)-site%z(jj+2)+cz,WWMH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi3 = fxi3 + forcex
       fyi3 = fyi3 + forcey
       fzi3 = fzi3 + forcez
       fxj2 = fxj2 + forcex
       fyj2 = fyj2 + forcey
       fzj2 = fzj2 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+3)-site%x(jj+3)+cx,site&
            & %y(ii+3)-site%y(jj+3)+cy,site%z(ii+3)-site%z(jj+3)+cz&
            & ,WWHH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+3)-site%x(jj+3)+cx,site%y(ii+3)&
            & -site%y(jj+3)+cy,site%z(ii+3)-site%z(jj+3)+cz,WWHH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi3 = fxi3 + forcex
       fyi3 = fyi3 + forcey
       fzi3 = fzi3 + forcez
       fxj3 = fxj3 + forcex
       fyj3 = fyj3 + forcey
       fzj3 = fzj3 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+3)-site%x(jj+4)+cx,site&
            & %y(ii+3)-site%y(jj+4)+cy,site%z(ii+3)-site%z(jj+4)+cz&
            & ,WWHH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+3)-site%x(jj+4)+cx,site%y(ii+3)&
            & -site%y(jj+4)+cy,site%z(ii+3)-site%z(jj+4)+cz,WWHH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi3 = fxi3 + forcex
       fyi3 = fyi3 + forcey
       fzi3 = fzi3 + forcez
       fxj4 = fxj4 + forcex
       fyj4 = fyj4 + forcey
       fzj4 = fzj4 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+4)-site%x(jj+2)+cx,site&
            & %y(ii+4)-site%y(jj+2)+cy,site%z(ii+4)-site%z(jj+2)+cz&
            & ,WWMH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+4)-site%x(jj+2)+cx,site%y(ii+4)&
            & -site%y(jj+2)+cy,site%z(ii+4)-site%z(jj+2)+cz,WWMH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi4 = fxi4 + forcex
       fyi4 = fyi4 + forcey
       fzi4 = fzi4 + forcez
       fxj2 = fxj2 + forcex
       fyj2 = fyj2 + forcey
       fzj2 = fzj2 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+4)-site%x(jj+3)+cx,site&
            & %y(ii+4)-site%y(jj+3)+cy,site%z(ii+4)-site%z(jj+3)+cz&
            & ,WWHH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+4)-site%x(jj+3)+cx,site%y(ii+4)&
            & -site%y(jj+3)+cy,site%z(ii+4)-site%z(jj+3)+cz,WWHH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi4 = fxi4 + forcex
       fyi4 = fyi4 + forcey
       fzi4 = fzi4 + forcez
       fxj3 = fxj3 + forcex
       fyj3 = fyj3 + forcey
       fzj3 = fzj3 + forcez
#ifdef EWALD
       call Coulomb_Force(ewald,site%x(ii+4)-site%x(jj+4)+cx,site&
            & %y(ii+4)-site%y(jj+4)+cy,site%z(ii+4)-site%z(jj+4)+cz&
            & ,WWHH,e0,forcex,forcey,forcez)
#else
       call Coulomb_Force(site%x(ii+4)-site%x(jj+4)+cx,site%y(ii+4)&
            & -site%y(jj+4)+cy,site%z(ii+4)-site%z(jj+4)+cz,WWHH,e0&
            & ,forcex,forcey,forcez)
#endif
       ep=ep+e0
       fxi4 = fxi4 + forcex
       fyi4 = fyi4 + forcey
       fzi4 = fzi4 + forcez
       fxj4 = fxj4 + forcex
       fyj4 = fyj4 + forcey
       fzj4 = fzj4 + forcez
       call LJ_Force(site%x(ii+1)-site%x(jj+1)+cx,site%y(ii+1)-site&
            & %y(jj+1)+cy,site%z(ii+1)-site%z(jj+1)+cz,AA,BB,e0&
            & ,forcex,forcey,forcez)
       ep = ep + e0
       fxi1 = fxi1 + forcex
       fyi1 = fyi1 + forcey
       fzi1 = fzi1 + forcez
       fxj1 = fxj1 + forcex
       fyj1 = fyj1 + forcey
       fzj1 = fzj1 + forcez
#else /*TEST12*/
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
       e0 = WWMM*radiusi
       f0 = WWMM*radiusi*dd
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
       e0 = WWMH*radiusi
       f0 = WWMH*radiusi*dd
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
       e0 = WWMH*radiusi
       f0 = WWMH*radiusi*dd
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
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = WWMH*radiusi
       f0 = WWMH*radiusi*dd
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
       e0 = WWMH*radiusi
       f0 = WWMH*radiusi*dd
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
       e0 = WWHH*radiusi
       f0 = WWHH*radiusi*dd
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
       e0 = WWHH*radiusi
       f0 = WWHH*radiusi*dd
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
       xj = site%x(jj+3)
       yj = site%y(jj+3)
       zj = site%z(jj+3)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = WWHH*radiusi
       f0 = WWHH*radiusi*dd
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
       e0 = WWHH*radiusi
       f0 = WWHH*radiusi*dd
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
#endif /*TEST12*/

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
       dhx=dhx*ep
       dhy=dhy*ep
       dhz=dhz*ep
       fxi5=fxi5+dhx
       fyi5=fyi5+dhy
       fzi5=fzi5+dhz
       fxj5=fxj5+dhx
       fyj5=fyj5+dhy
       fzj5=fzj5+dhz
       epsum=epsum+ep*hh
       !vrx = dx00*(fxi1+fxi2+fxi3+fxi4+fxi5)
       !vry = dy00*(fyi1+fyi2+fyi3+fyi4+fyi5)
       !vrz = dz00*(fzi1+fzi2+fzi3+fzi4+fzi5)
       !vrsum=vrsum+vrx+vry+vrz
       vrx = (fxi1+fxi2+fxi3+fxi4+fxi5)
       vry = (fyi1+fyi2+fyi3+fyi4+fyi5)
       vrz = (fzi1+fzi2+fzi3+fzi4+fzi5)
       vrmat(1,1) = dx00 * vrx
       vrmat(1,2) = dx00 * vry
       vrmat(1,3) = dx00 * vrz
       vrmat(2,1) = dy00 * vrx
       vrmat(2,2) = dy00 * vry
       vrmat(2,3) = dy00 * vrz
       vrmat(3,1) = dz00 * vrx
       vrmat(3,2) = dz00 * vry
       vrmat(3,3) = dz00 * vrz
       !if( i.eq.1 .and. j.lt.10 )then
       !   write(STDOUT,*) i,j,vrx,dx00,(site%x(ii+kk),kk=1,5)
       !endif
#ifdef SHOWFORCE
       !見苦しいが、対力はここで表示させざるをえない。
       !TIP4P以外は未実装
       !対力は反発の時正、引力のとき負。
       fxs=fxi1+fxi2+fxi3+fxi4+fxi5
       fys=fyi1+fyi2+fyi3+fyi4+fyi5
       fzs=fzi1+fzi2+fzi3+fzi4+fzi5
       fpair=dsqrt(fxs**2+fys**2+fzs**2)
       if(dx00*fxs+dy00*fys+dz00*fzs.lt.0) fpair=-fpair
       !出力は原則として配列を0から数える。
       !write(STDOUT,*) i-1,j-1,dsqrt(dx00**2+dy00**2+dz00**2),fpair
#endif
       !write(6,*) i,j,k,ep
       !
       ! 非vector計算機の場合には、この場で力の集計をした方がてっとり
       !     早い。
#ifdef VPOPTIMIZE
       !     ......計算結果をlist vectorに格納する。
       fxi(k*TIP4PSITE+0)=fxi1
       fxi(k*TIP4PSITE+1)=fxi2
       fxi(k*TIP4PSITE+2)=fxi3
       fxi(k*TIP4PSITE+3)=fxi4
       fxi(k*TIP4PSITE+4)=fxi5
       fyi(k*TIP4PSITE+0)=fyi1
       fyi(k*TIP4PSITE+1)=fyi2
       fyi(k*TIP4PSITE+2)=fyi3
       fyi(k*TIP4PSITE+3)=fyi4
       fyi(k*TIP4PSITE+4)=fyi5
       fzi(k*TIP4PSITE+0)=fzi1
       fzi(k*TIP4PSITE+1)=fzi2
       fzi(k*TIP4PSITE+2)=fzi3
       fzi(k*TIP4PSITE+3)=fzi4
       fzi(k*TIP4PSITE+4)=fzi5
       
       fxj(k*TIP4PSITE+0)=-fxj1
       fxj(k*TIP4PSITE+1)=-fxj2
       fxj(k*TIP4PSITE+2)=-fxj3
       fxj(k*TIP4PSITE+3)=-fxj4
       fxj(k*TIP4PSITE+4)=-fxj5
       fyj(k*TIP4PSITE+0)=-fyj1
       fyj(k*TIP4PSITE+1)=-fyj2
       fyj(k*TIP4PSITE+2)=-fyj3
       fyj(k*TIP4PSITE+3)=-fyj4
       fyj(k*TIP4PSITE+4)=-fyj5
       fzj(k*TIP4PSITE+0)=-fzj1
       fzj(k*TIP4PSITE+1)=-fzj2
       fzj(k*TIP4PSITE+2)=-fzj3
       fzj(k*TIP4PSITE+3)=-fzj4
       fzj(k*TIP4PSITE+4)=-fzj5
#else 
       site%fx(ii+1)=site%fx(ii+1)+fxi1
       site%fx(ii+2)=site%fx(ii+2)+fxi2
       site%fx(ii+3)=site%fx(ii+3)+fxi3
       site%fx(ii+4)=site%fx(ii+4)+fxi4
       site%fx(ii+5)=site%fx(ii+5)+fxi5
       site%fx(jj+1)=site%fx(jj+1)-fxj1
       site%fx(jj+2)=site%fx(jj+2)-fxj2
       site%fx(jj+3)=site%fx(jj+3)-fxj3
       site%fx(jj+4)=site%fx(jj+4)-fxj4
       site%fx(jj+5)=site%fx(jj+5)-fxj5
       
       site%fy(ii+1)=site%fy(ii+1)+fyi1
       site%fy(ii+2)=site%fy(ii+2)+fyi2
       site%fy(ii+3)=site%fy(ii+3)+fyi3
       site%fy(ii+4)=site%fy(ii+4)+fyi4
       site%fy(ii+5)=site%fy(ii+5)+fyi5
       site%fy(jj+1)=site%fy(jj+1)-fyj1
       site%fy(jj+2)=site%fy(jj+2)-fyj2
       site%fy(jj+3)=site%fy(jj+3)-fyj3
       site%fy(jj+4)=site%fy(jj+4)-fyj4
       site%fy(jj+5)=site%fy(jj+5)-fyj5

       site%fz(ii+1)=site%fz(ii+1)+fzi1
       site%fz(ii+2)=site%fz(ii+2)+fzi2
       site%fz(ii+3)=site%fz(ii+3)+fzi3
       site%fz(ii+4)=site%fz(ii+4)+fzi4
       site%fz(ii+5)=site%fz(ii+5)+fzi5
       site%fz(jj+1)=site%fz(jj+1)-fzj1
       site%fz(jj+2)=site%fz(jj+2)-fzj2
       site%fz(jj+3)=site%fz(jj+3)-fzj3
       site%fz(jj+4)=site%fz(jj+4)-fzj4
       site%fz(jj+5)=site%fz(jj+5)-fzj5
       !f90  -O6 2.724u 0.020s 0:03.05 89.8% 0+11k 6+13io 0pf+0w
       !kf90 -O6 2.751u 0.014s 0:03.12 88.4% 0+11k 7+11io 0pf+0w
       !この部分を除去すると
       !kf90 -O6 1.427u 0.021s 0:01.77 81.3% 0+11k 7+11io 0pf+0w
#endif /*VPOPTIMIZE*/
       !endif
    enddo
#ifdef VPOPTIMIZE
#ifdef TEST11
    fxm1(:,:)=0d0
    fym1(:,:)=0d0
    fzm1(:,:)=0d0
    fxm2(:,:)=0d0
    fym2(:,:)=0d0
    fzm2(:,:)=0d0
    fxm3(:,:)=0d0
    fym3(:,:)=0d0
    fzm3(:,:)=0d0
    fxm4(:,:)=0d0
    fym4(:,:)=0d0
    fzm4(:,:)=0d0
    fxm5(:,:)=0d0
    fym5(:,:)=0d0
    fzm5(:,:)=0d0
!OCL NOVREC
    do k=1,iv%npair
       i=iv%pair_i(k)
       j=iv%pair_j(k)
       fxm1(j,i)=fxj(k*5)
       fym1(j,i)=fyj(k*TIP4PSITE)
       fzm1(j,i)=fzj(k*TIP4PSITE)
       
       fxm2(j,i)=fxj(k*TIP4PSITE+1)
       fym2(j,i)=fyj(k*TIP4PSITE+1)
       fzm2(j,i)=fzj(k*TIP4PSITE+1)
       
       fxm3(j,i)=fxj(k*TIP4PSITE+2)
       fym3(j,i)=fyj(k*TIP4PSITE+2)
       fzm3(j,i)=fzj(k*TIP4PSITE+2)
       
       fxm4(j,i)=fxj(k*TIP4PSITE+3)
       fym4(j,i)=fyj(k*TIP4PSITE+3)
       fzm4(j,i)=fzj(k*TIP4PSITE+3)
       
       fxm5(j,i)=fxj(k*TIP4PSITE+4)
       fym5(j,i)=fyj(k*TIP4PSITE+4)
       fzm5(j,i)=fzj(k*TIP4PSITE+4)
    enddo
  !異種グループ間相互作用の場合は、一旦iに加わる力を集計する。
    if(.not.iv%isomol)then
       do i=1,mi%nmol
          !OCL NOVREC
          do j=1,mj%nmol
             site%fx((j-1)*TIP4PSITE+mj%offset+1)=site%fx((j-1)&
                  & *TIP4PSITE+mj%offset+1)+fxm1(j,i)
             site%fx((j-1)*TIP4PSITE+mj%offset+2)=site%fx((j-1)&
                  & *TIP4PSITE+mj%offset+2)+fxm2(j,i)
             site%fx((j-1)*TIP4PSITE+mj%offset+3)=site%fx((j-1)&
                  & *TIP4PSITE+mj%offset+3)+fxm3(j,i)
             site%fx((j-1)*TIP4PSITE+mj%offset+4)=site%fx((j-1)&
                  & *TIP4PSITE+mj%offset+4)+fxm4(j,i)
             site%fx((j-1)*TIP4PSITE+mj%offset+5)=site%fx((j-1)&
                  & *TIP4PSITE+mj%offset+5)+fxm5(j,i)
             site%fy((j-1)*TIP4PSITE+mj%offset+1)=site%fy((j-1)&
                  & *TIP4PSITE+mj%offset+1)+fym1(j,i)
             site%fy((j-1)*TIP4PSITE+mj%offset+2)=site%fy((j-1)&
                  & *TIP4PSITE+mj%offset+2)+fym2(j,i)
             site%fy((j-1)*TIP4PSITE+mj%offset+3)=site%fy((j-1)&
                  & *TIP4PSITE+mj%offset+3)+fym3(j,i)
             site%fy((j-1)*TIP4PSITE+mj%offset+4)=site%fy((j-1)&
                  & *TIP4PSITE+mj%offset+4)+fym4(j,i)
             site%fy((j-1)*TIP4PSITE+mj%offset+5)=site%fy((j-1)&
                  & *TIP4PSITE+mj%offset+5)+fym5(j,i)
             site%fz((j-1)*TIP4PSITE+mj%offset+1)=site%fz((j-1)&
                  & *TIP4PSITE+mj%offset+1)+fzm1(j,i)
             site%fz((j-1)*TIP4PSITE+mj%offset+2)=site%fz((j-1)&
                  & *TIP4PSITE+mj%offset+2)+fzm2(j,i)
             site%fz((j-1)*TIP4PSITE+mj%offset+3)=site%fz((j-1)&
                  & *TIP4PSITE+mj%offset+3)+fzm3(j,i)
             site%fz((j-1)*TIP4PSITE+mj%offset+4)=site%fz((j-1)&
                  & *TIP4PSITE+mj%offset+4)+fzm4(j,i)
             site%fz((j-1)*TIP4PSITE+mj%offset+5)=site%fz((j-1)&
                  & *TIP4PSITE+mj%offset+5)+fzm5(j,i)
          enddo
       enddo
       fxm1(:,:)=0d0
       fym1(:,:)=0d0
       fzm1(:,:)=0d0
       fxm2(:,:)=0d0
       fym2(:,:)=0d0
       fzm2(:,:)=0d0
       fxm3(:,:)=0d0
       fym3(:,:)=0d0
       fzm3(:,:)=0d0
       fxm4(:,:)=0d0
       fym4(:,:)=0d0
       fzm4(:,:)=0d0
       fxm5(:,:)=0d0
       fym5(:,:)=0d0
       fzm5(:,:)=0d0
    endif
  !OCL NOVREC
    do k=1,iv%npair
       i=iv%pair_i(k)
       j=iv%pair_j(k)
       fxm1(i,j)=fxi(k*TIP4PSITE)
       fym1(i,j)=fyi(k*TIP4PSITE)
       fzm1(i,j)=fzi(k*TIP4PSITE)
       
       fxm2(i,j)=fxi(k*TIP4PSITE+1)
       fym2(i,j)=fyi(k*TIP4PSITE+1)
       fzm2(i,j)=fzi(k*TIP4PSITE+1)
       
       fxm3(i,j)=fxi(k*TIP4PSITE+2)
       fym3(i,j)=fyi(k*TIP4PSITE+2)
       fzm3(i,j)=fzi(k*TIP4PSITE+2)
       
       fxm4(i,j)=fxi(k*TIP4PSITE+3)
       fym4(i,j)=fyi(k*TIP4PSITE+3)
       fzm4(i,j)=fzi(k*TIP4PSITE+3)
       
       fxm5(i,j)=fxi(k*TIP4PSITE+4)
       fym5(i,j)=fyi(k*TIP4PSITE+4)
       fzm5(i,j)=fzi(k*TIP4PSITE+4)
    enddo
    do j=1,mj%nmol
       !OCL NOVREC
       do i=1,mi%nmol
          site%fx((i-1)*TIP4PSITE+mi%offset+1)=site%fx((i-1)&
               & *TIP4PSITE+mi%offset+1)+fxm1(i,j)
          site%fx((i-1)*TIP4PSITE+mi%offset+2)=site%fx((i-1)&
               & *TIP4PSITE+mi%offset+2)+fxm2(i,j)
          site%fx((i-1)*TIP4PSITE+mi%offset+3)=site%fx((i-1)&
               & *TIP4PSITE+mi%offset+3)+fxm3(i,j)
          site%fx((i-1)*TIP4PSITE+mi%offset+4)=site%fx((i-1)&
               & *TIP4PSITE+mi%offset+4)+fxm4(i,j)
          site%fx((i-1)*TIP4PSITE+mi%offset+5)=site%fx((i-1)&
               & *TIP4PSITE+mi%offset+5)+fxm5(i,j)
          site%fy((i-1)*TIP4PSITE+mi%offset+1)=site%fy((i-1)&
               & *TIP4PSITE+mi%offset+1)+fym1(i,j)
          site%fy((i-1)*TIP4PSITE+mi%offset+2)=site%fy((i-1)&
               & *TIP4PSITE+mi%offset+2)+fym2(i,j)
          site%fy((i-1)*TIP4PSITE+mi%offset+3)=site%fy((i-1)&
               & *TIP4PSITE+mi%offset+3)+fym3(i,j)
          site%fy((i-1)*TIP4PSITE+mi%offset+4)=site%fy((i-1)&
               & *TIP4PSITE+mi%offset+4)+fym4(i,j)
          site%fy((i-1)*TIP4PSITE+mi%offset+5)=site%fy((i-1)&
               & *TIP4PSITE+mi%offset+5)+fym5(i,j)
          site%fz((i-1)*TIP4PSITE+mi%offset+1)=site%fz((i-1)&
               & *TIP4PSITE+mi%offset+1)+fzm1(i,j)
          site%fz((i-1)*TIP4PSITE+mi%offset+2)=site%fz((i-1)&
               & *TIP4PSITE+mi%offset+2)+fzm2(i,j)
          site%fz((i-1)*TIP4PSITE+mi%offset+3)=site%fz((i-1)&
               & *TIP4PSITE+mi%offset+3)+fzm3(i,j)
          site%fz((i-1)*TIP4PSITE+mi%offset+4)=site%fz((i-1)&
               & *TIP4PSITE+mi%offset+4)+fzm4(i,j)
          site%fz((i-1)*TIP4PSITE+mi%offset+5)=site%fz((i-1)&
               & *TIP4PSITE+mi%offset+5)+fzm5(i,j)
       enddo
    enddo
    deallocate(fxi)
    deallocate(fyi)
    deallocate(fzi)
    deallocate(fxj)
    deallocate(fyj)
    deallocate(fzj)
#else /*TEST11*/
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
          site%fx(ii+mi%offset)=site%fx(ii+mi%offset)+fxi(k*TIP4PSITE&
               & +s-1)
          site%fy(ii+mi%offset)=site%fy(ii+mi%offset)+fyi(k*TIP4PSITE&
               & +s-1)
          site%fz(ii+mi%offset)=site%fz(ii+mi%offset)+fzi(k*TIP4PSITE&
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
          site%fx(jj+mj%offset)=site%fx(jj+mj%offset)+fxj(k*TIP4PSITE&
               & +s-1)
          site%fy(jj+mj%offset)=site%fy(jj+mj%offset)+fyj(k*TIP4PSITE&
               & +s-1)
          site%fz(jj+mj%offset)=site%fz(jj+mj%offset)+fzj(k*TIP4PSITE&
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
#endif /*TEST11*/
#endif /*VPOPTIMIZE*/
    !write(6,*) kk
#ifdef VERBOSE
    write(STDERR,*) "INTERACTION PAIRS: ",count
#endif
  end subroutine energy_tip4p

  !
  !O-H距離で、結合を判定する。
  !
  subroutine tip4p_hb_by_oh_distance(iv,mi,mj,site,from,which,to)
    !use physconst_module
    use interaction_module
    use mol_module
    use site_module
    !use montecarlo_module
    implicit none
    type(sInteraction),intent(in),   target :: iv
    type(sMol)        ,intent(in)    :: mi,mj
    type(sSite)       ,intent(inout),target :: site
    integer,           intent(inout) :: from(*),which(*),to(*)
    real(kind=8) :: dx,dy,dz
    !real(kind=8) :: d6,d12,e0,radiusi
    !real(kind=8) :: dd
    !real(kind=8) :: qq,aa,bb,ljsig,ljeps
    integer :: i,j,k
    integer :: ii,jj
    integer :: isite,jsite
    real(kind=8) :: distanceMinimum,distance
    integer :: pairMinimum,pair
    do k=1,iv%Npair
       pair = 0
       pairMinimum = 0
       distanceMinimum = +1d10
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       !write(STDOUT,*) "i,j:",i,j
       isite = 1
       do jsite = 3,4
          !     write(STDERR,*) i,j
          ii = (i-1)*mi%nsite+mi%offset + isite
          jj = (j-1)*mj%nsite+mj%offset + jsite
          dx = site%x(ii) - site%x(jj)
          dy = site%y(ii) - site%y(jj)
          dz = site%z(ii) - site%z(jj)
          dx = dx - iv%ox(k)
          dy = dy - iv%oy(k)
          dz = dz - iv%oz(k)
          distance = (dx**2 + dy**2 + dz**2)
          pair = pair + 1
          if ( distance < distanceMinimum ) then
             distanceMinimum = distance
             pairMinimum     = pair
          endif
       enddo
       do isite = 3,4
          jsite = 1
          !     write(STDERR,*) i,j
          ii = (i-1)*mi%nsite+mi%offset + isite
          jj = (j-1)*mj%nsite+mj%offset + jsite
          dx = site%x(ii) - site%x(jj)
          dy = site%y(ii) - site%y(jj)
          dz = site%z(ii) - site%z(jj)
          dx = dx - iv%ox(k)
          dy = dy - iv%oy(k)
          dz = dz - iv%oz(k)
          distance = (dx**2 + dy**2 + dz**2)
          pair = pair + 1
          if ( distance < distanceMinimum ) then
             distanceMinimum = distance
             pairMinimum     = pair
          endif
       enddo
       !
       !determine by O-H distance
       !
       if ( distanceMinimum < 2.5d0**2 ) then
          if ( pairMinimum <= 2 ) then
             from( k )  = j
             which( k ) = pairMinimum - 1
             to( k )    = i
          else 
             from( k )  = i
             which( k ) = pairMinimum - 3
             to( k )    = j
          endif
       else
          from( k )  = 0
          which( k ) = 0
          to( k )    = 0
       endif
    enddo
  end subroutine tip4p_hb_by_oh_distance

  !
  !graph(WGPH,NGPH)を出力する。これらのフォーマットは混合物を想定して
  !いない。
  !
  subroutine tip4p_SaveGraph(iv,mi,site,file,outputtype)
    !use physconst_module
    use interaction_module
    use mol_module
    use site_module
    !use montecarlo_module
    implicit none
    type(sInteraction),intent(in),   target :: iv
    type(sMol)        ,intent(in)    :: mi
    type(sSite)       ,intent(inout),target :: site
    integer                          :: from(iv%npair),which(iv%npair),to(iv%npair)
    integer,intent(in) :: file,outputtype
    !
    ! local
    !
    integer :: i
    call interaction_hb_by_oh_distance(iv,mi,mi,site,from,which,to)
    if ( outputtype == 1 ) then
       !
       !WGPH
       !
       write(file,'("@WGPH")')
       write(file,*) mi%nmol
       do i=1,iv%npair
          if ( from(i) > 0 ) then
             write(file,'(99i5)') from(i)-1,to(i)-1,which(i)
          endif
       enddo
       write(file,'(99i5)') -1,-1,-1
    else if ( outputtype == 2 ) then
       !
       !NGPH
       !
       write(file,'("@NGPH")')
       write(file,*) mi%nmol
       do i=1,iv%npair
          if ( from(i) > 0 ) then
             write(file,'(99i5)') from(i)-1,to(i)-1
          endif
       enddo
       write(file,'(99i5)') -1,-1
    else
       write(STDERR,*) "unknown output type:", outputtype
    endif
  end subroutine tip4p_SaveGraph

  !
  !2本の水素結合の方向から、分子の配向を決定する。
  !Tools/WaterConfigでやってる方法が使えるはず。
  !しかし、もっと簡単にできないのか？？回転行列からQuaternionを逆算で
  !きればよいのだが・・・
  !
  subroutine tip4p_vectors_to_quat(x1,y1,z1,x2,y2,z2,qa,qb,qc,qd)
    implicit none
    real(kind=8) :: r
    real(kind=8) :: x1,y1,z1
    real(kind=8) :: x2,y2,z2
    real(kind=8) :: kx,ky,kz
    real(kind=8) :: jx,jy,jz
    real(kind=8) :: ix,iy,iz
    real(kind=8) :: xi,yi,zi
    real(kind=8) :: xj,yj,zj
    real(kind=8) :: ax,ay,az
    real(kind=8) :: qa,qb,qc,qd
    real(kind=8) :: ox,oy,oz
    real(kind=8) :: sinh,cosh
    real(kind=8) :: i0x,i0y,i0z,cosine,t
    real(kind=8) :: x0x,x0y,x0z
    !まず規格化する
    r=1.0/sqrt(x1*x1+y1*y1+z1*z1)
    x1 = x1 * r
    y1 = y1 * r
    z1 = z1 * r
    r=1.0/sqrt(x2*x2+y2*y2+z2*z2)
    x2 = x2 * r
    y2 = y2 * r
    z2 = z2 * r
    !2分ベクトル(z)
    kx=x1+x2
    ky=y1+y2
    kz=z1+z2
    r=1.0/sqrt(kx*kx+ky*ky+kz*kz)
    kx = kx * r
    ky = ky * r
    kz = kz * r
    !直交ベクトル(y)
    jx=x2-x1
    jy=y2-y1
    jz=z2-z1
    r=1.0/sqrt(jx*jx+jy*jy+jz*jz)
    jx = jx * r
    jy = jy * r
    jz = jz * r
    !法線ベクトル(x=y*z)
    ix=jy*kz-jz*ky
    iy=jz*kx-jx*kz
    iz=jx*ky-jy*kx
    !i軸をx軸に移す回転の軸は、iとxの2分面となる。
    !j軸をy軸に移す回転の軸は、jとyの2分面となる。
    !そして、それらを同時にみたす回転の軸は、それらの交線となる。
    !交線は、2つの面の法線のいずれとも直交する=外積である。
    xi=ix-1
    yi=iy
    zi=iz
    xj=jx
    yj=jy-1
    zj=jz
    ax=yi*zj-zi*yj
    ay=zi*xj-xi*zj
    az=xi*yj-yi*xj
    !回転軸aが求まった。
    x0x=1
    x0y=0
    x0z=0
    i0x=ix
    i0y=iy
    i0z=iz
    t=ax/(ax*ax+ay*ay+az*az)
    i0x = i0x - t*ax
    i0y = i0y - t*ay
    i0z = i0z - t*az
    x0x = x0x - t*ax
    x0y = x0y - t*ay
    x0z = x0z - t*az
    !check
    !/*printf("%f %f\n",i0x*ax+i0y*ay+i0z*az,x0x*ax+x0y*ay+x0z*az)*/
    r=1.0/sqrt(i0x*i0x+i0y*i0y+i0z*i0z)
    i0x = i0x * r
    i0y = i0y * r
    i0z = i0z * r
    r=1.0/sqrt(x0x*x0x+x0y*x0y+x0z*x0z)
    x0x = x0x * r
    x0y = x0y * r
    x0z = x0z * r
    /*inner product to determine angle*/
    cosine=i0x*x0x+i0y*x0y+i0z*x0z
    cosh=sqrt((1.0+cosine)*0.5)
    sinh=sqrt(1.0-cosh*cosh)
    /*outer product to determine direction*/
    ox=i0y*x0z-i0z*x0y
    oy=i0z*x0x-i0x*x0z
    oz=i0x*x0y-i0y*x0x
    if(ox*ax+oy*ay+oz*az<0)then
       sinh=-sinh
    endif
    r=1.0/sqrt(ax*ax+ay*ay+az*az)
    ax = ax * r
    ay = ay * r
    az = az * r
    qa=cosh
    qb=-sinh*ax
    qc=+sinh*ay
    qd=-sinh*az
  end subroutine tip4p_vectors_to_quat
end module tip4p_module

! -*- f90 -*-
#define AAMO ( 4d0*dsqrt(EPSM*EPSO)*AAMOS)
#define BBMO (-4d0*dsqrt(EPSM*EPSO)*BBMOS)
#define AAOMX ( 4d0*dsqrt(EPSWO*EPSM)*AAOMXS)
#define BBOMX (-4d0*dsqrt(EPSWO*EPSM)*BBOMXS)
#define AAOOX ( 4d0*dsqrt(EPSWO*EPSO)*AAOOXS)
#define BBOOX (-4d0*dsqrt(EPSWO*EPSO)*BBOOXS)
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

module oplsmeoh_module
  use common_module
  use physconst_module
  use interaction_module
  implicit none
  !OPLS MeOHでは、Meの電荷は0.265を使用
  real(kind=8), parameter :: QM =0.265d0
  real(kind=8), parameter :: QO =(-0.700d0)
  real(kind=8), parameter :: QH =0.435d0
  !real(kind=8), parameter :: QM =0d0
  !real(kind=8), parameter :: QO =0d0
  !real(kind=8), parameter :: QH =0d0
  !OPLS Methanol J. Chem. Phys. 105,11199(1996)
  real(kind=8), parameter :: EPSM =(0.207d0*CA)
  real(kind=8), parameter :: SIGM =3.775d0
  real(kind=8), parameter :: EPSO =(0.170d0*CA)
  real(kind=8), parameter :: SIGO =3.070d0
  !OPLS Mephanol, See J. Chem. Phys. 99,9428(1993)
  !real(kind=8), parameter :: EPSM =0.799d0
  !real(kind=8), parameter :: SIGM =3.840d0
  !real(kind=8), parameter :: EPSO =0.711d0
  !real(kind=8), parameter :: SIGO =3.070d0
  real(kind=8), parameter :: AAMM =( 4d0*EPSM*SIGM**12*100d0)
  real(kind=8), parameter :: BBMM =(-4d0*EPSM*SIGM** 6*100d0)
  real(kind=8), parameter :: AAOO =( 4d0*EPSO*SIGO**12*100d0)
  real(kind=8), parameter :: BBOO =(-4d0*EPSO*SIGO** 6*100d0)
  real(kind=8), parameter :: AAMOS =(((SIGM+SIGO)*0.5d0)**12*100d0)
  real(kind=8), parameter :: BBMOS =(((SIGM+SIGO)*0.5d0)** 6*100d0)
  real(kind=8), parameter :: QQMM =(QM*QM*COEFF * 100d0)
  real(kind=8), parameter :: QQMO =(QM*QO*COEFF * 100d0)
  real(kind=8), parameter :: QQMH =(QM*QH*COEFF * 100d0)
  real(kind=8), parameter :: QQOO =(QO*QO*COEFF * 100d0)
  real(kind=8), parameter :: QQHO =(QH*QO*COEFF * 100d0)
  real(kind=8), parameter :: QQHH =(QH*QH*COEFF * 100d0)
  
  !重心を含めて4点
  integer, parameter :: MEOHSITE = 4
  !
  character(len=8),parameter     :: meohName(MEOHSITE) = (/"Me", "O ", "H ", "  "/)

contains
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine oplsmeoh_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, MEOHSITE, LJ_COULOMB)
    !charges [Q]
    si%param(1,1) = QM
    si%param(1,2) = QO
    si%param(1,3) = QH
    si%param(1,4) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = EPSM
    si%param(2,2) = EPSO
    si%param(2,3) = 0d0
    si%param(2,4) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = SIGM
    si%param(3,2) = SIGO
    si%param(3,3) = 0d0
    si%param(3,4) = 0d0
  end subroutine oplsmeoh_setinteraction

  subroutine interaction_force_oplsmeoh(iv,mi,mj,site,ep1,vrmat)
    use mol_module
    use site_module
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    real(kind=8),intent(out) :: vrmat(3,3)
    real(kind=8), intent(INOUT) :: ep1(*)
    real(kind=8) :: fxi1,fyi1,fzi1
    real(kind=8) :: fxi2,fyi2,fzi2
    real(kind=8) :: fxi3,fyi3,fzi3
    real(kind=8) :: fxi4,fyi4,fzi4
    real(kind=8) :: fxj1,fyj1,fzj1
    real(kind=8) :: fxj2,fyj2,fzj2
    real(kind=8) :: fxj3,fyj3,fzj3
    real(kind=8) :: fxj4,fyj4,fzj4
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
#ifdef VERBOSE
    integer :: count
#endif
    !integer :: kk
#ifdef VPOPTIMIZE
    !0番目の要素は、loop内のifをへらすためのトリック
    real(kind=8),dimension(:),allocatable :: fxi,fyi,fzi,fxj,fyj,fzj
#ifdef TEST11
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm1,fym1,fzm1
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm2,fym2,fzm2
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm3,fym3,fzm3
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm4,fym4,fzm4
#endif
    real(kind=8) :: r_r0,r_r1,r_r0s,r_r1s,dr
    allocate(fxi(0:MEOHSITE*(iv%npair+1)))
    allocate(fyi(0:MEOHSITE*(iv%npair+1)))
    allocate(fzi(0:MEOHSITE*(iv%npair+1)))
    allocate(fxj(0:MEOHSITE*(iv%npair+1)))
    allocate(fyj(0:MEOHSITE*(iv%npair+1)))
    allocate(fzj(0:MEOHSITE*(iv%npair+1)))
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
    !kk=0
  !VPP用のベクトル化制御行
  !OCL VECTOR,REPEAT(1000000)
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
  !     write(STDERR,*) i,j
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+MEOHSITE)
       dy00 = site%y(ii+MEOHSITE)
       dz00 = site%z(ii+MEOHSITE)
       dx00 = dx00 - site%x(jj+MEOHSITE)
       dy00 = dy00 - site%y(jj+MEOHSITE)
       dz00 = dz00 - site%z(jj+MEOHSITE)
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
  
          ep=0d0
  
          !Me-Me
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
          radiusi = dsqrt(dd)
          e0 = QQMM*radiusi
          f0 = QQMM*radiusi*dd
          d6 = dd*dd*dd
          d12= d6*d6
          e0 = e0 + AAMM*d12 + BBMM*d6
          f0 = f0 + dd*(AAMM*12d0*d12 + BBMM*6d0*d6)
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
  
          !Me-O
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
          radiusi = dsqrt(dd)
          e0 = QQMO*radiusi
          f0 = QQMO*radiusi*dd
          d6 = dd*dd*dd
          d12= d6*d6
          e0 = e0 + AAMO*d12 + BBMO*d6
          f0 = f0 + dd*(AAMO*12d0*d12 + BBMO*6d0*d6)
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
  
          !Me-H
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
          radiusi = dsqrt(dd)
          e0 = QQMH*radiusi
          f0 = QQMH*radiusi*dd
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
  
          !O-Me
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
          radiusi = dsqrt(dd)
          e0 = QQMO*radiusi
          f0 = QQMO*radiusi*dd
          d6 = dd*dd*dd
          d12= d6*d6
          e0 = e0 + AAMO*d12 + BBMO*d6
          f0 = f0 + dd*(AAMO*12d0*d12 + BBMO*6d0*d6)
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
  
          !O-O
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
          e0 = QQOO*radiusi
          f0 = QQOO*radiusi*dd
          d6 = dd*dd*dd
          d12= d6*d6
          e0 = e0 + AAOO*d12 + BBOO*d6
          f0 = f0 + dd*(AAOO*12d0*d12 + BBOO*6d0*d6)
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
  
          !O-H
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
          e0 = QQHO*radiusi
          f0 = QQHO*radiusi*dd
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
  
          !H-Me
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
          fxj1 = fxj1 + forcex
          fyj1 = fyj1 + forcey
          fzj1 = fzj1 + forcez
  
          !H-O
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
          e0 = QQHO*radiusi
          f0 = QQHO*radiusi*dd
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
  
          !H-H
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
          dhx=dhx*ep
          dhy=dhy*ep
          dhz=dhz*ep
          fxi4=fxi4+dhx
          fyi4=fyi4+dhy
          fzi4=fzi4+dhz
          fxj4=fxj4+dhx
          fyj4=fyj4+dhy
          fzj4=fzj4+dhz
          ep1(k) = ep * hh
          !vrx = dx00*(fxi1+fxi2+fxi3+fxi4)
          !vry = dy00*(fyi1+fyi2+fyi3+fyi4)
          !vrz = dz00*(fzi1+fzi2+fzi3+fzi4)
          !vrsum=vrsum+vrx+vry+vrz
          vrx = (fxi1+fxi2+fxi3+fxi4)
          vry = (fyi1+fyi2+fyi3+fyi4)
          vrz = (fzi1+fzi2+fzi3+fzi4)
          vrmat(1,1) = dx00 * vrx
          vrmat(1,2) = dx00 * vry
          vrmat(1,3) = dx00 * vrz
          vrmat(2,1) = dy00 * vrx
          vrmat(2,2) = dy00 * vry
          vrmat(2,3) = dy00 * vrz
          vrmat(3,1) = dz00 * vrx
          vrmat(3,2) = dz00 * vry
          vrmat(3,3) = dz00 * vrz
          !write(6,*) i,j,k,ep
          !     非vector計算機の場合には、この場で力の集計をした方がてっとり
          !     早い。
#ifdef VPOPTIMIZE
          !     ......計算結果をlist vectorに格納する。
          fxi(k*MEOHSITE+0)=fxi1
          fxi(k*MEOHSITE+1)=fxi2
          fxi(k*MEOHSITE+2)=fxi3
          fxi(k*MEOHSITE+3)=fxi4
          fyi(k*MEOHSITE+0)=fyi1
          fyi(k*MEOHSITE+1)=fyi2
          fyi(k*MEOHSITE+2)=fyi3
          fyi(k*MEOHSITE+3)=fyi4
          fzi(k*MEOHSITE+0)=fzi1
          fzi(k*MEOHSITE+1)=fzi2
          fzi(k*MEOHSITE+2)=fzi3
          fzi(k*MEOHSITE+3)=fzi4
  
          fxj(k*MEOHSITE+0)=-fxj1
          fxj(k*MEOHSITE+1)=-fxj2
          fxj(k*MEOHSITE+2)=-fxj3
          fxj(k*MEOHSITE+3)=-fxj4
          fyj(k*MEOHSITE+0)=-fyj1
          fyj(k*MEOHSITE+1)=-fyj2
          fyj(k*MEOHSITE+2)=-fyj3
          fyj(k*MEOHSITE+3)=-fyj4
          fzj(k*MEOHSITE+0)=-fzj1
          fzj(k*MEOHSITE+1)=-fzj2
          fzj(k*MEOHSITE+2)=-fzj3
          fzj(k*MEOHSITE+3)=-fzj4
#else 
          site%fx(ii+1)=site%fx(ii+1)+fxi1
          site%fx(ii+2)=site%fx(ii+2)+fxi2
          site%fx(ii+3)=site%fx(ii+3)+fxi3
          site%fx(ii+4)=site%fx(ii+4)+fxi4
          site%fx(jj+1)=site%fx(jj+1)-fxj1
          site%fx(jj+2)=site%fx(jj+2)-fxj2
          site%fx(jj+3)=site%fx(jj+3)-fxj3
          site%fx(jj+4)=site%fx(jj+4)-fxj4
  
          site%fy(ii+1)=site%fy(ii+1)+fyi1
          site%fy(ii+2)=site%fy(ii+2)+fyi2
          site%fy(ii+3)=site%fy(ii+3)+fyi3
          site%fy(ii+4)=site%fy(ii+4)+fyi4
          site%fy(jj+1)=site%fy(jj+1)-fyj1
          site%fy(jj+2)=site%fy(jj+2)-fyj2
          site%fy(jj+3)=site%fy(jj+3)-fyj3
          site%fy(jj+4)=site%fy(jj+4)-fyj4
  
          site%fz(ii+1)=site%fz(ii+1)+fzi1
          site%fz(ii+2)=site%fz(ii+2)+fzi2
          site%fz(ii+3)=site%fz(ii+3)+fzi3
          site%fz(ii+4)=site%fz(ii+4)+fzi4
          site%fz(jj+1)=site%fz(jj+1)-fzj1
          site%fz(jj+2)=site%fz(jj+2)-fzj2
          site%fz(jj+3)=site%fz(jj+3)-fzj3
          site%fz(jj+4)=site%fz(jj+4)-fzj4
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
  !OCL NOVREC
    do k=1,iv%npair
       i=iv%pair_i(k)
       j=iv%pair_j(k)
       fxm1(j,i)=fxj(k*MEOHSITE)
       fym1(j,i)=fyj(k*MEOHSITE)
       fzm1(j,i)=fzj(k*MEOHSITE)
  
       fxm2(j,i)=fxj(k*MEOHSITE+1)
       fym2(j,i)=fyj(k*MEOHSITE+1)
       fzm2(j,i)=fzj(k*MEOHSITE+1)
  
       fxm3(j,i)=fxj(k*MEOHSITE+2)
       fym3(j,i)=fyj(k*MEOHSITE+2)
       fzm3(j,i)=fzj(k*MEOHSITE+2)
  
       fxm4(j,i)=fxj(k*MEOHSITE+3)
       fym4(j,i)=fyj(k*MEOHSITE+3)
       fzm4(j,i)=fzj(k*MEOHSITE+3)
    enddo
    !異種グループ間相互作用の場合は、一旦iに加わる力を集計する。
    if(.not.iv%isomol)then
       do i=1,mi%nmol
  !OCL NOVREC
          do j=1,mj%nmol
             site%fx((j-1)*MEOHSITE+mj%offset+1)=site%fx((j-1)*MEOHSITE+mj%offset+1)+fxm1(j,i)
             site%fx((j-1)*MEOHSITE+mj%offset+2)=site%fx((j-1)*MEOHSITE+mj%offset+2)+fxm2(j,i)
             site%fx((j-1)*MEOHSITE+mj%offset+3)=site%fx((j-1)*MEOHSITE+mj%offset+3)+fxm3(j,i)
             site%fx((j-1)*MEOHSITE+mj%offset+4)=site%fx((j-1)*MEOHSITE+mj%offset+4)+fxm4(j,i)
             site%fy((j-1)*MEOHSITE+mj%offset+1)=site%fy((j-1)*MEOHSITE+mj%offset+1)+fym1(j,i)
             site%fy((j-1)*MEOHSITE+mj%offset+2)=site%fy((j-1)*MEOHSITE+mj%offset+2)+fym2(j,i)
             site%fy((j-1)*MEOHSITE+mj%offset+3)=site%fy((j-1)*MEOHSITE+mj%offset+3)+fym3(j,i)
             site%fy((j-1)*MEOHSITE+mj%offset+4)=site%fy((j-1)*MEOHSITE+mj%offset+4)+fym4(j,i)
             site%fz((j-1)*MEOHSITE+mj%offset+1)=site%fz((j-1)*MEOHSITE+mj%offset+1)+fzm1(j,i)
             site%fz((j-1)*MEOHSITE+mj%offset+2)=site%fz((j-1)*MEOHSITE+mj%offset+2)+fzm2(j,i)
             site%fz((j-1)*MEOHSITE+mj%offset+3)=site%fz((j-1)*MEOHSITE+mj%offset+3)+fzm3(j,i)
             site%fz((j-1)*MEOHSITE+mj%offset+4)=site%fz((j-1)*MEOHSITE+mj%offset+4)+fzm4(j,i)
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
    endif
  !OCL NOVREC
    do k=1,iv%npair
       i=iv%pair_i(k)
       j=iv%pair_j(k)
       fxm1(i,j)=fxi(k*MEOHSITE)
       fym1(i,j)=fyi(k*MEOHSITE)
       fzm1(i,j)=fzi(k*MEOHSITE)
  
       fxm2(i,j)=fxi(k*MEOHSITE+1)
       fym2(i,j)=fyi(k*MEOHSITE+1)
       fzm2(i,j)=fzi(k*MEOHSITE+1)
  
       fxm3(i,j)=fxi(k*MEOHSITE+2)
       fym3(i,j)=fyi(k*MEOHSITE+2)
       fzm3(i,j)=fzi(k*MEOHSITE+2)
  
       fxm4(i,j)=fxi(k*MEOHSITE+3)
       fym4(i,j)=fyi(k*MEOHSITE+3)
       fzm4(i,j)=fzi(k*MEOHSITE+3)
    enddo
    do j=1,mj%nmol
  !OCL NOVREC
       do i=1,mi%nmol
          site%fx((i-1)*MEOHSITE+mi%offset+1)=site%fx((i-1)*MEOHSITE+mi%offset+1)+fxm1(i,j)
          site%fx((i-1)*MEOHSITE+mi%offset+2)=site%fx((i-1)*MEOHSITE+mi%offset+2)+fxm2(i,j)
          site%fx((i-1)*MEOHSITE+mi%offset+3)=site%fx((i-1)*MEOHSITE+mi%offset+3)+fxm3(i,j)
          site%fx((i-1)*MEOHSITE+mi%offset+4)=site%fx((i-1)*MEOHSITE+mi%offset+4)+fxm4(i,j)
          site%fy((i-1)*MEOHSITE+mi%offset+1)=site%fy((i-1)*MEOHSITE+mi%offset+1)+fym1(i,j)
          site%fy((i-1)*MEOHSITE+mi%offset+2)=site%fy((i-1)*MEOHSITE+mi%offset+2)+fym2(i,j)
          site%fy((i-1)*MEOHSITE+mi%offset+3)=site%fy((i-1)*MEOHSITE+mi%offset+3)+fym3(i,j)
          site%fy((i-1)*MEOHSITE+mi%offset+4)=site%fy((i-1)*MEOHSITE+mi%offset+4)+fym4(i,j)
          site%fz((i-1)*MEOHSITE+mi%offset+1)=site%fz((i-1)*MEOHSITE+mi%offset+1)+fzm1(i,j)
          site%fz((i-1)*MEOHSITE+mi%offset+2)=site%fz((i-1)*MEOHSITE+mi%offset+2)+fzm2(i,j)
          site%fz((i-1)*MEOHSITE+mi%offset+3)=site%fz((i-1)*MEOHSITE+mi%offset+3)+fzm3(i,j)
          site%fz((i-1)*MEOHSITE+mi%offset+4)=site%fz((i-1)*MEOHSITE+mi%offset+4)+fzm4(i,j)
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
          i=mi%lvmol(ii)
          s=mi%lvsite(ii)
          k=iv%partner_i(i,l)
          site%fx(ii+mi%offset)=site%fx(ii+mi%offset)+fxi(k*MEOHSITE+s-1)
          site%fy(ii+mi%offset)=site%fy(ii+mi%offset)+fyi(k*MEOHSITE+s-1)
          site%fz(ii+mi%offset)=site%fz(ii+mi%offset)+fzi(k*MEOHSITE+s-1)
       enddo
    enddo
    do l=1,iv%maxpartner_j
       do jj=1,mj%nsite*mj%nmol
          j=mj%lvmol(jj)
          s=mj%lvsite(jj)
          k=iv%partner_j(j,l)
          site%fx(jj+mj%offset)=site%fx(jj+mj%offset)+fxj(k*MEOHSITE+s-1)
          site%fy(jj+mj%offset)=site%fy(jj+mj%offset)+fyj(k*MEOHSITE+s-1)
          site%fz(jj+mj%offset)=site%fz(jj+mj%offset)+fzj(k*MEOHSITE+s-1)
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
  end subroutine interaction_force_oplsmeoh
  
  !
  !計算順序を変更。OPLS汎用に改良しやすいかも。
  !平成15年1月10日(金)動作確認
  !速度はまだチェックしていない。
  !集計部分はVP用に書かれていない。
  !
  subroutine force2_oplsmeoh(iv,mi,mj,site,ep,vrmat)
    use mol_module
    use site_module
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    real(kind=8),intent(out) :: vrmat(3,3)
    real(kind=8),intent(inout) :: ep(*)
    real(kind=8) :: dx,dy,dz,dx00,dy00,dz00
    real(kind=8) :: dhx,dhy,dhz
    real(kind=8) :: d6,d12,e0,f0,radiusi
    real(kind=8) :: xi,yi,zi,xj,yj,zj
    real(kind=8) forcex,forcey,forcez
    real(kind=8) ox,oy,oz
    real(kind=8) dd,hh,dh
    real(kind=8) :: vrx,vry,vrz
    integer :: i,j,k,l
    integer :: ii,jj
    !integer :: kk
    !0番目の要素は、loop内のifをへらすためのトリック
    real(kind=8),dimension(:,:),allocatable :: fxi,fyi,fzi,fxj,fyj,fzj
    allocate(fxi(iv%npair,MEOHSITE))
    allocate(fyi(iv%npair,MEOHSITE))
    allocate(fzi(iv%npair,MEOHSITE))
    allocate(fxj(iv%npair,MEOHSITE))
    allocate(fyj(iv%npair,MEOHSITE))
    allocate(fzj(iv%npair,MEOHSITE))
    fxi(:,:) = 0d0
    fyi(:,:) = 0d0
    fzi(:,:) = 0d0
    fxj(:,:) = 0d0
    fyj(:,:) = 0d0
    fzj(:,:) = 0d0
    ep(1:iv%npair)   = 0d0
    !kk=0
    !VPP用のベクトル化制御行
    !OCL VECTOR,REPEAT(1000000)
    !Me-Me
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
  !     write(STDERR,*) i,j
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       xi = site%x(ii+1)
       yi = site%y(ii+1)
       zi = site%z(ii+1)
       xj = site%x(jj+1)
       yj = site%y(jj+1)
       zj = site%z(jj+1)
       dx = xi-xj-ox
       dy = yi-yj-oy
       dz = zi-zj-oz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQMM*radiusi
       f0 = QQMM*radiusi*dd
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = e0 + AAMM*d12 + BBMM*d6
       f0 = f0 + dd*(AAMM*12d0*d12 + BBMM*6d0*d6)
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi(k,1) = fxi(k,1) + forcex
       fyi(k,1) = fyi(k,1) + forcey
       fzi(k,1) = fzi(k,1) + forcez
       fxj(k,1) = fxj(k,1) + forcex
       fyj(k,1) = fyj(k,1) + forcey
       fzj(k,1) = fzj(k,1) + forcez
       ep(k)    = ep(k) + e0
    enddo
    !Me-O
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
  !     write(STDERR,*) i,j
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       xi = site%x(ii+1)
       yi = site%y(ii+1)
       zi = site%z(ii+1)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj-ox
       dy = yi-yj-oy
       dz = zi-zj-oz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQMO*radiusi
       f0 = QQMO*radiusi*dd
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = e0 + AAMO*d12 + BBMO*d6
       f0 = f0 + dd*(AAMO*12d0*d12 + BBMO*6d0*d6)
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi(k,1) = fxi(k,1) + forcex
       fyi(k,1) = fyi(k,1) + forcey
       fzi(k,1) = fzi(k,1) + forcez
       fxj(k,2) = fxj(k,2) + forcex
       fyj(k,2) = fyj(k,2) + forcey
       fzj(k,2) = fzj(k,2) + forcez
       ep(k)               = ep(k) + e0
    enddo
    !Me-H
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
  !     write(STDERR,*) i,j
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       xi = site%x(ii+1)
       yi = site%y(ii+1)
       zi = site%z(ii+1)
       xj = site%x(jj+3)
       yj = site%y(jj+3)
       zj = site%z(jj+3)
       dx = xi-xj-ox
       dy = yi-yj-oy
       dz = zi-zj-oz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQMH*radiusi
       f0 = QQMH*radiusi*dd
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi(k,1) = fxi(k,1) + forcex
       fyi(k,1) = fyi(k,1) + forcey
       fzi(k,1) = fzi(k,1) + forcez
       fxj(k,3) = fxj(k,3) + forcex
       fyj(k,3) = fyj(k,3) + forcey
       fzj(k,3) = fzj(k,3) + forcez
       ep(k)               = ep(k) + e0
    enddo
    !O-Me
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
  !     write(STDERR,*) i,j
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+1)
       yj = site%y(jj+1)
       zj = site%z(jj+1)
       dx = xi-xj-ox
       dy = yi-yj-oy
       dz = zi-zj-oz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQMO*radiusi
       f0 = QQMO*radiusi*dd
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = e0 + AAMO*d12 + BBMO*d6
       f0 = f0 + dd*(AAMO*12d0*d12 + BBMO*6d0*d6)
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi(k,2) = fxi(k,2) + forcex
       fyi(k,2) = fyi(k,2) + forcey
       fzi(k,2) = fzi(k,2) + forcez
       fxj(k,1) = fxj(k,1) + forcex
       fyj(k,1) = fyj(k,1) + forcey
       fzj(k,1) = fzj(k,1) + forcez
       ep(k)               = ep(k) + e0
    enddo
    !O-O
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
  !     write(STDERR,*) i,j
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj-ox
       dy = yi-yj-oy
       dz = zi-zj-oz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQOO*radiusi
       f0 = QQOO*radiusi*dd
       d6 = dd*dd*dd
       d12= d6*d6
       e0 = e0 + AAOO*d12 + BBOO*d6
       f0 = f0 + dd*(AAOO*12d0*d12 + BBOO*6d0*d6)
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi(k,2) = fxi(k,2) + forcex
       fyi(k,2) = fyi(k,2) + forcey
       fzi(k,2) = fzi(k,2) + forcez
       fxj(k,2) = fxj(k,2) + forcex
       fyj(k,2) = fyj(k,2) + forcey
       fzj(k,2) = fzj(k,2) + forcez
       ep(k)               = ep(k) + e0
    enddo
    !O-H
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
  !     write(STDERR,*) i,j
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+3)
       yj = site%y(jj+3)
       zj = site%z(jj+3)
       dx = xi-xj-ox
       dy = yi-yj-oy
       dz = zi-zj-oz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQHO*radiusi
       f0 = QQHO*radiusi*dd
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi(k,2) = fxi(k,2) + forcex
       fyi(k,2) = fyi(k,2) + forcey
       fzi(k,2) = fzi(k,2) + forcez
       fxj(k,3) = fxj(k,3) + forcex
       fyj(k,3) = fyj(k,3) + forcey
       fzj(k,3) = fzj(k,3) + forcez
       ep(k)               = ep(k) + e0
    enddo
    !H-Me
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
  !     write(STDERR,*) i,j
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       xi = site%x(ii+3)
       yi = site%y(ii+3)
       zi = site%z(ii+3)
       xj = site%x(jj+1)
       yj = site%y(jj+1)
       zj = site%z(jj+1)
       dx = xi-xj-ox
       dy = yi-yj-oy
       dz = zi-zj-oz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQMH*radiusi
       f0 = QQMH*radiusi*dd
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi(k,3) = fxi(k,3) + forcex
       fyi(k,3) = fyi(k,3) + forcey
       fzi(k,3) = fzi(k,3) + forcez
       fxj(k,1) = fxj(k,1) + forcex
       fyj(k,1) = fyj(k,1) + forcey
       fzj(k,1) = fzj(k,1) + forcez
       ep(k)               = ep(k) + e0
    enddo
    !H-O
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
  !     write(STDERR,*) i,j
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       xi = site%x(ii+3)
       yi = site%y(ii+3)
       zi = site%z(ii+3)
       xj = site%x(jj+2)
       yj = site%y(jj+2)
       zj = site%z(jj+2)
       dx = xi-xj-ox
       dy = yi-yj-oy
       dz = zi-zj-oz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQHO*radiusi
       f0 = QQHO*radiusi*dd
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi(k,3) = fxi(k,3) + forcex
       fyi(k,3) = fyi(k,3) + forcey
       fzi(k,3) = fzi(k,3) + forcez
       fxj(k,2) = fxj(k,2) + forcex
       fyj(k,2) = fyj(k,2) + forcey
       fzj(k,2) = fzj(k,2) + forcez
       ep(k)               = ep(k) + e0
    enddo
    !H-H
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
  !     write(STDERR,*) i,j
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       xi = site%x(ii+3)
       yi = site%y(ii+3)
       zi = site%z(ii+3)
       xj = site%x(jj+3)
       yj = site%y(jj+3)
       zj = site%z(jj+3)
       dx = xi-xj-ox
       dy = yi-yj-oy
       dz = zi-zj-oz
       dd = 1d0/(dx**2 + dy**2 + dz**2)
       radiusi = dsqrt(dd)
       e0 = QQHH*radiusi
       f0 = QQHH*radiusi*dd
       forcex = dx*f0
       forcey = dy*f0
       forcez = dz*f0
       fxi(k,3) = fxi(k,3) + forcex
       fyi(k,3) = fyi(k,3) + forcey
       fzi(k,3) = fzi(k,3) + forcez
       fxj(k,3) = fxj(k,3) + forcex
       fyj(k,3) = fyj(k,3) + forcey
       fzj(k,3) = fzj(k,3) + forcez
       ep(k)               = ep(k) + e0
    enddo
    !
    !apply smoothing function
    !
    do l=1,MEOHSITE
       do k=1,iv%Npair
          hh                  = iv%eratio( k )
          fxi(k,l) = fxi(k,l) * hh
          fyi(k,l) = fyi(k,l) * hh
          fzi(k,l) = fzi(k,l) * hh
       enddo
    enddo
    do l=1,MEOHSITE
       do k=1,iv%Npair
          hh                  = iv%eratio( k )
          fxj(k,l) = fxj(k,l) * hh
          fyj(k,l) = fyj(k,l) * hh
          fzj(k,l) = fzj(k,l) * hh
       enddo
    enddo
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
  !     write(STDERR,*) i,j
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+MEOHSITE) - site%x(jj+MEOHSITE) - iv%ox(k)
       dy00 = site%y(ii+MEOHSITE) - site%y(jj+MEOHSITE) - iv%oy(k)
       dz00 = site%z(ii+MEOHSITE) - site%z(jj+MEOHSITE) - iv%oz(k)
       dh   = iv%fratio( k )
       dhx  = dh*dx00*ep(k)
       dhy  = dh*dy00*ep(k)
       dhz  = dh*dz00*ep(k)
       fxi(k,MEOHSITE) = fxi(k,MEOHSITE) + dhx
       fyi(k,MEOHSITE) = fyi(k,MEOHSITE) + dhy
       fzi(k,MEOHSITE) = fzi(k,MEOHSITE) + dhz
       fxj(k,MEOHSITE) = fxj(k,MEOHSITE) + dhx
       fyj(k,MEOHSITE) = fyj(k,MEOHSITE) + dhy
       fzj(k,MEOHSITE) = fzj(k,MEOHSITE) + dhz
       do l=1,MEOHSITE
          site%fx(ii+l) = site%fx(ii+l) + fxi(k,l)
          site%fy(ii+l) = site%fy(ii+l) + fyi(k,l)
          site%fz(ii+l) = site%fz(ii+l) + fzi(k,l)
       enddo
       do l=1,MEOHSITE
          site%fx(jj+l) = site%fx(jj+l) - fxj(k,l)
          site%fy(jj+l) = site%fy(jj+l) - fyj(k,l)
          site%fz(jj+l) = site%fz(jj+l) - fzj(k,l)
       enddo
    enddo
    do k=1,iv%Npair
       ep(k) = ep(k) * iv%eratio( k )
    enddo
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+MEOHSITE) - site%x(jj+MEOHSITE) - iv%ox(k)
       dy00 = site%y(ii+MEOHSITE) - site%y(jj+MEOHSITE) - iv%oy(k)
       dz00 = site%z(ii+MEOHSITE) - site%z(jj+MEOHSITE) - iv%oz(k)
       vrx = 0d0
       vry = 0d0
       vrz = 0d0
       do l=1,MEOHSITE
          vrx = vrx + fxi(k,l)
          vry = vry + fyi(k,l)
          vrz = vrz + fzi(k,l)
       enddo
       vrmat(1,1) = dx00 * vrx
       vrmat(1,2) = dx00 * vry
       vrmat(1,3) = dx00 * vrz
       vrmat(2,1) = dy00 * vrx
       vrmat(2,2) = dy00 * vry
       vrmat(2,3) = dy00 * vrz
       vrmat(3,1) = dz00 * vrx
       vrmat(3,2) = dz00 * vry
       vrmat(3,3) = dz00 * vrz
    enddo
    deallocate(fxi)
    deallocate(fyi)
    deallocate(fzi)
    deallocate(fxj)
    deallocate(fyj)
    deallocate(fzj)
  end subroutine force2_oplsmeoh

end module oplsmeoh_module


module oplsmeoh_tip4p_module
  use tip4p_constant_module
  use oplsmeoh_module
!TIP4P-OPLSMeOH interaction
  real(kind=8), parameter :: MMX= (WM*QM*COEFF * 100d0)
  real(kind=8), parameter :: MOX= (WM*QO*COEFF * 100d0)
  real(kind=8), parameter :: MHX= (WM*QH*COEFF * 100d0)
  real(kind=8), parameter :: HMX= (WH*QM*COEFF * 100d0)
  real(kind=8), parameter :: HOX= (WH*QO*COEFF * 100d0)
  real(kind=8), parameter :: HHX= (WH*QH*COEFF * 100d0)
  real(kind=8), parameter :: AAOMXS =(((SIGWO+SIGM)*0.5d0)**12*100d0)
  real(kind=8), parameter :: BBOMXS =(((SIGWO+SIGM)*0.5d0)** 6*100d0)
  real(kind=8), parameter :: AAOOXS =(((SIGWO+SIGO)*0.5d0)**12*100d0)
  real(kind=8), parameter :: BBOOXS =(((SIGWO+SIGO)*0.5d0)** 6*100d0)

contains
  subroutine interaction_force_tip4p_oplsmeoh(iv,wa,me,site,ep1,vrmat)
    use mol_module
    use site_module
    use tip4p_constant_module
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: wa,me
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
#ifdef VERBOSE
    integer :: count
#endif
    !integer :: kk
#ifdef VPOPTIMIZE
    !0番目の要素は、loop内のifをへらすためのトリック
    real(kind=8),dimension(:),allocatable :: fxi,fyi,fzi,fxj,fyj,fzj
    real(kind=8) :: r_r0,r_r1,r_r0s,r_r1s,dr
    allocate(fxi(0:TIP4PSITE*(iv%npair+1)))
    allocate(fyi(0:TIP4PSITE*(iv%npair+1)))
    allocate(fzi(0:TIP4PSITE*(iv%npair+1)))
    allocate(fxj(0:MEOHSITE*(iv%npair+1)))
    allocate(fyj(0:MEOHSITE*(iv%npair+1)))
    allocate(fzj(0:MEOHSITE*(iv%npair+1)))
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
    !kk=0
  !VPP用のベクトル化制御行
  !OCL VECTOR,REPEAT(1000000)
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       !write(STDERR,*) i,j,iv%npair0,iv%npair
       ii = (i-1)*wa%nsite+wa%offset
       jj = (j-1)*me%nsite+me%offset
       !重心の座標
       dx00 = site%x(ii+TIP4PSITE)
       dy00 = site%y(ii+TIP4PSITE)
       dz00 = site%z(ii+TIP4PSITE)
       !重心の座標
       dx00 = dx00 - site%x(jj+MEOHSITE)
       dy00 = dy00 - site%y(jj+MEOHSITE)
       dz00 = dz00 - site%z(jj+MEOHSITE)
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

          ep=0d0
  
          !O-Me
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
          !radiusi = dsqrt(dd)
          d6 = dd*dd*dd
          d12= d6*d6
          e0 = AAOMX*d12 + BBOMX*d6
          f0 = dd*(AAOMX*12d0*d12 + BBOMX*6d0*d6)
          ep = ep + e0
          !write(6,*) "om",e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi1 = fxi1 + forcex
          fyi1 = fyi1 + forcey
          fzi1 = fzi1 + forcez
          fxj1 = fxj1 + forcex
          fyj1 = fyj1 + forcey
          fzj1 = fzj1 + forcez
  
          !O-O
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
          !radiusi = dsqrt(dd)
          d6 = dd*dd*dd
          d12= d6*d6
          e0 = AAOOX*d12 + BBOOX*d6
          f0 = dd*(AAOOX*12d0*d12 + BBOOX*6d0*d6)
          !write(6,*) "oo",e0,AAOOX,BBOOX,AAOMX,BBOMX,AA,BB
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
  
          !O-H
          
          !M-Me
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
          radiusi = dsqrt(dd)
          e0 = MMX*radiusi
          f0 = MMX*radiusi*dd
          !write(6,*) "mm",e0
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
  
          !M-O
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
          e0 = MOX*radiusi
          f0 = MOX*radiusi*dd
          !write(6,*) "mo",e0
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
  
          !M-H
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
          e0 = MHX*radiusi
          f0 = MHX*radiusi*dd
          !write(6,*) "mh",e0
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
  
          !H1-Me
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
          radiusi = dsqrt(dd)
          e0 = HMX*radiusi
          f0 = HMX*radiusi*dd
          !!write(6,*) "hm",e0
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
  
          !H1-O
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
          e0 = HOX*radiusi
          f0 = HOX*radiusi*dd
          !write(6,*) "ho",e0
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
  
          !H1-H
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
          e0 = HHX*radiusi
          f0 = HHX*radiusi*dd
          !write(6,*) "hh",e0
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
  
          !H2-Me
          xi = site%x(ii+4)
          yi = site%y(ii+4)
          zi = site%z(ii+4)
          xj = site%x(jj+1)
          yj = site%y(jj+1)
          zj = site%z(jj+1)
          dx = xi-xj+cx
          dy = yi-yj+cy
          dz = zi-zj+cz
          dd = 1d0/(dx**2 + dy**2 + dz**2)
          radiusi = dsqrt(dd)
          e0 = HMX*radiusi
          f0 = HMX*radiusi*dd
          !write(6,*) "hm",e0
          ep = ep + e0
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi4 = fxi4 + forcex
          fyi4 = fyi4 + forcey
          fzi4 = fzi4 + forcez
          fxj1 = fxj1 + forcex
          fyj1 = fyj1 + forcey
          fzj1 = fzj1 + forcez
  
          !H2-O
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
          e0 = HOX*radiusi
          f0 = HOX*radiusi*dd
          !write(6,*) "ho",e0
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
  
          !H2-H
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
          e0 = HHX*radiusi
          f0 = HHX*radiusi*dd
          !write(6,*) "hh",e0
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
          dhx=dhx*ep
          dhy=dhy*ep
          dhz=dhz*ep
          fxi5=fxi5+dhx
          fyi5=fyi5+dhy
          fzi5=fzi5+dhz
          fxj4=fxj4+dhx
          fyj4=fyj4+dhy
          fzj4=fzj4+dhz
          ep1(k) = ep * hh
          !混合物の場合もこれで良いのか？？？
          !平成14年9月13日(金)圧力が正しくない原因ではないか
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
          !write(6,*) i,j,k,ep
          !     非vector計算機の場合には、この場で力の集計をした方がてっとり
          !     早い。
#ifdef VPOPTIMIZE
          !     ......計算結果をlist vectorに格納する。
          fxi(k*MEOHSITE+0)=fxi1
          fxi(k*MEOHSITE+1)=fxi2
          fxi(k*MEOHSITE+2)=fxi3
          fxi(k*MEOHSITE+3)=fxi4
          fxi(k*MEOHSITE+4)=fxi5
          fyi(k*MEOHSITE+0)=fyi1
          fyi(k*MEOHSITE+1)=fyi2
          fyi(k*MEOHSITE+2)=fyi3
          fyi(k*MEOHSITE+3)=fyi4
          fyi(k*MEOHSITE+4)=fyi5
          fzi(k*MEOHSITE+0)=fzi1
          fzi(k*MEOHSITE+1)=fzi2
          fzi(k*MEOHSITE+2)=fzi3
          fzi(k*MEOHSITE+3)=fzi4
          fzi(k*MEOHSITE+4)=fzi5
  
          fxj(k*MEOHSITE+0)=-fxj1
          fxj(k*MEOHSITE+1)=-fxj2
          fxj(k*MEOHSITE+2)=-fxj3
          fxj(k*MEOHSITE+3)=-fxj4
          fyj(k*MEOHSITE+0)=-fyj1
          fyj(k*MEOHSITE+1)=-fyj2
          fyj(k*MEOHSITE+2)=-fyj3
          fyj(k*MEOHSITE+3)=-fyj4
          fzj(k*MEOHSITE+0)=-fzj1
          fzj(k*MEOHSITE+1)=-fzj2
          fzj(k*MEOHSITE+2)=-fzj3
          fzj(k*MEOHSITE+3)=-fzj4
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
  
          site%fy(ii+1)=site%fy(ii+1)+fyi1
          site%fy(ii+2)=site%fy(ii+2)+fyi2
          site%fy(ii+3)=site%fy(ii+3)+fyi3
          site%fy(ii+4)=site%fy(ii+4)+fyi4
          site%fy(ii+5)=site%fy(ii+5)+fyi5
          site%fy(jj+1)=site%fy(jj+1)-fyj1
          site%fy(jj+2)=site%fy(jj+2)-fyj2
          site%fy(jj+3)=site%fy(jj+3)-fyj3
          site%fy(jj+4)=site%fy(jj+4)-fyj4
  
          site%fz(ii+1)=site%fz(ii+1)+fzi1
          site%fz(ii+2)=site%fz(ii+2)+fzi2
          site%fz(ii+3)=site%fz(ii+3)+fzi3
          site%fz(ii+4)=site%fz(ii+4)+fzi4
          site%fz(ii+5)=site%fz(ii+5)+fzi5
          site%fz(jj+1)=site%fz(jj+1)-fzj1
          site%fz(jj+2)=site%fz(jj+2)-fzj2
          site%fz(jj+3)=site%fz(jj+3)-fzj3
          site%fz(jj+4)=site%fz(jj+4)-fzj4
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
       do ii=1,wa%nsite*wa%nmol
          i=wa%lvmol(ii)
          s=wa%lvsite(ii)
          k=iv%partner_i(i,l)
          site%fx(ii+wa%offset)=site%fx(ii+wa%offset)+fxi(k*TIP4PSITE+s-1)
          site%fy(ii+wa%offset)=site%fy(ii+wa%offset)+fyi(k*TIP4PSITE+s-1)
          site%fz(ii+wa%offset)=site%fz(ii+wa%offset)+fzi(k*TIP4PSITE+s-1)
       enddo
    enddo
    do l=1,iv%maxpartner_j
       do jj=1,me%nsite*me%nmol
          j=me%lvmol(jj)
          s=me%lvsite(jj)
          k=iv%partner_j(j,l)
          site%fx(jj+me%offset)=site%fx(jj+me%offset)+fxj(k*MEOHSITE+s-1)
          site%fy(jj+me%offset)=site%fy(jj+me%offset)+fyj(k*MEOHSITE+s-1)
          site%fz(jj+me%offset)=site%fz(jj+me%offset)+fzj(k*MEOHSITE+s-1)
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
  end subroutine interaction_force_tip4p_oplsmeoh

end module oplsmeoh_tip4p_module

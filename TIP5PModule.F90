! -*- f90 -*-

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

module tip5p_module
  use common_module
  use physconst_module
  use interaction_module
  implicit none
  !TIP5P
  real(kind=8), parameter :: AA5P = ( 4d0*0.6694d0*3.12d0**12 * 100d0)
  real(kind=8), parameter :: BB5P = (-4d0*0.6694d0*3.12d0** 6 * 100d0)
  real(kind=8), parameter :: CC5P = 0.2410d0
  !comparison test with TIP4P
  !real(kind=8), parameter :: AA5P = ( 600000d0*CA * 100d0)
  !real(kind=8), parameter :: BB5P = (-610d0*CA * 100d0)
  !real(kind=8), parameter :: CC5P = 0.52d0
  real(kind=8), parameter :: MM5P = ( CC5P*CC5P*COEFF * 100d0)
  real(kind=8), parameter :: HM5P = (-CC5P*CC5P*COEFF * 100d0)
  real(kind=8), parameter :: HH5P = ( CC5P*CC5P*COEFF * 100d0)
  !重心を含めて6点
  integer, parameter :: TIP5PSITE = 6
  character(len=8),parameter :: tip5pName(TIP5PSITE)=(/"O", "H", "H", " ", " ", " "/)

contains
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine tip5p_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    call si_allocate(si, TIP5PSITE, LJ_COULOMB)
    !site order: OHHMM&CoM
    !charges [Q]
    si%param(1,1) = 0d0
    si%param(1,2) = CC5P
    si%param(1,3) = CC5P
    si%param(1,4) =-CC5P
    si%param(1,5) =-CC5P
    si%param(1,6) = 0d0
    !LJ-eps [kJ/mol]
    si%param(2,1) = 0.6694d0
    si%param(2,2) = 0d0
    si%param(2,3) = 0d0
    si%param(2,4) = 0d0
    si%param(2,5) = 0d0
    si%param(2,6) = 0d0
    !LJ-sig [AA]
    si%param(3,1) = 3.12d0
    si%param(3,2) = 0d0
    si%param(3,3) = 0d0
    si%param(3,4) = 0d0
    si%param(3,5) = 0d0
    si%param(3,6) = 0d0
  end subroutine tip5p_setinteraction

#ifdef EWALD
  subroutine interaction_force_tip5p(ewald,iv,mi,mj,site,epsum,vrsum)
#else
  subroutine interaction_force_tip5p(iv,mi,mj,site,epsum,vrsum)
#endif
    use mol_module
    use site_module
#ifdef EWALD
    use ewald_module
    type(sEwald),intent(in) :: ewald
#endif
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    real(kind=8),intent(out) :: epsum,vrsum
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
#ifdef VERBOSE
    integer :: count
#endif
    !integer :: kk
#ifdef VPOPTIMIZE
    !0番目の要素は、loop内のifをへらすためのトリック
#ifdef TEST10
    real(kind=8),dimension(0:MAXpair,TIP5PSITE) :: fxiq,fyiq,fziq,fxjq,fyjq,fzjq
#else /*TEST10*/
    real(kind=8),dimension(:),allocatable :: fxi,fyi,fzi,fxj,fyj,fzj
#ifdef TEST11
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm1,fym1,fzm1
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm2,fym2,fzm2
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm3,fym3,fzm3
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm4,fym4,fzm4
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm5,fym5,fzm5
    real(kind=8),dimension(MAXmol,MAXmol) :: fxm6,fym6,fzm6
#endif
#endif /*TEST10*/
    real(kind=8) :: r_r0,r_r1,r_r0s,r_r1s,dr
    allocate(fxi(0:TIP5PSITE*(iv%npair+1)))
    allocate(fyi(0:TIP5PSITE*(iv%npair+1)))
    allocate(fzi(0:TIP5PSITE*(iv%npair+1)))
    allocate(fxj(0:TIP5PSITE*(iv%npair+1)))
    allocate(fyj(0:TIP5PSITE*(iv%npair+1)))
    allocate(fzj(0:TIP5PSITE*(iv%npair+1)))
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
    vrsum=0d0
    !kk=0
  !VPP用のベクトル化制御行
  !OCL VECTOR,REPEAT(1000000)
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       dx00 = site%x(ii+TIP5PSITE)
       dy00 = site%y(ii+TIP5PSITE)
       dz00 = site%z(ii+TIP5PSITE)
       dx00 = dx00 - site%x(jj+TIP5PSITE)
       dy00 = dy00 - site%y(jj+TIP5PSITE)
       dz00 = dz00 - site%z(jj+TIP5PSITE)
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
#ifdef TEST12
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+2)-site%x(jj+2)+cx,site%y(ii+2)-site%y(jj+2)+cy,site%z(ii+2)-site%z(jj+2)+cz,HH5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+2)-site%x(jj+2)+cx,site%y(ii+2)-site%y(jj+2)+cy,site%z(ii+2)-site%z(jj+2)+cz,HH5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi2 = fxi2 + forcex
          fyi2 = fyi2 + forcey
          fzi2 = fzi2 + forcez
          fxj2 = fxj2 + forcex
          fyj2 = fyj2 + forcey
          fzj2 = fzj2 + forcez
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+2)-site%x(jj+3)+cx,site%y(ii+2)-site%y(jj+3)+cy,site%z(ii+2)-site%z(jj+3)+cz,HH5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+2)-site%x(jj+3)+cx,site%y(ii+2)-site%y(jj+3)+cy,site%z(ii+2)-site%z(jj+3)+cz,HH5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi2 = fxi2 + forcex
          fyi2 = fyi2 + forcey
          fzi2 = fzi2 + forcez
          fxj3 = fxj3 + forcex
          fyj3 = fyj3 + forcey
          fzj3 = fzj3 + forcez
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+3)-site%x(jj+2)+cx,site%y(ii+3)-site%y(jj+2)+cy,site%z(ii+3)-site%z(jj+2)+cz,HH5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+3)-site%x(jj+2)+cx,site%y(ii+3)-site%y(jj+2)+cy,site%z(ii+3)-site%z(jj+2)+cz,HH5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi3 = fxi3 + forcex
          fyi3 = fyi3 + forcey
          fzi3 = fzi3 + forcez
          fxj2 = fxj2 + forcex
          fyj2 = fyj2 + forcey
          fzj2 = fzj2 + forcez
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+3)-site%x(jj+3)+cx,site%y(ii+3)-site%y(jj+3)+cy,site%z(ii+3)-site%z(jj+3)+cz,HH5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+3)-site%x(jj+3)+cx,site%y(ii+3)-site%y(jj+3)+cy,site%z(ii+3)-site%z(jj+3)+cz,HH5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi3 = fxi3 + forcex
          fyi3 = fyi3 + forcey
          fzi3 = fzi3 + forcez
          fxj3 = fxj3 + forcex
          fyj3 = fyj3 + forcey
          fzj3 = fzj3 + forcez
  
  
  
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+4)-site%x(jj+4)+cx,site%y(ii+4)-site%y(jj+4)+cy,site%z(ii+4)-site%z(jj+4)+cz,MM5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+4)-site%x(jj+4)+cx,site%y(ii+4)-site%y(jj+4)+cy,site%z(ii+4)-site%z(jj+4)+cz,MM5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi4 = fxi4 + forcex
          fyi4 = fyi4 + forcey
          fzi4 = fzi4 + forcez
          fxj4 = fxj4 + forcex
          fyj4 = fyj4 + forcey
          fzj4 = fzj4 + forcez
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+4)-site%x(jj+5)+cx,site%y(ii+4)-site%y(jj+5)+cy,site%z(ii+4)-site%z(jj+5)+cz,MM5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+4)-site%x(jj+5)+cx,site%y(ii&
               & +4)-site%y(jj+5)+cy,site%z(ii+4)-site%z(jj+5)+cz&
               & ,MM5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi4 = fxi4 + forcex
          fyi4 = fyi4 + forcey
          fzi4 = fzi4 + forcez
          fxj5 = fxj5 + forcex
          fyj5 = fyj5 + forcey
          fzj5 = fzj5 + forcez
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+5)-site%x(jj+4)+cx,site%y(ii+5)-site%y(jj+4)+cy,site%z(ii+5)-site%z(jj+4)+cz,MM5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+5)-site%x(jj+4)+cx,site%y(ii+5)-site%y(jj+4)+cy,site%z(ii+5)-site%z(jj+4)+cz,MM5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi5 = fxi5 + forcex
          fyi5 = fyi5 + forcey
          fzi5 = fzi5 + forcez
          fxj4 = fxj4 + forcex
          fyj4 = fyj4 + forcey
          fzj4 = fzj4 + forcez
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+5)-site%x(jj+5)+cx,site%y(ii+5)-site%y(jj+5)+cy,site%z(ii+5)-site%z(jj+5)+cz,MM5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+5)-site%x(jj+5)+cx,site%y(ii+5)-site%y(jj+5)+cy,site%z(ii+5)-site%z(jj+5)+cz,MM5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi5 = fxi5 + forcex
          fyi5 = fyi5 + forcey
          fzi5 = fzi5 + forcez
          fxj5 = fxj5 + forcex
          fyj5 = fyj5 + forcey
          fzj5 = fzj5 + forcez
  
  
  
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+2)-site%x(jj+4)+cx,site%y(ii+2)-site%y(jj+4)+cy,site%z(ii+2)-site%z(jj+4)+cz,HM5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+2)-site%x(jj+4)+cx,site%y(ii+2)-site%y(jj+4)+cy,site%z(ii+2)-site%z(jj+4)+cz,HM5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi2 = fxi2 + forcex
          fyi2 = fyi2 + forcey
          fzi2 = fzi2 + forcez
          fxj4 = fxj4 + forcex
          fyj4 = fyj4 + forcey
          fzj4 = fzj4 + forcez
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+2)-site%x(jj+5)+cx,site%y(ii+2)-site%y(jj+5)+cy,site%z(ii+2)-site%z(jj+5)+cz,HM5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+2)-site%x(jj+5)+cx,site%y(ii+2)-site%y(jj+5)+cy,site%z(ii+2)-site%z(jj+5)+cz,HM5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi2 = fxi2 + forcex
          fyi2 = fyi2 + forcey
          fzi2 = fzi2 + forcez
          fxj5 = fxj5 + forcex
          fyj5 = fyj5 + forcey
          fzj5 = fzj5 + forcez
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+3)-site%x(jj+4)+cx,site%y(ii+3)-site%y(jj+4)+cy,site%z(ii+3)-site%z(jj+4)+cz,HM5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+3)-site%x(jj+4)+cx,site%y(ii+3)-site%y(jj+4)+cy,site%z(ii+3)-site%z(jj+4)+cz,HM5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi3 = fxi3 + forcex
          fyi3 = fyi3 + forcey
          fzi3 = fzi3 + forcez
          fxj4 = fxj4 + forcex
          fyj4 = fyj4 + forcey
          fzj4 = fzj4 + forcez
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+3)-site%x(jj+5)+cx,site%y(ii+3)-site%y(jj+5)+cy,site%z(ii+3)-site%z(jj+5)+cz,HM5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+3)-site%x(jj+5)+cx,site%y(ii+3)-site%y(jj+5)+cy,site%z(ii+3)-site%z(jj+5)+cz,HM5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi3 = fxi3 + forcex
          fyi3 = fyi3 + forcey
          fzi3 = fzi3 + forcez
          fxj5 = fxj5 + forcex
          fyj5 = fyj5 + forcey
          fzj5 = fzj5 + forcez
  
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+4)-site%x(jj+2)+cx,site%y(ii+4)-site%y(jj+2)+cy,site%z(ii+4)-site%z(jj+2)+cz,HM5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+4)-site%x(jj+2)+cx,site%y(ii+4)-site%y(jj+2)+cy,site%z(ii+4)-site%z(jj+2)+cz,HM5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi4 = fxi4 + forcex
          fyi4 = fyi4 + forcey
          fzi4 = fzi4 + forcez
          fxj2 = fxj2 + forcex
          fyj2 = fyj2 + forcey
          fzj2 = fzj2 + forcez
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+4)-site%x(jj+3)+cx,site%y(ii+4)-site%y(jj+3)+cy,site%z(ii+4)-site%z(jj+3)+cz,HM5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+4)-site%x(jj+3)+cx,site%y(ii+4)-site%y(jj+3)+cy,site%z(ii+4)-site%z(jj+3)+cz,HM5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi4 = fxi4 + forcex
          fyi4 = fyi4 + forcey
          fzi4 = fzi4 + forcez
          fxj3 = fxj3 + forcex
          fyj3 = fyj3 + forcey
          fzj3 = fzj3 + forcez
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+5)-site%x(jj+2)+cx,site%y(ii+5)-site%y(jj+2)+cy,site%z(ii+5)-site%z(jj+2)+cz,HM5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+5)-site%x(jj+2)+cx,site%y(ii+5)-site%y(jj+2)+cy,site%z(ii+5)-site%z(jj+2)+cz,HM5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi5 = fxi5 + forcex
          fyi5 = fyi5 + forcey
          fzi5 = fzi5 + forcez
          fxj2 = fxj2 + forcex
          fyj2 = fyj2 + forcey
          fzj2 = fzj2 + forcez
#ifdef EWALD
          call Coulomb_Force(ewald,site%x(ii+5)-site%x(jj+3)+cx,site%y(ii+5)-site%y(jj+3)+cy,site%z(ii+5)-site%z(jj+3)+cz,HM5P,e0,forcex,forcey,forcez)
#else
          call Coulomb_Force(site%x(ii+5)-site%x(jj+3)+cx,site%y(ii+5)-site%y(jj+3)+cy,site%z(ii+5)-site%z(jj+3)+cz,HM5P,e0,forcex,forcey,forcez)
#endif
          ep=ep+e0
          fxi5 = fxi5 + forcex
          fyi5 = fyi5 + forcey
          fzi5 = fzi5 + forcez
          fxj3 = fxj3 + forcex
          fyj3 = fyj3 + forcey
          fzj3 = fzj3 + forcez
          call LJ_Force(site%x(ii+1)-site%x(jj+1)+cx,site%y(ii+1)-site%y(jj+1)+cy,site%z(ii+1)-site%z(jj+1)+cz,AA5P,BB5P,e0,forcex,forcey,forcez)
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
          e0 = HH5P*radiusi
          f0 = HH5P*radiusi*dd
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
          e0 = HH5P*radiusi
          f0 = HH5P*radiusi*dd
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
          e0 = HH5P*radiusi
          f0 = HH5P*radiusi*dd
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
          e0 = HH5P*radiusi
          f0 = HH5P*radiusi*dd
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
          e0 = MM5P*radiusi
          f0 = MM5P*radiusi*dd
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
          e0 = MM5P*radiusi
          f0 = MM5P*radiusi*dd
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
          e0 = MM5P*radiusi
          f0 = MM5P*radiusi*dd
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
          e0 = MM5P*radiusi
          f0 = MM5P*radiusi*dd
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
          e0 = HM5P*radiusi
          f0 = HM5P*radiusi*dd
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
          e0 = HM5P*radiusi
          f0 = HM5P*radiusi*dd
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
          e0 = HM5P*radiusi
          f0 = HM5P*radiusi*dd
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
          e0 = HM5P*radiusi
          f0 = HM5P*radiusi*dd
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
          e0 = HM5P*radiusi
          f0 = HM5P*radiusi*dd
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
          e0 = HM5P*radiusi
          f0 = HM5P*radiusi*dd
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
          e0 = HM5P*radiusi
          f0 = HM5P*radiusi*dd
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
          e0 = HM5P*radiusi
          f0 = HM5P*radiusi*dd
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
          e0 = AA5P*d12 + BB5P*d6
          f0 = dd*(AA5P*12d0*d12 + BB5P*6d0*d6)
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
          fxi6=fxi6*hh
          fyi6=fyi6*hh
          fzi6=fzi6*hh
  
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
          dhx=dhx*ep
          dhy=dhy*ep
          dhz=dhz*ep
          fxi6=fxi6+dhx
          fyi6=fyi6+dhy
          fzi6=fzi6+dhz
          fxj6=fxj6+dhx
          fyj6=fyj6+dhy
          fzj6=fzj6+dhz
          epsum=epsum+ep*hh
          vrx = dx00*(fxi1+fxi2+fxi3+fxi4+fxi5+fxi6)
          vry = dy00*(fyi1+fyi2+fyi3+fyi4+fyi5+fyi6)
          vrz = dz00*(fzi1+fzi2+fzi3+fzi4+fzi5+fzi6)
          vrsum=vrsum+vrx+vry+vrz
          !write(6,*) i,j,k,ep
          !     非vector計算機の場合には、この場で力の集計をした方がてっとり
          !     早い。
#ifdef VPOPTIMIZE
          !     ......計算結果をlist vectorに格納する。
          fxi(k*TIP5PSITE+0)=fxi1
          fxi(k*TIP5PSITE+1)=fxi2
          fxi(k*TIP5PSITE+2)=fxi3
          fxi(k*TIP5PSITE+3)=fxi4
          fxi(k*TIP5PSITE+4)=fxi5
          fxi(k*TIP5PSITE+5)=fxi6
          fyi(k*TIP5PSITE+0)=fyi1
          fyi(k*TIP5PSITE+1)=fyi2
          fyi(k*TIP5PSITE+2)=fyi3
          fyi(k*TIP5PSITE+3)=fyi4
          fyi(k*TIP5PSITE+4)=fyi5
          fyi(k*TIP5PSITE+5)=fyi6
          fzi(k*TIP5PSITE+0)=fzi1
          fzi(k*TIP5PSITE+1)=fzi2
          fzi(k*TIP5PSITE+2)=fzi3
          fzi(k*TIP5PSITE+3)=fzi4
          fzi(k*TIP5PSITE+4)=fzi5
          fzi(k*TIP5PSITE+5)=fzi6
  
          fxj(k*TIP5PSITE+0)=-fxj1
          fxj(k*TIP5PSITE+1)=-fxj2
          fxj(k*TIP5PSITE+2)=-fxj3
          fxj(k*TIP5PSITE+3)=-fxj4
          fxj(k*TIP5PSITE+4)=-fxj5
          fxj(k*TIP5PSITE+5)=-fxj6
          fyj(k*TIP5PSITE+0)=-fyj1
          fyj(k*TIP5PSITE+1)=-fyj2
          fyj(k*TIP5PSITE+2)=-fyj3
          fyj(k*TIP5PSITE+3)=-fyj4
          fyj(k*TIP5PSITE+4)=-fyj5
          fyj(k*TIP5PSITE+5)=-fyj6
          fzj(k*TIP5PSITE+0)=-fzj1
          fzj(k*TIP5PSITE+1)=-fzj2
          fzj(k*TIP5PSITE+2)=-fzj3
          fzj(k*TIP5PSITE+3)=-fzj4
          fzj(k*TIP5PSITE+4)=-fzj5
          fzj(k*TIP5PSITE+5)=-fzj6
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
    fxm6(:,:)=0d0
    fym6(:,:)=0d0
    fzm6(:,:)=0d0
  !OCL NOVREC
    do k=1,iv%npair
       i=iv%pair_i(k)
       j=iv%pair_j(k)
       fxm1(j,i)=fxj(k*TIP5PSITE)
       fym1(j,i)=fyj(k*TIP5PSITE)
       fzm1(j,i)=fzj(k*TIP5PSITE)
  
       fxm2(j,i)=fxj(k*TIP5PSITE+1)
       fym2(j,i)=fyj(k*TIP5PSITE+1)
       fzm2(j,i)=fzj(k*TIP5PSITE+1)
  
       fxm3(j,i)=fxj(k*TIP5PSITE+2)
       fym3(j,i)=fyj(k*TIP5PSITE+2)
       fzm3(j,i)=fzj(k*TIP5PSITE+2)
  
       fxm4(j,i)=fxj(k*TIP5PSITE+3)
       fym4(j,i)=fyj(k*TIP5PSITE+3)
       fzm4(j,i)=fzj(k*TIP5PSITE+3)
  
       fxm5(j,i)=fxj(k*TIP5PSITE+4)
       fym5(j,i)=fyj(k*TIP5PSITE+4)
       fzm5(j,i)=fzj(k*TIP5PSITE+4)
  
       fxm6(j,i)=fxj(k*TIP5PSITE+5)
       fym6(j,i)=fyj(k*TIP5PSITE+5)
       fzm6(j,i)=fzj(k*TIP5PSITE+5)
    enddo
    !異種グループ間相互作用の場合は、一旦iに加わる力を集計する。
    if(.not.iv%isomol)then
       do i=1,mi%nmol
  !OCL NOVREC
          do j=1,mj%nmol
             site%fx((j-1)*TIP5PSITE+mj%offset+1)=site%fx((j-1)*TIP5PSITE+mj%offset+1)+fxm1(j,i)
             site%fx((j-1)*TIP5PSITE+mj%offset+2)=site%fx((j-1)*TIP5PSITE+mj%offset+2)+fxm2(j,i)
             site%fx((j-1)*TIP5PSITE+mj%offset+3)=site%fx((j-1)*TIP5PSITE+mj%offset+3)+fxm3(j,i)
             site%fx((j-1)*TIP5PSITE+mj%offset+4)=site%fx((j-1)*TIP5PSITE+mj%offset+4)+fxm4(j,i)
             site%fx((j-1)*TIP5PSITE+mj%offset+5)=site%fx((j-1)*TIP5PSITE+mj%offset+5)+fxm5(j,i)
             site%fx((j-1)*TIP5PSITE+mj%offset+6)=site%fx((j-1)*TIP5PSITE+mj%offset+6)+fxm6(j,i)
             site%fy((j-1)*TIP5PSITE+mj%offset+1)=site%fy((j-1)*TIP5PSITE+mj%offset+1)+fym1(j,i)
             site%fy((j-1)*TIP5PSITE+mj%offset+2)=site%fy((j-1)*TIP5PSITE+mj%offset+2)+fym2(j,i)
             site%fy((j-1)*TIP5PSITE+mj%offset+3)=site%fy((j-1)*TIP5PSITE+mj%offset+3)+fym3(j,i)
             site%fy((j-1)*TIP5PSITE+mj%offset+4)=site%fy((j-1)*TIP5PSITE+mj%offset+4)+fym4(j,i)
             site%fy((j-1)*TIP5PSITE+mj%offset+5)=site%fy((j-1)*TIP5PSITE+mj%offset+5)+fym5(j,i)
             site%fy((j-1)*TIP5PSITE+mj%offset+6)=site%fy((j-1)*TIP5PSITE+mj%offset+6)+fym6(j,i)
             site%fz((j-1)*TIP5PSITE+mj%offset+1)=site%fz((j-1)*TIP5PSITE+mj%offset+1)+fzm1(j,i)
             site%fz((j-1)*TIP5PSITE+mj%offset+2)=site%fz((j-1)*TIP5PSITE+mj%offset+2)+fzm2(j,i)
             site%fz((j-1)*TIP5PSITE+mj%offset+3)=site%fz((j-1)*TIP5PSITE+mj%offset+3)+fzm3(j,i)
             site%fz((j-1)*TIP5PSITE+mj%offset+4)=site%fz((j-1)*TIP5PSITE+mj%offset+4)+fzm4(j,i)
             site%fz((j-1)*TIP5PSITE+mj%offset+5)=site%fz((j-1)*TIP5PSITE+mj%offset+5)+fzm5(j,i)
             site%fz((j-1)*TIP5PSITE+mj%offset+6)=site%fz((j-1)*TIP5PSITE+mj%offset+6)+fzm6(j,i)
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
       fxm6(:,:)=0d0
       fym6(:,:)=0d0
       fzm6(:,:)=0d0
    endif
  !OCL NOVREC
    do k=1,iv%npair
       i=iv%pair_i(k)
       j=iv%pair_j(k)
       fxm1(i,j)=fxi(k*TIP5PSITE)
       fym1(i,j)=fyi(k*TIP5PSITE)
       fzm1(i,j)=fzi(k*TIP5PSITE)
  
       fxm2(i,j)=fxi(k*TIP5PSITE+1)
       fym2(i,j)=fyi(k*TIP5PSITE+1)
       fzm2(i,j)=fzi(k*TIP5PSITE+1)
  
       fxm3(i,j)=fxi(k*TIP5PSITE+2)
       fym3(i,j)=fyi(k*TIP5PSITE+2)
       fzm3(i,j)=fzi(k*TIP5PSITE+2)
  
       fxm4(i,j)=fxi(k*TIP5PSITE+3)
       fym4(i,j)=fyi(k*TIP5PSITE+3)
       fzm4(i,j)=fzi(k*TIP5PSITE+3)
  
       fxm5(i,j)=fxi(k*TIP5PSITE+4)
       fym5(i,j)=fyi(k*TIP5PSITE+4)
       fzm5(i,j)=fzi(k*TIP5PSITE+4)
  
       fxm6(i,j)=fxi(k*TIP5PSITE+5)
       fym6(i,j)=fyi(k*TIP5PSITE+5)
       fzm6(i,j)=fzi(k*TIP5PSITE+5)
    enddo
    do j=1,mj%nmol
  !OCL NOVREC
       do i=1,mi%nmol
          site%fx((i-1)*TIP5PSITE+mi%offset+1)=site%fx((i-1)*TIP5PSITE+mi%offset+1)+fxm1(i,j)
          site%fx((i-1)*TIP5PSITE+mi%offset+2)=site%fx((i-1)*TIP5PSITE+mi%offset+2)+fxm2(i,j)
          site%fx((i-1)*TIP5PSITE+mi%offset+3)=site%fx((i-1)*TIP5PSITE+mi%offset+3)+fxm3(i,j)
          site%fx((i-1)*TIP5PSITE+mi%offset+4)=site%fx((i-1)*TIP5PSITE+mi%offset+4)+fxm4(i,j)
          site%fx((i-1)*TIP5PSITE+mi%offset+5)=site%fx((i-1)*TIP5PSITE+mi%offset+5)+fxm5(i,j)
          site%fx((i-1)*TIP5PSITE+mi%offset+6)=site%fx((i-1)*TIP5PSITE+mi%offset+6)+fxm6(i,j)
  
          site%fy((i-1)*TIP5PSITE+mi%offset+1)=site%fy((i-1)*TIP5PSITE+mi%offset+1)+fym1(i,j)
          site%fy((i-1)*TIP5PSITE+mi%offset+2)=site%fy((i-1)*TIP5PSITE+mi%offset+2)+fym2(i,j)
          site%fy((i-1)*TIP5PSITE+mi%offset+3)=site%fy((i-1)*TIP5PSITE+mi%offset+3)+fym3(i,j)
          site%fy((i-1)*TIP5PSITE+mi%offset+4)=site%fy((i-1)*TIP5PSITE+mi%offset+4)+fym4(i,j)
          site%fy((i-1)*TIP5PSITE+mi%offset+5)=site%fy((i-1)*TIP5PSITE+mi%offset+5)+fym5(i,j)
          site%fy((i-1)*TIP5PSITE+mi%offset+6)=site%fy((i-1)*TIP5PSITE+mi%offset+6)+fym6(i,j)
  
          site%fz((i-1)*TIP5PSITE+mi%offset+1)=site%fz((i-1)*TIP5PSITE+mi%offset+1)+fzm1(i,j)
          site%fz((i-1)*TIP5PSITE+mi%offset+2)=site%fz((i-1)*TIP5PSITE+mi%offset+2)+fzm2(i,j)
          site%fz((i-1)*TIP5PSITE+mi%offset+3)=site%fz((i-1)*TIP5PSITE+mi%offset+3)+fzm3(i,j)
          site%fz((i-1)*TIP5PSITE+mi%offset+4)=site%fz((i-1)*TIP5PSITE+mi%offset+4)+fzm4(i,j)
          site%fz((i-1)*TIP5PSITE+mi%offset+5)=site%fz((i-1)*TIP5PSITE+mi%offset+5)+fzm5(i,j)
          site%fz((i-1)*TIP5PSITE+mi%offset+6)=site%fz((i-1)*TIP5PSITE+mi%offset+6)+fzm6(i,j)
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
          site%fx(ii+mi%offset)=site%fx(ii+mi%offset)+fxi(k*TIP5PSITE+s-1)
          site%fy(ii+mi%offset)=site%fy(ii+mi%offset)+fyi(k*TIP5PSITE+s-1)
          site%fz(ii+mi%offset)=site%fz(ii+mi%offset)+fzi(k*TIP5PSITE+s-1)
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
          site%fx(jj+mj%offset)=site%fx(jj+mj%offset)+fxj(k*TIP5PSITE+s-1)
          site%fy(jj+mj%offset)=site%fy(jj+mj%offset)+fyj(k*TIP5PSITE+s-1)
          site%fz(jj+mj%offset)=site%fz(jj+mj%offset)+fzj(k*TIP5PSITE+s-1)
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
    vrsum=+vrsum/3d0
    return
  end subroutine interaction_force_tip5p
end module tip5p_module

! -*- f90 -*-
!
!一時配列を使用する。すこし速くなるが、やはり配列のアクセスが遅い。
!
#undef TEST20
!
!aggressiveに座標配列を展開する。遅い。
!
#undef TEST21
!
!2重ループの順序を変えた。これは効いた。
!
#define TEST22
!
!xi,yi,ziに座標を一旦入れる。Genesis5での処理を真似ただけ。効果なし。
!
#undef TEST23
!
!force2で、水に特化して配列を完全に展開する。ここまでやってやっと少し
!速くなった。(53sec->41sec by e1/Mixture/input512)
!ifcは変化僅か(25.1sec->23.4sec)この程度であればあえて展開する必要は感
!じない。
!undefで固定平成15年4月10日(木)
#undef TEST24



!
!Standard Interaction: combination of LJ and Coulomb
!
module standard_interaction_module
  implicit none
  integer, parameter :: NONE=0, LJ_COULOMB=1
  type sStdInt
     integer :: numberOfSites
     integer :: mode
     integer :: numberOfParams
     real(kind=8),pointer :: param(:,:)
  end type sStdInt

contains
  subroutine si_allocate(si, nsite, mode)
    type(sStdInt),intent(inout) :: si
    integer,   intent(IN)       :: nsite
    integer, intent(IN)  :: mode
    si%numberOfSites  = nsite
    si%mode           = mode
    si%numberOfParams = 0
    if ( mode == LJ_COULOMB )then
       si%numberOfParams = 3
    endif
    allocate(si%param(si%numberOfParams,nsite))
  end subroutine si_allocate

  subroutine si_deallocate(si)
    type(sStdInt),intent(inout) :: si
    si%numberOfSites  = 0
    si%mode           = 0
    si%numberOfParams = 0
    deallocate(si%param)
  end subroutine si_deallocate
  
  subroutine StdForce(iv,mi,mj,si,sj,site,ep,vrmat)
    use physconst_module
    use interaction_module
    use mol_module
    use site_module
    implicit none
    type(sInteraction),intent(in),   target :: iv
    type(sMol)        ,intent(in)    :: mi,mj
    type(sStdInt)     ,intent(in)    :: si,sj
    type(sSite)       ,intent(inout),target :: site
    real(kind=8),intent(out) :: vrmat(3,3)
    real(kind=8),intent(inout) :: ep(*)
    real(kind=8) :: dx,dy,dz
    real(kind=8) :: dhx,dhy,dhz
    real(kind=8) :: d6,d12,e0,f0,radiusi
    real(kind=8) forcex,forcey,forcez
    real(kind=8) dd,hh,dh
    real(kind=8) :: vr(3), d00(3), dx00, dy00, dz00
    real(kind=8) :: qq,aa,bb,ljsig,ljeps
    integer :: i,j,k,l
    integer :: ii,jj
    integer :: isite,jsite
    !integer :: kk
    real(kind=8),dimension(:,:),allocatable :: fxi,fyi,fzi,fxj,fyj,fzj
#ifdef TEST20
    !
    !もっと積極的に、一旦一時配列に構造体要素をコピーして使用する。
    !
    real(kind=8), allocatable :: six(:),siy(:),siz(:)
    real(kind=8), allocatable :: sjx(:),sjy(:),sjz(:)
#endif
#ifdef TEST21
    !
    !もっと積極的に、一旦一時配列に構造体要素をコピーして使用する。
    !
    real(kind=8), allocatable :: six(:,:),siy(:,:),siz(:,:)
    real(kind=8), allocatable :: sjx(:,:),sjy(:,:),sjz(:,:)
    real(kind=8), allocatable :: dxv(:),dyv(:),dzv(:)
    !固定配列にしたところで速くはならない。
    !real(kind=8)              :: dxv(MAXPAIR)
#endif
#ifdef TEST23
    real(kind=8) :: xi,yi,zi,xj,yj,zj
#endif
#undef DEBUG0
#ifdef DEBUG0
    real(kind=8) :: esum
#endif
    !
    !構造体の値を一時変数に入れておくと明らかに速くなる場合がある。
    !というより、構造体の処理が遅すぎる。
    !
    !平成15年1月16日(木)これでも速くならない。
    !
    integer :: mi_nsite,mj_nsite
    integer :: mi_offset,mj_offset
    real(kind=8), pointer :: site_x(:),iv_ox(:)
    !real(kind=8)          :: iv_ox(MAXPAIR)
    !write(*,*) "maxpair",MAXPAIR
    mi_nsite = mi%nsite
    mj_nsite = mj%nsite
    mi_offset = mi%offset
    mj_offset = mj%offset
    site_x => site%x
    iv_ox  => iv%ox
    !iv_ox(1:iv%npair) = iv%ox(1:iv%npair)
    !write(STDERR,*) site%x(1),site%y(1),site%z(1)
    !write(STDERR,*) site%x(2),site%y(2),site%z(2)
#ifdef TEST20
    allocate(six(mi_nsite*mi%nmol))
    allocate(siy(mi_nsite*mi%nmol))
    allocate(siz(mi_nsite*mi%nmol))
    allocate(sjx(mj_nsite*mj%nmol))
    allocate(sjy(mj_nsite*mj%nmol))
    allocate(sjz(mj_nsite*mj%nmol))
    do i=1,mi_nsite*mi%nmol
       six(i) = site%x(i + mi_offset)
       siy(i) = site%y(i + mi_offset)
       siz(i) = site%z(i + mi_offset)
    enddo
    do i=1,mj_nsite*mj%nmol
       sjx(i) = site%x(i + mj_offset)
       sjy(i) = site%y(i + mj_offset)
       sjz(i) = site%z(i + mj_offset)
    enddo
#endif
#ifdef TEST21
    allocate(six(iv%npair,si%numberOfSites))
    allocate(siy(iv%npair,si%numberOfSites))
    allocate(siz(iv%npair,si%numberOfSites))
    allocate(sjx(iv%npair,sj%numberOfSites))
    allocate(sjy(iv%npair,sj%numberOfSites))
    allocate(sjz(iv%npair,sj%numberOfSites))
    allocate(dxv(iv%npair))
    allocate(dyv(iv%npair))
    allocate(dzv(iv%npair))
    do isite=1,si%numberOfSites
       do k=1,iv%npair
          i = iv%pair_i(k)
          ii = (i-1)*mi_nsite+mi_offset + isite
          six(k,isite) = site_x(ii)
          siy(k,isite) = site%y(ii)
          siz(k,isite) = site%z(ii)
       enddo
    enddo
    do jsite=1,sj%numberOfSites
       do k=1,iv%npair
          j = iv%pair_j(k)
          jj = (j-1)*mj_nsite+mj_offset + jsite
          sjx(k,jsite) = site_x(jj)
          sjy(k,jsite) = site%y(jj)
          sjz(k,jsite) = site%z(jj)
       enddo
    enddo
#endif
    !
    !fxiの引数の順番はこれでよい。平成15年1月15日(水)
    !
    allocate(fxi(iv%npair,si%numberOfSites))
    allocate(fyi(iv%npair,si%numberOfSites))
    allocate(fzi(iv%npair,si%numberOfSites))
    allocate(fxj(iv%npair,sj%numberOfSites))
    allocate(fyj(iv%npair,sj%numberOfSites))
    allocate(fzj(iv%npair,sj%numberOfSites))
    !
    !初期化に時間がかかっている。
    !
    fxi(:,:) = 0d0
    fyi(:,:) = 0d0
    fzi(:,:) = 0d0
    fxj(:,:) = 0d0
    fyj(:,:) = 0d0
    fzj(:,:) = 0d0
    ep(1:iv%npair)    = 0d0
    !kk=0
    !VPP用のベクトル化制御行
    !OCL VECTOR,REPEAT(1000000)
    !write(STDERR,*) si%numberOfSites, sj%numberOfSites
    do isite = 1, si%numberOfSites
       do jsite = 1, sj%numberOfSites
          !
          !Product of charges
          !
          qq = si%param(1,isite) * sj%param(1,jsite) * COEFF * KX2X * J2I * iv%scale
          !
          !Product of eps
          !
          aa = si%param(2,isite) * sj%param(2,jsite)
          !write(STDERR,*) isite,jsite,qq,aa
          if ( qq /= 0d0 ) then
             if ( aa /= 0d0 ) then
                !
                !LJ + Coulomb
                !
                ljeps = sqrt( aa )
                ljsig = ( si%param(3,isite) + sj%param(3,jsite) ) * 0.5d0
                ljsig = ljsig**6
                aa    = 4d0 * ljeps * ljsig * ljsig * KX2X * J2I * iv%scale
                bb    =-4d0 * ljeps * ljsig         * KX2X * J2I * iv%scale
#ifdef TEST21
                !dxv(1:iv%npair) = six(1:iv%npair,isite)
                dxv(:) = six(:,isite)
                dyv(:) = siy(:,isite)
                dzv(:) = siz(:,isite)
                !dxv(1:iv%npair) = dxv(1:iv%npair) - sjx(1:iv%npair,jsite)
                dxv(:) = dxv(:) - sjx(:,jsite)
                dyv(:) = dyv(:) - sjy(:,jsite)
                dzv(:) = dzv(:) - sjz(:,jsite)
                !dxv(1:iv%npair) = dxv(1:iv%npair) - iv_ox(1:iv%npair)
                dxv(:) = dxv(:) - iv%ox(:)
                dyv(:) = dyv(:) - iv%oy(:)
                dzv(:) = dzv(:) - iv%oz(:)
#endif
                do k=1,iv%Npair
                   i = iv%pair_i(k)
                   j = iv%pair_j(k)
                   !     write(STDERR,*) i,j
#ifdef TEST20
                   ii = (i-1)*mi_nsite+isite
                   jj = (j-1)*mj_nsite+jsite
                   dx = six(ii) - sjx(jj)
                   dy = siy(ii) - sjy(jj)
                   dz = siz(ii) - sjz(jj)
#elif defined (TEST21)
                   !dx = six(k,isite) - sjx(k,jsite)
                   !dy = siy(k,isite) - sjy(k,jsite)
                   !dz = siz(k,isite) - sjz(k,jsite)
#else
                   ii = (i-1)*mi_nsite+mi_offset + isite
                   jj = (j-1)*mj_nsite+mj_offset + jsite
                   dx = site_x(ii) - site_x(jj)
                   dy = site%y(ii) - site%y(jj)
                   dz = site%z(ii) - site%z(jj)
#endif
#ifdef TEST21
                   dd = 1d0/(dxv(k)**2 + dyv(k)**2 + dzv(k)**2)
#else
                   dx = dx - iv_ox(k)
                   dy = dy - iv%oy(k)
                   dz = dz - iv%oz(k)
                   dd = 1d0/(dx**2 + dy**2 + dz**2)
#endif
                   radiusi = dsqrt(dd)
                   e0 = qq*radiusi
                   f0 = qq*radiusi*dd
                   d6 = dd*dd*dd
                   d12= d6*d6
                   e0 = e0 + aa*d12 + bb*d6
                   f0 = f0 + dd*(aa*12d0*d12 + bb*6d0*d6)
#ifdef TEST21
                   forcex = dxv(k)*f0
                   forcey = dyv(k)*f0
                   forcez = dzv(k)*f0
#else
                   forcex = dx*f0
                   forcey = dy*f0
                   forcez = dz*f0
#endif
                   fxi(k,isite) = fxi(k,isite) + forcex
                   fyi(k,isite) = fyi(k,isite) + forcey
                   fzi(k,isite) = fzi(k,isite) + forcez
                   fxj(k,jsite) = fxj(k,jsite) + forcex
                   fyj(k,jsite) = fyj(k,jsite) + forcey
                   fzj(k,jsite) = fzj(k,jsite) + forcez
                   ep(k)        = ep(k) + e0
                enddo
             else
                !
                !Coulomb Only
                !
#ifdef TEST21
                !dxv(1:iv%npair) = six(1:iv%npair,isite)
                dxv(:) = six(:,isite)
                dyv(:) = siy(:,isite)
                dzv(:) = siz(:,isite)
                !dxv(1:iv%npair) = dxv(1:iv%npair) - sjx(1:iv%npair,jsite)
                dxv(:) = dxv(:) - sjx(:,jsite)
                dyv(:) = dyv(:) - sjy(:,jsite)
                dzv(:) = dzv(:) - sjz(:,jsite)
                !dxv(1:iv%npair) = dxv(1:iv%npair) - iv_ox(1:iv%npair)
                dxv(:) = dxv(:) - iv%ox(:)
                dyv(:) = dyv(:) - iv%oy(:)
                dzv(:) = dzv(:) - iv%oz(:)
#endif
#ifdef DEBUG0
                esum=0d0
#endif
                do k=1,iv%Npair
                   i = iv%pair_i(k)
                   j = iv%pair_j(k)
                   !write(STDERR,*) k,i,j
#ifdef TEST20
                   ii = (i-1)*mi_nsite+isite
                   jj = (j-1)*mj_nsite+jsite
                   dx = six(ii) - sjx(jj)
                   dy = siy(ii) - sjy(jj)
                   dz = siz(ii) - sjz(jj)
#elif defined (TEST21)
                   !dx = six(k,isite) - sjx(k,jsite)
                   !dy = siy(k,isite) - sjy(k,jsite)
                   !dz = siz(k,isite) - sjz(k,jsite)
#else
                   ii = (i-1)*mi_nsite+mi_offset + isite
                   jj = (j-1)*mj_nsite+mj_offset + jsite
                   !
                   !!!この、何でもない差を計算する部分が致命的に遅い。
                   !Genesis5でいうと、この部分はxi,xjの差を求める部分。
                   !おそらく、サイトの組みあわせのループを外にもっていっ
                   !たせいで、大きな配列site%x()が離散的にアクセスされ
                   !るようになり、キャッシュから溢れて遅くなっているの
                   !だと思う。巧妙な方法でループを入れかえることはでき
                   !るけど、それでもGenesis5の場合のように相互作用が
                   !hardcodingされている場合ほどの速度には至らないかも
                   !しれない。平成15年4月10日(木)一応、巧妙なバージョ
                   !ンを試作してみよう。
                   !
                   !
#ifdef TEST23
                   xi = site%x(ii)
                   yi = site%y(ii)
                   zi = site%z(ii)
                   xj = site%x(jj)
                   yj = site%y(jj)
                   zj = site%z(jj)
                   dx = xi - xj
                   dy = yi - yj
                   dz = zi - zj
#else /*TEST23*/
                   dx = site_x(ii) - site_x(jj)
                   dy = site%y(ii) - site%y(jj)
                   dz = site%z(ii) - site%z(jj)
#endif /*TEST23*/
#endif
#ifdef TEST21
                   dd = 1d0/(dxv(k)**2 + dyv(k)**2 + dzv(k)**2)
#else
                   dx = dx - iv_ox(k)
                   dy = dy - iv%oy(k)
                   dz = dz - iv%oz(k)
                   dd = 1d0/(dx**2 + dy**2 + dz**2)
#endif
                   radiusi = dsqrt(dd)
                   e0 = qq*radiusi
                   f0 = qq*radiusi*dd
#ifdef TEST21
                   forcex = dxv(k)*f0
                   forcey = dyv(k)*f0
                   forcez = dzv(k)*f0
#else
                   forcex = dx*f0
                   forcey = dy*f0
                   forcez = dz*f0
#endif
                   fxi(k,isite) = fxi(k,isite) + forcex
                   fyi(k,isite) = fyi(k,isite) + forcey
                   fzi(k,isite) = fzi(k,isite) + forcez
                   fxj(k,jsite) = fxj(k,jsite) + forcex
                   fyj(k,jsite) = fyj(k,jsite) + forcey
                   fzj(k,jsite) = fzj(k,jsite) + forcez
                   ep(k)        = ep(k) + e0
#ifdef DEBUG0
                   esum         = esum + e0
#endif
                enddo
#ifdef DEBUG0
                write(STDERR,*) isite,jsite,esum
#endif
             endif
          else
             if ( aa /= 0d0 ) then
                !
                !LJ only
                !
                ljeps = sqrt( aa )
                ljsig = ( si%param(3,isite) + sj%param(3,jsite) ) * 0.5d0
                ljsig = ljsig**6
                aa    = 4d0 * ljeps * ljsig * ljsig * KX2X * J2I * iv%scale
                bb    =-4d0 * ljeps * ljsig         * KX2X * J2I * iv%scale
#ifdef TEST21
                !dxv(1:iv%npair) = six(1:iv%npair,isite)
                dxv(:) = six(:,isite)
                dyv(:) = siy(:,isite)
                dzv(:) = siz(:,isite)
                !dxv(1:iv%npair) = dxv(1:iv%npair) - sjx(1:iv%npair,jsite)
                dxv(:) = dxv(:) - sjx(:,jsite)
                dyv(:) = dyv(:) - sjy(:,jsite)
                dzv(:) = dzv(:) - sjz(:,jsite)
                !dxv(1:iv%npair) = dxv(1:iv%npair) - iv_ox(1:iv%npair)
                dxv(:) = dxv(:) - iv%ox(:)
                dyv(:) = dyv(:) - iv%oy(:)
                dzv(:) = dzv(:) - iv%oz(:)
#endif
#ifdef DEBUG0
                esum=0d0
#endif
                do k=1,iv%Npair
                   i = iv%pair_i(k)
                   j = iv%pair_j(k)
                   !write(STDERR,*) k,i,j
#ifdef TEST20
                   ii = (i-1)*mi_nsite+isite
                   jj = (j-1)*mj_nsite+jsite
                   dx = six(ii) - sjx(jj)
                   dy = siy(ii) - sjy(jj)
                   dz = siz(ii) - sjz(jj)
#elif defined (TEST21)
                   !dx = six(k,isite) - sjx(k,jsite)
                   !dy = siy(k,isite) - sjy(k,jsite)
                   !dz = siz(k,isite) - sjz(k,jsite)
#else
                   ii = (i-1)*mi_nsite+mi_offset + isite
                   jj = (j-1)*mj_nsite+mj_offset + jsite
                   dx = site_x(ii) - site_x(jj)
                   dy = site%y(ii) - site%y(jj)
                   dz = site%z(ii) - site%z(jj)
#endif
#ifdef TEST21
                   dd = 1d0/(dxv(k)**2 + dyv(k)**2 + dzv(k)**2)
#else
                   dx = dx - iv_ox(k)
                   dy = dy - iv%oy(k)
                   dz = dz - iv%oz(k)
                   dd = 1d0/(dx**2 + dy**2 + dz**2)
#endif
                   d6 = dd*dd*dd
                   d12= d6*d6
                   e0 = aa*d12 + bb*d6
                   !write(6,*) k
                   !write(6,*) i,isite,ii,mi_offset,mi_nsite
                   !write(6,*) j,jsite,jj,mj_offset,mj_nsite
                   !write(6,*) e0
                   f0 = dd*(aa*12d0*d12 + bb*6d0*d6)
#ifdef TEST21
                   forcex = dxv(k)*f0
                   forcey = dyv(k)*f0
                   forcez = dzv(k)*f0
#else
                   forcex = dx*f0
                   forcey = dy*f0
                   forcez = dz*f0
#endif
                   fxi(k,isite) = fxi(k,isite) + forcex
                   fyi(k,isite) = fyi(k,isite) + forcey
                   fzi(k,isite) = fzi(k,isite) + forcez
                   fxj(k,jsite) = fxj(k,jsite) + forcex
                   fyj(k,jsite) = fyj(k,jsite) + forcey
                   fzj(k,jsite) = fzj(k,jsite) + forcez
                   ep(k)        = ep(k) + e0
#ifdef DEBUG0
                   esum         = esum + e0
#endif
                enddo
#ifdef DEBUG0
                write(STDERR,*) isite,jsite,esum
#endif
             endif
          endif
       enddo
    enddo
    !write( STDOUT,* ) (fxi(1,i), i=1,7)
    !write( STDOUT,* ) (fyi(1,i), i=1,7)
    !write( STDOUT,* ) (fzi(1,i), i=1,7)
    !write( STDOUT,* ) (fxj(1,i), i=1,7)
    !write( STDOUT,* ) (fyj(1,i), i=1,7)
    !write( STDOUT,* ) (fzj(1,i), i=1,7)
    !
    !apply smoothing function
    !
    isite = si%numberOfSites
    jsite = sj%numberOfSites
    do l=1,isite
       do k=1,iv%Npair
          !
          !hhはtruncation減衰係数
          !
          hh                  = iv%eratio( k )
          fxi(k,l) = fxi(k,l) * hh
          fyi(k,l) = fyi(k,l) * hh
          fzi(k,l) = fzi(k,l) * hh
       enddo
    enddo
    do l=1,jsite
       do k=1,iv%Npair
          !
          !hhはtruncation減衰係数
          !
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
       dx00 = site_x(ii+isite) - site_x(jj+jsite) - iv_ox(k)
       dy00 = site%y(ii+isite) - site%y(jj+jsite) - iv%oy(k)
       dz00 = site%z(ii+isite) - site%z(jj+jsite) - iv%oz(k)
       dh   = iv%fratio( k )
       dhx  = dh*dx00*ep(k)
       dhy  = dh*dy00*ep(k)
       dhz  = dh*dz00*ep(k)
       !
       !重心サイトに、カットオフによって派生する力を加える。
       !
       fxi(k,isite) = fxi(k,isite) + dhx
       fyi(k,isite) = fyi(k,isite) + dhy
       fzi(k,isite) = fzi(k,isite) + dhz
       fxj(k,jsite) = fxj(k,jsite) + dhx
       fyj(k,jsite) = fyj(k,jsite) + dhy
       fzj(k,jsite) = fzj(k,jsite) + dhz
#ifdef TEST22
    enddo
    do l=1,isite
       do k=1,iv%Npair
          i = iv%pair_i(k)
          j = iv%pair_j(k)
          ii = (i-1)*mi%nsite+mi%offset
          jj = (j-1)*mj%nsite+mj%offset
          site%fx(ii+l) = site%fx(ii+l) + fxi(k,l)
          site%fy(ii+l) = site%fy(ii+l) + fyi(k,l)
          site%fz(ii+l) = site%fz(ii+l) + fzi(k,l)
       enddo
    enddo
    do l=1,jsite
       do k=1,iv%Npair
          i = iv%pair_i(k)
          j = iv%pair_j(k)
          ii = (i-1)*mi%nsite+mi%offset
          jj = (j-1)*mj%nsite+mj%offset
          site%fx(jj+l) = site%fx(jj+l) - fxj(k,l)
          site%fy(jj+l) = site%fy(jj+l) - fyj(k,l)
          site%fz(jj+l) = site%fz(jj+l) - fzj(k,l)
       enddo
    enddo
#else
       do l=1,isite
          site%fx(ii+l) = site%fx(ii+l) + fxi(k,l)
          site%fy(ii+l) = site%fy(ii+l) + fyi(k,l)
          site%fz(ii+l) = site%fz(ii+l) + fzi(k,l)
       enddo
       do l=1,jsite
          site%fx(jj+l) = site%fx(jj+l) - fxj(k,l)
          site%fy(jj+l) = site%fy(jj+l) - fyj(k,l)
          site%fz(jj+l) = site%fz(jj+l) - fzj(k,l)
       enddo
    enddo
#endif
    do k=1,iv%Npair
       ep(k) = ep(k) * iv%eratio( k )
    enddo
    !do i=1,5
    !   ii = (i-1)*mi%nsite+mi%offset
    !   write(6,*) site%fx(ii+isite), site%fy(ii+isite), site%fz(ii+isite)
    !enddo

    vrmat(:,:) = 0d0
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       ii = (i-1)*mi%nsite+mi%offset
       jj = (j-1)*mj%nsite+mj%offset
       d00(1) = site_x(ii+isite) - site_x(jj+jsite) - iv_ox(k)
       d00(2) = site%y(ii+isite) - site%y(jj+jsite) - iv%oy(k)
       d00(3) = site%z(ii+isite) - site%z(jj+jsite) - iv%oz(k)
       vr(1) = 0d0
       vr(2) = 0d0
       vr(3) = 0d0
       do l=1,isite
          vr(1) = vr(1) + fxi(k,l)
          vr(2) = vr(2) + fyi(k,l)
          vr(3) = vr(3) + fzi(k,l)
       enddo
       do i=1,3
          do j=1,3
             vrmat(j,i) = vrmat(j,i) + d00(j) * vr(i)
          enddo
       enddo
    enddo
    deallocate(fxi)
    deallocate(fyi)
    deallocate(fzi)
    deallocate(fxj)
    deallocate(fyj)
    deallocate(fzj)
#ifdef TEST20
    deallocate(six)
    deallocate(siy)
    deallocate(siz)
    deallocate(sjx)
    deallocate(sjy)
    deallocate(sjz)
#endif
#ifdef TEST21
    deallocate(six)
    deallocate(siy)
    deallocate(siz)
    deallocate(sjx)
    deallocate(sjy)
    deallocate(sjz)
    !deallocate(dxv)
    deallocate(dyv)
    deallocate(dzv)
#endif
  end subroutine StdForce


  subroutine StdPairPotential(iv,mi,mj,si,sj,site,ep )
    use physconst_module
    use interaction_module
    use mol_module
    use site_module
    use montecarlo_module
    implicit none
    type(sInteraction),intent(in),   target :: iv
    type(sMol)        ,intent(in)    :: mi,mj
    type(sStdInt)     ,intent(in)    :: si,sj
    type(sSite)       ,intent(inout),target :: site
    real(kind=8),      intent(out)   :: ep(*)
    real(kind=8) :: dx,dy,dz
    real(kind=8) :: d6,d12,e0,radiusi
    real(kind=8) :: dd
    real(kind=8) :: qq,aa,bb,ljsig,ljeps
    integer :: i,j,k
    integer :: ii,jj
    integer :: isite,jsite
    !logical, intent(IN):: verbose
    ep (1:iv%npair)   = 0d0
    do isite = 1, si%numberOfSites
       do jsite = 1, sj%numberOfSites
          !
          !Product of charges
          !
          qq = si%param(1,isite) * sj%param(1,jsite) * COEFF * KX2X * J2I * iv%scale
          !
          !Product of eps
          !
          aa = si%param(2,isite) * sj%param(2,jsite)
          !write(STDERR,*) isite,jsite,qq,aa
          if ( qq /= 0d0 ) then
             if ( aa /= 0d0 ) then
                !
                !LJ + Coulomb
                !
                ljeps = sqrt( aa )
                ljsig = ( si%param(3,isite) + sj%param(3,jsite) ) * 0.5d0
                ljsig = ljsig**6
                aa    = 4d0 * ljeps * ljsig * ljsig * KX2X * J2I * iv%scale
                bb    =-4d0 * ljeps * ljsig         * KX2X * J2I * iv%scale
                do k=1,iv%Npair
                   i = iv%pair_i(k)
                   j = iv%pair_j(k)
                   !     write(STDERR,*) i,j
                   ii = (i-1)*mi%nsite+mi%offset + isite
                   jj = (j-1)*mj%nsite+mj%offset + jsite
                   dx = site%x(ii) - site%x(jj)
                   dy = site%y(ii) - site%y(jj)
                   dz = site%z(ii) - site%z(jj)
                   dx = dx - iv%ox(k)
                   dy = dy - iv%oy(k)
                   dz = dz - iv%oz(k)
                   dd = 1d0/(dx**2 + dy**2 + dz**2)
                   radiusi = dsqrt(dd)
                   e0 = qq*radiusi
                   d6 = dd*dd*dd
                   d12= d6*d6
                   e0 = e0 + aa*d12 + bb*d6
                   ep(k)        = ep(k) + e0
                   !if ( verbose ) write(6,*) "ep", k, ep(k), e0
                enddo
             else
                !
                !Coulomb Only
                !
                do k=1,iv%Npair
                   i = iv%pair_i(k)
                   j = iv%pair_j(k)
                   !     write(STDERR,*) i,j
                   ii = (i-1)*mi%nsite+mi%offset + isite
                   jj = (j-1)*mj%nsite+mj%offset + jsite
                   dx = site%x(ii) - site%x(jj)
                   dy = site%y(ii) - site%y(jj)
                   dz = site%z(ii) - site%z(jj)
                   dx = dx - iv%ox(k)
                   dy = dy - iv%oy(k)
                   dz = dz - iv%oz(k)
                   dd = 1d0/(dx**2 + dy**2 + dz**2)
                   radiusi = dsqrt(dd)
                   e0 = qq*radiusi
                   ep(k)        = ep(k) + e0
                   !if ( verbose ) write(6,*) "ep", k, ep(k), e0
                enddo
             endif
          else
             if ( aa /= 0d0 ) then
                !
                !LJ only
                !
                ljeps = sqrt( aa )
                ljsig = ( si%param(3,isite) + sj%param(3,jsite) ) * 0.5d0
                ljsig = ljsig**6
                aa    = 4d0 * ljeps * ljsig * ljsig * KX2X * J2I * iv%scale
                bb    =-4d0 * ljeps * ljsig         * KX2X * J2I * iv%scale
                do k=1,iv%Npair
                   i = iv%pair_i(k)
                   j = iv%pair_j(k)
                   !     write(STDERR,*) i,j
                   ii = (i-1)*mi%nsite+mi%offset + isite
                   jj = (j-1)*mj%nsite+mj%offset + jsite
                   dx = site%x(ii) - site%x(jj)
                   dy = site%y(ii) - site%y(jj)
                   dz = site%z(ii) - site%z(jj)
                   dx = dx - iv%ox(k)
                   dy = dy - iv%oy(k)
                   dz = dz - iv%oz(k)
                   dd = 1d0/(dx**2 + dy**2 + dz**2)
                   d6 = dd*dd*dd
                   d12= d6*d6
                   e0 = aa*d12 + bb*d6
                   ep(k)        = ep(k) + e0
                   !if ( verbose ) then
                   !   write(6,*) "ep", k,  site%y(ii),site%y(jj)
                   !endif
                enddo
             endif
          endif
       enddo
    enddo
    !
    !apply smoothing function
    !
    do k=1,iv%Npair
       ep(k) = ep(k) * iv%eratio( k )
       !if ( verbose ) write(6,*) "ep", k, ep(k), e0
    enddo
    !write(STDERR,*) "std",ep(17)
  end subroutine StdPairPotential


  subroutine StdOHLength(iv,mi,mj,si,sj,site,ep )
    use physconst_module
    use interaction_module
    use mol_module
    use site_module
    use montecarlo_module
    implicit none
    type(sInteraction),intent(in),   target :: iv
    type(sMol)        ,intent(in)    :: mi,mj
    type(sStdInt)     ,intent(in)    :: si,sj
    type(sSite)       ,intent(inout),target :: site
    real(kind=8),      intent(out)   :: ep(*)
    real(kind=8) :: dx,dy,dz
    real(kind=8) :: dd
    integer :: i,j,k
    integer :: ii,jj
    integer :: isite,jsite

    ep (1:iv%npair)   = 1000000000d0
    do isite = 1, si%numberOfSites
       if ( mi%name(isite)(1:1) .eq. "O" ) then
          do jsite = 1, sj%numberOfSites
             if ( mj%name(jsite)(1:1) .eq. "H" ) then
                do k=1,iv%Npair
                   i = iv%pair_i(k)
                   j = iv%pair_j(k)
                   !     write(STDERR,*) i,j
                   ii = (i-1)*mi%nsite+mi%offset + isite
                   jj = (j-1)*mj%nsite+mj%offset + jsite
                   dx = site%x(ii) - site%x(jj)
                   dy = site%y(ii) - site%y(jj)
                   dz = site%z(ii) - site%z(jj)
                   dx = dx - iv%ox(k)
                   dy = dy - iv%oy(k)
                   dz = dz - iv%oz(k)
                   dd = sqrt(dx**2 + dy**2 + dz**2)
                   if ( dd < ep(k) ) ep(k) = dd
                enddo
             endif
          enddo
       endif
    enddo
    do isite = 1, si%numberOfSites
       if ( mi%name(isite)(1:1) .eq. "H" ) then
          do jsite = 1, sj%numberOfSites
             if ( mj%name(jsite)(1:1) .eq. "O" ) then
                do k=1,iv%Npair
                   i = iv%pair_i(k)
                   j = iv%pair_j(k)
                   !     write(STDERR,*) i,j
                   ii = (i-1)*mi%nsite+mi%offset + isite
                   jj = (j-1)*mj%nsite+mj%offset + jsite
                   dx = site%x(ii) - site%x(jj)
                   dy = site%y(ii) - site%y(jj)
                   dz = site%z(ii) - site%z(jj)
                   dx = dx - iv%ox(k)
                   dy = dy - iv%oy(k)
                   dz = dz - iv%oz(k)
                   dd = sqrt(dx**2 + dy**2 + dz**2)
                   if ( dd < ep(k) ) ep(k) = dd
                enddo
             endif
          enddo
       endif
    enddo
  end subroutine StdOHLength



#ifdef SOLVATIONTEST
  subroutine force_sol_one( mi,mj,site,joint,i,j,ox,oy,oz,ep)
    use site_module
    use mol_module
    use cutoff_module
    use external_module
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    type(sJointHB), intent(IN)  :: joint
    real(kind=8) :: fxi2,fyi2,fzi2
    real(kind=8) :: fxi3,fyi3,fzi3
    real(kind=8) :: fxj1,fyj1,fzj1
    real(kind=8) :: dx,dy,dz,dx00,dy00,dz00
    real(kind=8) :: dhx,dhy,dhz
    real(kind=8) :: d6,d12,e0,f0,radiusi
    real(kind=8) :: xi,yi,zi,xj,yj,zj
    real(kind=8) forcex,forcey,forcez
    real(kind=8), intent(IN)  :: ox,oy,oz
    real(kind=8) cx,cy,cz
    real(kind=8), intent(INOUT) :: ep
    real(kind=8) dd,d
    integer :: i,j,k
    integer :: ii,jj
    logical :: cond1, cond2
    ii = (i-1)*mi%nsite+mi%offset
    jj = (j-1)*mj%nsite+mj%offset
    !dx00 = site%x(ii+mi%nsite)
    !dy00 = site%y(ii+mi%nsite)
    !dz00 = site%z(ii+mi%nsite)
    !dx00 = dx00 - site%x(jj+mj%nsite)
    !dy00 = dy00 - site%y(jj+mj%nsite)
    !dz00 = dz00 - site%z(jj+mj%nsite)
    !dx00 = dx00-ox
    !dy00 = dy00-oy
    !dz00 = dz00-oz
    !write(STDERR,*) i,j,dx00,dy00,dz00
    cx = -ox
    cy = -oy
    cz = -oz
    fxi2 = 0d0
    fyi2 = 0d0
    fzi2 = 0d0
    fxi3 = 0d0
    fyi3 = 0d0
    fzi3 = 0d0
    
    fxj1 = 0d0
    fyj1 = 0d0
    fzj1 = 0d0
    
    !
    !まず、水素1に注目する。
    !結合相手が指定されている場合にのみ引力や斥力が働く。
    !
    ep=0d0
    if ( joint%partner( 1, i ) .ne. 0 ) then
       xi = site%x(ii+2)
       yi = site%y(ii+2)
       zi = site%z(ii+2)
       xj = site%x(jj+1)
       yj = site%y(jj+1)
       zj = site%z(jj+1)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = (dx**2 + dy**2 + dz**2)
       !引力が働く条件
       ! 1. jが結合相手として指定されていて
       ! 2. 距離が遠い
       cond1 = ( joint%partner( 1, i ) .eq. j ) .and. ( joint%balance**2 < dd )
       !斥力が働く条件
       ! 1. j以外の結合相手が指定されていて
       ! 2. 距離が近い
       cond2 = ( joint%partner( 1, i ) .ne. j ) .and. ( dd < joint%balance**2 )
       if ( cond1 .or. cond2 ) then
          d  = sqrt(dd)
          ep = ep + joint%a * (d - joint%balance)**joint%b
          f0 = -joint%a * joint%b * ( d - joint%balance )**( joint%b - 1 ) / d
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi2 = fxi2 + forcex
          fyi2 = fyi2 + forcey
          fzi2 = fzi2 + forcez
          fxj1 = fxj1 + forcex
          fyj1 = fyj1 + forcey
          fzj1 = fzj1 + forcez
          !write(STDERR,*) i,j,"(1)", cond1
       endif
    endif
    if ( joint%partner( 2, i ) .ne. 0 ) then
       xi = site%x(ii+3)
       yi = site%y(ii+3)
       zi = site%z(ii+3)
       xj = site%x(jj+1)
       yj = site%y(jj+1)
       zj = site%z(jj+1)
       dx = xi-xj+cx
       dy = yi-yj+cy
       dz = zi-zj+cz
       dd = (dx**2 + dy**2 + dz**2)
       !引力が働く条件
       ! 1. jが結合相手として指定されていて
       ! 2. 距離が遠い
       cond1 = ( joint%partner( 2, i ) .eq. j ) .and. ( joint%balance**2 < dd )
       !斥力が働く条件
       ! 1. j以外の結合相手が指定されていて
       ! 2. 距離が近い
       cond2 = ( joint%partner( 2, i ) .ne. j ) .and. ( dd < joint%balance**2 )
       if ( cond1 .or. cond2 ) then
          d  = sqrt(dd)
          ep = ep + joint%a * (d - joint%balance)**joint%b
          f0 = -joint%a * joint%b * ( d - joint%balance )**( joint%b - 1 ) / d
          forcex = dx*f0
          forcey = dy*f0
          forcez = dz*f0
          fxi3 = fxi3 + forcex
          fyi3 = fyi3 + forcey
          fzi3 = fzi3 + forcez
          fxj1 = fxj1 + forcex
          fyj1 = fyj1 + forcey
          fzj1 = fzj1 + forcez
          !write(STDERR,*) i,j,"(2)", cond1
       endif
    endif
    site%fx(ii+2)=site%fx(ii+2)+fxi2
    site%fx(ii+3)=site%fx(ii+3)+fxi3
    
    site%fy(ii+2)=site%fy(ii+2)+fyi2
    site%fy(ii+3)=site%fy(ii+3)+fyi3
    
    site%fz(ii+2)=site%fz(ii+2)+fzi2
    site%fz(ii+3)=site%fz(ii+3)+fzi3
    
    site%fx(jj+1)=site%fx(jj+1)-fxj1
    site%fy(jj+1)=site%fy(jj+1)-fyj1
    site%fz(jj+1)=site%fz(jj+1)-fzj1
  end subroutine force_sol_one
       


  subroutine force_sol(iv,mi,mj,site,joint,ep1)
    use site_module
    use mol_module
    use cutoff_module
    use external_module
    use interaction_module
    type(sInteraction),intent(in) :: iv
    type(sMol),intent(in) :: mi,mj
    type(sSite),intent(inout) :: site
    real(kind=8), intent(INOUT) :: ep1(*)
    type(sJointHB), intent(IN)  :: joint
    real(kind=8) ox,oy,oz
    integer :: i,j,k
    ep1(1:iv%npair) = 0d0
    !kk=0
    !VPP用のベクトル化制御行
    !OCL VECTOR,REPEAT(1000000)
    do k=1,iv%Npair
       !平成16年5月19日(水)write文を入れるとAddress Errorを回避できる。なぜ？
       ox = iv%ox(k)
       oy = iv%oy(k)
       oz = iv%oz(k)
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       !
       !引力を加算するのは、結合i-->jが指定されている場合。
       !
       if ( joint%members( i ) .ne. 0 ) then
       !斥力を加算するのは、iがメンバーで、引力が働かない場合。
       !要するに、iがメンバーでなければ、全く考慮の必要なし。
       !write(STDERR,*) i,j,k,iv%Npair,mi%offset,mj%offset
          call force_sol_one(mi,mj,site,joint,i,j,ox,oy,oz,ep1(k))
       endif
       !
       !対表にはi-->j片方しか載っていないので反対向きも計算。
       !(!!!本当はmi == mjの場合のみ)
       !
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       ox = -iv%ox(k)
       oy = -iv%oy(k)
       oz = -iv%oz(k)
       if ( joint%members( j ) .ne. 0 ) then
          call force_sol_one(mj,mi,site,joint,j,i,ox,oy,oz,ep1(k))
       endif
!       write(STDERR,*)
    enddo
  end subroutine force_sol
#endif



end module standard_interaction_module

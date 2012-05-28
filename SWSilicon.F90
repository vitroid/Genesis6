! -*- f90 -*-

!
!Stillinger-Weber's silicon model archive #00055
!

!
! 2004-08-31 Lattice energy is checked.
! 2004-08-31 Energy conservation is checked.
!

module swsilicon_module
  use physconst_module
  use interaction_module
  implicit none
  real(kind=8), parameter :: SWSILICON_AA = 7.044556277d0
  real(kind=8), parameter :: SWSILICON_BB = 0.6022245584d0
  real(kind=8), parameter :: SWSILICON_A =  1.80d0
  real(kind=8), parameter :: SWSILICON_GAMMA =  1.20d0
  real(kind=8), parameter :: SWSILICON_LAMBDA = 21.0d0
  real(kind=8), parameter,private :: EPSILON = 50d3 * CA * J2I
  real(kind=8), parameter,private :: SIGMA = 2.0951d0
  real(kind=8), parameter :: SWSILICON_A_REAL = SWSILICON_A * SIGMA
  !重心を含めて2点
  integer, parameter :: SWSILICONSITE = 2
  !
  character(len=8),parameter :: swsiliconName(SWSILICONSITE) = (/"Si", "  "/)
  !
  !相互作用行列。将来はStdInteractionに組みこむ。
  !
  integer, parameter :: INTR_TYPE_NONE=0, INTR_TYPE_LJCOULOMB=1, INTR_TYPE_LJ=2, INTR_TYPE_COULOMB=3, INTR_TYPE_SJ=4, INTR_TYPE_SW=5
  type intr
     integer :: type
     real(kind=8) :: param( 10 )
  end type intr
  type(intr),private :: intrMatrix(SWSILICONSITE, SWSILICONSITE)
  character(len=8), parameter :: swsilicon_id08 = "SWSILICO"

contains

  subroutine SWSilicon_ScaleLambda( value )
    real(kind=8), intent(IN) :: value
    intrMatrix(1,1)%param(5) = value * SWSILICON_LAMBDA
  end subroutine SWSilicon_ScaleLambda

  subroutine Swsilicon_GetMass( mass )
    implicit none
    real(kind=8) :: mass(SWSILICONSITE)
    mass(1) = 28.09d0
    mass(2) = 0d0
  end subroutine Swsilicon_GetMass
  !
  !StdInteractionを使う場合の初期化
  !
  subroutine swsilicon_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    integer :: i,j
    call si_allocate(si, SWSILICONSITE, LJ_COULOMB)
    intrMatrix(1:1, 1:1)%type = INTR_TYPE_SW
    intrMatrix(1:1, 2)%type   = INTR_TYPE_NONE
    intrMatrix(2, 1:1)%type   = INTR_TYPE_NONE
    intrMatrix(2, 2)%type     = INTR_TYPE_NONE

    intrMatrix(1:1, 1:1)%param(1) = swsilicon_aa * EPSILON
    intrMatrix(1:1, 1:1)%param(2) = swsilicon_bb
    intrMatrix(1:1, 1:1)%param(3) = swsilicon_a
    !Three-body interaction parameters; should not be here.
    intrMatrix(1:1, 1:1)%param(4) = swsilicon_gamma
    intrMatrix(1:1, 1:1)%param(5) = swsilicon_lambda
  end subroutine swsilicon_setinteraction

  subroutine force_swsilicon(iv,mi,mj,site,ep,vrmat)
    use physconst_module
    use interaction_module
    use mol_module
    use site_module
    implicit none
    type(sInteraction),intent(in),   target :: iv
    type(sMol)        ,intent(in)    :: mi,mj
    type(sSite)       ,intent(inout),target :: site
    real(kind=8),intent(out) :: vrmat(3,3)
    real(kind=8),intent(inout) :: ep(*)
    real(kind=8) :: dx,dy,dz
    real(kind=8) :: dhx,dhy,dhz
    real(kind=8) :: d6,d12,e0,f0,radiusi, radius
    real(kind=8) forcex,forcey,forcez
    real(kind=8) dd,hh,dh
    real(kind=8) :: vr(3), d00(3), dx00, dy00, dz00
    real(kind=8) :: aa,bb,r0
    integer :: i,j,k,l
    integer :: ii,jj
    integer :: isite,jsite
    !integer :: kk
    real(kind=8),dimension(:,:),allocatable :: fxi,fyi,fzi,fxj,fyj,fzj
#ifdef TEST23
    real(kind=8) :: xi,yi,zi,xj,yj,zj
#endif
#undef DEBUG0
#ifdef DEBUG0
    real(kind=8) :: esum
#endif
    integer :: mi_nsite,mj_nsite
    integer :: mi_offset,mj_offset
    !real(kind=8), pointer :: site_x(:),iv_ox(:)
    !real(kind=8)          :: iv_ox(MAXPAIR)
    !write(*,*) "maxpair",MAXPAIR
    mi_nsite = mi%nsite
    mj_nsite = mj%nsite
    mi_offset = mi%offset
    mj_offset = mj%offset
    !site_x => site%x
    !iv_ox  => iv%ox
    !
    !fxiの引数の順番はこれでよい。平成15年1月15日(水)
    !
    !write(STDERR,*) iv%npair,SWSILICONSITE
    allocate(fxi(iv%npair,SWSILICONSITE))
    allocate(fyi(iv%npair,SWSILICONSITE))
    allocate(fzi(iv%npair,SWSILICONSITE))
    allocate(fxj(iv%npair,SWSILICONSITE))
    allocate(fyj(iv%npair,SWSILICONSITE))
    allocate(fzj(iv%npair,SWSILICONSITE))
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
    do isite = 1, 1
       do jsite = 1, 1
          !if ( intrMatrix(isite, jsite)%type .eq. INTR_TYPE_SW ) then
          aa = intrMatrix(isite, jsite)%param(1)
          bb = intrMatrix(isite, jsite)%param(2)
          r0 = intrMatrix(isite, jsite)%param(3)

          !ljeps = sqrt( aa )
          !ljsig = ( si%param(3,isite) + sj%param(3,jsite) ) * 0.5d0
          !ljsig = ljsig**6
          !aa    = 4d0 * ljeps * ljsig * ljsig * 100d0 * iv%scale
          !bb    =-4d0 * ljeps * ljsig         * 100d0 * iv%scale
          do k=1,iv%Npair
             i = iv%pair_i(k)
             j = iv%pair_j(k)
             !     write(STDERR,*) i,j
             ii = (i-1)*mi_nsite+mi_offset + isite
             jj = (j-1)*mj_nsite+mj_offset + jsite
             dx = site%x(ii) - site%x(jj)
             dy = site%y(ii) - site%y(jj)
             dz = site%z(ii) - site%z(jj)
             dx = dx - iv%ox(k)
             dy = dy - iv%oy(k)
             dz = dz - iv%oz(k)
             dx = dx / SIGMA
             dy = dy / SIGMA
             dz = dz / SIGMA
             !
             !SW-Silicon pair potential has singular point at cutoff.
             !Be careful.
             !
             dd = (dx**2 + dy**2 + dz**2)
             radius = dsqrt(dd)
             radiusi = 1d0 / radius
             dd = 1d0 / dd
             d6 = dd*dd*dd
             !if ( radius .lt. 0.75d0 ) then
             !   write(STDERR,*) i,j,radius
             !endif
             !d12= d6*d6
             !e0 = e0 + aa*d12 + bb*d6
             !f0 = f0 + dd*(aa*12d0*d12 + bb*6d0*d6)
             e0 = aa * ( bb * radius**(-4) - 1d0 ) * exp( 1d0/( radius - r0 ) )
             !write(STDERR,*) "PAIR", radius, e0 / EPSILON
             f0 = 1d0 / SIGMA * aa * exp( 1d0/(radius - r0) ) * ( 4d0*bb*radius**(-6) +(BB*radius**(-4) - 1d0)/( radius * ( radius - r0 )**2 ) )
             !write(STDERR,*) "PAIR", radius, e0 / EPSILON, f0 / EPSILON * SIGMA**2
             !write(6,*) radius, e0,f0
             forcex = dx*f0
             forcey = dy*f0
             forcez = dz*f0
             fxi(k,isite) = fxi(k,isite) + forcex
             fyi(k,isite) = fyi(k,isite) + forcey
             fzi(k,isite) = fzi(k,isite) + forcez
             fxj(k,jsite) = fxj(k,jsite) + forcex
             fyj(k,jsite) = fyj(k,jsite) + forcey
             fzj(k,jsite) = fzj(k,jsite) + forcez
             ep(k)        = ep(k) + e0
          enddo
       enddo
    enddo
    !
    !apply smoothing function
    !
    isite = SWSILICONSITE
    jsite = SWSILICONSITE
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
       dx00 = site%x(ii+isite) - site%x(jj+jsite) - iv%ox(k)
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
       d00(1) = site%x(ii+isite) - site%x(jj+jsite) - iv%ox(k)
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
  end subroutine force_swsilicon

  subroutine force_swsilicon_triplet(iv,triplet,mi,site,ep,vrmat)
    use triplet_module
    use physconst_module
    use interaction_module
    use mol_module
    use site_module
    implicit none
    type(sInteraction),intent(in),   target :: iv
    type(sTriplet),intent(in),       target :: triplet
    type(sMol)        ,intent(in)           :: mi
    type(sSite)       ,intent(inout),target :: site
    real(kind=8),intent(out)                :: vrmat(3,3)
    real(kind=8),intent(inout)              :: ep(*)
    real(kind=8) :: fxi( triplet%num )
    real(kind=8) :: fyi( triplet%num )
    real(kind=8) :: fzi( triplet%num )
    real(kind=8) :: fxj( triplet%num )
    real(kind=8) :: fyj( triplet%num )
    real(kind=8) :: fzj( triplet%num )
    real(kind=8) :: fxk( triplet%num )
    real(kind=8) :: fyk( triplet%num )
    real(kind=8) :: fzk( triplet%num )

    integer :: tri

    integer :: ii,jj,kk
    integer :: center, left, right
    integer :: clpair, crpair
    real(kind=8) :: e0
    real(kind=8) :: xij, yij, zij
    real(kind=8) :: xik, yik, zik
    real(kind=8) :: aij, aik, a, b, c, d
    real(kind=8) :: rij, rrij
    real(kind=8) :: rik, rrik
    real(kind=8) :: dadxi, dadyi, dadzi
    real(kind=8) :: dbdxi, dbdyi, dbdzi
    real(kind=8) :: dcdxi, dcdyi, dcdzi
    real(kind=8) :: dddxi, dddyi, dddzi
    real(kind=8) :: dadxj, dadyj, dadzj
    real(kind=8) :: dbdxj, dbdyj, dbdzj
    real(kind=8) :: dcdxj, dcdyj, dcdzj
    real(kind=8) :: dddxj, dddyj, dddzj
    real(kind=8) :: dadxk, dadyk, dadzk
    real(kind=8) :: dbdxk, dbdyk, dbdzk
    real(kind=8) :: dcdxk, dcdyk, dcdzk
    real(kind=8) :: dddxk, dddyk, dddzk

    real(kind=8) :: ega, b32, b31, gamma, lambda, r0

    !Pressure calculation
    vrmat(:,:) = 0d0
    !Interaction parameters
    r0     = intrMatrix(1,1)%param(3)
    gamma  = intrMatrix(1,1)%param(4)
    lambda = intrMatrix(1,1)%param(5)
    do tri = 1, triplet%num
       !write(STDERR,*) 1,tri, triplet%center(1)
       center = triplet%center( tri )
       left   = triplet%left( tri )
       right  = triplet%right( tri )
       clpair = triplet%centerleft( tri )
       crpair = triplet%centerright( tri )
       ! center-left pair
       ii = ( center-1 )*mi%nsite+mi%offset + 1
       jj = ( left  -1 )*mi%nsite+mi%offset + 1
       xij = site%x( jj ) - site%x( ii )
       yij = site%y( jj ) - site%y( ii )
       zij = site%z( jj ) - site%z( ii )
       if ( 0 < clpair ) then
          xij = xij + iv%ox( clpair )
          yij = yij + iv%oy( clpair )
          zij = zij + iv%oz( clpair )
       else
          xij = xij - iv%ox( -clpair )
          yij = yij - iv%oy( -clpair )
          zij = zij - iv%oz( -clpair )
       endif
       ! center-right pair
       ii = ( center-1 )*mi%nsite+mi%offset + 1
       kk = ( right -1 )*mi%nsite+mi%offset + 1
       xik = site%x( kk ) - site%x( ii )
       yik = site%y( kk ) - site%y( ii )
       zik = site%z( kk ) - site%z( ii )
       if ( 0 < crpair ) then
          xik = xik + iv%ox( crpair )
          yik = yik + iv%oy( crpair )
          zik = zik + iv%oz( crpair )
       else
          xik = xik - iv%ox( -crpair )
          yik = yik - iv%oy( -crpair )
          zik = zik - iv%oz( -crpair )
       endif
       xij = xij / SIGMA
       yij = yij / SIGMA
       zij = zij / SIGMA
       xik = xik / SIGMA
       yik = yik / SIGMA
       zik = zik / SIGMA
       ! potential energy
       rrij = ( xij**2 + yij**2 + zij**2 )
       rrik = ( xik**2 + yik**2 + zik**2 )
       rij  = sqrt( rrij )
       rik  = sqrt( rrik )
       aij  = 1d0 / ( rij - r0 )
       aik  = 1d0 / ( rik - r0 )
       a    = aij + aik
       c    = ( xij*xik + yij*yik + zij*zik )
       d    = rij * rik
       b    = c / d
       ega  = exp( gamma * a )
       b31  = ( b + 1d0/3d0 )
       b32  = ( b + 1d0/3d0 ) ** 2
       !e0   = epsilon * lambda * exp( gamma * a ) * ( b + 1d0/3d0 ) ** 2
       e0   = epsilon * lambda * ega * b32
       ep(tri) = e0
       ! force. see note
       !write(STDERR,*) 2,tri, triplet%center(1)
       dadxj   = -xij / rij * aij**2
       dadyj   = -yij / rij * aij**2
       dadzj   = -zij / rij * aij**2
       dcdxj   = xik
       dcdyj   = yik
       dcdzj   = zik
       dddxj   = xij * rik / rij
       dddyj   = yij * rik / rij
       dddzj   = zij * rik / rij
       dbdxj   = dcdxj / d - c * dddxj / d**2
       dbdyj   = dcdyj / d - c * dddyj / d**2
       dbdzj   = dcdzj / d - c * dddzj / d**2
       !fxj(tri)= -epsilon * lambda / sigma * exp( gamma * a ) * ( gamma * dadxj * ( b + 1d0/3d0 )**2 + 2 * dbdxj * ( b + 1d0/3d0 ) )
       !fyj(tri)= -epsilon * lambda / sigma * exp( gamma * a ) * ( gamma * dadyj * ( b + 1d0/3d0 )**2 + 2 * dbdyj * ( b + 1d0/3d0 ) )
       !fzj(tri)= -epsilon * lambda / sigma * exp( gamma * a ) * ( gamma * dadzj * ( b + 1d0/3d0 )**2 + 2 * dbdzj * ( b + 1d0/3d0 ) )
       fxj(tri)= -epsilon * lambda / sigma * ega * ( gamma * dadxj * b32 + 2 * dbdxj * b31 )
       fyj(tri)= -epsilon * lambda / sigma * ega * ( gamma * dadyj * b32 + 2 * dbdyj * b31 )
       fzj(tri)= -epsilon * lambda / sigma * ega * ( gamma * dadzj * b32 + 2 * dbdzj * b31 )

       !write(STDERR,*) 3,tri, triplet%center(1)
       dadxk   = -xik / rik * aik**2
       dadyk   = -yik / rik * aik**2
       dadzk   = -zik / rik * aik**2
       dcdxk   = xij
       dcdyk   = yij
       dcdzk   = zij
       dddxk   = xik * rij / rik
       dddyk   = yik * rij / rik
       dddzk   = zik * rij / rik
       dbdxk   = dcdxk / d - c * dddxk / d**2
       dbdyk   = dcdyk / d - c * dddyk / d**2
       dbdzk   = dcdzk / d - c * dddzk / d**2
       fxk(tri)= -epsilon * lambda / sigma * ega * ( gamma * dadxk * b32 + 2 * dbdxk * b31 )
       fyk(tri)= -epsilon * lambda / sigma * ega * ( gamma * dadyk * b32 + 2 * dbdyk * b31 )
       fzk(tri)= -epsilon * lambda / sigma * ega * ( gamma * dadzk * b32 + 2 * dbdzk * b31 )

       !write(STDERR,*) 4,tri, triplet%center(1)
       dadxi   = -dadxj - dadxk
       dadyi   = -dadyj - dadyk
       dadzi   = -dadzj - dadzk
       dcdxi   = -xij - xik
       dcdyi   = -yij - yik
       dcdzi   = -zij - zik
       dddxi   = -xij * rik / rij - xik * rij / rik
       dddyi   = -yij * rik / rij - yik * rij / rik
       dddzi   = -zij * rik / rij - zik * rij / rik
       dbdxi   = dcdxi / d - c * dddxi / d**2
       dbdyi   = dcdyi / d - c * dddyi / d**2
       dbdzi   = dcdzi / d - c * dddzi / d**2
       fxi(tri)= -epsilon * lambda / sigma * ega * ( gamma * dadxi * b32 + 2 * dbdxi * b31 )
       fyi(tri)= -epsilon * lambda / sigma * ega * ( gamma * dadyi * b32 + 2 * dbdyi * b31 )
       fzi(tri)= -epsilon * lambda / sigma * ega * ( gamma * dadzi * b32 + 2 * dbdzi * b31 )

       !write(STDERR,*) 5,tri, triplet%center(1)
       vrmat(1,1) = vrmat(1,1) + xij * fxj( tri ) + xik * fxk( tri )
       vrmat(2,1) = vrmat(2,1) + yij * fxj( tri ) + yik * fxk( tri )
       vrmat(3,1) = vrmat(3,1) + zij * fxj( tri ) + zik * fxk( tri )
       vrmat(1,2) = vrmat(1,2) + xij * fyj( tri ) + xik * fyk( tri )
       vrmat(2,2) = vrmat(2,2) + yij * fyj( tri ) + yik * fyk( tri )
       vrmat(3,2) = vrmat(3,2) + zij * fyj( tri ) + zik * fyk( tri )
       vrmat(1,3) = vrmat(1,3) + xij * fzj( tri ) + xik * fzk( tri )
       vrmat(2,3) = vrmat(2,3) + yij * fzj( tri ) + yik * fzk( tri )
       vrmat(3,3) = vrmat(3,3) + zij * fzj( tri ) + zik * fzk( tri )
       !write(STDERR,*) 6,tri, triplet%center(1)

    enddo
    ! unscaling for xij etc.
    vrmat(:,:) = vrmat(:,:) * SIGMA

    do tri=1,triplet%num
       center = triplet%center( tri )
       left   = triplet%left( tri )
       right  = triplet%right( tri )
       !write(STDERR,*) tri,center, left, right
       ii = ( center-1 )*mi%nsite+mi%offset + 1
       jj = ( left  -1 )*mi%nsite+mi%offset + 1
       kk = ( right -1 )*mi%nsite+mi%offset + 1
       site%fx(ii) = site%fx(ii) + fxi(tri)
       site%fy(ii) = site%fy(ii) + fyi(tri)
       site%fz(ii) = site%fz(ii) + fzi(tri)
       site%fx(jj) = site%fx(jj) + fxj(tri)
       site%fy(jj) = site%fy(jj) + fyj(tri)
       site%fz(jj) = site%fz(jj) + fzj(tri)
       site%fx(kk) = site%fx(kk) + fxk(tri)
       site%fy(kk) = site%fy(kk) + fyk(tri)
       site%fz(kk) = site%fz(kk) + fzk(tri)
    enddo
    
  end subroutine force_swsilicon_triplet

end module swsilicon_module

! -*- f90 -*-

!See archive #00060
!まだなんにも変えてない。!
!StdInteractionには乗らないので、forceを書かなければならない。


module sjbenzene_module
  use common_module
  use physconst_module
  use interaction_module
  implicit none
  real(kind=8), parameter :: sjbenzene_acc = 78998d0  !Kcal/mol
  real(kind=8), parameter :: sjbenzene_bcc = 3.6d0    ! /A
  real(kind=8), parameter :: sjbenzene_ccc = 519.3d0  !Kcal A^6/mol
  real(kind=8), parameter :: sjbenzene_ahh = 2384.6d0
  real(kind=8), parameter :: sjbenzene_bhh = 3.74d0
  real(kind=8), parameter :: sjbenzene_chh = 24.624d0
  real(kind=8), parameter :: sjbenzene_ach = 3888.0d0
  real(kind=8), parameter :: sjbenzene_bch = 3.415d0
  real(kind=8), parameter :: sjbenzene_cch = 124.416d0
  real(kind=8), parameter :: sjbenzene_qc  = -0.11d0
  !real(kind=8), parameter :: sjbenzene_qc  = -0.00d0
  real(kind=8), parameter :: sjbenzene_cc  =  1.395d0 !A  value from archive#00062
  real(kind=8), parameter :: sjbenzene_ch  =  1.084d0 !A  value from archive#00062
  real(kind=8), parameter :: sjbenzene_hmass=1d0, sjbenzene_cmass=12d0
  !real(kind=8), parameter, private :: epsoh = sqrt( sjbenzene_epso * sjbenzene_epsh )
  !重心を含めて13点
  integer, parameter :: SJBENZENESITE = 13
  character(len=8),parameter :: sjbenzeneName(SJBENZENESITE)=(/"C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H", " "/)
  !
  !相互作用行列。将来はStdInteractionに組みこむ。
  !
  integer, parameter :: INTR_TYPE_NONE=0, INTR_TYPE_LJCOULOMB=1, INTR_TYPE_LJ=2, INTR_TYPE_COULOMB=3, INTR_TYPE_SJ=4
  type intr
     integer :: type
     real(kind=8) :: param( 10 )
  end type intr
  type(intr),private :: intrMatrix(SJBENZENESITE, SJBENZENESITE)
  character(len=8), parameter :: sjbenzene_id08 = "SJBENZEN"

contains

  subroutine Rigid_Sjbenzene_Constructor(r)
    use rigid_module
    type(sRigid),intent(INOUT) :: r
    real(kind=8) :: comx,comy,comz
    integer      :: i
    real(kind=8) :: theta
    real(kind=8) :: mass(SJBENZENESITE)
    real(kind=8) :: rh, rc
    rc = sjbenzene_cc
    rh = sjbenzene_cc + sjbenzene_ch
    ! 13 sites; CCCCCCHHHHHH and center of mass
    do i = 1 , 6
       theta = dble(i-1) * pi / 3d0
       r%molx( i ) = rc * cos( theta )
       r%moly( i ) = rc * sin( theta )
       r%molz( i ) = 0d0
       r%molx( i + 6 ) = rh * cos( theta )
       r%moly( i + 6 ) = rh * sin( theta )
       r%molz( i + 6 ) = 0d0
       mass( i )       = sjbenzene_cmass
       mass( i + 6 )   = sjbenzene_hmass
    enddo
    r%molx(SJBENZENESITE) = 0d0
    r%moly(SJBENZENESITE) = 0d0
    r%molz(SJBENZENESITE) = 0d0

    r%mass = 6d0 * ( sjbenzene_hmass + sjbenzene_cmass )
    r%massi= 1d0/r%mass

    r%Ixx = 0d0
    r%Iyy = 0d0
    r%Izz = 0d0
    do i=1,SJBENZENESITE - 1
       r%Ixx = r%Ixx + mass(i)*(r%moly(i)**2 + r%molz(i)**2)
       r%Iyy = r%Iyy + mass(i)*(r%molz(i)**2 + r%molx(i)**2)
       r%Izz = r%Izz + mass(i)*(r%molx(i)**2 + r%moly(i)**2)
    enddo
    r%ixxi = 1d0/r%ixx
    r%iyyi = 1d0/r%iyy
    r%izzi = 1d0/r%izz
  end subroutine Rigid_Sjbenzene_Constructor



  subroutine sjbenzene_setinteraction(si)
    use standard_interaction_module
    type(sStdInt),intent(inout) :: si
    integer :: i,j
    call si_allocate(si, SJBENZENESITE, LJ_COULOMB)
    intrMatrix(1:12, 1:12)%type = INTR_TYPE_SJ
    intrMatrix(1:12, 13)%type   = INTR_TYPE_NONE
    intrMatrix(13, 1:12)%type   = INTR_TYPE_NONE
    intrMatrix(13, 13)%type     = INTR_TYPE_NONE

    intrMatrix(1:6, 1:6)%param(1) = sjbenzene_acc
    intrMatrix(1:6, 1:6)%param(2) = sjbenzene_bcc
    intrMatrix(1:6, 1:6)%param(3) = sjbenzene_ccc
    intrMatrix(1:6, 1:6)%param(4) = sjbenzene_qc * sjbenzene_qc

    intrMatrix(1:6, 7:12)%param(1) = sjbenzene_ach
    intrMatrix(1:6, 7:12)%param(2) = sjbenzene_bch
    intrMatrix(1:6, 7:12)%param(3) = sjbenzene_cch
    intrMatrix(1:6, 7:12)%param(4) = - sjbenzene_qc * sjbenzene_qc

    intrMatrix(7:12, 1:6)%param(1) = sjbenzene_ach
    intrMatrix(7:12, 1:6)%param(2) = sjbenzene_bch
    intrMatrix(7:12, 1:6)%param(3) = sjbenzene_cch
    intrMatrix(7:12, 1:6)%param(4) = - sjbenzene_qc * sjbenzene_qc

    intrMatrix(7:12, 7:12)%param(1) = sjbenzene_ahh
    intrMatrix(7:12, 7:12)%param(2) = sjbenzene_bhh
    intrMatrix(7:12, 7:12)%param(3) = sjbenzene_chh
    intrMatrix(7:12, 7:12)%param(4) = sjbenzene_qc * sjbenzene_qc

    intrMatrix(1:12, 1:12)%param(1) = intrMatrix(1:12, 1:12)%param(1) * CA * 1d3 * J2I
    intrMatrix(1:12, 1:12)%param(3) = intrMatrix(1:12, 1:12)%param(3) * CA * 1d3 * J2I
    intrMatrix(1:12, 1:12)%param(4) = intrMatrix(1:12, 1:12)%param(4) * COEFF * 1d3 * J2I ! * iv%scale
    
  end subroutine sjbenzene_setinteraction



  subroutine force_sjbenzene(iv,mi,mj,site,ep,vrmat)
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
    real(kind=8) :: qq,aa,bb,cc !,ljsig,ljeps
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
    !
    !構造体の値を一時変数に入れておくと明らかに速くなる場合がある。
    !というより、構造体の処理が遅すぎる。
    !
    !平成15年1月16日(木)これでも速くならない。
    !
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
    allocate(fxi(iv%npair,SJBENZENESITE))
    allocate(fyi(iv%npair,SJBENZENESITE))
    allocate(fzi(iv%npair,SJBENZENESITE))
    allocate(fxj(iv%npair,SJBENZENESITE))
    allocate(fyj(iv%npair,SJBENZENESITE))
    allocate(fzj(iv%npair,SJBENZENESITE))
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
    do isite = 1, 12
       do jsite = 1, 12
          !if ( intrMatrix(isite, jsite)%type .eq. INTR_TYPE_SJ ) then
          aa = intrMatrix(isite, jsite)%param(1)
          bb = intrMatrix(isite, jsite)%param(2)
          cc = intrMatrix(isite, jsite)%param(3)
          qq = intrMatrix(isite, jsite)%param(4)

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
             !
             !exp-6 potential is strongly attractive when the distance is 
             ! too short. Be careful.
             !
             dd = (dx**2 + dy**2 + dz**2)
             radius = dsqrt(dd)
             radiusi = 1d0 / radius
             dd = 1d0 / dd
             e0 = qq * radiusi
             f0 = qq * dd * radiusi
             d6 = dd*dd*dd
             !if ( radius .lt. 0.75d0 ) then
             !   write(STDERR,*) i,j,radius
             !endif
             !d12= d6*d6
             !e0 = e0 + aa*d12 + bb*d6
             !f0 = f0 + dd*(aa*12d0*d12 + bb*6d0*d6)
             e0 = e0 +      aa * exp( -bb * radius ) -       cc * d6
             f0 = f0 + aa * bb * exp( -bb * radius ) * radiusi - 6d0 * cc * d6 * dd
             !write(STDERR,*) aa,bb,cc,dd,e0,f0
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
    isite = SJBENZENESITE
    jsite = SJBENZENESITE
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
  end subroutine force_sjbenzene

end module sjbenzene_module

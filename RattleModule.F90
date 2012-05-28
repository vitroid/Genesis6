! -*- f90 -*-
!
!Structure内に、massや速度も保持する。また、サイト番号も付けかえる。そうすることで、このモジュールの設計を独立にできる。
!

module rattle_module
  implicit none
  real(kind=8),parameter  :: TOL=1d-12,TOL2=TOL*2d0
  type sRattle
     !
     !束縛を受けるサイト対
     !ii,jjは番号を付けなおしたあとのサイト番号
     !
     integer, dimension(:), pointer :: pair_i,pair_j
     integer, dimension(:), pointer :: pair_ii,pair_jj
     !
     !サイト対の数
     !
     integer :: Npair,maxpair
     !
     !束縛されるサイトのみで番号を付けなおす。(必要か？)
     !
     integer, pointer :: site(:)
     integer :: Nsite
     !
     !サイトの質量。indexは付けなおし後の番号
     !
     real(kind=8), pointer :: mass(:),massi(:)
     !
     !サイトの速度。
     !rattle_positionでは、位置補正によって生じる速度補正量が返り値として入る。
     !rattle_velocityでは与えられた速度が補正されて返り値となる。
     !
     real(kind=8), pointer :: vx(:),vy(:),vz(:)
     !
     !あるサイトが束縛されている相手のリスト。
     !
     !integer, dimension(MAXrattle_site,MAXrattle_partner)::partner
     !integer, dimension(MAXrattle_site) :: Npartner
     !
     !定数
     !
     real(kind=8), pointer, dimension(:) :: reducedMass,rr,reducedMassrri
     !
     !ワークエリア
     !
     real(kind=8), pointer, dimension(:) :: deltax,deltay,deltaz
  end type sRattle

contains

  !
  !一意な番号(内部番号)を付けなおす。(Vector化のため)
  !
  function rattle_register_site(rattle,i,mass)
    implicit none
    integer, intent(in) :: i
    real(kind=8), intent(in) :: mass
    type(sRattle), intent(INOUT) :: rattle
    integer :: k,rattle_register_site
    do k=1,rattle%nsite
       if(rattle%site(k).eq.i)exit
    enddo
    if(k.gt.rattle%nsite)then
       rattle%nsite    = k
       rattle%site(k)  = i
       rattle%mass(k)  = mass
       rattle%massi(k) = 1d0 / mass
    endif
    rattle_register_site=k
  end function rattle_register_site

  !
  !内部番号の逆引き。効率がわるい。
  !
  function query_site(rattle,i)
    implicit none
    integer, intent(in) :: i
    type(sRattle), intent(INOUT) :: rattle
    integer :: k,query_site
    do k=1,rattle%nsite
       if(rattle%site(k).eq.i)exit
    enddo
    query_site = k
  end function query_site

  subroutine rattle_register_pair(rattle,i,j,length)
    use site_module
    use mol_module
    implicit none
    type(sRattle) :: rattle
    real(kind=8), intent(IN) :: length
    integer,      intent(IN) :: i,j
    integer :: k,ii,jj
    rattle%Npair = rattle%Npair + 1
    k = rattle%Npair
    rattle%pair_i(k) = i
    rattle%pair_j(k) = j
    rattle%rr(k)     = length**2
    ii = query_site(rattle, i)
    jj = query_site(rattle, j)
    rattle%reducedMass(k)    = 1d0/( rattle%massi(ii) + rattle%massi(jj))
    rattle%reducedMassrri(k) = rattle%reducedMass(k) / rattle%rr(k)
    rattle%pair_ii(k) = ii
    rattle%pair_jj(k) = jj
    !rattle%Npartner(ii)      = rattle%Npartner(ii)+1
    !rattle%partner(ii,rattle%Npartner(ii))=k
    !rattle%Npartner(jj)      = rattle%Npartner(jj)+1
    !rattle%partner(jj,rattle%Npartner(jj))=-k
  end subroutine rattle_register_pair

  subroutine rattle_constructor(rattle, maxsite, maxpair)
    implicit none
    type(sRattle),intent(out) :: rattle
    integer,          intent(in)  :: maxpair, maxsite
    integer i,j,ii,jj,k
    rattle%Npair=0
    rattle%Nsite = 0
    rattle%maxpair = maxpair
    !do i=1,MAXrattle_site
    !   rattle%Npartner(i)=0
    !enddo
    allocate(rattle%deltax(maxpair))
    allocate(rattle%deltay(maxpair))
    allocate(rattle%deltaz(maxpair))
    allocate(rattle%rr(maxpair))
    allocate(rattle%reducedMass(maxpair))
    allocate(rattle%reducedMassrri(maxpair))
    allocate(rattle%pair_i(maxpair))
    allocate(rattle%pair_j(maxpair))
    allocate(rattle%pair_ii(maxpair))
    allocate(rattle%pair_jj(maxpair))
    allocate(rattle%mass(maxsite))
    allocate(rattle%massi(maxsite))
    allocate(rattle%vx(maxsite))
    allocate(rattle%vy(maxsite))
    allocate(rattle%vz(maxsite))
    allocate(rattle%site(maxsite))
  end subroutine rattle_constructor
  
  subroutine rattle_prepare(rattle,site)
    use site_module
    implicit none
    type(sRattle) :: rattle
    type(sSite) :: site
    integer :: i,j,ij
    do ij=1,rattle%Npair
       i = rattle%pair_i(ij)
       j = rattle%pair_j(ij)
       rattle%deltax(ij) = site%x(i)-site%x(j)
       rattle%deltay(ij) = site%y(i)-site%y(j)
       rattle%deltaz(ij) = site%z(i)-site%z(j)
    enddo
  end subroutine rattle_prepare
  
  !
  !aka "SHAKE"
  !
  subroutine rattle_position(rattle,site)
    use site_module
    implicit none
    type(sRattle) :: rattle
    type(sSite) :: site
    !     local variables
    real(kind=8) ::  prod,gab,dd,gabi,gabj
    real(kind=8) ::  deltax,deltay,deltaz
    integer k,l,si,ij,ii,jj
    integer ncorrect
    real(kind=8), dimension(rattle%maxpair) :: fxi,fyi,fzi,fxj,fyj,fzj
    real(kind=8), dimension(rattle%nsite)   :: x,y,z
    !
    ! 3分子だけでも収束まで50回ぐらいループしなければいけないようだ。
    !     なんとか高速化できないか？行列演算した方が速いかもしれない。
    !
    call rattle_prepare( rattle, site )

    !     ...まず総てのサイトの新座標を準備する。
    do k=1,rattle%nsite
       si = rattle%site(k)
       x(k) = site%x(si)
       y(k) = site%y(si)
       z(k) = site%z(si)
    enddo
    do
       ncorrect=0
       do ij=1,rattle%Npair
          ii = rattle%pair_ii(ij)
          jj = rattle%pair_jj(ij)
          deltax = x(ii) - x(jj)
          deltay = y(ii) - y(jj)
          deltaz = z(ii) - z(jj)
          dd = deltax**2 + deltay**2 + deltaz**2
          dd = rattle%rr(ij) - dd
          !write(6,*) "Diff: ",ij,dd
          if( rattle%rr(ij) * TOL2 < dabs(dd) )then
             prod = deltax*rattle%deltax(ij) + deltay*rattle%deltay(ij) + deltaz*rattle%deltaz(ij)
             if( prod < rattle%rr(ij) * 1d-6 )then
                write(STDERR,*) "Rattle Failure."
                stop
             endif
             gab = dd*rattle%reducedMass(ij)/(2d0*prod)
             gabi = gab*rattle%massi(ii)
             gabj = gab*rattle%massi(jj)
             fxi(ij) = rattle%deltax(ij) * gabi
             fyi(ij) = rattle%deltay(ij) * gabi
             fzi(ij) = rattle%deltaz(ij) * gabi
             fxj(ij) =-rattle%deltax(ij) * gabj
             fyj(ij) =-rattle%deltay(ij) * gabj
             fzj(ij) =-rattle%deltaz(ij) * gabj
             ncorrect = ncorrect+1
          else
             fxi(ij) = 0d0
             fyi(ij) = 0d0
             fzi(ij) = 0d0
             fxj(ij) = 0d0
             fyj(ij) = 0d0
             fzj(ij) = 0d0
          endif
       enddo
       if(ncorrect.gt.0)then
          !     count=count+1
#ifdef VECTOR
          do l=1,MAXrattle_partner
             do k=1,Nsite
                !i = rattle%site(j)
                !k=rattle%partner(i,l)
                ij=rattle%partner(k,l)
                if(ij.gt.0)then
                   x(k) = x(k) + fxi(ij)
                   y(k) = y(k) + fyi(ij)
                   z(k) = z(k) + fzi(ij)
                endif
                if(ij.lt.0)then
                   ij=-ij
                   x(k) = x(k) + fxj(ij)
                   y(k) = y(k) + fyj(ij)
                   z(k) = z(k) + fzj(ij)
                endif
             enddo
          enddo
#else 
          do ij=1,rattle%Npair
             ii = rattle%pair_ii(ij)
             jj = rattle%pair_jj(ij)
             x(ii) = x(ii) + fxi(ij)
             y(ii) = y(ii) + fyi(ij)
             z(ii) = z(ii) + fzi(ij)
             x(jj) = x(jj) + fxj(ij)
             y(jj) = y(jj) + fyj(ij)
             z(jj) = z(jj) + fzj(ij)
          enddo
#endif
        !     do k=1,Npair
        !     i = pair_i(k)
        !     j = pair_j(k)
        !     deltax = site%x(i)-site%x(j)
        !     deltay = site%y(i)-site%y(j)
        !     deltaz = site%z(i)-site%z(j)
        !     write(6,*) k,deltax**2+deltay**2+deltaz**2
        !     enddo

        else
           exit
        endif
     enddo
     !     write(7,*) count
     !     ...位置の補正にともなう速度の修正はとりあえず行わない。(Verletなら
     !     いいのだが、Gearの場合には何次微分まで修正する必要があるだろう
     do k=1,rattle%Nsite
        si = rattle%site(k)
        !
        !速度の補正量をvxに返す。それらをflexに反映させる作業はここでは行わない。
        !
        rattle%vx(k) = x(k) - site%x(si)
        rattle%vy(k) = y(k) - site%y(si)
        rattle%vz(k) = z(k) - site%z(si)
        site%x(si) = x(k)
        site%y(si) = y(k)
        site%z(si) = z(k)
     enddo
  end subroutine rattle_position

  subroutine rattle_velocity( rattle, site )
    use site_module
    implicit none
    type(sRattle) :: rattle
    type(sSite) :: site
    !     args
    !     local variables
    real(kind=8) ::  prod,gab
    real(kind=8) ::  deltax,deltay,deltaz
    integer k,l,ij,ii,jj
    integer ncorrect
    real(kind=8), dimension(rattle%maxpair) :: fxi,fyi,fzi,fxj,fyj,fzj
    
    call rattle_prepare(rattle,site)
    !
    !rattle%vxにはあらかじめ各サイトの速度を入れておく。補正して返す。
    !
    do
       ncorrect=0
       do ij=1,rattle%Npair
          ii = rattle%pair_ii(ij)
          jj = rattle%pair_jj(ij)
          deltax = rattle%vx(ii) - rattle%vx(jj)
          deltay = rattle%vy(ii) - rattle%vy(jj)
          deltaz = rattle%vz(ii) - rattle%vz(jj)
          prod = deltax*rattle%deltax(ij) + deltay*rattle%deltay(ij) + deltaz*rattle%deltaz(ij)
          gab = -prod * rattle%reducedMassrri(ij)
          !write(STDOUT,*) ncorrect, gab
          if( TOL < dabs(gab))then
             deltax = rattle%deltax(ij)*gab
             deltay = rattle%deltay(ij)*gab
             deltaz = rattle%deltaz(ij)*gab
             fxi(ij) = deltax
             fyi(ij) = deltay
             fzi(ij) = deltaz
             fxj(ij) =-deltax
             fyj(ij) =-deltay
             fzj(ij) =-deltaz
             ncorrect = ncorrect+1
          else
             fxi(ij) = 0d0
             fyi(ij) = 0d0
             fzi(ij) = 0d0
             fxj(ij) = 0d0
             fyj(ij) = 0d0
             fzj(ij) = 0d0
          endif
       enddo
       if(ncorrect.gt.0)then
#ifdef VECTOR
          do l=1,MAXrattle_partner
             do k=1,rattle%Nsite
                ij=rattle%partner(k,l)
                if(ij.gt.0)then
                   rattle%vx(k) = rattle%vx(k) + fxi(ij) * rattle%massi(k)
                   rattle%vy(k) = rattle%vy(k) + fyi(ij) * rattle%massi(k)
                   rattle%vz(k) = rattle%vz(k) + fzi(ij) * rattle%massi(k)
                endif
                if(ij.lt.0)then
                   ij=-ij
                   rattle%vx(k) = rattle%vx(k) + fxj(ij)*rattle%massi(k)
                   rattle%vy(k) = rattle%vy(k) + fyj(ij)*rattle%massi(k)
                   rattle%vz(k) = rattle%vz(k) + fzj(ij)*rattle%massi(k)
                endif
             enddo
          enddo
#else 
          do ij=1,rattle%Npair
             ii = rattle%pair_ii(ij)
             jj = rattle%pair_jj(ij)
             rattle%vx(ii) = rattle%vx(ii) + fxi(ij) * rattle%massi(ii)
             rattle%vy(ii) = rattle%vy(ii) + fyi(ij) * rattle%massi(ii)
             rattle%vz(ii) = rattle%vz(ii) + fzi(ij) * rattle%massi(ii)
             rattle%vx(jj) = rattle%vx(jj) + fxj(ij) * rattle%massi(jj)
             rattle%vy(jj) = rattle%vy(jj) + fyj(ij) * rattle%massi(jj)
             rattle%vz(jj) = rattle%vz(jj) + fzj(ij) * rattle%massi(jj)
          enddo
#endif
        else
           exit
        endif
     enddo
  end subroutine rattle_velocity
end module rattle_module

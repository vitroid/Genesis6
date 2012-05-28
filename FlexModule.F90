! -*- f90 -*-
!
!RigidModule/FlexModuleの役割は、サイト単位で受けた力を、分子単位に束
!ねることにある。
!
#undef FLEXDEBUG1
!
!柔軟分子の運動を扱う。Argonなどの単原子分子も含む。
!
!相互作用点(サイト)に関する説明はSiteModule.F90のコメントを参照。
!
!柔軟分子の相互作用のTruncationは、重心間距離で計算する。
!慣習上、分子の最後のサイトをカットオフの基準サイ
!トとする。CO2のように重心位置に相互作用点がある場合は、サイトの順序を
!O、O、Cとする。また、もっと柔軟なethanol分子の場合は、重心仮想サイトを一
!つ追加して7サイト模型とする。単原子分子であっても2サイト(相互作用点と
!しての1サイトと、重心としての1サイト)として取り扱う。
!
module flex_module
  use common_module
  use mol_module
  use box_module
  use site_module
  use vector_module
  use error_module
  implicit none
  integer, parameter :: FLEX_AR3A=0, FLEX_APC5=1, FLEX_MDVW=2, FLEX_FL3A=3
  type sFlex1
     sequence
     real(kind=8),dimension(5) :: dx,dy,dz
     integer                   :: rattlesite
  end type sFlex1
  
  type sFlex
     sequence
     !
     !site数だけmassを確保しておく。
     !
     real(kind=8),pointer :: mass(:)
     !
     !extra information for sites
     !こちらもサイト数分確保
     !
     type(sFlex1),pointer :: xinfo(:)
  end type sFlex
  
  interface save
     module procedure flex_save
  end interface

  interface save_binary
     module procedure flex_savebinary
  end interface

contains

  subroutine flex_allocate(r,m)
    use mol_module
    type(sFlex),intent(inout) :: r
    type(sMol),intent(in)     :: m
#ifdef VERBOSE
    write(STDERR,*) "Flex_allocate:", m%nmol * m%nsite
#endif
    allocate(r%mass(m%nmol * m%nsite))
    allocate(r%xinfo(m%nmol * m%nsite))
  end subroutine flex_allocate
  
!
!まったくチェックしていない
!
#if defined(PVM) || defined(MPI)
  subroutine flex_collectforce(NPROCS,MYRANK,r,m,s)
    use mol_module
    use site_module
    implicit none
#ifdef MPI
    include "mpif.h"
    integer,intent(IN) :: NPROCS,MYRANK
    integer :: IERR
#endif
#ifdef PVM
    include "fpvm3.h"
    integer,intent(IN) :: NPROCS,MYRANK
    integer :: IERR,tid,msgtag
#endif
    type(sFlex) :: r
    type(sSite) :: s
    type(sMol) :: m
    integer :: j,k
    real(kind=8),dimension(m%nmol*m%nsite) :: forcex,forcey,forcez
#ifdef MPI
    call MPI_ALLREDUCE(s%fx(m%offset+1),forcex,m%nmol*m%nsite
    & ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(s%fy(m%offset+1),forcey,m%nmol*m%nsite
    & ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(s%fz(m%offset+1),forcez,m%nmol*m%nsite
    & ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
    !
    !再分配
    !
    do j=1,m%nsite * m%nmol
       k=j + m%offset
       s%fx(k) = forcex(j)
       s%fy(k) = forcey(j)
       s%fz(k) = forcez(j)
    enddo
#else /*MPI*/
    call die( 0, "Flex 1" )
#endif
    return
  end subroutine flex_collectforce
#endif

  subroutine flex_predict(r,m,s)
    type(sMol), intent(in)    :: m
    type(sFlex),intent(INOUT) :: r
    type(sSite),intent(INOUT) :: s
    integer i,j,mol,site
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          j = i + m%offset
          s%x(j) = s%x(j) + r%xinfo(i)%dx(1)+r%xinfo(i)%dx(2)+r&
               & %xinfo(i)%dx(3)+r%xinfo(i)%dx(4)+r%xinfo(i)%dx(5)
          r%xinfo(i)%dx(1) = r%xinfo(i)%dx(1)+2d0*r%xinfo(i)%dx(2)+3d0*r&
               & %xinfo(i)%dx(3)+4d0*r%xinfo(i)%dx(4)+5d0*r%xinfo(i)%dx(5)
          r%xinfo(i)%dx(2) = r%xinfo(i)%dx(2)+3d0*r%xinfo(i)%dx(3)+6d0*r&
               & %xinfo(i)%dx(4)+10d0*r%xinfo(i)%dx(5)
          r%xinfo(i)%dx(3) = r%xinfo(i)%dx(3)+4d0*r%xinfo(i)%dx(4)+10d0*r&
               & %xinfo(i)%dx(5)
          r%xinfo(i)%dx(4) = r%xinfo(i)%dx(4)+5d0*r%xinfo(i)%dx(5)
          s%y(j) = s%y(j) + r%xinfo(i)%dy(1)+r%xinfo(i)%dy(2)+r&
               & %xinfo(i)%dy(3)+r%xinfo(i)%dy(4)+r%xinfo(i)%dy(5)
          r%xinfo(i)%dy(1) = r%xinfo(i)%dy(1)+2d0*r%xinfo(i)%dy(2)+3d0*r&
               & %xinfo(i)%dy(3)+4d0*r%xinfo(i)%dy(4)+5d0*r%xinfo(i)%dy(5)
          r%xinfo(i)%dy(2) = r%xinfo(i)%dy(2)+3d0*r%xinfo(i)%dy(3)+6d0*r&
               & %xinfo(i)%dy(4)+10d0*r%xinfo(i)%dy(5)
          r%xinfo(i)%dy(3) = r%xinfo(i)%dy(3)+4d0*r%xinfo(i)%dy(4)+10d0*r&
               & %xinfo(i)%dy(5)
          r%xinfo(i)%dy(4) = r%xinfo(i)%dy(4)+5d0*r%xinfo(i)%dy(5)
          s%z(j) = s%z(j) + r%xinfo(i)%dz(1)+r%xinfo(i)%dz(2)+r&
               & %xinfo(i)%dz(3)+r%xinfo(i)%dz(4)+r%xinfo(i)%dz(5)
          r%xinfo(i)%dz(1) = r%xinfo(i)%dz(1)+2d0*r%xinfo(i)%dz(2)+3d0*r&
               & %xinfo(i)%dz(3)+4d0*r%xinfo(i)%dz(4)+5d0*r%xinfo(i)%dz(5)
          r%xinfo(i)%dz(2) = r%xinfo(i)%dz(2)+3d0*r%xinfo(i)%dz(3)+6d0*r&
               & %xinfo(i)%dz(4)+10d0*r%xinfo(i)%dz(5)
          r%xinfo(i)%dz(3) = r%xinfo(i)%dz(3)+4d0*r%xinfo(i)%dz(4)+10d0*r&
               & %xinfo(i)%dz(5)
          r%xinfo(i)%dz(4) = r%xinfo(i)%dz(4)+5d0*r%xinfo(i)%dz(5)
       enddo
    enddo
  end subroutine flex_predict

  subroutine Flex_Ek(r,m,t,ekt)
    use time_module
    use mol_module
    implicit none
    type(sMol),intent(IN) :: m
    type(sFlex),intent(IN) :: r
    type(sTime),intent(IN) :: t
    real(kind=8),intent(out) :: ekt
    integer i,site,mol
    ekt=0
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          ekt = ekt + r%mass(i) * (r%xinfo(i)%dx(1)**2+r%xinfo(i)%dy(1)**2+r%xinfo(i)%dz(1)**2)
       enddo
    enddo
    ekt = ekt * 0.5d0*(1d0/(t%Dt**2))
  end subroutine Flex_Ek

  subroutine Flex_PressureTensor(r,m,t,vrmat)
    use time_module
    type(sMol),intent(IN) :: m
    type(sFlex),intent(IN) :: r
    type(sTime),intent(IN) :: t
    real(kind=8),intent(out) :: vrmat(3,3)
    integer :: i,site,mol
    vrmat(:,:)=0d0
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          vrmat(1,1) = vrmat(1,1) + r%mass(i) * (r%xinfo(i)%dx(1) *r%xinfo(i)%dx(1))
          vrmat(1,2) = vrmat(1,2) + r%mass(i) * (r%xinfo(i)%dx(1) *r%xinfo(i)%dy(1))
          vrmat(1,3) = vrmat(1,3) + r%mass(i) * (r%xinfo(i)%dx(1) *r%xinfo(i)%dz(1))
          vrmat(2,1) = vrmat(2,1) + r%mass(i) * (r%xinfo(i)%dy(1) *r%xinfo(i)%dx(1))
          vrmat(2,2) = vrmat(2,2) + r%mass(i) * (r%xinfo(i)%dy(1) *r%xinfo(i)%dy(1))
          vrmat(2,3) = vrmat(2,3) + r%mass(i) * (r%xinfo(i)%dy(1) *r%xinfo(i)%dz(1))
          vrmat(3,1) = vrmat(3,1) + r%mass(i) * (r%xinfo(i)%dz(1) *r%xinfo(i)%dx(1))
          vrmat(3,2) = vrmat(3,2) + r%mass(i) * (r%xinfo(i)%dz(1) *r%xinfo(i)%dy(1))
          vrmat(3,3) = vrmat(3,3) + r%mass(i) * (r%xinfo(i)%dz(1) *r%xinfo(i)%dz(1))
       enddo
    enddo
    vrmat(:,:) = vrmat(:,:) / (t%Dt**2)
  end subroutine Flex_PressureTensor

  subroutine Flex_ScaleVelocity(r,m,ratio)
    use mol_module
    implicit none
    type(sMol),intent(IN) :: m
    type(sFlex),intent(INOUT) :: r
    real(kind=8),intent(IN) :: ratio
    integer i,mol,site
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          r%xinfo( i )%dx(1) = r%xinfo( i )%dx(1) * ratio
          r%xinfo( i )%dy(1) = r%xinfo( i )%dy(1) * ratio
          r%xinfo( i )%dz(1) = r%xinfo( i )%dz(1) * ratio
       enddo
    enddo
  end subroutine Flex_ScaleVelocity
  
  subroutine Flex_Version
    write(STDERR,*) "$Id: Flex.F90,v 1.22 2002/11/01 07:22:49 matto&
         & Exp $"
    return
  end subroutine Flex_Version

  subroutine flex_vreset2(r,m,t,temp,random)
    use physconst_module
    use random_module
    use time_module
    use mol_module
    implicit none
    type(sFlex),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    type(sTime),intent(IN) :: t
    !temperature
    real(kind=8),intent(in) :: temp
    type(pRandom), pointer :: random
    integer :: i, site, mol
    real(kind=8) :: r1,r2
    !velocity distribution is normal distribution (avg=0, var=sig**2
    ! =kT/m)
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          r1=Random_GetNext( random )
          r2=Random_GetNext( random )
          r%xinfo(i)%dx(1)=Random_Normal(0d0,temp*rgas*J2I/r%mass(i),r1,r2)&
               & *t%Dt
          r1=Random_GetNext( random )
          r2=Random_GetNext( random )
          r%xinfo(i)%dy(1)=Random_Normal(0d0,temp*rgas*J2I/r%mass(i),r1,r2)&
               & *t%Dt
          r1=Random_GetNext( random )
          r2=Random_GetNext( random )
          r%xinfo(i)%dz(1)=Random_Normal(0d0,temp*rgas*J2I/r%mass(i),r1,r2)&
               & *t%Dt
       enddo
    enddo
  end subroutine flex_vreset2

! Get total momenta of the group r
  subroutine flex_gettotalmomenta(r,m,t,momentx,momenty,momentz,mass)
    use time_module
    use mol_module
    implicit none
    type(sFlex),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    type(sTime),intent(IN) :: t
    real(kind=8),intent(out) :: momentx,momenty,momentz,mass
    integer :: i,mol,site
    momentx=0d0
    momenty=0d0
    momentz=0d0
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          momentx = momentx + r%xinfo(i)%dx(1) * r%mass(i)
          momenty = momenty + r%xinfo(i)%dy(1) * r%mass(i)
          momentz = momentz + r%xinfo(i)%dz(1) * r%mass(i)
       enddo
    enddo
    momentx = momentx / t%dt
    momenty = momenty / t%dt
    momentz = momentz / t%dt
  end subroutine flex_gettotalmomenta

  subroutine flex_getmomentum(r,m,moment,mass)
    use time_module
    use mol_module
    implicit none
    type(sFlex),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    type(vector3),intent(OUT)   :: moment
    real(kind=8), intent(OUT)   :: mass
    integer :: i,mol,site
    moment%vec(:) = 0d0
    mass          = 0d0
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          mass    = mass + r%mass(i)
          moment%vec(1) = moment%vec(1) + r%xinfo(i)%dx(1) * r%mass(i)
          moment%vec(2) = moment%vec(2) + r%xinfo(i)%dy(1) * r%mass(i)
          moment%vec(3) = moment%vec(3) + r%xinfo(i)%dz(1) * r%mass(i)
       enddo
    enddo
  end subroutine flex_getmomentum

  subroutine flex_addmomenta(flex,m,moment)
    type(sFlex),intent(INOUT) :: flex
    type(sMol),intent(IN) :: m
    type(vector3),intent(IN) :: moment
    real(kind=8) :: velx,vely,velz
    real(kind=8) :: nmol
    integer      :: i,mol,site
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          flex%xinfo(i)%dx(1) = flex%xinfo(i)%dx(1) + moment%vec(1)
          flex%xinfo(i)%dy(1) = flex%xinfo(i)%dy(1) + moment%vec(2)
          flex%xinfo(i)%dz(1) = flex%xinfo(i)%dz(1) + moment%vec(3)
       enddo
    enddo
  end subroutine flex_addmomenta
  
! Add total momenta to the group r
  subroutine flex_addmomenta_obsolete(r,m,t,momentx,momenty,momentz)
    use time_module
    use mol_module
    implicit none
    type(sFlex),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    type(sTime),intent(IN) :: t
    real(kind=8),intent(IN) :: momentx,momenty,momentz
    real(kind=8) :: velx,vely,velz
    real(kind=8) :: nmol
    real(kind=8) :: mass
    integer      :: i,mol,site
    mass=0d0
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          mass = mass + r%mass(i)
       enddo
    enddo
    velx = momentx * t%dt / mass
    vely = momenty * t%dt / mass
    velz = momentz * t%dt / mass
    r%xinfo(1:m%nmol)%dx(1) = r%xinfo(1:m%nmol)%dx(1) + velx
    r%xinfo(1:m%nmol)%dy(1) = r%xinfo(1:m%nmol)%dy(1) + vely
    r%xinfo(1:m%nmol)%dz(1) = r%xinfo(1:m%nmol)%dz(1) + velz
  end subroutine flex_addmomenta_obsolete

  subroutine flex_correct(r,m,s,t,fNose,nose,fAndersen,andersen)
    use physconst_module
    use mol_module
    use time_module
    use nose_module
    use andersen_module
    implicit none
    type(sMol),intent(IN) :: m
    type(sFlex),intent(INOUT) :: r
    type(sSite),intent(INOUT) :: s
    type(sTime),intent(IN) :: t
    type(sNose),intent(IN) :: nose
    type(sAndersen),intent(IN) :: andersen
    logical,intent(IN) :: fNose,fAndersen
    integer i,j,mol,site
    real(kind=8) deltax,deltay,deltaz
    real(kind=8) newfx,newfy,newfz
    real(kind=8) :: nosecoeff,acoeff
    real(kind=8) :: extendtx, extendty,extendtz
    real(kind=8) :: extendrx, extendry,extendrz
    extendtx = 0d0
    extendty = 0d0
    extendtz = 0d0
    extendrx = 0d0
    extendry = 0d0
    extendrz = 0d0
    if(fNose)then
       nosecoeff = nose%zeta0*t%dt
       extendtx = extendtx + nosecoeff * 0.5d0
       extendty = extendty + nosecoeff * 0.5d0
       extendtz = extendtz + nosecoeff * 0.5d0
       extendrx = extendrx + nosecoeff
       extendry = extendry + nosecoeff
       extendrz = extendrz + nosecoeff
    endif
    if(fAndersen) then
       !acoeffは上田本(9.17)のdV / ( 3 V dt )
       if ( andersen%mode .eq. orz ) then
          acoeff = andersen%v1 / (2d0*andersen%v0)
          extendtz = extendtz + acoeff
       else if ( andersen%mode .eq. orthorhombic ) then
          !Ueda's formula. nph.cでは問題なく動いていたのでこれでいいはず。
          acoeff = andersen%v1 / (2d0*3d0*andersen%v0)
          !Okazaki's fomula
          !acoeff = andersen%v1 / (3d0*andersen%v0)
          extendtx = extendtx + acoeff
          extendty = extendty + acoeff
          extendtz = extendtz + acoeff
       endif
    endif
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          j = i + m%offset
          newfx = s%fx(j) * (t%dt*t%dt*0.5d0/r%Mass(i))
          newfy = s%fy(j) * (t%dt*t%dt*0.5d0/r%Mass(i))
          newfz = s%fz(j) * (t%dt*t%dt*0.5d0/r%Mass(i))
          newfx = newfx - r%xinfo(i)%dx(1) * extendtx
          newfy = newfy - r%xinfo(i)%dy(1) * extendty
          newfz = newfz - r%xinfo(i)%dz(1) * extendtz
          deltax = newfx - r%xinfo(i)%dx(2)
          deltay = newfy - r%xinfo(i)%dy(2)
          deltaz = newfz - r%xinfo(i)%dz(2)
          s%x(j) = s%x(j) + Gear5_0*deltax
          s%y(j) = s%y(j) + Gear5_0*deltay
          s%z(j) = s%z(j) + Gear5_0*deltaz
          r%xinfo(i)%dx(1)=r%xinfo(i)%dx(1) + Gear5_1*deltax
          r%xinfo(i)%dy(1)=r%xinfo(i)%dy(1) + Gear5_1*deltay
          r%xinfo(i)%dz(1)=r%xinfo(i)%dz(1) + Gear5_1*deltaz
          r%xinfo(i)%dx(2)=newfx
          r%xinfo(i)%dy(2)=newfy
          r%xinfo(i)%dz(2)=newfz
          r%xinfo(i)%dx(3)=r%xinfo(i)%dx(3) + Gear5_3*deltax
          r%xinfo(i)%dy(3)=r%xinfo(i)%dy(3) + Gear5_3*deltay
          r%xinfo(i)%dz(3)=r%xinfo(i)%dz(3) + Gear5_3*deltaz
          r%xinfo(i)%dx(4)=r%xinfo(i)%dx(4) + Gear5_4*deltax
          r%xinfo(i)%dy(4)=r%xinfo(i)%dy(4) + Gear5_4*deltay
          r%xinfo(i)%dz(4)=r%xinfo(i)%dz(4) + Gear5_4*deltaz
          r%xinfo(i)%dx(5)=r%xinfo(i)%dx(5) + Gear5_5*deltax
          r%xinfo(i)%dy(5)=r%xinfo(i)%dy(5) + Gear5_5*deltay
          r%xinfo(i)%dz(5)=r%xinfo(i)%dz(5) + Gear5_5*deltaz
       enddo
    enddo
  end subroutine flex_correct

!後方互換
  subroutine flex_correct2(r,m,s,t,nose,andersen)
    use mol_module
    use andersen_module
    use nose_module
    use time_module
    implicit none
    type(sMol),intent(IN) :: m
    type(sFlex),intent(INOUT) :: r
    type(sSite),intent(INOUT) :: s
    type(sTime),intent(IN) :: t
    type(sNose),intent(IN) :: nose
    type(sAndersen),intent(IN),optional :: andersen
    if(present(andersen))then
       call flex_correct(r,m,s,t,nose%active,nose,andersen%mode.ne.noandersen,andersen)
    else
       call flex_correct(r,m,s,t,nose%active,nose,.false.,andersen)
    endif
  end subroutine flex_correct2

  function flex_serialize_position(flex,s,m,p,offset)
    use quat_module
    integer :: flex_serialize_position
    type(sFlex),intent(in)  :: flex
    type(sSite),intent(in)  :: s
    type(sMol),intent(in)    :: m
    real(kind=8),intent(inout) :: p(*)
    integer,intent(inout)    :: offset
    integer                  :: i,j,mol,site
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          j = i + m%offset
          p(offset + i + m%nmol*(m%nsite-1)*0)=s%x(j)
          p(offset + i + m%nmol*(m%nsite-1)*1)=s%y(j)
          p(offset + i + m%nmol*(m%nsite-1)*2)=s%z(j)
       enddo
    enddo
    flex_serialize_position = offset + m%nmol*(m%nsite-1) * 3
  end function flex_serialize_position

  integer function flex_unserialize_position(flex,s,m,p,offset)
    use quat_module
    use mol_module
    implicit none
    type(sFlex),intent(inout) :: flex
    type(sSite),intent(inout) :: s
    type(sMol),intent(in)      :: m
    real(kind=8),intent(in)    :: p(*)
    integer,intent(inout)      :: offset
    integer                    :: i,j,site,mol
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          j = i + m%offset
          s%x(j) = p(offset + i + m%nmol*(m%nsite-1)*0)
          s%y(j) = p(offset + i + m%nmol*(m%nsite-1)*1)
          s%z(j) = p(offset + i + m%nmol*(m%nsite-1)*2)
       enddo
    enddo
    flex_unserialize_position = offset + m%nmol*(m%nsite-1) * 3
  end function flex_unserialize_position

  subroutine Flex_SetMass( flex, m, mass )
    type(sFlex),intent(inout) :: flex
    type(sMol),intent(in)     :: m
    real(kind=8),intent(in)   :: mass(*)
    !
    !local
    !
    integer :: i,mol,site
    allocate(flex%mass(m%nmol*m%nsite))
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          flex%mass(i) = mass(site)
       enddo
    enddo
  end subroutine Flex_SetMass
  !
  !単原子分子のみ
  !
  subroutine Flex_ReadAR3A(r,m,si,dreset,file)
    use site_module
    type(sFlex),intent(INOUT) :: r
    type(sMol),intent(INOUT) :: m
    type(sSite),intent(INOUT) :: si
    logical,intent(in) :: dreset
    integer,intent(IN) :: file
    integer :: n,i,s,mol,site,k
    if ( m%nsite /= 2 ) then
       write(STDERR,*) "AR3A is for monatomic molecule. Use FL3A instead.", m%nsite
       call die( 0, "Flex 2" )
    endif
    read(file,*) n
    call Mol_Allocate(m,si,n,FLEX_MODE)
    call Flex_Allocate(r,m)
    if(dreset)then
       do i=1,n
          do s=1,5
             r%xinfo(i)%dx(s)=0
             r%xinfo(i)%dy(s)=0
             r%xinfo(i)%dz(s)=0
          enddo
       enddo
    endif
    do mol = 1, m%nmol
       do site = 1, m%nsite-1
          k = (mol-1)*m%nsite + site + m%offset
          read(file,*) si%x(k),si%y(k),si%z(k)
          !write(STDERR,*) i,si%x(i),si%y(i),si%z(i)
       enddo
    enddo
  end subroutine Flex_ReadAR3A

  subroutine Flex_ReadFL3A(r,m,si,dreset,file)
    use site_module
    type(sFlex),intent(INOUT) :: r
    type(sMol),intent(INOUT) :: m
    type(sSite),intent(INOUT) :: si
    logical,intent(in) :: dreset
    integer,intent(IN) :: file
    integer :: site,mol,nsite,nmol,i,s,k
    read(file,*) nsite, nmol
    if ( nsite+1 /= m%nsite ) then
       call die( error_different_num_of_molecules, "Flex 3" )
    endif
    call Mol_Allocate(m,si,nmol,FLEX_MODE)
    call Flex_Allocate(r,m)
    if(dreset)then
       do mol=1,nmol
          do site=1,nsite
             k = (mol-1)*m%nsite + site + m%offset
             do s=1,5
                r%xinfo(k)%dx(s)=0
                r%xinfo(k)%dy(s)=0
                r%xinfo(k)%dz(s)=0
             enddo
          enddo
       enddo
    endif
    do mol = 1, m%nmol
       do site = 1, nsite
          k = (mol-1)*m%nsite + site + m%offset
          read(file,*) si%x(k),si%y(k),si%z(k)
          !write(STDERR,*) i,si%x(i),si%y(i),si%z(i)
       enddo
    enddo
  end subroutine Flex_ReadFL3A

  subroutine flex_toInternalUnit( flex, m, time )
    use time_module
    type(sMol),intent(in)     :: m
    type(sFlex),intent(INOUT) :: flex
    type(sTime),intent(IN)    :: time
    integer :: i,s,mol,site
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          do s=1,5
             flex%xinfo(i)%dx(s) = flex%xinfo(i)%dx(s) / time%dt
             flex%xinfo(i)%dy(s) = flex%xinfo(i)%dy(s) / time%dt
             flex%xinfo(i)%dz(s) = flex%xinfo(i)%dz(s) / time%dt
          enddo
       enddo
    enddo
  end subroutine flex_toInternalUnit

  subroutine flex_toExternalUnit( flex, m, time )
    use time_module
    type(sMol),intent(in)     :: m
    type(sFlex),intent(INOUT) :: flex
    type(sTime),intent(IN)    :: time
    integer :: i,s,mol,site
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          do s=1,5
             flex%xinfo(i)%dx(s) = flex%xinfo(i)%dx(s) * time%dt
             flex%xinfo(i)%dy(s) = flex%xinfo(i)%dy(s) * time%dt
             flex%xinfo(i)%dz(s) = flex%xinfo(i)%dz(s) * time%dt
          enddo
       enddo
    enddo
  end subroutine flex_toExternalUnit

  subroutine Flex_GetOffset(flex,m,s,moment,mass)
    type(sFlex),  intent(IN)    :: flex
    type(sMol),   intent(IN)    :: m
    type(sSite),  intent(IN)    :: s
    type(vector3),intent(OUT)   :: moment
    real(kind=8), intent(OUT)   :: mass
    !
    !local
    !
    real(kind=8) :: x,y,z
    integer :: i,j,mol,site
    moment%vec(:) =0d0
    mass = 0d0
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          j = i + m%offset
          mass          = mass + flex%mass(i)
          moment%vec(1) = moment%vec(1) + s%x(j) * flex%mass(i)
          moment%vec(2) = moment%vec(2) + s%y(j) * flex%mass(i)
          moment%vec(3) = moment%vec(3) + s%z(j) * flex%mass(i)
       enddo
    enddo
  end subroutine Flex_GetOffset

  !
  ! Get total momenta of the group r 
  ! momentの値にはdtがかかっている。
  !
  subroutine Flex_GetAngularMomenta(r,m,s,moment)
    type(sFlex),  intent(IN)    :: r
    type(sMol),   intent(IN)    :: m
    type(sSite),  intent(IN)    :: s
    type(vector3),intent(OUT)   :: moment
    !
    !local
    !
    real(kind=8) :: x,y,z,vx,vy,vz
    integer :: i,j,mol,site
    moment%vec(:) =0d0
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          j = i + m%offset
          x = s%x(j)
          y = s%y(j)
          z = s%z(j)
          vx = r%xinfo(i)%dx(1)
          vy = r%xinfo(i)%dy(1)
          vz = r%xinfo(i)%dz(1)
          moment%vec(1) = moment%vec(1) + ( y*vz - z*vy ) * r%mass(j)
          moment%vec(2) = moment%vec(2) + ( z*vx - x*vz ) * r%mass(j)
          moment%vec(3) = moment%vec(3) + ( x*vy - y*vx ) * r%mass(j)
       enddo
    enddo
  end subroutine Flex_GetAngularMomenta

  subroutine flex_getinertiatensor(flex,m,s,t)
    type(sFlex),  intent(IN) :: flex
    type(sMol),  intent(IN)  :: m
    type(sSite), intent(IN)  :: s
    real(kind=8),intent(OUT) :: t(3,3)
    !
    !Local variables
    !
    real(kind=8) :: x,y,z,mass
    integer :: i,j,mol,site
    t(:,:) = 0d0
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          j = i + m%offset
          x = s%x(j)
          y = s%y(j)
          z = s%z(j)
          mass = flex%mass(i)
          t(1,1) = t(1,1) + mass*( y**2 + z**2 )
          t(1,2) = t(1,2) - mass*( x*y )
          t(1,3) = t(1,3) - mass*( x*z )
          t(2,1) = t(2,1) - mass*( y*x )
          t(2,2) = t(2,2) + mass*( x**2 + z**2 )
          t(2,3) = t(2,3) - mass*( y*z )
          t(3,1) = t(3,1) - mass*( z*x )
          t(3,2) = t(3,2) - mass*( z*y )
          t(3,3) = t(3,3) + mass*( x**2 + y**2 )
       enddo
    enddo
  end subroutine flex_getinertiatensor

! Add total momenta to the group r
  subroutine flex_AddAngularVelocity(r,m,s,w)
    type(sFlex),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    type(sSite),  intent(IN) :: s
    type(vector3),intent(in) :: w
    !
    !Local variables
    !
    real(kind=8) :: velx,vely,velz
    real(kind=8) :: x,y,z
    integer      :: i,j,mol,site
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          j = i + m%offset
          !
          !goldstein (5-1)
          !
          x = s%x(j)
          y = s%y(j)
          z = s%z(j)
          velx = w%vec(2) * z - w%vec(3) * y
          vely = w%vec(3) * x - w%vec(1) * z
          velz = w%vec(1) * y - w%vec(2) * x
          r%xinfo(i)%dx(1) = r%xinfo(i)%dx(1) + velx
          r%xinfo(i)%dy(1) = r%xinfo(i)%dy(1) + vely
          r%xinfo(i)%dz(1) = r%xinfo(i)%dz(1) + velz
       enddo
    enddo
  end subroutine flex_AddAngularVelocity

  subroutine flex_savebinary( flex, mol, site, file, mode )
    type(sFlex), intent(in) :: flex
    type(sMol), intent(in)   :: mol
    type(sSite), intent(IN)  :: site
    integer, intent(in)      :: mode, file
    write(file) "@FIXC"
    write(file) mol%isFixed
    if ( mode .eq. FLEX_APC5 ) then
       call Flex_WriteBinaryAPC5( flex, mol, site, file )
    endif
  end subroutine flex_savebinary

  subroutine Flex_WriteBinaryAPC5(r,m,s,file)
    type(sFlex),intent(IN) :: r
    type(sMol),intent(IN) :: m
    type(sSite),intent(IN) :: s
    integer,intent(IN) :: file
    integer :: n,i,mol,site,j
    character(len=8) :: id
    id = m%id(1:8)
    n=m%nmol
    write(file) "@ID08"
    write(file) id
    write(file)"@APC5"
    write(file) n
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          j = i + m%offset
          write(file) s%x(j),s%y(j),s%z(j),r%mass(i)
       enddo
    enddo
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          call flex1_writebinary(r%xinfo(i),file)
       enddo
    enddo
  end subroutine Flex_WriteBinaryAPC5

  subroutine Flex_ReadBinaryAPC5(r,m,s,file)
    type(sFlex),intent(INOUT) :: r
    type(sMol),intent(INOUT)  :: m
    type(sSite),intent(INOUT) :: s
    integer,intent(IN) :: file
    integer :: i,mol,site,j,n
    read(file) n
    !
    !本当はここでやるべきことではない。program中で明示的に確保すべき
    !
    call Mol_Allocate(m,s,n,FLEX_MODE)
    call Flex_Allocate(r,m)
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          j = i + m%offset
          read(file) s%x(j),s%y(j),s%z(j),r%mass(i)
       enddo
    enddo
    do site=1, m%nsite-1
       do mol=1, m%nmol
          i = (mol-1)*m%nsite + site
          call flex1_readbinary(r%xinfo(i),file)
       enddo
    enddo
  end subroutine Flex_ReadBinaryAPC5

  subroutine flex1_writebinary(flex1, file)
    type(sFlex1), intent(IN) :: flex1
    integer, intent(IN)      :: file
    ! local
    integer :: i
    write(file) (flex1%dx(i),flex1%dy(i),flex1%dz(i),i=1,5)
  end subroutine flex1_writebinary

  subroutine flex1_readbinary(flex1, file)
    type(sFlex1), intent(INOUT) :: flex1
    integer, intent(IN)         :: file
    ! local
    integer :: i
    read(file) (flex1%dx(i),flex1%dy(i),flex1%dz(i),i=1,5)
  end subroutine flex1_readbinary

  subroutine Flex_WriteAR3A(r,m,s,file)
    type(sFlex),intent(IN) :: r
    type(sMol),intent(IN) :: m
    type(sSite), intent(IN) :: s
    integer,intent(IN) :: file
    integer :: n,i,mol,site,j
    n=m%nmol
    if ( m%nsite .ne. 2 ) then
       write(STDERR,*) "Warning: Monatomic molecule only."
    endif
    write(file,'("@ID08")')
    write(file,'(a8)') m%id
    write(file,'("@AR3A")')
    write(file,*) n
    do mol=1, m%nmol
       do site=1, m%nsite-1
          i = (mol-1)*m%nsite + site
          j = i + m%offset
          write(file,'(3(e17.10,1x))') s%x(j),s%y(j),s%z(j)
       enddo
    enddo
  end subroutine Flex_WriteAR3A

  subroutine Flex_WriteMDVW(r,m,s,file)
    type(sFlex),intent(IN) :: r
    type(sMol),intent(IN) :: m
    type(sSite),intent(IN) :: s
    integer,intent(IN) :: file
    integer :: n,i,mol,site,j,k

    real(kind=8) :: dx,dy,dz,dd
    do mol=1, m%nmol
       do site=1, m%nsite-1
          i = (mol-1)*m%nsite + site
          j = i + m%offset
          if ( m%name(site)(1:1) .ne. " " ) then
             write(file,'(a8,3(1x,e17.9E3))') m%name(site),s%x(j),s%y(j),s%z(j)
          endif
       enddo
       !
       !for debug.
       !
       k = (mol-1)*m%nsite + m%offset
       do i=1,m%nsite-1
          do j=i+1,m%nsite-1
             dx = s%x(k+i) - s%x(k+j)
             dy = s%y(k+i) - s%y(k+j)
             dz = s%z(k+i) - s%z(k+j)
             write(file,*) mol,i,j,dx**2+dy**2+dz**2
          enddo
       enddo
    enddo
  end subroutine Flex_WriteMDVW

  subroutine flex_save( flex, mol, site, file, mode )
    type(sFlex), intent(in) :: flex
    type(sMol), intent(in)   :: mol
    type(sSite), intent(in)  :: site
    integer, intent(in)      :: mode, file
    if ( mode .eq. FLEX_AR3A ) then
       write(file,'("@FIXC")')
       write(file,*) mol%isFixed
       call Flex_WriteAR3A( flex, mol, site, file )
    else if ( mode .eq. FLEX_MDVW ) then
       call Flex_WriteMDVW( flex, mol, site, file )
    endif
  end subroutine flex_save

  subroutine flex_setcom( flex, mol, site )
    use site_module
    type(sFlex), intent(IN) :: flex
    type(sSite), intent(INOUT) :: site
    type(sMol),  intent(IN) :: mol
    !
    !local
    !
    integer :: nmol, i,j,k
    real(kind=8) :: totalmass, mass1, comx,comy,comz
    nmol = mol%nmol
    k = 0
    do i=1,nmol
       totalmass = 0d0
       comx      = 0d0
       comy      = 0d0
       comz      = 0d0
       do j=1,mol%nsite - 1
          k = k + 1
          mass1 = flex%mass(k)
          comx = comx + site%x(k + mol%offset) * mass1
          comy = comy + site%y(k + mol%offset) * mass1
          comz = comz + site%z(k + mol%offset) * mass1
          totalmass = totalmass + mass1
       enddo
       k = k + 1
       !
       ! last site is Center of Mass
       !
       site%x(k + mol%offset) = comx / totalmass
       site%y(k + mol%offset) = comy / totalmass
       site%z(k + mol%offset) = comz / totalmass
    enddo
  end subroutine flex_setcom
  
  subroutine flex_distributeforce( flex, mol, site )
    use site_module
    type(sFlex), intent(IN) :: flex
    type(sMol),  intent(IN) :: mol
    type(sSite), intent(INOUT) :: site
    !
    !local
    !
    integer :: nmol, i,j,k
    real(kind=8) :: totalmass, mass1, comfx, comfy, comfz
    nmol = mol%nmol
    do i=1,nmol
       totalmass = 0d0
       !
       !まず各分子の質量を求める。
       !
       do j=1,mol%nsite - 1
          !k = mol%offset + (i-1)*mol%nsite + j
          k = (i-1)*mol%nsite + j
          totalmass = totalmass + flex%mass(k)
       enddo
       k = mol%offset + i*mol%nsite
       comfx = site%fx(k)
       comfy = site%fy(k)
       comfz = site%fz(k)
       !
       !カットオフ作用点に加わった力を、質量比で分配する。
       !
       do j=1,mol%nsite - 1
          !k = mol%offset + (i-1)*mol%nsite + j
          k = (i-1)*mol%nsite + j
          mass1 = flex%mass(k) / totalmass
          site%fx( k + mol%offset ) = site%fx( k + mol%offset ) + comfx * mass1
          site%fy( k + mol%offset ) = site%fy( k + mol%offset ) + comfy * mass1
          site%fz( k + mol%offset ) = site%fz( k + mol%offset ) + comfz * mass1
       enddo
    enddo
  end subroutine flex_distributeforce

end module flex_module

! -*- f90 -*-
#undef RIGIDDEBUG1
!
!分子を構成するサイトのうち、最後のサイト(通常は重心)を、カットオフの
!基準とする。順序に注意
!
module rigid_module
  use common_module
  use mol_module
  use box_module
  use vector_module
  implicit none
  integer, parameter :: RIGID_NX4A=0, RIGID_WTG2=1, RIGID_MDVW=2, RIGID_WTG3=3
  type sRigid1
     sequence
     !real(kind=8) :: comx,comy,com%vec(3)
     type(vector3) :: com
     real(kind=8),dimension(5) :: dx,dy,dz
     real(kind=8),dimension(5) :: wx,wy,wz
     type(vector4),dimension(5) :: quat
     real(kind=8) :: t11,t12,t13,t21,t22,t23,t31,t32,t33
     real(kind=8) :: forcex,forcey,forcez
     real(kind=8) :: torquex,torquey,torquez
     real(kind=8),dimension(MAXsite) :: intrax,intray,intraz
  end type sRigid1
  
  type sRigid
     sequence
     real(kind=8) Ixxi,Iyyi,Izzi,massi
     real(kind=8) Ixx,Iyy,Izz,mass
     real(kind=8),dimension(MAXsite) ::  Molx,Moly,Molz
     type(sRigid1),dimension(MAXMOL) :: mol
  end type sRigid
  
  interface save
     module procedure rigid_save
  end interface

  interface save_binary
     module procedure rigid_savebinary
  end interface

contains

  subroutine Rigid_relocate(r,m,b)
    type(sRigid),intent(inout) :: r
    type(sMol),intent(in) :: m
    type(sBox),intent(in) :: b
    type(vector3)         :: cell
    !real(kind=8)          :: ox,oy,oz
    integer :: i
    do i=1,m%nmol
       call box_renormalize(b,r%mol(i)%com,cell)
    enddo
  end subroutine Rigid_relocate
  
  subroutine rigid_allocate(r,n)
    type(sRigid),intent(inout) :: r
    integer,intent(in) :: n
  end subroutine rigid_allocate
  
  subroutine rigid_setmatrix(r,m)
    type(sRigid),intent(inout) :: r
    type(sMol),intent(in) :: m
    real(kind=8) :: qa,qb,qc,qd
    integer i
    !OCL NOVREC
    do i=1,m%nmol
       qa = r%mol(i)%quat(1)%vec(1)
       qb = r%mol(i)%quat(1)%vec(2)
       qc = r%mol(i)%quat(1)%vec(3)
       qd = r%mol(i)%quat(1)%vec(4)
       r%mol(i)%t11=(qa*qa+qb*qb-(qc*qc+qd*qd))
       r%mol(i)%t12=-2.0d0*(qa*qd+qb*qc)
       r%mol(i)%t13=2.0d0*(qb*qd-qa*qc)
       r%mol(i)%t21=2.0d0*(qa*qd-qb*qc)
       r%mol(i)%t22=qa*qa+qc*qc-(qb*qb+qd*qd)
       r%mol(i)%t23=-2.0d0*(qa*qb+qc*qd)
       r%mol(i)%t31=2.0d0*(qa*qc+qb*qd)
       r%mol(i)%t32=2.0d0*(qa*qb-qc*qd)
       r%mol(i)%t33=qa*qa+qd*qd-(qb*qb+qc*qc)
     enddo
  end subroutine rigid_setmatrix
  
  !
  !実際には力は全然関係ない。rigid_setsitepositionなどとすべき。
  !
  subroutine rigid_setsiteposition0(r,m,s)
    use site_module
    type(sRigid),intent(inout) :: r
    type(sSite),intent(out) :: s
    type(sMol),intent(in) :: m
    integer i,j,k!,k0
#ifdef VPOPTIMIZE
    !OCL NOVREC
    do k=1,m%nmol*m%nsite
       j=m%lvsite(k)
       i=m%lvmol(k)
       r%mol(i)%intrax(j)=r%mol(i)%t11*r%molx(j)+r%mol(i)%t12*r&
            & %moly(j)+r%mol(i)%t13*r%molz(j)
       r%mol(i)%intray(j)=r%mol(i)%t21*r%molx(j)+r%mol(i)%t22*r&
            & %moly(j)+r%mol(i)%t23*r%molz(j)
       r%mol(i)%intraz(j)=r%mol(i)%t31*r%molx(j)+r%mol(i)%t32*r&
            & %moly(j)+r%mol(i)%t33*r%molz(j)
       s%x(k+m%offset)=r%mol(i)%com%vec(1)+r%mol(i)%intrax(j)
       s%y(k+m%offset)=r%mol(i)%com%vec(2)+r%mol(i)%intray(j)
       s%z(k+m%offset)=r%mol(i)%com%vec(3)+r%mol(i)%intraz(j)
    enddo
#else
    do j=1,m%nsite
       do i=1,m%nmol
          r%mol(i)%intrax(j)=r%mol(i)%t11*r%molx(j)+r%mol(i)%t12*r&
               & %moly(j)+r%mol(i)%t13*r%molz(j)
          r%mol(i)%intray(j)=r%mol(i)%t21*r%molx(j)+r%mol(i)%t22*r&
               & %moly(j)+r%mol(i)%t23*r%molz(j)
          r%mol(i)%intraz(j)=r%mol(i)%t31*r%molx(j)+r%mol(i)%t32*r&
               & %moly(j)+r%mol(i)%t33*r%molz(j)
       enddo
    enddo
    k=m%offset
    do i=1,m%nmol
       do j=1,m%nsite
          k=k+1
          s%x(k)=r%mol(i)%com%vec(1)+r%mol(i)%intrax(j)
          s%y(k)=r%mol(i)%com%vec(2)+r%mol(i)%intray(j)
          s%z(k)=r%mol(i)%com%vec(3)+r%mol(i)%intraz(j)
       enddo
    enddo
#endif
    !test output for yaplot
    !write(10,6)
    !k=m%offset
    !do i=1,m%nmol
    !   k=k+1
    !   k0=k
    !   do j=2,m%nsite
    !      k=k+1
    !      write(10,1) s%x(k0),s%y(k0),s%z(k0),s%x(k),s%y(k),s%z(k)
    !   enddo
    !enddo
    !1 format("l",6(f10.5))
    !6 format("@ 6")
    return
  end subroutine rigid_setsiteposition0

  subroutine rigid_setsiteposition(r,m,s)
    use site_module
    type(sRigid),intent(inout) :: r
    type(sSite),intent(out) :: s
    type(sMol),intent(in) :: m
    call rigid_setmatrix(r,m)
    call rigid_setsiteposition0(r,m,s)
  end subroutine rigid_setsiteposition
  
  
#if defined(PVM) || defined(MPI)
  subroutine rigid_collectforce(NPROCS,MYRANK,r,m,s)
#else
  subroutine rigid_collectforce(r,m,s)
#endif
    use site_module
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
    type(sRigid) :: r
    type(sSite) :: s
    type(sMol) :: m
    integer :: i,j,k
#if defined(VPOPTIMIZE) || defined(MPI) || defined(PVM)
    real(kind=8),dimension(MAXmol) :: forcex,forcey,forcez,torquex&
         & ,torquey,torquez
#ifdef MPI
    real(kind=8),dimension(MAXmol) :: tx,ty,tz
#endif
    forcex(:)=0d0
    forcey(:)=0d0
    forcez(:)=0d0
    torquex(:)=0d0
    torquey(:)=0d0
    torquez(:)=0d0
    do j=1,m%nsite
       !OCL VECTOR,NOVREC
       do i=1,m%nmol
          k=(i-1)*m%nsite+j+m%offset
          forcex(i) = forcex(i)+s%fx(k)
          forcey(i) = forcey(i)+s%fy(k)
          forcez(i) = forcez(i)+s%fz(k)
          torquex(i) = torquex(i) + r%mol(i)%intray(j)*s%fz(k) - r&
               & %mol(i)%intraz(j)*s%fy(k) 
          torquey(i) = torquey(i) + r%mol(i)%intraz(j)*s%fx(k) - r&
               & %mol(i)%intrax(j)*s%fz(k) 
          torquez(i) = torquez(i) + r%mol(i)%intrax(j)*s%fy(k) - r&
               & %mol(i)%intray(j)*s%fx(k) 
       enddo
    enddo
#ifdef MPI
    do i=1,m%nmol
       tx(i)=r%mol(i)%t11*torquex(i)+r%mol(i)%t21*torquey(i)+r&
            & %mol(i)%t31*torquez(i)
       ty(i)=r%mol(i)%t12*torquex(i)+r%mol(i)%t22*torquey(i)+r&
            & %mol(i)%t32*torquez(i)
       tz(i)=r%mol(i)%t13*torquex(i)+r%mol(i)%t23*torquey(i)+r&
            & %mol(i)%t33*torquez(i)
    enddo
    call MPI_ALLREDUCE(forcex,r%mol(i)%forcex,m%nmol&
         & ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(forcey,r%mol(i)%forcey,m%nmol&
         & ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(forcez,r%mol(i)%forcez,m%nmol&
         & ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(tx,r%mol(i)%torquex,m%nmol&
         & ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(ty,r%mol(i)%torquey,m%nmol&
         & ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
    call MPI_ALLREDUCE(tz,r%mol(i)%torquez,m%nmol&
         & ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
#else /*MPI*/
    !OCL VECTOR,NOVREC
    do i=1,m%nmol
       r%mol(i)%forcex=forcex(i)
       r%mol(i)%forcey=forcey(i)
       r%mol(i)%forcez=forcez(i)
       r%mol(i)%torquex = r%mol(i)%t11*torquex(i)+r%mol(i)%t21&
            & *torquey(i)+r%mol(i)%t31*torquez(i)
       r%mol(i)%torquey = r%mol(i)%t12*torquex(i)+r%mol(i)%t22&
            & *torquey(i)+r%mol(i)%t32*torquez(i)
       r%mol(i)%torquez = r%mol(i)%t13*torquex(i)+r%mol(i)%t23&
            & *torquey(i)+r%mol(i)%t33*torquez(i)
       !if(i==1)write(6,*) "TX",forcex,torquex,torquey,torquez,r&
       !     & %t11,r%mol(i)%t21,r%mol(i)%t31
    enddo
#ifdef PVM
    call PVMFreduce(PvmSum,r%mol(i)%forcex,m%nmol,PVM_DOUBLE,0&
         & ,"WaterGrp",0,IERR)
    call PVMFreduce(PvmSum,r%mol(i)%forcey,m%nmol,PVM_DOUBLE,0&
         & ,"WaterGrp",0,IERR)
    call PVMFreduce(PvmSum,r%mol(i)%forcez,m%nmol,PVM_DOUBLE,0&
         & ,"WaterGrp",0,IERR)
    call PVMFreduce(PvmSum,r%mol(i)%torquex,m%nmol,PVM_DOUBLE,0&
         & ,"WaterGrp",0,IERR)
    call PVMFreduce(PvmSum,r%mol(i)%torquey,m%nmol,PVM_DOUBLE,0&
         & ,"WaterGrp",0,IERR)
    call PVMFreduce(PvmSum,r%mol(i)%torquez,m%nmol,PVM_DOUBLE,0&
         & ,"WaterGrp",0,IERR)
    msgtag=5
    if(MYRANK==0)then
       call PVMFparent(tid)
       call PVMFinitsend(PVMDEFAULT,IERR)
       call PVMFpack(REAL8,r%mol(i)%forcex,m%nmol,1,IERR)
       call PVMFpack(REAL8,r%mol(i)%forcey,m%nmol,1,IERR)
       call PVMFpack(REAL8,r%mol(i)%forcez,m%nmol,1,IERR)
       call PVMFpack(REAL8,r%mol(i)%torquex,m%nmol,1,IERR)
       call PVMFpack(REAL8,r%mol(i)%torquey,m%nmol,1,IERR)
       call PVMFpack(REAL8,r%mol(i)%torquez,m%nmol,1,IERR)
       call PVMFbcast("WaterGrp",msgtag,IERR)
    else
       call PVMFparent(tid)
       call PVMFrecv(tid,msgtag,IERR)
       call PVMFunpack(REAL8,r%mol(i)%forcex,m%nmol,1,IERR)
       call PVMFunpack(REAL8,r%mol(i)%forcey,m%nmol,1,IERR)
       call PVMFunpack(REAL8,r%mol(i)%forcez,m%nmol,1,IERR)
       call PVMFunpack(REAL8,r%mol(i)%torquex,m%nmol,1,IERR)
       call PVMFunpack(REAL8,r%mol(i)%torquey,m%nmol,1,IERR)
       call PVMFunpack(REAL8,r%mol(i)%torquez,m%nmol,1,IERR)
    endif
#endif  /*PVM*/
#endif  /*MPI*/
#else
    real(kind=8) forcex,forcey,forcez
    real(kind=8) torquex,torquey,torquez
    k=m%offset
    do i=1,m%nmol
       torquex=0d0
       torquey=0d0
       torquez=0d0
       forcex=0d0
       forcey=0d0
       forcez=0d0
       do j=1,m%nsite
          k=k+1
          forcex = forcex+s%fx(k)
          forcey = forcey+s%fy(k)
          forcez = forcez+s%fz(k)
          torquex = torquex + r%mol(i)%intray(j)*s%fz(k) - r%mol(i)&
               & %intraz(j)*s%fy(k) 
          torquey = torquey + r%mol(i)%intraz(j)*s%fx(k) - r%mol(i)&
               & %intrax(j)*s%fz(k) 
          torquez = torquez + r%mol(i)%intrax(j)*s%fy(k) - r%mol(i)&
               & %intray(j)*s%fx(k) 
       enddo
       r%mol(i)%forcex=forcex
       r%mol(i)%forcey=forcey
       r%mol(i)%forcez=forcez
       r%mol(i)%torquex = r%mol(i)%t11*torquex+r%mol(i)%t21*torquey&
            & +r%mol(i)%t31*torquez
       r%mol(i)%torquey = r%mol(i)%t12*torquex+r%mol(i)%t22*torquey&
            & +r%mol(i)%t32*torquez
       r%mol(i)%torquez = r%mol(i)%t13*torquex+r%mol(i)%t23*torquey&
            & +r%mol(i)%t33*torquez
       !debug output for yaplot
       !write(10,5)
       !write(10,1) r%mol(i)%comx,r%mol(i)%comy,r%mol(i)%com%vec(3),r
       ! %mol(i)%comx+r%mol(i)%forcex*1d-4,r%mol(i)%comy+r%mol(i)
       ! %forcey*1d-4,r%mol(i)%com%vec(3)+r%mol(i)%forcez*1d-4
       !5    format("@ 5")
       !1    format("l",6(f10.5))
       !if(i==1)write(6,*) "TX",forcex,torquex,torquey,torquez,r&
       !     & %t11,r%mol(i)%t21,r%mol(i)%t31
    enddo
#endif
    return
  end subroutine rigid_collectforce


  subroutine Rigid_TIP4P_Constructor2(r,hmass)
    use physconst_module
    type(sRigid),intent(INOUT) :: r
    real(kind=8),intent(IN) :: hmass
    real(kind=8) :: omass,ohlen,angle,oz,hy,hz,slide,omlen
    real(kind=8),dimension(5) :: mass
    integer :: i
    !molecular weight
    omass=16d0
    !angstrom
    ohlen=0.9572d0
    omlen=0.15d0
    !degree
    angle=104.52d0
    
    !coordinate relative to oxygen position
    hy = ohlen * dsin( angle * PI / 360d0)
    hz = ohlen * dcos( angle * PI / 360d0)
    
    !slide to the com
    slide = hz * (hmass * 2d0) / (hmass * 2d0 + omass)
    hz    = hz - slide
    oz    = -slide
    omlen = omlen - slide
    
    !oxygen
    mass(1)  = omass
    r%molx(1)=0d0
    r%moly(1)=0d0
    r%molz(1)=oz
    !charge site
    mass(2)  = 0d0
    r%molx(2)=0d0
    r%moly(2)=0d0
    r%molz(2)=omlen
    !hydrogen 1
    mass(3)  = hmass
    r%molx(3)=0d0
    r%moly(3)=hy
    r%molz(3)=hz
    !hydrogen 2
    mass(4)  = hmass
    r%molx(4)=0d0
    r%moly(4)=-hy
    r%molz(4)=hz
    !center of mass
    mass(5)  = 0d0
    r%molx(5)=0d0
    r%moly(5)=0d0
    r%molz(5)=0d0
    
    r%Ixx = 0d0
    r%Iyy = 0d0
    r%Izz = 0d0
    r%mass = 0d0
    do i=1,5
       r%Ixx = r%Ixx + mass(i)*(r%moly(i)**2 + r%molz(i)**2)
       r%Iyy = r%Iyy + mass(i)*(r%molz(i)**2 + r%molx(i)**2)
       r%Izz = r%Izz + mass(i)*(r%molx(i)**2 + r%moly(i)**2)
       r%mass = r%mass + mass(i)
    enddo
    
    r%ixxi = 1d0/r%ixx
    r%iyyi = 1d0/r%iyy
    r%izzi = 1d0/r%izz
    r%massi= 1d0/r%mass
    return
  end subroutine Rigid_TIP4P_Constructor2

  subroutine Rigid_OPLSMeOH_Constructor(r)
    type(sRigid),intent(INOUT) :: r
    real(kind=8),dimension(4) :: mass
    integer :: i
    mass(1)=15d0
    mass(2)=16d0
    mass(3)=1d0
    mass(4)=0d0
    r%molx(1)=0.0133d0
    r%molx(2)=-0.0634d0
    r%molx(3)=0.8152d0
    r%molx(4)=0d0
    r%moly(1)=0.7666d0
    r%moly(2)=-0.6559d0
    r%moly(3)=-1.0042d0
    r%moly(4)=0d0
    r%molz(1)=0d0
    r%molz(2)=0d0
    r%molz(3)=0d0
    r%molz(4)=0d0
    r%Ixx = 0d0
    r%Iyy = 0d0
    r%Izz = 0d0
    r%mass = 0d0
    do i=1,4
       r%Ixx = r%Ixx + mass(i)*(r%moly(i)**2 + r%molz(i)**2)
       r%Iyy = r%Iyy + mass(i)*(r%molz(i)**2 + r%molx(i)**2)
       r%Izz = r%Izz + mass(i)*(r%molx(i)**2 + r%moly(i)**2)
       r%mass = r%mass + mass(i)
    enddo
    r%ixxi = 1d0/r%ixx
    r%iyyi = 1d0/r%iyy
    r%izzi = 1d0/r%izz
    r%massi= 1d0/r%mass
    return
  end subroutine Rigid_OPLSMeOH_Constructor


  subroutine Rigid_TIP5P_Constructor(r)
    use physconst_module
    type(sRigid),intent(INOUT) :: r
    real(kind=8),parameter :: theta=104.52d0*PI/180d0,phi=109.47d0*PI/180d0
    real(kind=8),parameter :: l1=0.9572d0,l2=0.70d0
    !comparison test with TIP4P
    !real(kind=8),parameter :: theta=104.52*PIright/180d0,phi=0.
    ! *PIright/180d0
    !real(kind=8),parameter :: l1=0.9572,l2=-0.15
    real(kind=8),parameter :: hmass=1d0,omass=16d0
    real(kind=8) :: zshift
    integer :: i
    r%molx(1)=0d0
    r%moly(1)=0d0
    r%molz(1)=0d0
    r%molx(2)=0d0
    r%moly(2)=l1*dsin(0.5d0*theta)
    r%molz(2)=l1*dcos(0.5d0*theta)
    r%molx(3)=0d0
    r%moly(3)=-l1*dsin(0.5d0*theta)
    r%molz(3)=l1*dcos(0.5d0*theta)
    r%molx(4)=l2*dsin(0.5d0*phi)
    r%moly(4)=0d0
    r%molz(4)=-l2*dcos(0.5d0*phi)
    r%molx(5)=-l2*dsin(0.5d0*phi)
    r%moly(5)=0d0
    r%molz(5)=-l2*dcos(0.5d0*phi)
    zshift=r%molz(2)*2d0*hmass/(2d0*hmass+omass)
    do i=1,5
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
    r%ixxi = 0d0
    r%iyyi = 0d0
    r%izzi = 0d0
    if ( r%ixx /= 0d0 ) r%ixxi = 1d0/r%ixx
    if ( r%iyy /= 0d0 ) r%iyyi = 1d0/r%iyy
    if ( r%izz /= 0d0 ) r%izzi = 1d0/r%izz
    r%mass = 2d0*hmass + omass
    r%massi= 1d0/r%mass
  end subroutine Rigid_TIP5P_Constructor


  subroutine rigid_predict(r,m)
    type(sMol),intent(in) :: m
    type(sRigid),intent(INOUT) :: r
    integer i
    real(kind=8) re
    do i=1,m%nmol
       r%mol(i)%com%vec(1) = r%mol(i)%com%vec(1)+r%mol(i)%dx(1)+r%mol(i)%dx(2)+r&
            & %mol(i)%dx(3)+r%mol(i)%dx(4)+r%mol(i)%dx(5)
       r%mol(i)%dx(1) = r%mol(i)%dx(1)+2d0*r%mol(i)%dx(2)+3d0*r&
            & %mol(i)%dx(3)+4d0*r%mol(i)%dx(4)+5d0*r%mol(i)%dx(5)
       r%mol(i)%dx(2) = r%mol(i)%dx(2)+3d0*r%mol(i)%dx(3)+6d0*r&
            & %mol(i)%dx(4)+10d0*r%mol(i)%dx(5)
       r%mol(i)%dx(3) = r%mol(i)%dx(3)+4d0*r%mol(i)%dx(4)+10d0*r&
            & %mol(i)%dx(5)
       r%mol(i)%dx(4) = r%mol(i)%dx(4)+5d0*r%mol(i)%dx(5)
       r%mol(i)%com%vec(2) = r%mol(i)%com%vec(2)+r%mol(i)%dy(1)+r%mol(i)%dy(2)+r&
            & %mol(i)%dy(3)+r%mol(i)%dy(4)+r%mol(i)%dy(5)
       r%mol(i)%dy(1) = r%mol(i)%dy(1)+2d0*r%mol(i)%dy(2)+3d0*r&
            & %mol(i)%dy(3)+4d0*r%mol(i)%dy(4)+5d0*r%mol(i)%dy(5)
       r%mol(i)%dy(2) = r%mol(i)%dy(2)+3d0*r%mol(i)%dy(3)+6d0*r&
            & %mol(i)%dy(4)+10d0*r%mol(i)%dy(5)
       r%mol(i)%dy(3) = r%mol(i)%dy(3)+4d0*r%mol(i)%dy(4)+10d0*r&
            & %mol(i)%dy(5)
       r%mol(i)%dy(4) = r%mol(i)%dy(4)+5d0*r%mol(i)%dy(5)
       r%mol(i)%com%vec(3) = r%mol(i)%com%vec(3)+r%mol(i)%dz(1)+r%mol(i)%dz(2)+r&
            & %mol(i)%dz(3)+r%mol(i)%dz(4)+r%mol(i)%dz(5)
       r%mol(i)%dz(1) = r%mol(i)%dz(1)+2d0*r%mol(i)%dz(2)+3d0*r&
            & %mol(i)%dz(3)+4d0*r%mol(i)%dz(4)+5d0*r%mol(i)%dz(5)
       r%mol(i)%dz(2) = r%mol(i)%dz(2)+3d0*r%mol(i)%dz(3)+6d0*r&
            & %mol(i)%dz(4)+10d0*r%mol(i)%dz(5)
       r%mol(i)%dz(3) = r%mol(i)%dz(3)+4d0*r%mol(i)%dz(4)+10d0*r&
            & %mol(i)%dz(5)
       r%mol(i)%dz(4) = r%mol(i)%dz(4)+5d0*r%mol(i)%dz(5)
       r%mol(i)%wx(1) = r%mol(i)%wx(1)+r%mol(i)%wx(2)+r%mol(i)%wx(3)&
            & +r%mol(i)%wx(4)+r%mol(i)%wx(5)
       r%mol(i)%wx(2) = r%mol(i)%wx(2)+2d0*r%mol(i)%wx(3)+3d0*r&
            & %mol(i)%wx(4)+4d0*r%mol(i)%wx(5)
       r%mol(i)%wx(3) = r%mol(i)%wx(3)+3d0*r%mol(i)%wx(4)+6d0*r&
            & %mol(i)%wx(5)
       r%mol(i)%wx(4) = r%mol(i)%wx(4)+4d0*r%mol(i)%wx(5)
       r%mol(i)%wy(1) = r%mol(i)%wy(1)+r%mol(i)%wy(2)+r%mol(i)%wy(3)&
            & +r%mol(i)%wy(4)+r%mol(i)%wy(5)
       r%mol(i)%wy(2) = r%mol(i)%wy(2)+2d0*r%mol(i)%wy(3)+3d0*r&
            & %mol(i)%wy(4)+4d0*r%mol(i)%wy(5)
       r%mol(i)%wy(3) = r%mol(i)%wy(3)+3d0*r%mol(i)%wy(4)+6d0*r&
            & %mol(i)%wy(5)
       r%mol(i)%wy(4) = r%mol(i)%wy(4)+4d0*r%mol(i)%wy(5)
       r%mol(i)%wz(1) = r%mol(i)%wz(1)+r%mol(i)%wz(2)+r%mol(i)%wz(3)&
            & +r%mol(i)%wz(4)+r%mol(i)%wz(5)
       r%mol(i)%wz(2) = r%mol(i)%wz(2)+2d0*r%mol(i)%wz(3)+3d0*r&
            & %mol(i)%wz(4)+4d0*r%mol(i)%wz(5)
       r%mol(i)%wz(3) = r%mol(i)%wz(3)+3d0*r%mol(i)%wz(4)+6d0*r&
            & %mol(i)%wz(5)
       r%mol(i)%wz(4) = r%mol(i)%wz(4)+4d0*r%mol(i)%wz(5)
       r%mol(i)%quat(1)%vec(1) = r%mol(i)%quat(1)%vec(1)+r%mol(i)%quat(2)%vec(1)+r%mol(i)%quat(3)%vec(1)&
            & +r%mol(i)%quat(4)%vec(1)+r%mol(i)%quat(5)%vec(1)
       r%mol(i)%quat(2)%vec(1) = r%mol(i)%quat(2)%vec(1)+2d0*r%mol(i)%quat(3)%vec(1)+3d0*r&
            & %mol(i)%quat(4)%vec(1)+4d0*r%mol(i)%quat(5)%vec(1)
       r%mol(i)%quat(1)%vec(2) = r%mol(i)%quat(1)%vec(2)+r%mol(i)%quat(2)%vec(2)+r%mol(i)%quat(3)%vec(2)&
            & +r%mol(i)%quat(4)%vec(2)+r%mol(i)%quat(5)%vec(2)
       r%mol(i)%quat(2)%vec(2) = r%mol(i)%quat(2)%vec(2)+2d0*r%mol(i)%quat(3)%vec(2)+3d0*r&
            & %mol(i)%quat(4)%vec(2)+4d0*r%mol(i)%quat(5)%vec(2)
       r%mol(i)%quat(1)%vec(3) = r%mol(i)%quat(1)%vec(3)+r%mol(i)%quat(2)%vec(3)+r%mol(i)%quat(3)%vec(3)&
            & +r%mol(i)%quat(4)%vec(3)+r%mol(i)%quat(5)%vec(3)
       r%mol(i)%quat(2)%vec(3) = r%mol(i)%quat(2)%vec(3)+2d0*r%mol(i)%quat(3)%vec(3)+3d0*r&
            & %mol(i)%quat(4)%vec(3)+4d0*r%mol(i)%quat(5)%vec(3)
       r%mol(i)%quat(1)%vec(4) = r%mol(i)%quat(1)%vec(4)+r%mol(i)%quat(2)%vec(4)+r%mol(i)%quat(3)%vec(4)&
            & +r%mol(i)%quat(4)%vec(4)+r%mol(i)%quat(5)%vec(4)
       r%mol(i)%quat(2)%vec(4) = r%mol(i)%quat(2)%vec(4)+2d0*r%mol(i)%quat(3)%vec(4)+3d0*r&
            & %mol(i)%quat(4)%vec(4)+4d0*r%mol(i)%quat(5)%vec(4)
       re =1d0/dsqrt(r%mol(i)%quat(1)%vec(1)**2+r%mol(i)%quat(1)%vec(2)**2+r%mol(i)&
            & %quat(1)%vec(3)**2+r%mol(i)%quat(1)%vec(4)**2)
       r%mol(i)%quat(1)%vec(1)=r%mol(i)%quat(1)%vec(1)*re
       r%mol(i)%quat(1)%vec(2)=r%mol(i)%quat(1)%vec(2)*re
       r%mol(i)%quat(1)%vec(3)=r%mol(i)%quat(1)%vec(3)*re
       r%mol(i)%quat(1)%vec(4)=r%mol(i)%quat(1)%vec(4)*re
    enddo
  !write(7,*) r%mol(i)%com%vec(1)(1),r%mol(i)%com%vec(2)(1),r%mol(i)%com%vec(3)(1)
    return
  end subroutine rigid_predict

  subroutine rigid_propagatemomenta(r,m,slice)
    type(sMol),intent(in) :: m
    type(sRigid),intent(INOUT) :: r
    real(kind=8) :: slice
    integer i
    do i=1,m%nmol
       r%mol(i)%dx(1) = r%mol(i)%dx(1) + 2d0 * r%mol(i)%dx(2) * slice
       r%mol(i)%dy(1) = r%mol(i)%dy(1) + 2d0 * r%mol(i)%dy(2) * slice
       r%mol(i)%dz(1) = r%mol(i)%dz(1) + 2d0 * r%mol(i)%dz(2) * slice
       r%mol(i)%wx(1) = r%mol(i)%wx(1) + r%mol(i)%wx(2) * slice
       r%mol(i)%wy(1) = r%mol(i)%wy(1) + r%mol(i)%wy(2) * slice
       r%mol(i)%wz(1) = r%mol(i)%wz(1) + r%mol(i)%wz(2) * slice
    enddo
  end subroutine rigid_propagatemomenta

  subroutine Rigid_Ek(r,m,t,ekt,ekr)
    use time_module
    type(sMol),intent(IN) :: m
    type(sRigid),intent(IN) :: r
    type(sTime),intent(IN) :: t
    real(kind=8),intent(out) :: ekt,ekr
    real(kind=8) ekrx,ekry,ekrz
    integer i
    ekrx=0d0
    ekry=0d0
    ekrz=0d0
    ekt=0
    do i=1,m%Nmol
       ekt = ekt + r%mol(i)%dx(1)**2+r%mol(i)%dy(1)**2+r%mol(i)%dz(1)**2
       ekrx = ekrx + r%mol(i)%wx(1)**2
       ekry = ekry + r%mol(i)%wy(1)**2
       ekrz = ekrz + r%mol(i)%wz(1)**2
    enddo
    ekr = (ekrx*r%Ixx+ekry*r%Iyy+ekrz*r%Izz)*0.5d0
    ekt = ekt * r%Mass*0.5d0*(1d0/(t%Dt**2))
    return
  end subroutine Rigid_Ek

  subroutine Rigid_PressureTensor(r,m,t,vrmat)
    use time_module
    type(sMol),intent(IN) :: m
    type(sRigid),intent(IN) :: r
    type(sTime),intent(IN) :: t
    real(kind=8),intent(out) :: vrmat(3,3)
    integer :: i,j,k
    vrmat(:,:)=0d0
    do i=1,m%Nmol
       vrmat(1,1) = vrmat(1,1) + r%mol(i)%dx(1) * r%mol(i)%dx(1)
       vrmat(1,2) = vrmat(1,2) + r%mol(i)%dx(1) * r%mol(i)%dy(1)
       vrmat(1,3) = vrmat(1,3) + r%mol(i)%dx(1) * r%mol(i)%dz(1)
       vrmat(2,1) = vrmat(2,1) + r%mol(i)%dy(1) * r%mol(i)%dx(1)
       vrmat(2,2) = vrmat(2,2) + r%mol(i)%dy(1) * r%mol(i)%dy(1)
       vrmat(2,3) = vrmat(2,3) + r%mol(i)%dy(1) * r%mol(i)%dz(1)
       vrmat(3,1) = vrmat(3,1) + r%mol(i)%dz(1) * r%mol(i)%dx(1)
       vrmat(3,2) = vrmat(3,2) + r%mol(i)%dz(1) * r%mol(i)%dy(1)
       vrmat(3,3) = vrmat(3,3) + r%mol(i)%dz(1) * r%mol(i)%dz(1)
    enddo
    vrmat(:,:) = vrmat(:,:) * r%Mass*(1d0/(t%Dt**2))
  end subroutine Rigid_PressureTensor

  subroutine Rigid_ScaleVelocity(r,m,ratio)
    use time_module
    type(sMol),intent(IN) :: m
    type(sRigid),intent(INOUT) :: r
    real(kind=8),intent(IN) :: ratio
    r%mol(1:m%nmol)%dx(1)= r%mol(1:m%nmol)%dx(1)*ratio
    r%mol(1:m%nmol)%dy(1)= r%mol(1:m%nmol)%dy(1)*ratio
    r%mol(1:m%nmol)%dz(1)= r%mol(1:m%nmol)%dz(1)*ratio
    r%mol(1:m%nmol)%wx(1)= r%mol(1:m%nmol)%wx(1)*ratio
    r%mol(1:m%nmol)%wy(1)= r%mol(1:m%nmol)%wy(1)*ratio
    r%mol(1:m%nmol)%wz(1)= r%mol(1:m%nmol)%wz(1)*ratio
    return
  end subroutine Rigid_ScaleVelocity
  
  subroutine Rigid_ReadBinaryWTG2(r,m,si,file)
    use site_module
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(INOUT) :: m
    type(sSite),intent(INOUT) :: si
    integer,intent(IN) :: file
    integer :: n,i,s
    read(file) n
    !
    !本当はここでやるべきことではない。program中で明示的に確保すべき
    !
    call Mol_Allocate(m,si,n,RIGID_MODE)
    call Rigid_Allocate(r,n)
    read(file) &
         (r%mol(i)%com%vec(1),i=1,n)&
         ,(r%mol(i)%com%vec(2),i=1,n)&
         ,(r%mol(i)%com%vec(3),i=1,n)&
         ,((r%mol(i)%dx(s),i=1,n),s=1,5)&
         ,((r%mol(i)%dy(s),i=1,n),s=1,5)&
         ,((r%mol(i)%dz(s),i=1,n),s=1,5)&
         ,((r%mol(i)%quat(s)%vec(1),i=1,n),s=1,5)&
         ,((r%mol(i)%quat(s)%vec(2),i=1,n),s=1,5)&
         ,((r%mol(i)%quat(s)%vec(3),i=1,n),s=1,5)&
         ,((r%mol(i)%quat(s)%vec(4),i=1,n),s=1,5)&
         ,((r%mol(i)%wx(s),i=1,n),s=1,5)&
         ,((r%mol(i)%wy(s),i=1,n),s=1,5)&
         ,((r%mol(i)%wz(s),i=1,n),s=1,5)
    return
  end subroutine Rigid_ReadBinaryWTG2
  
  subroutine Rigid_ReadNX4A(r,m,si,file)
    use site_module
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(INOUT) :: m
    type(sSite),intent(INOUT) :: si
    integer,intent(IN) :: file
    call Rigid_ReadNX4A_2(r,m,si,.true.,file)
  end subroutine Rigid_ReadNX4A
  
  subroutine Rigid_ResetDerivatives( rigid, n )
    type(sRigid),intent(INOUT) :: rigid
    integer, intent(IN) :: n
    integer :: i,s
    do i=1,n
       do s=1,5
          rigid%mol(i)%dx(s)=0
          rigid%mol(i)%dy(s)=0
          rigid%mol(i)%dz(s)=0
          rigid%mol(i)%wx(s)=0
          rigid%mol(i)%wy(s)=0
          rigid%mol(i)%wz(s)=0
          rigid%mol(i)%quat(s)%vec(1)=0
          rigid%mol(i)%quat(s)%vec(2)=0
          rigid%mol(i)%quat(s)%vec(3)=0
          rigid%mol(i)%quat(s)%vec(4)=0
       enddo
    enddo
  end subroutine Rigid_ResetDerivatives

  subroutine Rigid_ReadNX4A_2(r,m,si,dreset,file)
    use site_module
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(INOUT) :: m
    type(sSite),intent(INOUT) :: si
    logical,intent(in) :: dreset
    integer,intent(IN) :: file
    integer :: n,i,s
    read(file,*) n
    call Mol_Allocate(m,si,n,RIGID_MODE)
    call Rigid_Allocate(r,n)
    r%mol(:)%com%vec(1)=0
    r%mol(:)%com%vec(2)=0
    r%mol(:)%com%vec(3)=0
    if(dreset)then
       call Rigid_ResetDerivatives( r, n )
    endif
    call Rigid_ReadOverNX4A( r, n, file )
  end subroutine Rigid_ReadNX4A_2

  subroutine Rigid_ReadOverNX4A(r,n,file)
    type(sRigid),intent(INOUT) :: r
    integer,intent(IN) :: file
    integer :: n,i
    !
    !読みこむのはやめた。平成16年7月6日(火)
    !
    !read(file,*) n
    do i=1,n
    read(file,*) r%mol(i)%com%vec(1),r%mol(i)%com%vec(2),r%mol(i)%com%vec(3),r%mol(i)&
         %quat(1)%vec(1),r%mol(i)%quat(1)%vec(2),r%mol(i)%quat(1)&
         & %vec(3),r%mol(i)%quat(1)%vec(4)
    enddo
  end subroutine Rigid_ReadOverNX4A
  
  subroutine Rigid_ReadOverBinaryNX4A(r,n,file)
    type(sRigid),intent(INOUT) :: r
    integer,intent(IN) :: file
    integer :: n,i
    !
    !読みこむのはやめた。平成16年7月6日(火)
    !
    !read(file,*) n
    read(file) (r%mol(i)%com%vec(1),r%mol(i)%com%vec(2),r%mol(i)%com%vec(3),r%mol(i)&
         %quat(1)%vec(1),r%mol(i)%quat(1)%vec(2),r%mol(i)%quat(1)&
         & %vec(3),r%mol(i)%quat(1)%vec(4), i=1,n )
  end subroutine Rigid_ReadOverBinaryNX4A
  
  subroutine Rigid_ReadNX3A(r,m,si,file)
    use quat_module
    use site_module
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(INOUT) :: m
    type(sSite),intent(INOUT) :: si
    integer,intent(IN) :: file
    integer :: n,i,s
    real(kind=8) :: ea,eb,ec
    read(file,*) n
    call Mol_Allocate(m,si,n, RIGID_MODE )
    call Rigid_Allocate(r,n)
    r%mol(:)%com%vec(1)=0
    r%mol(:)%com%vec(2)=0
    r%mol(:)%com%vec(3)=0
    call Rigid_ResetDerivatives( r, n )
    do i=1,n
       read(file,*) r%mol(i)%com%vec(1),r%mol(i)%com%vec(2),r%mol(i)%com%vec(3),ea,eb,ec
       call abc2abcd(ea,eb,ec,r%mol(i)%quat(1)%vec(1),r%mol(i)%quat(1)%vec(2),r%mol(i)&
            & %quat(1)%vec(3),r%mol(i)%quat(1)%vec(4))
    enddo
  end subroutine Rigid_ReadNX3A

!Give compat for NPT2
  subroutine Rigid_ReadBinaryNTNK(r,m,si,file)
    use site_module
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(INOUT) :: m
    type(sSite),intent(INOUT) :: si
    integer,intent(IN) :: file
    real(kind=8) :: dum
    integer :: n,i,s
    call Rigid_Allocate(r,MAXmol)
    READ(file) N,(r%mol(i)%com%vec(1),r%mol(i)%com%vec(2),r%mol(i)%com%vec(3),r%mol(i)%dx(1),r&
       & %mol(i)%dy(1),r%mol(i)%dz(1),r%mol(i)%dx(2),r%mol(i)%dy(2),r&
       & %mol(i)%dz(2),r%mol(i)%dx(3),r%mol(i)%dy(3),r%mol(i)%dz(3),r&
       & %mol(i)%dx(4),r%mol(i)%dy(4),r%mol(i)%dz(4),r%mol(i)%dx(5),r&
       & %mol(i)%dy(5),r%mol(i)%dz(5),r%mol(i)%quat(1)%vec(1),r%mol(i)%quat(1)%vec(2),r&
       & %mol(i)%quat(1)%vec(3),r%mol(i)%quat(1)%vec(4),r%mol(i)%quat(2)%vec(1),r%mol(i)%quat(2)%vec(2),r&
       & %mol(i)%quat(2)%vec(3),r%mol(i)%quat(2)%vec(4),r%mol(i)%quat(3)%vec(1),r%mol(i)%quat(3)%vec(2),r&
       & %mol(i)%quat(3)%vec(3),r%mol(i)%quat(3)%vec(4),r%mol(i)%quat(4)%vec(1),r%mol(i)%quat(4)%vec(2),r&
       & %mol(i)%quat(4)%vec(3),r%mol(i)%quat(4)%vec(4),r%mol(i)%quat(5)%vec(1),r%mol(i)%quat(5)%vec(2),r&
       & %mol(i)%quat(5)%vec(3),r%mol(i)%quat(5)%vec(4),r%mol(i)%wx(1),r%mol(i)%wy(1),r&
       & %mol(i)%wz(1),r%mol(i)%wx(2),r%mol(i)%wy(2),r%mol(i)%wz(2),r&
       & %mol(i)%wx(3),r%mol(i)%wy(3),r%mol(i)%wz(3),r%mol(i)%wx(4),r&
       & %mol(i)%wy(4),r%mol(i)%wz(4),r%mol(i)%wx(5),r%mol(i)%wy(5),r&
       & %mol(i)%wz(5),dum,dum,dum,dum,dum,dum,I=1,N)
    do i=1,n
       do s=1,5
          r%mol(i)%wx(s)=r%mol(i)%wx(s)/0021.8852721616497249235d0
          r%mol(i)%wy(s)=r%mol(i)%wy(s)/0021.8852721616497249235d0
          r%mol(i)%wz(s)=r%mol(i)%wz(s)/0021.8852721616497249235d0
       enddo
    enddo
    !write(6,*) n,r%mol(i)%wx(1,1)
    call Mol_Allocate(m,si,n, RIGID_MODE )
    return
  end subroutine Rigid_ReadBinaryNTNK
  
  subroutine Rigid_WriteBinaryWTG2(r,m,file)
    type(sRigid),intent(IN) :: r
    type(sMol),intent(IN) :: m
    integer,intent(IN) :: file
    integer :: n,i,s
    character(len=8) :: id
    id = m%id(1:8)
    n=m%nmol
    write(file) "@ID08"
    write(file) id
    write(file)"@WTG2"
    write(file) n
    write(file) &
         (r%mol(i)%com%vec(1),i=1,n)&
         ,(r%mol(i)%com%vec(2),i=1,n)&
         ,(r%mol(i)%com%vec(3),i=1,n)&
         ,((r%mol(i)%dx(s),i=1,n),s=1,5)&
         ,((r%mol(i)%dy(s),i=1,n),s=1,5)&
         ,((r%mol(i)%dz(s),i=1,n),s=1,5)&
         ,((r%mol(i)%quat(s)%vec(1),i=1,n),s=1,5)&
         ,((r%mol(i)%quat(s)%vec(2),i=1,n),s=1,5)&
         ,((r%mol(i)%quat(s)%vec(3),i=1,n),s=1,5)&
         ,((r%mol(i)%quat(s)%vec(4),i=1,n),s=1,5)&
         ,((r%mol(i)%wx(s),i=1,n),s=1,5)&
         ,((r%mol(i)%wy(s),i=1,n),s=1,5)&
         ,((r%mol(i)%wz(s),i=1,n),s=1,5)
    return
  end subroutine Rigid_WriteBinaryWTG2
  
  subroutine Rigid_WriteWTG3(r,m,t,file)
    use time_module
    type(sRigid),intent(IN) :: r
    type(sMol),intent(IN) :: m
    type(sTime),intent(IN) :: t
    integer,intent(IN) :: file
    integer :: n,i,s
    character(len=8) :: id
    id = m%id(1:8)
    n=m%nmol
    write(file,'("@ID08")')
    write(file,'(a8)') id
    write(file,'("@WTG3")')
    write(file,*) n
    do i=1,n
       call Rigid1_WriteWTG3( r%mol(i), t, file )
    enddo
  end subroutine Rigid_WriteWTG3
  
  subroutine Rigid1_WriteWTG3( r, t, file )
    use time_module
    type(sRigid1), intent(IN) :: r
    type(sTime),intent(IN) :: t
    integer,intent(IN) :: file
    !
    !速度の単位はA/ps
    !角速度の単位はrad/ps
    !
    write(file,'(13(e17.10,1x))') r%com%vec(1),r%com%vec(2),r%com%vec(3),&
         r%quat(1)%vec(1),r%quat(1)%vec(2),r%quat(1)%vec(3),r%quat(1)%vec(4)&
         ,r%dx(1)/t%dt,r%dy(1)/t%dt,r%dz(1)/t%dt,&
         r%wx(1),r%wy(1),r%wz(1)
    return
  end subroutine Rigid1_WriteWTG3

  subroutine Rigid_WriteNX4A(r,m,file)
    type(sRigid),intent(IN) :: r
    type(sMol),intent(IN) :: m
    integer,intent(IN) :: file
    integer :: n,i
    n=m%nmol
    write(file,'("@ID08")')
    write(file,'(a8)') m%id
    write(file,'("@NX4A")')
    write(file,*) n
    write(file,'(7(e17.10,1x))') &
         (r%mol(i)%com%vec(1),&
         r%mol(i)%com%vec(2),&
         r%mol(i)%com%vec(3),&
         r%mol(i)%quat(1)%vec(1),&
         r%mol(i)%quat(1)%vec(2),&
         r%mol(i)%quat(1)%vec(3),&
         r%mol(i)%quat(1)%vec(4),i=1,n)
    return
  end subroutine Rigid_WriteNX4A

  subroutine Rigid_WriteBinaryNX4A(r,n,file)
    type(sRigid),intent(IN) :: r
    integer,intent(IN) :: file
    integer :: n,i
    write(file) "@NX4A"
    write(file) n
    write(file) &
         (r%mol(i)%com%vec(1),&
         r%mol(i)%com%vec(2),&
         r%mol(i)%com%vec(3),&
         r%mol(i)%quat(1)%vec(1),&
         r%mol(i)%quat(1)%vec(2),&
         r%mol(i)%quat(1)%vec(3),&
         r%mol(i)%quat(1)%vec(4),i=1,n)
  end subroutine Rigid_WriteBinaryNX4A

  subroutine Rigid_ReadBinaryNX4A(r,m,si,file)
    use site_module
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(INOUT)   :: m
    type(sSite),intent(INOUT)  :: si
    integer,intent(IN) :: file
    integer :: n
    integer :: i
    read(file) n
    call Mol_Allocate(m,si,n,RIGID_MODE)
    call Rigid_Allocate(r,n)
    r%mol(:)%com%vec(1)=0
    r%mol(:)%com%vec(2)=0
    r%mol(:)%com%vec(3)=0
    call Rigid_ResetDerivatives( r, n )
    read(file) &
         (r%mol(i)%com%vec(1),&
         r%mol(i)%com%vec(2),&
         r%mol(i)%com%vec(3),&
         r%mol(i)%quat(1)%vec(1),&
         r%mol(i)%quat(1)%vec(2),&
         r%mol(i)%quat(1)%vec(3),&
         r%mol(i)%quat(1)%vec(4),i=1,n)
  end subroutine Rigid_ReadBinaryNX4A

  subroutine rigid_water_writeyaplot(r,m,file,layer,fBox,box)
    type(sRigid),intent(inout) :: r
    type(sMol),intent(in) :: m
    logical,intent(in) :: fBox
    type(sBox),intent(in) :: box
    integer,intent(in) :: file,layer
    integer i
    type(vector3)             :: xyz,cell
    write(file,3) layer
3   format("y ",i4)
    do i=1,m%nmol
       xyz%vec(:) = r%mol(i)%com%vec(:)
       if(fBox)then
          call box_renormalize(box,xyz,cell)
       endif
       write(file,1) "l ",xyz%vec(1)+r%mol(i)%intrax(3),xyz%vec(2)+r%mol(i)%intray(3),xyz%vec(3)&
            & +r&
            & %mol(i)%intraz(3),xyz%vec(1)+r%mol(i)%intrax(2),xyz%vec(2)+r%mol(i)&
            & %intray(2),xyz%vec(3)+r%mol(i)%intraz(2)
       write(file,1) "l ",xyz%vec(1)+r%mol(i)%intrax(4),xyz%vec(2)+r%mol(i)%intray(4),xyz%vec(3)&
            & +r&
            & %mol(i)%intraz(4),xyz%vec(1)+r%mol(i)%intrax(2),xyz%vec(2)+r%mol(i)&
            & %intray(2),xyz%vec(3)+r%mol(i)%intraz(2)
1      format(a2,6(1x,f7.2))
    enddo
  end subroutine rigid_water_writeyaplot
  
  subroutine yaplot_nextpage(file)
    integer,intent(in) :: file
    write(file,2)
2   format("")
  end subroutine yaplot_nextpage

!(x1+x2+x3+x4+x5)/5
!=(x1+(x1+d1)+(x1+d1+d2)+(x1+d1+d2+d3)+(x1+d1+d2+d3+d4))/5
!=x1+d1*4/5+d2*3/5+d3*2/5+d4*1/5
  
  subroutine Rigid_Version
    write(STDERR,*) "$Id: Rigid.F90,v 1.22 2002/11/01 07:22:49 matto&
         & Exp $"
    return
  end subroutine Rigid_Version

  
  subroutine rigid_vreset_obsolete(r,m,t,temp,seed)
    use physconst_module
    use random_module
    use time_module
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    type(sTime),intent(IN) :: t
    !temperature
    real(kind=8),intent(in) :: temp
    integer,intent(in) :: seed
    integer :: i
    real(kind=8) :: r1,r2
    real(kind=8),dimension(12*MAXMOL) :: random
    !velocity distribution is normal distribution (avg=0, var=sig**2
    ! =kT/m)
    !
    !mtを直接呼ぶべきではない。
    !
    call mtrngv(seed,m%nmol*12)
    call mtrndv(random,m%nmol*12)
    do i=1,m%nmol
       r1=random(0*m%nmol+i)
       r2=random(1*m%nmol+i)
       r%mol(i)%dx(1)=Random_Normal(0d0,temp*rgas*J2I*r%massi,r1,r2)&
            & *t%Dt
       r1=random(2*m%nmol+i)
       r2=random(3*m%nmol+i)
       r%mol(i)%dy(1)=Random_Normal(0d0,temp*rgas*J2I*r%massi,r1,r2)&
            & *t%Dt
       r1=random(4*m%nmol+i)
       r2=random(5*m%nmol+i)
       r%mol(i)%dz(1)=Random_Normal(0d0,temp*rgas*J2I*r%massi,r1,r2)&
            & *t%Dt
       r1=random(6*m%nmol+i)
       r2=random(7*m%nmol+i)
       r%mol(i)%wx(1)=Random_Normal(0d0,temp*rgas*J2I*r%ixxi,r1,r2)
       r1=random(8*m%nmol+i)
       r2=random(9*m%nmol+i)
       r%mol(i)%wy(1)=Random_Normal(0d0,temp*rgas*J2I*r%iyyi,r1,r2)
       r1=random(10*m%nmol+i)
       r2=random(11*m%nmol+i)
       r%mol(i)%wz(1)=Random_Normal(0d0,temp*rgas*J2I*r%izzi,r1,r2)
    enddo
  end subroutine rigid_vreset_obsolete

  subroutine rigid_vreset2(r,m,t,temp,rand)
    use physconst_module
    use random_module
    use time_module
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    type(sTime),intent(IN) :: t
    !temperature
    real(kind=8),intent(in) :: temp
    type(pRandom), pointer :: rand
    integer :: i
    real(kind=8) :: r1,r2
    !velocity distribution is normal distribution (avg=0, var=sig**2
    ! =kT/m)
    do i=1,m%nmol
       r1=Random_GetNext( rand )
       r2=Random_GetNext( rand )
       r%mol(i)%dx(1)=Random_Normal(0d0,temp*rgas*J2I*r%massi,r1,r2)&
            & *t%Dt
       r1=Random_GetNext( rand )
       r2=Random_GetNext( rand )
       r%mol(i)%dy(1)=Random_Normal(0d0,temp*rgas*J2I*r%massi,r1,r2)&
            & *t%Dt
       r1=Random_GetNext( rand )
       r2=Random_GetNext( rand )
       r%mol(i)%dz(1)=Random_Normal(0d0,temp*rgas*J2I*r%massi,r1,r2)&
            & *t%Dt
       r1=Random_GetNext( rand )
       r2=Random_GetNext( rand )
       r%mol(i)%wx(1)=Random_Normal(0d0,temp*rgas*J2I*r%ixxi,r1,r2)
       r1=Random_GetNext( rand )
       r2=Random_GetNext( rand )
       r%mol(i)%wy(1)=Random_Normal(0d0,temp*rgas*J2I*r%iyyi,r1,r2)
       r1=Random_GetNext( rand )
       r2=Random_GetNext( rand )
       r%mol(i)%wz(1)=Random_Normal(0d0,temp*rgas*J2I*r%izzi,r1,r2)
    enddo
  end subroutine rigid_vreset2

! Get total momenta of the group r
  subroutine rigid_gettotalmomenta_obsolete(r,m,t,momentx,momenty,momentz,mass)
    use time_module
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    type(sTime),intent(IN) :: t
    real(kind=8),intent(out) :: momentx,momenty,momentz,mass
    real(kind=8) :: dti
    integer :: i
    momentx=0d0
    momenty=0d0
    momentz=0d0
    mass = r%mass * m%nmol
    dti=1d0/t%dt
    do i=1,m%nmol
       momentx = momentx + r%mol(i)%dx(1)
       momenty = momenty + r%mol(i)%dy(1)
       momentz = momentz + r%mol(i)%dz(1)
    enddo
    momentx = momentx * r%mass * dti
    momenty = momenty * r%mass * dti
    momentz = momentz * r%mass * dti
  end subroutine rigid_gettotalmomenta_obsolete
  
! Add total momenta to the group r
  subroutine rigid_addmomenta_obsolete(r,m,t,momentx,momenty,momentz)
    use time_module
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    type(sTime),intent(IN) :: t
    real(kind=8),intent(IN) :: momentx,momenty,momentz
    real(kind=8) :: velx,vely,velz
    real(kind=8) :: nmol
    nmol=m%nmol
    !
    !たぶん間違っている。角速度は一律だが速度は一律ではない。
    !
    velx = momentx * t%dt / (nmol*r%mass)
    vely = momenty * t%dt / (nmol*r%mass)
    velz = momentz * t%dt / (nmol*r%mass)
    r%mol(1:m%nmol)%dx(1) = r%mol(1:m%nmol)%dx(1) + velx
    r%mol(1:m%nmol)%dy(1) = r%mol(1:m%nmol)%dy(1) + vely
    r%mol(1:m%nmol)%dz(1) = r%mol(1:m%nmol)%dz(1) + velz
  end subroutine rigid_addmomenta_obsolete

  subroutine Rigid_GetOffset(r,m,moment,mass)
    type(sRigid), intent(IN)    :: r
    type(sMol),   intent(IN)    :: m
    type(vector3),intent(OUT)   :: moment
    real(kind=8), intent(OUT)   :: mass
    !
    !local
    !
    real(kind=8) :: x,y,z
    integer :: i
    moment%vec(:) =0d0
    mass=r%mass * m%nmol
    do i=1,m%nmol
       x = r%mol(i)%com%vec(1)
       y = r%mol(i)%com%vec(2)
       z = r%mol(i)%com%vec(3)
       moment%vec(1) = moment%vec(1) + x
       moment%vec(2) = moment%vec(2) + y
       moment%vec(3) = moment%vec(3) + z
    enddo
    moment%vec(:) = moment%vec(:) * r%mass
  end subroutine Rigid_GetOffset

  subroutine Rigid_GetMomentum(r,m,moment,mass)
    type(sRigid), intent(IN)    :: r
    type(sMol),   intent(IN)    :: m
    type(vector3),intent(OUT)   :: moment
    real(kind=8), intent(OUT)   :: mass
    !
    !local
    !
    real(kind=8) :: vx,vy,vz
    integer :: i
    moment%vec(:) = 0d0
    mass          = r%mass * m%nmol
    do i=1,m%nmol
       vx = r%mol(i)%dx(1)
       vy = r%mol(i)%dy(1)
       vz = r%mol(i)%dz(1)
       moment%vec(1) = moment%vec(1) + vx
       moment%vec(2) = moment%vec(2) + vy
       moment%vec(3) = moment%vec(3) + vz
    enddo
    moment%vec(:) = moment%vec(:) * r%mass
  end subroutine Rigid_GetMomentum

  subroutine Rigid_AddMomenta(r,m,w)
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    type(vector3),intent(in) :: w
    !
    !Local variables
    !
    integer      :: i
    do i=1,m%nmol
       r%mol(i)%dx(1) = r%mol(i)%dx(1) + w%vec(1)
       r%mol(i)%dy(1) = r%mol(i)%dy(1) + w%vec(2)
       r%mol(i)%dz(1) = r%mol(i)%dz(1) + w%vec(3)
    enddo
  end subroutine Rigid_AddMomenta

  subroutine Rigid_ShiftPosition(r,m,w)
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    type(vector3),intent(in) :: w
    !
    !Local variables
    !
    integer      :: i
    do i=1,m%nmol
       r%mol(i)%com%vec(:) = r%mol(i)%com%vec(:) + w%vec(:)
    enddo
  end subroutine Rigid_ShiftPosition
  !
  ! Get total momenta of the group r 
  ! momentの値にはdtがかかっている。
  !
  subroutine Rigid_GetAngularMomenta(r,m,moment)
    type(sRigid), intent(IN)    :: r
    type(sMol),   intent(IN)    :: m
    type(vector3),intent(OUT)   :: moment
    !
    !local
    !
    real(kind=8) :: x,y,z,vx,vy,vz
    integer :: i
    moment%vec(:) =0d0
    do i=1,m%nmol
       x = r%mol(i)%com%vec(1)
       y = r%mol(i)%com%vec(2)
       z = r%mol(i)%com%vec(3)
       vx = r%mol(i)%dx(1)
       vy = r%mol(i)%dy(1)
       vz = r%mol(i)%dz(1)
       moment%vec(1) = moment%vec(1) + y*vz - z*vy
       moment%vec(2) = moment%vec(2) + z*vx - x*vz
       moment%vec(3) = moment%vec(3) + x*vy - y*vx
    enddo
    moment%vec(:) = moment%vec(:) * r%mass
  end subroutine Rigid_GetAngularMomenta
  
  subroutine rigid_getinertiatensor(r,m,t)
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    real(kind=8),intent(OUT) :: t(3,3)
    !
    !Local variables
    !
    real(kind=8) :: x,y,z
    integer :: i
    t(:,:) = 0d0
    do i=1,m%nmol
       x = r%mol(i)%com%vec(1)
       y = r%mol(i)%com%vec(2)
       z = r%mol(i)%com%vec(3)
       t(1,1) = t(1,1) + y**2 + z**2
       t(1,2) = t(1,2) - x*y
       t(1,3) = t(1,3) - x*z
       t(2,1) = t(2,1) - y*x
       t(2,2) = t(2,2) + x**2 + z**2
       t(2,3) = t(2,3) - y*z
       t(3,1) = t(3,1) - z*x
       t(3,2) = t(3,2) - z*y
       t(3,3) = t(3,3) + x**2 + y**2
    enddo
    t(:,:) = t(:,:) * r%mass
  end subroutine rigid_getinertiatensor

! Add total momenta to the group r
  subroutine rigid_AddAngularVelocity(r,m,w)
    type(sRigid),intent(INOUT) :: r
    type(sMol),intent(IN) :: m
    type(vector3),intent(in) :: w
    !
    !Local variables
    !
    real(kind=8) :: velx,vely,velz
    real(kind=8) :: x,y,z
    integer      :: i
    do i=1,m%nmol
       !
       !goldstein (5-1)
       !
       x = r%mol(i)%com%vec(1)
       y = r%mol(i)%com%vec(2)
       z = r%mol(i)%com%vec(3)
       velx = w%vec(2) * z - w%vec(3) * y
       vely = w%vec(3) * x - w%vec(1) * z
       velz = w%vec(1) * y - w%vec(2) * x
       r%mol(i)%dx(1) = r%mol(i)%dx(1) + velx
       r%mol(i)%dy(1) = r%mol(i)%dy(1) + vely
       r%mol(i)%dz(1) = r%mol(i)%dz(1) + velz
    enddo
  end subroutine rigid_AddAngularVelocity

subroutine rigid_correct(r,m,t,fNose,nose,fAndersen,andersen)
  use physconst_module
  use mol_module
  use time_module
  use nose_module
  use andersen_module
  type(sMol),intent(IN) :: m
  type(sRigid),intent(INOUT) :: r
  type(sTime),intent(IN) :: t
  type(sNose),intent(IN) :: nose
  type(sAndersen),intent(IN) :: andersen
  logical,intent(IN) :: fNose,fAndersen
  integer i
  real(kind=8) newwa,newwb,newwc,newwd
  real(kind=8) deltawa,deltawb,deltawc,deltawd
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
        ! 上田本
        acoeff = andersen%v1 / (2d0*3d0*andersen%v0)
        ! 岡崎本P.115 (7.47)
        ! %v1 == v-dot dt, %v0 == Volume
        ! ただし、次の計算のためにあらかじめ1/2が掛けてある。
        !acoeff = andersen%v1 / (3d0*andersen%v0)
        extendtx = extendtx + acoeff
        extendty = extendty + acoeff
        extendtz = extendtz + acoeff
     endif
  endif
  do i=1,m%nmol
     newfx = r%mol(i)%forcex*(t%dt*t%dt*0.5d0*r%Massi)
     newfy = r%mol(i)%forcey*(t%dt*t%dt*0.5d0*r%Massi)
     newfz = r%mol(i)%forcez*(t%dt*t%dt*0.5d0*r%Massi)
     newfx=newfx - r%mol(i)%dx(1) * extendtx
     newfy=newfy - r%mol(i)%dy(1) * extendty
     newfz=newfz - r%mol(i)%dz(1) * extendtz
     deltax = newfx - r%mol(i)%dx(2)
     deltay = newfy - r%mol(i)%dy(2)
     deltaz = newfz - r%mol(i)%dz(2)
     r%mol(i)%com%vec(1) = r%mol(i)%com%vec(1) + Gear5_0*deltax
     r%mol(i)%com%vec(2) = r%mol(i)%com%vec(2) + Gear5_0*deltay
     r%mol(i)%com%vec(3) = r%mol(i)%com%vec(3) + Gear5_0*deltaz
     r%mol(i)%dx(1)=r%mol(i)%dx(1) + Gear5_1*deltax
     r%mol(i)%dy(1)=r%mol(i)%dy(1) + Gear5_1*deltay
     r%mol(i)%dz(1)=r%mol(i)%dz(1) + Gear5_1*deltaz
     r%mol(i)%dx(2)=newfx
     r%mol(i)%dy(2)=newfy
     r%mol(i)%dz(2)=newfz
     r%mol(i)%dx(3)=r%mol(i)%dx(3) + Gear5_3*deltax
     r%mol(i)%dy(3)=r%mol(i)%dy(3) + Gear5_3*deltay
     r%mol(i)%dz(3)=r%mol(i)%dz(3) + Gear5_3*deltaz
     r%mol(i)%dx(4)=r%mol(i)%dx(4) + Gear5_4*deltax
     r%mol(i)%dy(4)=r%mol(i)%dy(4) + Gear5_4*deltay
     r%mol(i)%dz(4)=r%mol(i)%dz(4) + Gear5_4*deltaz
     r%mol(i)%dx(5)=r%mol(i)%dx(5) + Gear5_5*deltax
     r%mol(i)%dy(5)=r%mol(i)%dy(5) + Gear5_5*deltay
     r%mol(i)%dz(5)=r%mol(i)%dz(5) + Gear5_5*deltaz
     !if(i==1)write(6,*) "FX",newfx,deltax/newfx,r%mol(i)%com%vec(1)(i),r%mol(i)%dx(1)&
     !     & ,r%mol(i)%dx(2),r%mol(i)%dx(3),r%mol(i)%dx(4),r%mol(i)%dx(5)
     newwa = (-r%mol(i)%wx(1) * r%mol(i)%quat(1)%vec(2) &
              +r%mol(i)%wy(1) * r%mol(i)%quat(1)%vec(3) &
              -r%mol(i)%wz(1) * r%mol(i)%quat(1)%vec(4)) * t%dt * 0.5d0
     newwb = ( r%mol(i)%wx(1) * r%mol(i)%quat(1)%vec(1) &
              -r%mol(i)%wy(1) * r%mol(i)%quat(1)%vec(4) &
              -r%mol(i)%wz(1) * r%mol(i)%quat(1)%vec(3)) * t%dt * 0.5d0
     newwc = (-r%mol(i)%wx(1) * r%mol(i)%quat(1)%vec(4) &
              -r%mol(i)%wy(1) * r%mol(i)%quat(1)%vec(1) &
              +r%mol(i)%wz(1) * r%mol(i)%quat(1)%vec(2)) * t%dt * 0.5d0
     newwd = ( r%mol(i)%wx(1) * r%mol(i)%quat(1)%vec(3) &
              +r%mol(i)%wy(1) * r%mol(i)%quat(1)%vec(2) &
              +r%mol(i)%wz(1) * r%mol(i)%quat(1)%vec(1)) * t%dt * 0.5d0
     !if(i==1)write(6,*) "wx",r%mol(i)%wx(1)*t%dt,r%mol(i)%quat(1)%vec(1)
     !if(i.eq.1)write(STDOUT,*) newwa,newwb,newwc,newwd
     deltawa = newwa - r%mol(i)%quat(2)%vec(1)
     deltawb = newwb - r%mol(i)%quat(2)%vec(2)
     deltawc = newwc - r%mol(i)%quat(2)%vec(3)
     deltawd = newwd - r%mol(i)%quat(2)%vec(4)
     r%mol(i)%quat(1)%vec(1) = r%mol(i)%quat(1)%vec(1) + Gear4_0*deltawa
     r%mol(i)%quat(1)%vec(2) = r%mol(i)%quat(1)%vec(2) + Gear4_0*deltawb
     r%mol(i)%quat(1)%vec(3) = r%mol(i)%quat(1)%vec(3) + Gear4_0*deltawc
     r%mol(i)%quat(1)%vec(4) = r%mol(i)%quat(1)%vec(4) + Gear4_0*deltawd
     !if(i==1)write(6,*) "qa",newwa,deltawa/newwa,r%mol(i)%quat(1)%vec(1)
     !     Goldstein P.268
     newfx=(r%mol(i)%torquex+(r%Iyy-r%Izz)*r%mol(i)%wy(1)*r%mol(i)%wz(1))*t%dt*r%Ixxi
     newfy=(r%mol(i)%torquey+(r%Izz-r%Ixx)*r%mol(i)%wz(1)*r%mol(i)%wx(1))*t%dt*r%Iyyi
     newfz=(r%mol(i)%torquez+(r%Ixx-r%Iyy)*r%mol(i)%wx(1)*r%mol(i)%wy(1))*t%dt*r%Izzi
     !if(i==1)write(6,*) "torque",r%mol(i)%torquex(i)*t%dt**2/18d0,(r%Iyy&
     !& -r%Izz)*r%mol(i)%wy(1)*r&
     !& %wz(1)*t%dt**2/18d0,newfx*t%dt,(r%Iyy-r%Izz)/18d0
     newfx = newfx - r%mol(i)%wx(1) * extendrx
     newfy = newfy - r%mol(i)%wy(1) * extendry
     newfz = newfz - r%mol(i)%wz(1) * extendrz

     deltax = newfx - r%mol(i)%wx(2)
     !     write(6,*) "deltax",deltax,torquex(i)
     deltay = newfy - r%mol(i)%wy(2)
     deltaz = newfz - r%mol(i)%wz(2)
     r%mol(i)%wx(1) = r%mol(i)%wx(1) + Gear4_0*deltax
     r%mol(i)%wy(1) = r%mol(i)%wy(1) + Gear4_0*deltay
     r%mol(i)%wz(1) = r%mol(i)%wz(1) + Gear4_0*deltaz
     r%mol(i)%wx(2) = newfx
     r%mol(i)%wy(2) = newfy
     r%mol(i)%wz(2) = newfz
     r%mol(i)%wx(3) = r%mol(i)%wx(3) + Gear4_2*deltax
     r%mol(i)%wy(3) = r%mol(i)%wy(3) + Gear4_2*deltay
     r%mol(i)%wz(3) = r%mol(i)%wz(3) + Gear4_2*deltaz
     r%mol(i)%wx(4) = r%mol(i)%wx(4) + Gear4_3*deltax
     r%mol(i)%wy(4) = r%mol(i)%wy(4) + Gear4_3*deltay
     r%mol(i)%wz(4) = r%mol(i)%wz(4) + Gear4_3*deltaz
     r%mol(i)%wx(5) = r%mol(i)%wx(5) + Gear4_4*deltax
     r%mol(i)%wy(5) = r%mol(i)%wy(5) + Gear4_4*deltay
     r%mol(i)%wz(5) = r%mol(i)%wz(5) + Gear4_4*deltaz
     !if(i==1)write(6,*) "newfx",newfx,deltax,r%wx(1),r%wx(2),r&
     !     & %wx(3),r%wx(4),r%wx(5)

     r%mol(i)%quat(2)%vec(3)=&
          (-r%mol(i)%wx(1)*r%mol(i)%quat(1)%vec(4)-r%mol(i)%wy(1)*r&
          & %mol(i)%quat(1)%vec(1)+r%mol(i)%wz(1)*r%mol(i)%quat(1)&
          & %vec(2))*t%dt*0.5d0
     r%mol(i)%quat(2)%vec(2)=&
          ( r%mol(i)%wx(1)*r%mol(i)%quat(1)%vec(1)-r%mol(i)%wy(1)*r&
          & %mol(i)%quat(1)%vec(4)-r%mol(i)%wz(1)*r%mol(i)%quat(1)&
          & %vec(3))*t%dt*0.5d0
     r%mol(i)%quat(2)%vec(4)=&
          ( r%mol(i)%wx(1)*r%mol(i)%quat(1)%vec(3)+r%mol(i)%wy(1)*r&
          & %mol(i)%quat(1)%vec(2)+r%mol(i)%wz(1)*r%mol(i)%quat(1)&
          & %vec(1))*t%dt*0.5d0
     r%mol(i)%quat(2)%vec(1)=&
          (-r%mol(i)%wx(1)*r%mol(i)%quat(1)%vec(2)+r%mol(i)%wy(1)*r&
          & %mol(i)%quat(1)%vec(3)-r%mol(i)%wz(1)*r%mol(i)%quat(1)&
          & %vec(4))*t%dt*0.5d0
     
     r%mol(i)%quat(3)%vec(3)=(-r%mol(i)%wx(2)*r%mol(i)%quat(1)%vec(4)-r%mol(i)%wy(2)*r&
          & %mol(i)%quat(1)%vec(1)+r%mol(i)%wz(2)*r%mol(i)%quat(1)%vec(2)-r%mol(i)&
          & %wx(1)*r%mol(i)%quat(2)%vec(4)-r%mol(i)%wy(1)*r%mol(i)%quat(2)%vec(1)+r&
          & %mol(i)%wz(1)*r%mol(i)%quat(2)%vec(2))*t%dt*0.25d0
     r%mol(i)%quat(3)%vec(2)=( r%mol(i)%wx(2)*r%mol(i)%quat(1)%vec(1)-r%mol(i)%wy(2)*r&
          & %mol(i)%quat(1)%vec(4)-r%mol(i)%wz(2)*r%mol(i)%quat(1)%vec(3)+r%mol(i)&
          & %wx(1)*r%mol(i)%quat(2)%vec(1)-r%mol(i)%wy(1)*r%mol(i)%quat(2)%vec(4)-r&
          & %mol(i)%wz(1)*r%mol(i)%quat(2)%vec(3))*t%dt*0.25d0
     r%mol(i)%quat(3)%vec(4)=( r%mol(i)%wx(2)*r%mol(i)%quat(1)%vec(3)+r%mol(i)%wy(2)*r&
          & %mol(i)%quat(1)%vec(2)+r%mol(i)%wz(2)*r%mol(i)%quat(1)%vec(1)+r%mol(i)&
          & %wx(1)*r%mol(i)%quat(2)%vec(3)+r%mol(i)%wy(1)*r%mol(i)%quat(2)%vec(2)+r&
          & %mol(i)%wz(1)*r%mol(i)%quat(2)%vec(1))*t%dt*0.25d0
     r%mol(i)%quat(3)%vec(1)=(-r%mol(i)%wx(2)*r%mol(i)%quat(1)%vec(2)+r%mol(i)%wy(2)*r&
          & %mol(i)%quat(1)%vec(3)-r%mol(i)%wz(2)*r%mol(i)%quat(1)%vec(4)-r%mol(i)&
          & %wx(1)*r%mol(i)%quat(2)%vec(2)+r%mol(i)%wy(1)*r%mol(i)%quat(2)%vec(3)-r&
          & %mol(i)%wz(1)*r%mol(i)%quat(2)%vec(4))*t%dt*0.25d0
     
     r%mol(i)%quat(4)%vec(3)=(-r%mol(i)%wx(3)*r%mol(i)%quat(1)%vec(4)-r%mol(i)%wy(3)*r&
          & %mol(i)%quat(1)%vec(1)+r%mol(i)%wz(3)*r%mol(i)%quat(1)%vec(2)-r%mol(i)&
          & %wx(1)*r%mol(i)%quat(3)%vec(4)-r%mol(i)%wy(1)*r%mol(i)%quat(3)%vec(1)+r&
          & %mol(i)%wz(1)*r%mol(i)%quat(3)%vec(2)-r%mol(i)%wx(2)*r%mol(i)&
          & %quat(2)%vec(4)-r%mol(i)%wy(2)*r%mol(i)%quat(2)%vec(1)+r%mol(i)%wz(2)*r&
          & %mol(i)%quat(2)%vec(2))*t%dt*(1d0/6d0)
     r%mol(i)%quat(4)%vec(2)=( r%mol(i)%wx(3)*r%mol(i)%quat(1)%vec(1)-r%mol(i)%wy(3)*r&
          & %mol(i)%quat(1)%vec(4)-r%mol(i)%wz(3)*r%mol(i)%quat(1)%vec(3)+r%mol(i)&
          & %wx(1)*r%mol(i)%quat(3)%vec(1)-r%mol(i)%wy(1)*r%mol(i)%quat(3)%vec(4)-r&
          & %mol(i)%wz(1)*r%mol(i)%quat(3)%vec(3)+r%mol(i)%wx(2)*r%mol(i)&
          & %quat(2)%vec(1)-r%mol(i)%wy(2)*r%mol(i)%quat(2)%vec(4)-r%mol(i)%wz(2)*r&
          & %mol(i)%quat(2)%vec(3))*t%dt*(1d0/6d0)
     r%mol(i)%quat(4)%vec(4)=( r%mol(i)%wx(3)*r%mol(i)%quat(1)%vec(3)+r%mol(i)%wy(3)*r&
          & %mol(i)%quat(1)%vec(2)+r%mol(i)%wz(3)*r%mol(i)%quat(1)%vec(1)+r%mol(i)&
          & %wx(1)*r%mol(i)%quat(3)%vec(3)+r%mol(i)%wy(1)*r%mol(i)%quat(3)%vec(2)+r&
          & %mol(i)%wz(1)*r%mol(i)%quat(3)%vec(1)+r%mol(i)%wx(2)*r%mol(i)&
          & %quat(2)%vec(3)+r%mol(i)%wy(2)*r%mol(i)%quat(2)%vec(2)+r%mol(i)%wz(2)*r&
          & %mol(i)%quat(2)%vec(1))*t%dt*(1d0/6d0)
     r%mol(i)%quat(4)%vec(1)=(-r%mol(i)%wx(3)*r%mol(i)%quat(1)%vec(2)+r%mol(i)%wy(3)*r&
          & %mol(i)%quat(1)%vec(3)-r%mol(i)%wz(3)*r%mol(i)%quat(1)%vec(4)-r%mol(i)&
          & %wx(1)*r%mol(i)%quat(3)%vec(2)+r%mol(i)%wy(1)*r%mol(i)%quat(3)%vec(3)-r&
          & %mol(i)%wz(1)*r%mol(i)%quat(3)%vec(4)-r%mol(i)%wx(2)*r%mol(i)&
          & %quat(2)%vec(2)+r%mol(i)%wy(2)*r%mol(i)%quat(2)%vec(3)-r%mol(i)%wz(2)*r&
          & %mol(i)%quat(2)%vec(4))*t%dt*(1d0/6d0)

     r%mol(i)%quat(5)%vec(3)=(-r%mol(i)%wx(4)*r%mol(i)%quat(1)%vec(4)-r%mol(i)%wy(4)*r&
          & %mol(i)%quat(1)%vec(1)+r%mol(i)%wz(4)*r%mol(i)%quat(1)%vec(2)-r%mol(i)&
          & %wx(3)*r%mol(i)%quat(2)%vec(4)-r%mol(i)%wy(3)*r%mol(i)%quat(2)%vec(1)+r&
          & %mol(i)%wz(3)*r%mol(i)%quat(2)%vec(2)-r%mol(i)%wx(2)*r%mol(i)&
          & %quat(3)%vec(4)-r%mol(i)%wy(2)*r%mol(i)%quat(3)%vec(1)+r%mol(i)%wz(2)*r&
          & %mol(i)%quat(3)%vec(2)-r%mol(i)%wx(1)*r%mol(i)%quat(4)%vec(4)-r%mol(i)&
          & %wy(1)*r%mol(i)%quat(4)%vec(1)+r%mol(i)%wz(1)*r%mol(i)%quat(4)%vec(2))*t%dt&
          & *(1d0/8d0)
     r%mol(i)%quat(5)%vec(2)=( r%mol(i)%wx(4)*r%mol(i)%quat(1)%vec(1)-r%mol(i)%wy(4)*r&
          & %mol(i)%quat(1)%vec(4)-r%mol(i)%wz(4)*r%mol(i)%quat(1)%vec(3)+r%mol(i)&
          & %wx(3)*r%mol(i)%quat(2)%vec(1)-r%mol(i)%wy(3)*r%mol(i)%quat(2)%vec(4)-r&
          & %mol(i)%wz(3)*r%mol(i)%quat(2)%vec(3)+r%mol(i)%wx(2)*r%mol(i)&
          & %quat(3)%vec(1)-r%mol(i)%wy(2)*r%mol(i)%quat(3)%vec(4)-r%mol(i)%wz(2)*r&
          & %mol(i)%quat(3)%vec(3)+r%mol(i)%wx(1)*r%mol(i)%quat(4)%vec(1)-r%mol(i)&
          & %wy(1)*r%mol(i)%quat(4)%vec(4)-r%mol(i)%wz(1)*r%mol(i)%quat(4)%vec(3))*t%dt&
          & *(1d0/8d0)
     r%mol(i)%quat(5)%vec(4)=( r%mol(i)%wx(4)*r%mol(i)%quat(1)%vec(3)+r%mol(i)%wy(4)*r&
          & %mol(i)%quat(1)%vec(2)+r%mol(i)%wz(4)*r%mol(i)%quat(1)%vec(1)+r%mol(i)&
          & %wx(3)*r%mol(i)%quat(2)%vec(3)+r%mol(i)%wy(3)*r%mol(i)%quat(2)%vec(2)+r&
          & %mol(i)%wz(3)*r%mol(i)%quat(2)%vec(1)+r%mol(i)%wx(2)*r%mol(i)&
          & %quat(3)%vec(3)+r%mol(i)%wy(2)*r%mol(i)%quat(3)%vec(2)+r%mol(i)%wz(2)*r&
          & %mol(i)%quat(3)%vec(1)+r%mol(i)%wx(1)*r%mol(i)%quat(4)%vec(3)+r%mol(i)&
          & %wy(1)*r%mol(i)%quat(4)%vec(2)+r%mol(i)%wz(1)*r%mol(i)%quat(4)%vec(1))*t%dt&
          & *(1d0/8d0)
     r%mol(i)%quat(5)%vec(1)=(-r%mol(i)%wx(4)*r%mol(i)%quat(1)%vec(2)+r%mol(i)%wy(4)*r&
          & %mol(i)%quat(1)%vec(3)-r%mol(i)%wz(4)*r%mol(i)%quat(1)%vec(4)-r%mol(i)&
          & %wx(3)*r%mol(i)%quat(2)%vec(2)+r%mol(i)%wy(3)*r%mol(i)%quat(2)%vec(3)-r&
          & %mol(i)%wz(3)*r%mol(i)%quat(2)%vec(4)-r%mol(i)%wx(2)*r%mol(i)&
          & %quat(3)%vec(2)+r%mol(i)%wy(2)*r%mol(i)%quat(3)%vec(3)-r%mol(i)%wz(2)*r&
          & %mol(i)%quat(3)%vec(4)-r%mol(i)%wx(1)*r%mol(i)%quat(4)%vec(2)+r%mol(i)&
          & %wy(1)*r%mol(i)%quat(4)%vec(3)-r%mol(i)%wz(1)*r%mol(i)&
          & %quat(4)%vec(4))*t%dt&
          & *(1d0/8d0)
  enddo
  return
end subroutine rigid_correct

!後方互換
  subroutine rigid_correct2(r,m,t,nose,andersen)
    use mol_module
    use andersen_module
    use nose_module
    use time_module
    type(sMol),intent(IN) :: m
    type(sRigid),intent(INOUT) :: r
    type(sTime),intent(IN) :: t
    type(sNose),intent(IN) :: nose
    type(sAndersen),intent(IN),optional :: andersen
    if(present(andersen))then
       call rigid_correct(r,m,t,nose%active,nose,andersen%mode.ne.noandersen,andersen)
    else
       call rigid_correct(r,m,t,nose%active,nose,.false.,andersen)
    endif
  end subroutine rigid_correct2

  function rigid_serialize_position(rigid,mol,p,offset)
    use quat_module
    integer :: rigid_serialize_position
    type(sRigid),intent(in)  :: rigid
    type(sMol),intent(in)    :: mol
    real(kind=8),intent(inout) :: p(*)
    integer,intent(inout)    :: offset
    integer                  :: i
    real(kind=8)             :: ea,eb,ec
    do i=1,mol%nmol
       p(offset + i + mol%nmol*0)=rigid%mol(i)%com%vec(1)
       p(offset + i + mol%nmol*1)=rigid%mol(i)%com%vec(2)
       p(offset + i + mol%nmol*2)=rigid%mol(i)%com%vec(3)
       call abcd2abc(rigid%mol(i)%quat(1)%vec(1), &
            rigid%mol(i)%quat(1)%vec(2),          &
            rigid%mol(i)%quat(1)%vec(3),          &
            rigid%mol(i)%quat(1)%vec(4),ea,eb,ec)
       p(offset + i + mol%nmol*3)=ea
       p(offset + i + mol%nmol*4)=eb
       p(offset + i + mol%nmol*5)=ec
    enddo
    rigid_serialize_position = offset + mol%nmol * 6
  end function rigid_serialize_position

  integer function rigid_serialize_positionq(rigid,mol,p,offset)
    use quat_module
    type(sRigid),intent(in)  :: rigid
    type(sMol),intent(in)    :: mol
    real(kind=8),intent(inout) :: p(*)
    integer,intent(inout)    :: offset
    integer                  :: i
    do i=1,mol%nmol
       p(offset + i + mol%nmol*0 ) = rigid%mol(i)%com%vec(1)
       p(offset + i + mol%nmol*1 ) = rigid%mol(i)%com%vec(2)
       p(offset + i + mol%nmol*2 ) = rigid%mol(i)%com%vec(3)
       p(offset + i + mol%nmol*3 ) = rigid%mol(i)%quat(1)%vec(1)
       p(offset + i + mol%nmol*4 ) = rigid%mol(i)%quat(1)%vec(2)
       p(offset + i + mol%nmol*5 ) = rigid%mol(i)%quat(1)%vec(3)
       p(offset + i + mol%nmol*6 ) = rigid%mol(i)%quat(1)%vec(4)
    enddo
    rigid_serialize_positionq = offset + mol%nmol * 7
  end function rigid_serialize_positionq

  integer function rigid_unserialize_position(rigid,mol,p,offset)
    use quat_module
    type(sRigid),intent(inout) :: rigid
    type(sMol),intent(in)      :: mol
    real(kind=8),intent(in)    :: p(*)
    integer,intent(inout)      :: offset
    integer                    :: i
    real(kind=8)               :: ea,eb,ec,qa,qb,qc,qd
    do i=1,mol%nmol
       ea = p(offset + i + mol%nmol*3)
       eb = p(offset + i + mol%nmol*4)
       ec = p(offset + i + mol%nmol*5)
       call abc2abcd(ea,eb,ec,qa,qb,qc,qd)
       rigid%mol(i)%com%vec(1)     = p(offset + i + mol%nmol*0)
       rigid%mol(i)%com%vec(2)     = p(offset + i + mol%nmol*1)
       rigid%mol(i)%com%vec(3)     = p(offset + i + mol%nmol*2)
       rigid%mol(i)%quat(1)%vec(1) = qa
       rigid%mol(i)%quat(1)%vec(2) = qb
       rigid%mol(i)%quat(1)%vec(3) = qc
       rigid%mol(i)%quat(1)%vec(4) = qd
    enddo
    rigid_unserialize_position = offset + mol%nmol * 6
  end function rigid_unserialize_position

  integer function rigid_unserialize_positionq(rigid,mol,p,offset)
    use quat_module
    type(sRigid),intent(inout) :: rigid
    type(sMol),intent(in)      :: mol
    real(kind=8),intent(in)    :: p(*)
    integer,intent(inout)      :: offset
    integer                    :: i
    do i=1,mol%nmol
       rigid%mol(i)%com%vec(1)     = p(offset + i + mol%nmol*0)
       rigid%mol(i)%com%vec(2)     = p(offset + i + mol%nmol*1)
       rigid%mol(i)%com%vec(3)     = p(offset + i + mol%nmol*2)
       rigid%mol(i)%quat(1)%vec(1) = p(offset + i + mol%nmol*3)
       rigid%mol(i)%quat(1)%vec(2) = p(offset + i + mol%nmol*4)
       rigid%mol(i)%quat(1)%vec(3) = p(offset + i + mol%nmol*5)
       rigid%mol(i)%quat(1)%vec(4) = p(offset + i + mol%nmol*6)
    enddo
    rigid_unserialize_positionq = offset + mol%nmol * 7
  end function rigid_unserialize_positionq

  subroutine rigid_save( rigid, mol, file, mode, t )
    use time_module
    type(sRigid), intent(in) :: rigid
    type(sMol), intent(in)   :: mol
    integer, intent(in)      :: mode, file
    type(sTime), intent(in), optional  :: t
    if ( mode .eq. RIGID_NX4A ) then
       write(file,'("@FIXC")')
       write(file,*) mol%isFixed
       call Rigid_WriteNX4A( rigid, mol, file )
    else if ( mode .eq. RIGID_MDVW ) then
       call Rigid_WriteMDVW( rigid, mol, file )
    else if ( mode .eq. RIGID_WTG3 ) then
       call Rigid_WriteWTG3( rigid, mol, t, file )
    endif
  end subroutine rigid_save

  subroutine rigid_savebinary( rigid, mol, file, mode )
    type(sRigid), intent(in) :: rigid
    type(sMol), intent(in)   :: mol
    integer, intent(in)      :: mode, file
    write(file) "@FIXC"
    write(file) mol%isFixed
    if ( mode .eq. RIGID_WTG2 ) then
       call Rigid_WriteBinaryWTG2( rigid, mol, file )
    endif
  end subroutine rigid_savebinary

  subroutine Rigid_WriteMDVW(r,m,file)
    type(sRigid),intent(IN) :: r
    type(sMol),intent(IN) :: m
    integer,intent(IN) :: file
    integer :: n,mol,site,j
    !あらかじめrigid_setsitepositionで分子内座標が計算されているものとする。
    do mol=1, m%nmol
       do site=1, m%nsite-1
          if ( m%name(site)(1:1) .ne. " " ) then
             write(file,'(a8,3(1x,e17.10))') m%name(site) &
                  ,r%mol(mol)%com%vec(1)+r%mol(mol)%intrax(site) &
                  ,r%mol(mol)%com%vec(2)+r%mol(mol)%intray(site) &
                  ,r%mol(mol)%com%vec(3)+r%mol(mol)%intraz(site)
          endif
       enddo
    enddo
  end subroutine Rigid_WriteMDVW

  subroutine rigid_toInternalUnit( rigid, mol, ti )
    use time_module
    type(sRigid), intent(inout) :: rigid
    type(sMol),   intent(in)    :: mol
    type(sTime),  intent(in)    :: ti
    integer                     :: i,s
    do i=1,mol%nmol
       do s=1,5
          rigid%mol(i)%wx(s) = rigid%mol(i)%wx(s) / ti%dt
          rigid%mol(i)%wy(s) = rigid%mol(i)%wy(s) / ti%dt
          rigid%mol(i)%wz(s) = rigid%mol(i)%wz(s) / ti%dt
       enddo
    enddo
  end subroutine rigid_toInternalUnit

  subroutine rigid_toExternalUnit( rigid, mol, ti )
    use time_module
    type(sRigid), intent(inout) :: rigid
    type(sMol),   intent(in)    :: mol
    type(sTime),  intent(in)    :: ti
    integer                     :: i,s
    do i=1,mol%nmol
       do s=1,5
          rigid%mol(i)%wx(s) = rigid%mol(i)%wx(s) * ti%dt
          rigid%mol(i)%wy(s) = rigid%mol(i)%wy(s) * ti%dt
          rigid%mol(i)%wz(s) = rigid%mol(i)%wz(s) * ti%dt
       enddo
    enddo
  end subroutine rigid_toExternalUnit


end module rigid_module

! -*- f90 -*-
module interaction_module
#define TEST7
#ifdef TEST7
#define partner_i(x,y) partner_i_(y,x)
#define partner_j(x,y) partner_j_(y,x)
  !��������줫��������®���ʤä���
#endif

!
  !box_renormalize��vector����롣������®�����뤬���纹�Ϥʤ���
!undef�Ǹ���ʿ��15ǯ4��10��(��)
!
#undef TEST30
!
!inner_product()�������ѻߤ���ľ�˽񤤤Ƥߤ롣���̤���(�ܺ٤�Genesis5/input512��)
!define�Ǹ���ʿ��15ǯ4��10��(��)
!
#define TEST31


!
!box_renormalize��inlineŸ�����롣ʿ��15ǯ11��22��(��)
!������ơ������Ǥ�®�٤��Ѥ��ʤ����������ѡ�����ԥ塼�����ǰ㤦�Ϥ��ʤΤǤȤꤢ����define�ˤ��Ƥ�����ʿ��15ǯ11��22��(��)
!
#define TEST32
!
!DEC�Ǥϡ�MAXPAIR��parameterʸ���Ϥ��ȡ�internal error��Ҥ���������
!�����������������������ʤ��Τǡ������Ǻ�������롣
!
!ʿ��15ǯ6��19��(��)�ƤӺ�����⤦DEC�Ǥϼ¹Ԥ��ʤ���������
!#define MAXPAIR 1024*1024
  use common_module
  implicit none
  type sInteraction
     !sequence
     !���Х���٥��ȥ�
     real(kind=8),dimension(:),pointer :: ox,oy,oz
     !��ߺ��Ѥη������̾��1���ۼ�ʬ�Ҵ���ߺ��Ѥ�LB§���⶯�᤿���᤿�ꤹ��Τ˻��Ѥ��롣���ΤȤ�����StdInteraction��Ȥ��׻��ˤΤߺ��Ѥ��롣
     real(kind=8) :: scale
     
     logical :: isomol
     integer ::  nmol_i,nmol_j
     !����ɽ
     integer npair0
     integer,dimension(MAXPAIR) :: pair_i0,pair_j0
     !��ɽ
     integer :: Npair
     integer,dimension(MAXPAIR) :: pair_i,pair_j
     !
     !��ߺ��Ѹ��그�����̾��eratio=1d0, fratio=0d0������Ƥ�����
     !
     real(kind=8),dimension(MAXPAIR) :: eratio, fratio
     !������ݤΤ�����ܰ�
     integer :: maxpair
     !
     !partner������ˡ�neighbors�����Ѥ��Ƥߤ롣
     !
#ifdef VPOPTIMIZE
     !��ߺ��Ѥ�������ɽ�򤢤餫����������Ƥ�����
     !�Ϥ�׻��������ȡ�������¸�����Ȥޤǥ٥��ȥ벽���뤿�ᡣ
     !partner
     ! list�����������ˡ�ϡ��Ϥ����������ˤϤ��ޤ��Ψ���褯�ʤ���
     integer,dimension(MAXMol) :: npartner_i
     integer,dimension(MAXMol) :: npartner_j
#ifdef TEST7
     integer,dimension(MAXpartner,MAXmol) :: partner_i_,partner_j_
#else
     integer,dimension(MAXmol,MAXpartner) :: partner_i,partner_j
#endif /*TEST7*/
     integer :: maxpartner_i
     integer :: maxpartner_j
#endif /*VPOPTIMIZE*/
  end type sInteraction

contains
!�Ϥη׻���subroutine����kf90��Ȥ��Ф�����inlineŸ�����Ƥ����Τ�
!®�٤��㲼���ɤ��롣
!ʿ������ǯ�����(��)������ơ�������kap����ѥ����Ȥ�ʤ���
!inlineŸ���Ǥ��ʤ��Τ�����kap��Ȥ��ȷ�̤��礭�����ȤʤäƤ���褦����
!������α���롣
!vpp�Ǥ������ʤ���®�٤��㲼���ʤ���
!#define TEST12
#undef TEST12

!Ewald�Υ����ɤ�TEST12�����ˤ����������Ƥ��ʤ�(��ȴ��)
#ifdef EWALD
#define TEST12
#endif

#ifdef TEST12
!Must be inlined.
#ifdef EWALD
  subroutine Coulomb_Force(ewald,dx,dy,dz,a,e,fx,fy,fz)
#else
  subroutine Coulomb_Force(dx,dy,dz,a,e,fx,fy,fz)
#endif
# ifdef EWALD
    use ewald_module
    type(sEwald),intent(in) :: ewald
    real(kind=8) :: radius
#endif
    real(kind=8),intent(in) :: a
    real(kind=8),intent(out) :: e,fx,fy,fz
    real(kind=8),intent(in) :: dx,dy,dz
    real(kind=8) :: f0,dd,ddi,radiusi
    dd=dx**2+dy**2+dz**2
    ddi=1d0/dd
    radiusi=dsqrt(ddi)
#ifdef EWALD
    radius=radiusi*dd
    e=a*radiusi*erfc(ewald%kappa*radius)
    f0=(e+a*ewald%kpi2*dexp(-ewald%kapps * dd))*ddi
#else
    e = a*radiusi
    f0 = a*radiusi*ddi
#endif
    fx = dx*f0
    fy = dy*f0
    fz = dz*f0
  end subroutine Coulomb_Force
  
  subroutine LJ_Force(dx,dy,dz,a,b,e,fx,fy,fz)
    real(kind=8),intent(in) :: a,b
    real(kind=8),intent(out) :: e,fx,fy,fz
    real(kind=8),intent(in) :: dx,dy,dz
    real(kind=8) :: f0,dd,d6,d12
    dd = 1d0/(dx**2 + dy**2 + dz**2)
    d6 = dd*dd*dd
    d12= d6*d6
    e = a*d12 + b*d6
    f0 = dd*(a*12d0*d12 + b*6d0*d6)
    fx = dx*f0
    fy = dy*f0
    fz = dz*f0
  end subroutine LJ_Force
#endif /*TEST12*/
  
  !Register����Ѥ�����ߺ���ɽ��ƹ�����������ɬ���Ԥ�ɬ�פ�������
  subroutine Interaction_Reset(in)
    type(sInteraction),intent(INOUT) :: in
    in%npair=0
#ifdef VPOPTIMIZE
    in%maxpartner_i=0
    in%maxpartner_j=0
    in%npartner_j(:)=0
    in%npartner_i(:)=0
    !���ν������ճ��˽Ť���ɬ�׺Ǿ��¤ν�����ˤ���ˤϡ�
    in%partner_i(:,:)=0
    in%partner_j(:,:)=0
#endif /*VPOPTIMIZE*/
  end subroutine Interaction_Reset
  
#ifdef MPI
  subroutine interaction_alltoall(NPROCS,MYRANK,iv)
#else
  subroutine interaction_alltoall(iv)
#endif
#ifdef MPI
    integer,intent(IN) :: NPROCS,MYRANK
    integer :: node
#endif
    type(sInteraction),intent(inout) :: iv
    integer :: i,j,k,jmin,n
    if(iv%isomol)then
       n=iv%nmol_i*(iv%nmol_i-1)/2
    else
       n=iv%nmol_i*iv%nmol_j
    endif
    call Interaction_Reset(iv)
    !�ֲ��ץꥹ�ȥ٥��ȥ���������Ƥ�����
    k=0
    do i=1,iv%Nmol_i
       if(iv%isomol)then
          jmin=i+1
       else
          jmin=1
       endif
       do j=jmin,iv%Nmol_j
#if defined(MPI) || defined(PVM)
          !��롼�׳��μ�³���ʤΤǤ���������˽񤤤Ƥ�����
          if ( mod(node,NPROCS) == MYRANK ) then
             k=k+1
             iv%pair_i0(k)=i
             iv%pair_j0(k)=j
          endif
          node=node+1
#else /*MPI*/
          k=k+1
          iv%pair_i0(k)=i
          iv%pair_j0(k)=j
#endif /*MPI*/
       enddo
    enddo
    iv%npair0=k
  end subroutine interaction_alltoall
  
  
  subroutine interaction_initialize(iv,nmol_i,nmol_j)
    type(sInteraction),intent(out) :: iv
    integer,intent(in)             :: nmol_i
    integer,intent(in), optional   :: nmol_j
    iv%scale = 1d0
    iv%nmol_i = nmol_i
    !     if(present(nmol_j))then
    if(present(nmol_j) .and. nmol_j.ge.0)then
       iv%isomol=.false.
       iv%nmol_j = nmol_j
    else
       iv%isomol=.true.
       iv%nmol_j = nmol_i
    endif
    !�����list vector�Υ��������ǽ�Ū�˷��ꤷ�Ƥ�����ݤ��롣
    iv%maxpair=0
    if ( associated( iv%ox ) )then
       deallocate( iv%ox )
       deallocate( iv%oy )
       deallocate( iv%oz )
    endif
    allocate(iv%ox(MAXPAIR))
    allocate(iv%oy(MAXPAIR))
    allocate(iv%oz(MAXPAIR))
  end subroutine interaction_initialize
  
  !
  ! �Ϥη׻��ϡ��¤�rigid��flex������Ʊ��ǹ���ʤ���TIP4P�ξ��ʤ顢
  !     5�����ȥ�ǥ�ˤ��Ƥ��ޤ��Ф褤��    
  
  !     ��ߺ��Ѥε��Ҥ�ʬ�Ҥˤ�äưۤʤꡢ���̲��Ϻ���
  !     ����Ū���٤�
  !     �ü��®��
  !     ��ξ�������٤�褦�˥ץ���������Τ��褤�Ȼפ���
  
  !     �����Ǥϡ������ߺ��Ѥ�®��ͥ��ǵ��Ҥ��Ƥߤ롣
  
  !Homo��ߺ��Ѥξ��ϡ�����ȿ���Ѥ�Ʊ���ѿ��������Ф褤��
  !
  subroutine Interaction_Compress_obsolete(in,box,rr,xi,yi,zi,xj,yj&
       & ,zj,fNoCompression)
    use box_module
    use vector_module
!Grid_PR3��Grid_NeighborListHomo�θ�Ⱦ�򤽤Τޤޤ�äƤ���������Ǥ���
!�Τʤ鶦�̲����褦��
    type(sBox),intent(IN) :: box
    type(sInteraction),intent(INOUT) :: in
    logical,intent(in) :: fNoCompression
    real(kind=8),dimension(*),intent(in) :: xi,yi,zi,xj,yj,zj
    real(kind=8),intent(in) :: rr
    !integer,dimension(:),allocatable ::   pair_i,pair_j
    integer,dimension(:),allocatable :: pair_k
    integer,dimension(:),allocatable :: pair_k2
    integer :: i,j,k,newk
    real(kind=8) :: dd
#ifdef TEST30
    real(kind=8)              :: dx,dy,dz
#else
    type(vector3), dimension(:), allocatable :: cell
    type(vector3) :: delta
#endif
    call Interaction_Reset(in)
    !safely deallocated inside this routine.
    allocate(pair_k2(in%npair0))
    allocate(pair_k(in%npair0))
#ifdef TEST30
#else
    allocate(cell(in%npair0))
#endif
    !elseʸ��0�ˤ�������®���褦����
    !pair_k(1:npair)=0
#ifdef VERBOSE
    write(STDERR,*) "Uncompressed vector length:",in%npair0
#endif
    do k=1,in%npair0
       i=in%pair_i0(k)
       j=in%pair_j0(k)
#ifdef TEST30
       dx = xi(i)-xj(j)
       dy = yi(i)-yj(j)
       dz = zi(i)-zj(j)
#else
       delta%vec(1) = xi(i)-xj(j)
       delta%vec(2) = yi(i)-yj(j)
       delta%vec(3) = zi(i)-zj(j)
#endif
       !ʿ��13ǯ4��17��(��)���֥롼����
       !     ox(k) = dnint(dx*box%boxxi)*box%boxx
       !     oy(k) = dnint(dy*box%boxyi)*box%boxy
       !     oz(k) = dnint(dz*box%boxzi)*box%boxz
       !     dx = dx-ox(k)
       !     dy = dy-oy(k)
       !     dz = dz-oz(k)
#ifdef TEST30
       call box_renormalize0(box,dx,dy,dz,in%ox(k),in%oy(k),in%oz(k))
       dd = dx**2 + dy**2 + dz**2
#else
       call box_renormalize(box,delta,cell(k))
#ifdef TEST31
       dd = delta%vec(1)**2 + delta%vec(2)**2 + delta%vec(3)**2
#else
       dd = inner_product( delta, delta )
#endif
#endif
       if(dd < rr)then
          pair_k(k)=k
       else
          pair_k(k)=0
       endif
    enddo
    i=0
    !OCL VECTOR
    do k=1,in%npair0
       if(pair_k(k) /= 0)then
          i          = i + 1
          pair_k2(i) = pair_k(k)
       endif
    enddo
    deallocate(pair_k)
    in%npair=i
    !OCL VECTOR,NOVREC
    do newk=1,in%npair
       k=pair_k2(newk)
       in%pair_i(newk)=in%pair_i0(k)
       in%pair_j(newk)=in%pair_j0(k)
#ifdef TEST30
#else
       in%ox(newk)=cell(k)%vec(1)
       in%oy(newk)=cell(k)%vec(2)
       in%oz(newk)=cell(k)%vec(3)
#endif
    enddo
#ifdef VERBOSE
    write(STDERR,*) "Compressed vector length:",in%npair
#endif
    deallocate(pair_k2)  
#ifdef TEST30
#else
    deallocate(cell)
#endif
#ifdef VPOPTIMIZE
    !��������֤ˤĤ�Ƥ����ʤ��Ƥ⹽��ʤ����Ȥˤ����餺���֤�®���Ǥ�
    !�롣
    !OCL VECTOR
    if(fNoCompression)then
       do k=1,in%npair
          i=in%pair_i(k)
          j=in%pair_j(k)
          in%partner_i(i,j)=k
          in%partner_j(j,i)=k
       enddo
     !�ʲ��Ϥ����ǽ��������ɬ�פϤʤ����롼�פγ��ǽ��������Ф�����
       do i=1,in%nmol_i
          in%npartner_i(i)=in%nmol_j
       enddo
       do j=1,in%nmol_j
          in%npartner_j(j)=in%nmol_i
       enddo
       in%maxpartner_i=in%nmol_j
       in%maxpartner_j=in%nmol_i
    else
       !ʿ������ǯ�������(��)��ɤ�����ʬ���٤���
       do k=1,in%npair
          i=in%pair_i(k)
          in%npartner_i(i) = in%npartner_i(i)+1
          in%partner_i(i,in%npartner_i(i)) = k
       enddo
     !ʿ������ǯ�������(��)��ɤ�����ʬ���٤���
       do k=1,in%npair
          j=in%pair_j(k)
          in%npartner_j(j) = in%npartner_j(j)+1
          in%partner_j(j,in%npartner_j(j)) = k
       enddo
       in%maxpartner_i=0
       do i=1,in%nmol_i
          if(in%npartner_i(i).gt.in%maxpartner_i)then
             in%maxpartner_i=in%npartner_i(i)
          endif
       enddo
       in%maxpartner_j=0
       do j=1,in%nmol_j
          if(in%npartner_j(j).gt.in%maxpartner_j)then
             in%maxpartner_j=in%npartner_j(j)
          endif
       enddo
       !if(in%maxpartner_i.gt.MAXpartner)then
       !   write(STDERR,*) "Too many partners: ",in%maxpartner_i
       ! ,MAXpartner
       !   stop
       !endif
    endif
#endif
  end subroutine Interaction_Compress_obsolete

!��ߺ�����ɽ�Τ������ºݤ���ߺ��Ѥ��ʤ�(���åȥ��ճ���)�Ф�����
!��ߺ�����ɽ�򰵽̤��롣����򤹤뤳�Ȥǡ��Ǥ���֤Τ������Ϸ׻�����
!�η׻��̤�︺���뤳�Ȥ��Ǥ������Τ�®�٤����夹�롣
!An argument is added.
!subroutine Interaction_Compress(in,box,rr,xi,yi,zi,xj,yj,zj,fNoCompression)
  subroutine Interaction_Compress2(in,fbox,box,rr,comi,comj,fNoCompression)
    use vector_module
    use box_module
!Grid_PR3��Grid_NeighborListHomo�θ�Ⱦ�򤽤Τޤޤ�äƤ���������Ǥ���
!�Τʤ鶦�̲����褦��
    logical,intent(IN) :: fbox
    type(sBox),intent(IN) :: box
    type(sInteraction),intent(INOUT) :: in
    logical,intent(in) :: fNoCompression
    type(vector3),dimension(*),intent(in) :: comi,comj
    !rr<0�ξ��(̵�±���ߺ���)�Ͼ��ʬ�����ʤ���
    real(kind=8),intent(in) :: rr
    !integer,dimension(:),allocatable ::   pair_i,pair_j
    integer,dimension(:),allocatable :: pair_k
    integer,dimension(:),allocatable :: pair_k2
#ifdef TEST30
    real(kind=8)                     :: dx,dy,dz
#else
    type(vector3), dimension(:),allocatable :: cell
    type(vector3)                    :: delta
#endif
    integer :: i,j,k,newk
    real(kind=8)                     :: dd
    call Interaction_Reset(in)
    !safely deallocated inside this routine.
    allocate(pair_k2(in%npair0))
    allocate(pair_k(in%npair0))
#ifdef TEST30
#else
    allocate(cell(in%npair0))
#endif
    !elseʸ��0�ˤ�������®���褦����
    !pair_k(1:npair)=0
#ifdef VERBOSE
    write(STDERR,*) "Uncompressed vector length:",in%npair0
#endif
    if(rr.lt.0)then
       do k=1,in%npair0
          pair_k(k)=k
#ifdef TEST30
#else
          cell(k)%vec(:)=0d0
#endif
       enddo
    else
       if(fbox)then
          do k=1,in%npair0
             i=in%pair_i0(k)
             j=in%pair_j0(k)
#ifdef TEST30
             dx = comi(i)%vec(1) - comj(j)%vec(1)
             dy = comi(i)%vec(2) - comj(j)%vec(2)
             dz = comi(i)%vec(3) - comj(j)%vec(3)
             call box_renormalize0(box,dx,dy,dz,in%ox(k),in%oy(k),in%oz(k))
             dd = dx**2 + dy**2 + dz**2
#else
             delta%vec(:) = comi(i)%vec(:) - comj(j)%vec(:)
#ifdef TEST32
             cell(k)%vec(:)  = dnint( delta%vec(:) * box%invsize%vec(:) ) * box%size%vec(:)
             delta%vec(:) = delta%vec(:) - cell(k)%vec(:)
#else
             call box_renormalize(box,delta,cell(k))
#endif
#ifdef TEST31
             dd = delta%vec(1)**2 + delta%vec(2)**2 + delta%vec(3)**2
#else
             dd           = inner_product( delta, delta )
#endif
#endif
             !if ( i==1 .and. j==13 ) then
             !   write(STDERR,*) "PAIR(3)",i,j,delta%vec(1),cell(k)%vec(1)
             !endif
             if(dd < rr)then
                pair_k(k)=k
             else
                pair_k(k)=0
             endif
          enddo
       else
          do k=1,in%npair0
             i=in%pair_i0(k)
             j=in%pair_j0(k)
#ifdef TEST30
             dx = comi(i)%vec(1) - comj(j)%vec(1)
             dy = comi(i)%vec(2) - comj(j)%vec(2)
             dz = comi(i)%vec(3) - comj(j)%vec(3)
             dd = dx**2 + dy**2 + dz**2
#else
             delta%vec(:) = comi(i)%vec(:) - comj(j)%vec(:)
#ifdef TEST31
             dd = delta%vec(1)**2 + delta%vec(2)**2 + delta%vec(3)**2
#else
             dd           = inner_product( delta, delta )
#endif
#endif
             if(dd < rr)then
                pair_k(k)=k
             else
                pair_k(k)=0
             endif
#ifdef TEST30
             in%ox(k) = 0d0
             in%oy(k) = 0d0
             in%oz(k) = 0d0
#else
             cell(k)%vec(:) = 0d0
#endif
          enddo
       endif
    endif
    i=0
!OCL VECTOR
    do k=1,in%npair0
       if(pair_k(k).ne.0)then
          i=i+1
          pair_k2(i)=pair_k(k)
       endif
    enddo
    deallocate(pair_k)
    in%npair=i
!OCL VECTOR,NOVREC
    do newk=1,in%npair
       k=pair_k2(newk)
       in%pair_i(newk)=in%pair_i0(k)
       in%pair_j(newk)=in%pair_j0(k)
       !if(newk.le.15)write(STDERR,*) in%pair_i(newk),in%pair_j(newk)
#ifdef TEST30
       in%ox(newk)     = in%ox(k)
       in%oy(newk)     = in%oy(k)
       in%oz(newk)     = in%oz(k)
#else
       in%ox(newk)     = cell(k)%vec(1)
       in%oy(newk)     = cell(k)%vec(2)
       in%oz(newk)     = cell(k)%vec(3)
#endif
    enddo
#ifdef VERBOSE
    write(STDERR,*) "Compressed vector length:",in%npair
#endif
    deallocate(pair_k2)  
#ifdef TEST30
#else
    deallocate(cell)
#endif
#ifdef VPOPTIMIZE
    !
    !fNoCompression == false�ξ�硢����ߺ��Ѥ��Ǽ����2����������
    !�����ʤꤹ���ʤ��褦�ˡ����餬���ᰵ�̤��Ƥ�����512ʬ�Ҥ�ۤ����
    !���ʷϤǤϰ��̤��ʤ��ȥ����̵�Ѥ˿���������塢���׽������٤�
    !�ʤ롣(VP�ξ��)
    !
    !��������֤ˤĤ�Ƥ����ʤ��Ƥ⹽��ʤ����Ȥˤ����餺���֤�®���Ǥ�
    !�롣
!OCL VECTOR
    if(fNoCompression)then
       do k=1,in%npair
          i=in%pair_i(k)
          j=in%pair_j(k)
          in%partner_i(i,j)=k
          in%partner_j(j,i)=k
       enddo
     !�ʲ��Ϥ����ǽ��������ɬ�פϤʤ����롼�פγ��ǽ��������Ф�����
       do i=1,in%nmol_i
          in%npartner_i(i)=in%nmol_j
       enddo
       do j=1,in%nmol_j
          in%npartner_j(j)=in%nmol_i
       enddo
       in%maxpartner_i=in%nmol_j
       in%maxpartner_j=in%nmol_i
    else
       !ʿ������ǯ�������(��)��ɤ�����ʬ���٤���
       do k=1,in%npair
          i=in%pair_i(k)
          in%npartner_i(i) = in%npartner_i(i)+1
          in%partner_i(i,in%npartner_i(i)) = k
       enddo
       !ʿ������ǯ�������(��)��ɤ�����ʬ���٤���
       do k=1,in%npair
          j=in%pair_j(k)
          in%npartner_j(j) = in%npartner_j(j)+1
          in%partner_j(j,in%npartner_j(j)) = k
       enddo
       in%maxpartner_i=0
       do i=1,in%nmol_i
          if(in%npartner_i(i).gt.in%maxpartner_i)then
             in%maxpartner_i=in%npartner_i(i)
          endif
       enddo
       in%maxpartner_j=0
       do j=1,in%nmol_j
          if(in%npartner_j(j).gt.in%maxpartner_j)then
             in%maxpartner_j=in%npartner_j(j)
          endif
       enddo
       !if(in%maxpartner_i.gt.MAXpartner)then
       !   write(STDERR,*) "Too many partners: ",in%maxpartner_i
       ! ,MAXpartner
       !   stop
       !endif
    endif
#endif
    !
    !��ߺ��Ѹ���ؿ��˥ǥե�����ͤ��������롣
    !
    in%eratio( 1 : in%npair ) = 1d0
    in%fratio( 1 : in%npair ) = 0d0
  end subroutine Interaction_Compress2

!box��optional�ˤ��뤿��ˡ�Compress2�Ȥϰ����ν��֤����äƤ��뤳��
!������
  subroutine Interaction_Compress(in,comi,comj,fNoCompression,r,box)
    use vector_module
    use box_module
!Grid_PR3��Grid_NeighborListHomo�θ�Ⱦ�򤽤Τޤޤ�äƤ���������Ǥ���
!�Τʤ鶦�̲����褦��
    type(sBox),intent(IN),optional        :: box
    type(sInteraction),intent(INOUT)      :: in
    logical,intent(in)                    :: fNoCompression
    type(vector3),dimension(*),intent(in) :: comi,comj
    real(kind=8),intent(in)               :: r
    real(kind=8)                          :: rr
    if(r.le.0d0)then
       rr = -1d0
    else
       rr = r*r
    endif
    if(present(box))then
       call Interaction_Compress2(in,box%mode.ne.BOX_NONE,box,rr,comi,comj,fNoCompression)
    else
       call Interaction_Compress2(in,.false.,box,rr,comi,comj,fNoCompression)
    endif
  end subroutine Interaction_Compress

!Cluster�Ǥ����Ȥ�ʤ����������顢�Ƥ��Ȥ��˽񤯡��ϱշϤǤ��������ѹ�
!��ɬ�ס�
  function interaction_findisolated(iv,nmol,x,y,z)
    type(sInteraction),intent(in) :: iv
    real(kind=8),intent(in),dimension(*) :: x,y,z
    integer,intent(in) :: nmol
    real(kind=8) :: interaction_findisolated
    real(kind=8) :: dx,dy,dz,dd
    integer :: i,j,k
    !integer :: kk
    real(kind=8) :: r(MAXMOL),rmax
    !VPP�ѤΥ٥��ȥ벽�����
    !OCL VECTOR,REPEAT(1000000)
    do i=1,nmol
       r(i)=100d0
    enddo
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       dx = x(i)
       dy = y(i)
       dz = z(i)
       dx = dx - x(j)
       dy = dy - y(j)
       dz = dz - z(j)
       dd = dx**2+dy**2+dz**2
       if(r(i).gt.dd)r(i)=dd
       if(r(j).gt.dd)r(j)=dd
    enddo
    rmax=0d0
    do i=1,nmol
       if(rmax.lt.r(i))rmax=r(i)
    enddo
    interaction_findisolated=rmax
  end function interaction_findisolated

  subroutine Interaction_Show(iv)
    type(sInteraction),intent(in) :: iv
#ifdef VERBOSE
    write(6,*) iv%npair,"npair"
#ifdef VPOPTIMIZE
    write(6,*) iv%maxpartner_i,"maxpartner_i",iv%npartner_i(1)
    write(6,*) iv%maxpartner_j,"maxpartner_j",iv%npartner_j(1)
#endif
#endif
  end subroutine Interaction_Show

  !
  !O-H��Υ�ǡ�����Ƚ�ꤹ�롣
  !
  subroutine interaction_hb_by_oh_distance(iv,mi,mj,site,from,which,to)
    use mol_module
    use site_module
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
    integer :: pairMinimum,pair,whichH, whichMinimum
    do k=1,iv%Npair
       pair = 0
       pairMinimum = 0
       distanceMinimum = +1d10
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       !write(STDOUT,*) "i,j:",i,j
       do isite = 1, mi%nsite
          if ( mi%name(isite)(1:1) .eq. 'O' ) then
             whichH = 0
             do jsite = 1, mj%nsite
                if ( mj%name(jsite)(1:1) .eq. 'H' ) then
                   whichH = whichH + 1
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
                      whichMinimum    = whichH
                   endif
                endif
             enddo
          endif
       enddo
       do jsite = 1, mj%nsite
          if ( mj%name(jsite)(1:1) .eq. 'O' ) then
             whichH = 0
             do isite = 1, mi%nsite
                if ( mi%name(isite)(1:1) .eq. 'H' ) then
                   whichH = whichH + 1
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
                      whichMinimum    = whichH
                   endif
                endif
             enddo
          endif
       enddo
       !
       !determine by O-H distance
       !
       if ( distanceMinimum < 2.5d0**2 ) then
          if ( pairMinimum <= 2 ) then
             from( k )  = j
             which( k ) = whichMinimum
             to( k )    = i
          else 
             from( k )  = i
             which( k ) = whichMinimum
             to( k )    = j
          endif
       else
          from( k )  = 0
          which( k ) = 0
          to( k )    = 0
       endif
    enddo
  end subroutine interaction_hb_by_oh_distance

  !
  !graph(WGPH,NGPH)����Ϥ��롣�����Υե����ޥåȤϺ���ʪ�����ꤷ��
  !���ʤ���
  !
  subroutine interaction_SaveGraph(iv,mi,site,file,outputtype)
    use mol_module
    use site_module
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
  end subroutine interaction_SaveGraph
end module interaction_module
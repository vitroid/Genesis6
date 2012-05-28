! -*- f90 -*-
#define DEBUG
#undef DEBUG5
!
!ʿ��17ǯ3��10��(��)���������르�ꥺ��Ǥϥ���å��夹��ɬ�פϤۤȤ�ɤʤ���
!����å��夷�ʤ�����2/3���٤˷׻����֤�û�̤Ǥ��롣
!
#undef USE_Q6CACHE
!
!�ǽ��ǥХå�ʿ��17ǯ3��9��(��)���Υ�å�������ɽ��������ʤ�
!
#undef DEBUG20050309


!Q6��Umbrella MC��Ԥ�����λ�����
!�ҥ��ȥ�����Ȥ鷺�������ͤ򤽤Τޤ�OP�Ȥ��롣
!
!5A����Ρ����٤Ƥο��Ƿ��(��Υ3A�ʲ�)�����̤���Q6��׻����롣
!����ޤ�(ʿ��17ǯ2��23��(��))�ȤäƤ���Q6T�ϡ��ڤ�Ф�Ⱦ��6A�Ȥ��Ƥ�������
!�׻����֤�MC��ͭ���Թ�(��ʬ�׻��Τ���Υޡ�����ɬ�פˤʤ�)��û�̤�����
!
module q6local2
  !use set_module
  use neighbors_module
  !use distrib_module
  !use umbrella
  use common_module
  use set_module
#ifdef USE_Q6CACHE
#else
  use q6_module
#endif
  implicit none

  integer, parameter :: MAXBOND=10
  type sq6system
     type(sNeighbors) :: nei6(MAXMOL) ! �縵������ʬ��ɽ
     type(sNeighbors) :: nei3(MAXMOL) ! HB��Υ��nei6����������
     type(sNeighbors) :: nei5(MAXMOL) ! Q6�Υ�󥸡�nei6����������
     real(kind=8)     :: q6( MAXMOL ), sumq6

     !for 2
     complex(kind=8)  :: qlmsum( -6:6, MAXMOL )  ! foreach mol
     complex(kind=8)  :: qlm( -6:6, MAXBOND, MAXMOL ) ! foreach bond
     type(sSet)       :: neis( MAXMOL )          ! ����ʬ�Ҥ���5A�����ʬ��
     type(sSet)       :: neib( MAXBOND, MAXMOL )      ! �����礫��5A�����ʬ��
     integer          :: nbond( MAXMOL )         ! ����ʬ�Ҥ���5A����η��ο�
  end type sq6system

  type sQ6diff
     type(sNeighbors) :: newnei3(99)
     type(sNeighbors) :: newnei5(99)
     type(sNeighbors) :: newnei6(99)
     integer          :: newnei6idx(99)
     integer          :: nei6modified
     !for 2
     real(kind=8)     :: q6( MAXMOL ), sumq6
     complex(kind=8)  :: qlm( -6:6, MAXBOND ) ! foreach bond
     type(sSet)       :: neis            ! foreach mol
     type(sSet)       :: plus            ! foreach mol
     type(sSet)       :: minus           ! foreach mol
     type(sSet)       :: neib( MAXBOND )      ! foreach bond
  end type sQ6diff

  type sq6systemPtr
     type(sq6system), pointer :: p
  end type sq6systemPtr

  type sq6diffPtr
     type(sq6diff), pointer :: p
  end type sq6diffPtr

  interface new
     module procedure sq6system_new
  end interface

  interface done
     module procedure sq6system_done
  end interface

contains

  subroutine sq6system_new( tri )
    type( sq6system ), pointer :: tri
    allocate( tri )
    !
    !for qlmcache 10000�ǽ�ʬ��1000����­��ʤ����礭������Ȥޤ��٤��ʤ롣
    !10000����50000�����䤹�ȡ�isequal�θƽФ���������ΤǤ����ä��٤��ʤ뤫�⡣
    !10000����5000�˸��餷�Ƥ��Ϥ�isequal�θƽФ��������롣����migrate�Ǥ�����⡣
#ifdef USE_Q6CACHE
    call qlmcache_begin(5000)
#endif
  end subroutine sq6system_new

  subroutine sq6system_done( tri )
    type( sq6system ), pointer :: tri
    deallocate( tri )
  end subroutine sq6system_done


  subroutine sq6system_prepareall( q6s, nmol, ww, com, op )
    use interaction_module
    use vector_module
#ifdef USE_Q6CACHE
    use q6cache_module
#else
    use q6_module
#endif
    type(sq6system), pointer :: q6s
    integer, intent(IN) :: nmol
    type( sInteraction ), intent(IN) :: ww
    type(vector3), intent(IN) :: com( * )
    real(kind=8), intent(OUT), optional :: op

    complex(kind=8) :: dummy
    integer :: i,j,k,jj,kk,m
    real(kind=8) :: costh, newq6
    real(kind=8) :: pp, px, py, pz
    !
    !����ʬ�Ҥζ�˵ʬ�Ҥ�Ͽ���Ƥ�����
    !
    !call neighbors_all( nmol, q6s%nei3, ww, com, 0, 3d0, q6s%hbond )
    call neighbors_all2( nmol, q6s%nei6, ww, com, 0, ORDER_BY_LABEL, 6d0 )
    !
    !nei6�򸵤ˡ�ɬ�פ˱������������롣
    !
    do i=1, nmol
       call neighbor1_compress2( q6s%nei3(i), q6s%nei6(i), 0, ORDER_BY_LABEL, 3d0 )
       call neighbor1_compress2( q6s%nei5(i), q6s%nei6(i), 0, ORDER_BY_LABEL, 5d0 )
#ifdef DEBUG20050309
       write(6,*)"NNEAR", q6s%nei6(i)%nnear,q6s%nei5(i)%nnear,q6s%nei3(i)%nnear
#endif
    enddo
    !
    !���٤Ƥ�ʬ�ҤˤĤ���
    !
    do i=1,nmol
       !
       !����ʬ�Ҥζ�˵ʬ��(��ʬ��ޤ�)�ν�����롣
       !
       call new( q6s%neis(i), nmol, q6s%nei5(i)%nnear, q6s%nei5(i)%near, .false. )
       call join( q6s%neis(i), i )
    enddo
    do i=1,nmol
       !
       !ʬ��i�Τ��٤Ƥη��ˤĤ���
       !
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          !
          !ij�����˵�˴ޤ�ʬ�Ҥν�����롣
          !
          call intersection( q6s%neis(i), q6s%neis(j), q6s%neib( jj, i ) )
          !
          !ij����qlm�ͤ򤢤餫����׻����Ƥ�����
          !
          pp = 1d0 / sqrt( q6s%nei3(i)%dd(jj) )
          px = q6s%nei3(i)%dx(jj) * pp
          py = q6s%nei3(i)%dy(jj) * pp
          pz = q6s%nei3(i)%dz(jj) * pp
          do m=-6,6
#ifdef USE_Q6CACHE
             call qlm_cache( 6, m, px, py, pz, q6s%qlm( m, jj, i ) )
#else
             q6s%qlm( m, jj, i )  = qlm( 6, m, px, py, pz)
#endif
          enddo
          !write(6,*)q6s%qlm( -6, jj, i )
          !if ( i == 1 ) then
          !   write(STDERR,*) i,j,px,py,pz
          !endif
          !q6s%qlm( m, jj, i ) = qlm( 6, m, px, py, pz )
       enddo
    enddo

    !
    !qlmsum��ʬ�Ҥ��ȤΡ�qlm�ͤι�פ�׻����롣
    !
    q6s%qlmsum(:,:) = 0d0
    q6s%nbond(1:nmol) = 0
    !
    !�ºݤˤϡ��롼�פϷ�礴�Ȥ˲󤹡�
    !
    do i=1,nmol
       !
       !ʬ��i�Τ��٤Ƥη��ˤĤ���
       !
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          !
          !ij�����˵�˴ޤ�ʬ�ҤˤĤ���
          !
          do kk=1, q6s%neib( jj, i )%size
             k = q6s%neib( jj, i )%set(kk)
             q6s%nbond(k) = q6s%nbond(k) + 1
#ifdef DEBUG20050309
             if ( k==1 .and. i<j ) write(STDOUT,*) "NEIB",i,j !,q6s%qlm( -6, jj, i )
#endif
             do m=-6,6
                q6s%qlmsum(m,k) = q6s%qlmsum(m,k) + q6s%qlm( m, jj, i )
             enddo
          enddo
       enddo
    enddo
    !
    !Average all q6
    !
    q6s%sumq6 = 0d0
    do i=1, nmol
       newq6 = 0d0
       do m=-6,6
          newq6 = newq6 + abs(q6s%qlmsum(m,i))**2
       enddo
       newq6 = sqrt( newq6 * 4d0 * pi / ( 2d0 * 6 + 1d0 ) ) / q6s%nbond(i)
#ifdef DEBUG20050309
       write(6,*) "Q6:", q6s%nbond(i), newq6
#endif
       q6s%q6(i) = newq6
       q6s%sumq6 = q6s%sumq6 + q6s%q6(i)
    enddo
#ifdef USE_Q6CACHE
    call dict_stat()
#endif
    if ( present( op ) ) then
    !write(6,*) "q6:", q6s%sumq6 / nmol
       op = q6s%sumq6 / nmol
    endif
    !
    !�����ޤǤ�Q6Local.F90�ȴ�����Ʊ��ư����ǧ ʿ��17ǯ3��7��(��)
    !
  end subroutine sq6system_prepareall
    
  !
  !dirty�ʥ٥��ȥ�Τߤ�ꥹ�Ȥ�����ʬ�׻����ߤ롣ʿ��17ǯ3��2��(��)
  !
  subroutine sq6system_Differentiate( q6s, target, deltax, deltay, deltaz, nmol, q6d, op )
#ifdef USE_Q6CACHE
    use q6cache_module
#else
    use q6_module
#endif
    type( sq6system ) :: q6s
    type( sq6diff )   :: q6d
    integer, intent(IN) :: target
    real(kind=8), intent(IN) :: deltax,deltay,deltaz
    integer, intent(IN) :: nmol
    real(kind=8), intent(OUT) :: op
    
    integer :: i, j, k, ii, jj, kk, m
    real(kind=8) :: px, py, pz, pp, newq6

#ifdef DEBUG20050309
    call sq6system_write( q6s, nmol,STDOUT )
#endif
    !
    !��ư���ˡ�target�˿��Ƿ�礷�Ƥ�������
    !
#ifdef DEBUG20050309
    write(STDERR,*) "OLD NEAR", target, ":", ( q6s%nei3(target)%near(i), i=1, q6s%nei3(target)%nnear )
#endif
    !
    !�����η����˵�Ȥ���ʬ�ҤΡ�qlmsum����qlm�ͤ򺹤�����
    !
    do ii=1, q6s%nei3(target)%nnear
       i = q6s%nei3(target)%near(ii)
       !
       !target-i���ζ�˵��ʬ��k�ˤĤ���
       !
       do kk=1, q6s%neib( ii, target )%size
          k = q6s%neib( ii, target )%set(kk)
          !
          !nbond��qlmsum��q6s���ͤ�ľ���ѹ�����(reject���줿�������줬ɬ��)
          !
          q6s%nbond(k) = q6s%nbond(k) - 2
          do m=-6,6
             q6s%qlmsum(m,k) = q6s%qlmsum(m,k) - q6s%qlm( m, ii, target ) * 2
#ifdef DEBUG20050309
             if(  m==-6 ) then
                write(6,*) -target, i,k,q6s%qlm( m, ii, target )
             endif
#endif
          enddo
       enddo
    enddo
    !
    !��ư���target�ΰ��֤ǡ����ܥ٥��ȥ��������롣
    !
#ifdef DEBUG20050309
    write(6,*) "d",deltax,deltay,deltaz
#endif
    call neighbors_displacement( q6s%nei6, target, deltax, deltay, deltaz, q6d%nei6modified, q6d%newnei6, q6d%newnei6idx )
    do i=1, q6d%nei6modified
       call neighbor1_compress2( q6d%newnei3(i), q6d%newnei6(i), 0, ORDER_BY_LABEL, 3d0 )
       call neighbor1_compress2( q6d%newnei5(i), q6d%newnei6(i), 0, ORDER_BY_LABEL, 5d0 )
    enddo
    !
    !��ư��ˡ�target�˿��Ƿ�礷�Ƥ�������
    !
#ifdef DEBUG20050309
    write(STDERR,*) "NEW NEAR", target, ":", ( q6d%newnei3(1)%near(i), i=1, q6d%newnei3(1)%nnear )
#endif

    !
    !��ư��ˤ����뤳��
    !�ޤ�target����礹����꤬�ؤ�롣
    !�Ѳ������������ܤ���ʬ�ҤΥꥹ�Ȥ��ؤ��
    !����ʬ�Ҥ�qlmsum�ͤ��Ѳ����롣
    
    !target������ʬ�ҽ����ƹ������롣
    !�ۤ���ʬ�Ҥϡ���ư���Ƥ��ʤ��Τǡ��ۤ���ʬ�Ҥ�����ʬ�ҽ��������
    call new( q6d%neis, nmol, q6d%newnei5(1)%nnear, q6d%newnei5(1)%near, .false. )
    call join( q6d%neis, target )
    !
    !�����˶�˵�����äƤ���ʬ�Ҥν���
    !
    call difference( q6d%neis, q6s%neis(target), q6d%plus )
    !
    !��˵����ФƤ��ä�ʬ�Ҥν���
    !
    call difference( q6s%neis(target), q6d%neis, q6d%minus )

    !
    !�ޤ���target��ľ��ʬ�ҤȤδ֤ˤ����礬��ư���뤳�Ȥˤ���Ѳ��η׻�
    !
    !
    !1.target�ȡ�ľ������ʬ�Ҥ�����ʬ�ҽ���θ򺹤򹹿����롣
    do jj=1, q6d%newnei3(1)%nnear
       j = q6d%newnei3(1)%near(jj)
       call intersection( q6d%neis, q6s%neis(j), q6d%neib( jj ) )
       !
       !���η���qlm�ͤ򤢤餫����׻����Ƥ�����
       !
       pp = 1d0 / sqrt( q6d%newnei3(1)%dd(jj) )
       px = q6d%newnei3(1)%dx(jj) * pp
       py = q6d%newnei3(1)%dy(jj) * pp
       pz = q6d%newnei3(1)%dz(jj) * pp
       do m=-6,6
#ifdef USE_Q6CACHE
          call qlm_cache( 6, m, px, py, pz, q6d%qlm( m, jj ) )
#else
          q6d%qlm( m, jj )  = qlm( 6, m, px, py, pz )
#endif
#ifdef DEBUG20050309
          if ( m == -6 ) then
             write(6,*) target, j
             write(6,*) px,py,pz
             do k=1,q6s%nei3(target)%nnear
                if ( q6s%nei3(target)%near(k) == j ) then
                   write(6,*) q6s%nei3(target)%dx(k)/sqrt(q6s%nei3(target)%dd(k)),&
                        q6s%nei3(target)%dy(k)/sqrt(q6s%nei3(target)%dd(k)),&
                        q6s%nei3(target)%dz(k)/sqrt(q6s%nei3(target)%dd(k))
                endif
             enddo
          endif
#endif
       enddo
    enddo
    !
    !2.��ʬ�Ҥ�qlmsum�ˡ�������qlm�ͤ�û����롣
    !
    do jj=1, q6d%newnei3(1)%nnear
       !
       !j��target�ȿ��Ƿ�礷�Ƥ��롣
       !
       j = q6d%newnei3(1)%near(jj)
       !
       !target-j��礫��5A�����ʬ��k�ˤĤ���
       !
       do kk=1, q6d%neib( jj )%size
          k = q6d%neib( jj )%set(kk)
          q6s%nbond(k) = q6s%nbond(k) + 2
          do m=-6,6
             q6s%qlmsum( m, k ) = q6s%qlmsum( m, k ) + q6d%qlm( m, jj ) * 2
             !if( 32 == k .and. m==-6 ) then
#ifdef DEBUG20050309
             if(  m==-6 ) then
                write(6,*) target, j,k,q6d%qlm( m, jj )
             endif
#endif
          enddo
       enddo
    enddo
    !
    !��˵�˽и�����ʬ�Ҥϡ��⤷�����������������फ�⤷��ʤ���
    !����˵ʬ��i��{Q}�ˤĤ���
    do ii=1, q6d%plus%size
       i = q6d%plus%set(ii)
#ifdef DEBUG20050309
       write(6,*) "PLUS", i
#endif
       !
       !����ȷ�礹��ʬ��j�ˤĤ���
       !
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          !
          !�⤷��j��target��5A��˵�Ǥ���С�i-j���ζ�˵��target�����뤳�Ȥˤʤ롣
          !
          if ( set_exist( q6d%neis, j ) ) then
             !
             !���ξ��ϡ�target��q�ͤ˲û����롣
             !
             q6s%nbond(target) = q6s%nbond(target) + 2
             do m=-6,6
                q6s%qlmsum( m, target ) = q6s%qlmsum( m, target ) + q6s%qlm( m, jj, i ) * 2
             enddo
          endif
       enddo
    enddo
    !
    !��˵����ü�����ʬ�Ҥϡ�����μ¤˼�����
    !�ü�ʬ��i��{P}�ˤĤ���
    do ii=1, q6d%minus%size
       i = q6d%minus%set(ii)
#ifdef DEBUG20050309
       write(6,*) "MINUS", i
#endif
       !
       !����ȷ�礹��ʬ��j��{P,R}�ˤĤ���
       !
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          !
          !�⤷��j���Ȥ��target��5A��˵�Ǥ���С�i-j���ζ�˵��target�����ä����Ȥˤʤ롣
          !
          if ( exist( q6s%neis(target), j ) ) then
             !
             !���ξ��ϡ�target��q�ͤ��鸺�����롣
             !
#ifdef DEBUG20050309
             write(6,*) "MINUS-", j
#endif
             q6s%nbond(target) = q6s%nbond(target) - 2
             do m=-6,6
                q6s%qlmsum( m, target ) = q6s%qlmsum( m, target ) - q6s%qlm( m, jj, i ) * 2
             enddo
          endif
       enddo
    enddo
    !
    !Average all q6
    !
    q6d%sumq6 = 0d0
    do i=1, nmol
       newq6 = 0d0
       do m=-6,6
          newq6 = newq6 + abs(q6s%qlmsum(m,i))**2
       enddo
       newq6 = sqrt( newq6 * 4d0 * pi / ( 2d0 * 6 + 1d0 ) ) / q6s%nbond(i)
#ifdef DEBUG20050309
       write(6,*) i, "nbond", q6s%nbond(i),"Q6", newq6
#endif
       q6d%q6(i) = newq6
       q6d%sumq6 = q6d%sumq6 + q6d%q6(i)
    enddo
    
    op = q6d%sumq6 / nmol
    do ii=1, q6d%nei6modified
       i = q6d%newnei6idx(ii)
#ifdef DEBUG20050309
       write(STDOUT,'(a3,2i3,f8.5)') "Q6 ",target,i,q6d%q6(i)
#endif
    enddo

  end subroutine sq6system_Differentiate
  
  subroutine sq6system_Integrate( q6s, q6d, nmol, target )
    type( sq6system ) :: q6s
    type( sq6diff )   :: q6d
    integer, intent(IN) :: nmol, target
    
    type(sSet) :: oldnei3,newnei3,lost,kept,found
    integer :: i,j,k,ii,jj,kk
    complex(kind=8) :: qlmtemp( -6:6, MAXBOND ), qlmtemp2( -6:6 )
#ifdef DEBUG5
    type(sNeighbors) :: test(MAXMOL)
#endif    
    !
    !q6�򹹿�����
    !
    do i=1, nmol
       q6s%q6(i) = q6d%q6(i)
    enddo
    q6s%sumq6 = q6d%sumq6
    !
    !target���Ȥ�qlm����򹹿����롣
    !
    do jj=1, q6d%newnei3(1)%nnear
       !q6d%qlm( m, jj )��target��jj���ܤο��Ƿ���qlm��
       !
       !j��target�ȿ��Ƿ�礷�Ƥ��롣
       !
#ifdef DEBUG20050309
       write(6,*) "UPDATED", target, q6d%newnei3(1)%near(jj)
#endif
       q6s%qlm( -6:6, jj, target ) = q6d%qlm( -6:6, jj )
    enddo
    !
    !j���٤�ˤĤ��ơ�
    !target��ư���Τǡ���礬�ڤ줿�ꤢ�餿�˷�礷���ꤹ�롣
    !���˷�礷�Ƥ�����Τǡ���礬�ڤ줿����qlmɽ��Ĥ�ʤ���Ф����ʤ�����
    !�����˷�礷������qlmɽ���ɲä��ʤ���Ф����ʤ���
    !�ɤ���ξ��⡢newnei3��ν����nei3��ν���Ȥϲ���ط����ʤ��Τǡ�
    !�¤��ؤ��򤪤��ʤ�ʤ���Ф����ʤ���
    !
    !
    !neib������˹�������ΤϤ��ʤ����ѡ�neis����������褦��
    !
    !1.target���Ȥ�neis�ι���
    !
    q6s%neis(target) = q6d%neis
    !
    !2.target�ΰ�ư�ˤ�ꡢtarget����������ä�/�ФƤ��ä����ձ�����ʬ�Ҥ�neis���Ѳ����롣
    !
    do ii=1, q6d%plus%size
       !
       !i�Ͽ����˶�˵�ˤʤä�ʬ�ҡ�
       !
       i = q6d%plus%set(ii)
       !target��i�ζ�˵����Ͽ���롣
       call join( q6s%neis(i), target )
    enddo
    do ii=1, q6d%minus%size
       !
       !i�϶�˵����ž�Ф���ʬ�ҡ�
       !
       i = q6d%minus%set(ii)
       !target��i�ζ�˵���������
       call drop( q6s%neis(i), target )
    enddo
    !
    !neib�Ͻ���黻�����ʤΤǽ�����®����
    !�����ǺƷ׻����Ƥ��ޤ���
    !
    !target�˷�礷�Ƥ���������ư�ˤ���ڤ�Ƥ��ޤä�����neib
    !
    call new( oldnei3, MAXBOND, q6s%nei3(target)%nnear,q6s%nei3(target)%near,.false.)
    call new( newnei3, MAXBOND, q6d%newnei3(1)%nnear,  q6d%newnei3(1)%near,  .false.)
    call difference( oldnei3, newnei3, lost )
    call difference( newnei3, oldnei3, found )
    call intersection( oldnei3, newnei3, kept )
    !���Ȥη׻��ϡ�nei3�ι����Τ��ȤΤۤ��������ʤΤǤ������ǹԤ���

    !
    !�����˶�˵�����ä�/�Ф�ʬ�ҤΡ�neib�򹹿����롩
    !
    do ii=1, q6d%plus%size
       !
       !i�Ͽ����˶�˵�ˤʤä�ʬ�ҡ�
       !
       i = q6d%plus%set(ii)
       do jj=1, q6s%nei3(i)%nnear
          !
          !j��i�ȷ�礷�Ƥ���ʬ��
          !
          j = q6s%nei3(i)%near(jj)
          !
          !�⤷j��target�ζ�˵�Ǥ����
          !
          if ( exist( q6d%neis, j ) ) then
             !
             !���i-j��neib�ˤϡ�target��ޤࡣ
             !
             call join( q6s%neib( jj, i ), target )
             kk = sNeighbor1_Lookup( q6s%nei3(j), i )
             if ( 0 < kk ) then
                call join( q6s%neib( kk, j ), target )
             endif
          endif
       enddo
    enddo
    do ii=1, q6d%minus%size
       !
       !i�϶�˵����ž�Ф���ʬ��
       !
       i = q6d%minus%set(ii)
       do jj=1, q6s%nei3(i)%nnear
          !
          !j��i�ȷ�礷�Ƥ���ʬ��
          !
          j = q6s%nei3(i)%near(jj)
          !
          !j����target�ΰ�ư���ˤ϶�˵�ˤ��ä��ʤ�
          !
          if ( exist( q6s%neis(target), j ) ) then
             call drop( q6s%neib( jj, i ), target )
             kk = sNeighbor1_Lookup( q6s%nei3(j), i )
             if ( 0 < kk ) then
                call drop( q6s%neib( kk, j ), target )
             endif
          endif
       enddo
    enddo

    !
    !����ɽ�򹹿����롣
    !
    do i=1, q6d%nei6modified
       j = q6d%newnei6idx(i)
       q6s%nei6(j) = q6d%newnei6(i)
    enddo


    
    !
    !G��{lost}�Ϥ�Ȥ��target�ȷ�礷������ư�ˤ���礬�ڤ줿�����η���qlm�Ͻ�����G�λĤ�η���qlm��ƹ��ۤ���ɬ�פ����롣target�Ȥη��ʳ���ư���Ƥ��ʤ��Τǡ��ͤϵ줤�ͤ򤽤Τޤ޻Ȥ���Ϥ��������б��ط��������Ƥ���Τ��񵷡�
    !
    !
    !target������ʬ�Ҥ�qlm�ι�����nei3����٥��˥����Ȥ���Ƥ����ΤȤߤʤ��ƽ��������ˤ��롣
    !
    !�ޤ�����Ȥ�ȷ�礷�Ƥ���������礬�ڤ�Ƥ��ޤä���lost�Υ��С��ˤĤ���
    !
#ifdef DEBUG20050309
    write(6,*) "oldnei3"
    call sset_write(oldnei3, stdout )
    write(6,*) "newnei3"
    call sset_write(newnei3, stdout )
    write(6,*) lost%size, "lost, ", kept%size, "kept, ", found%size, " found."
#endif
    do ii=1, lost%size
       i = lost%set(ii)
       jj = 1
       !
       !nnear�ϰ�ư��������ʬ�ҿ���
       !

       do kk=1, q6s%nei3(i)%nnear
          k = q6s%nei3(i)%near(kk)
#ifdef DEBUG20050309
          write(6,*) "LOST", i,k
#endif
          if ( k /= target ) then
             q6s%qlm( -6:6, jj, i ) = q6s%qlm( -6:6, kk, i )
             jj = jj + 1
          endif
       enddo
    enddo
    !
    !��礬�ڤ��Ѥ���⤷�ʤ�����ʬ�ҤˤĤ���
    !
    do ii=1, kept%size
       i = kept%set(ii)
       !
       !target-i����qlm�򤢤餫����Ĵ�٤Ƥ�����
       !
       jj = sNeighbor1_Lookup( q6d%newnei3(1), i )
       if ( 0 < jj ) then
          qlmtemp2( -6:6 ) = q6d%qlm( -6:6, jj )
       endif
       kk = sNeighbor1_Lookup( q6s%nei3(i), target )
       if ( 0 < kk ) then
          k = q6s%nei3(i)%near(kk)
#ifdef DEBUG20050309
          write(6,*) "KEPT", i,k
#endif
          !
          !target�����Ѥ˷׻�����qlm�ͤ���Ѥ��롣
          !
          q6s%qlm( -6:6, kk, i ) = qlmtemp2( -6:6 )
       endif
    enddo
    !
    !�����˷�礬�Ǥ������ϡ�
    !qlm�ε줤�ͤ򥳥ԡ����Ĥġ�target�Ȥη��˴ؤ��Ƥ�q3d���ͤ�Ȥ���
    !
    !
    !i�ϡ�target��6A��˵��ʬ��
    !
    do ii=1, q6d%nei6modified
       i = q6d%newnei6idx(ii)
       !
       !�⤷i��{found}�ʤ�
       !
       if ( exist( found, i ) ) then
          !qlm������֤���Τǡ���������ƨ������
          !
          qlmtemp( -6:6, 1:q6d%newnei3(ii)%nnear-1 ) = q6s%qlm( -6:6, 1:q6d%newnei3(ii)%nnear-1, i )
          !
          !target-i����qlm�򤢤餫����Ĵ�٤Ƥ�����
          !
          jj = sNeighbor1_Lookup( q6d%newnei3(1), i )
          if ( 0 < jj ) then
             qlmtemp2( -6:6 ) = q6d%qlm( -6:6, jj )
          endif
          !
          !��ư������֤Ǥ�i������ʬ��{k}(target��ޤ�)�ˤĤ���
          !
          jj = 1
          do kk=1, q6d%newnei3(ii)%nnear
             k = q6d%newnei3(ii)%near(kk)
#ifdef DEBUG20050309
             write(6,*) "FOUND", i,k
#endif
             !
             !��Ȥ�Ȥ�target�����ܤ��Ƥ��ʤ��ä��Τǡ�qlm�ͤϷ׻�����Ƥ��ʤ�
             !
             if ( k /= target ) then
#ifdef DEBUG20050309
                write(6,*) "(FOUND)=", kk,i,jj
#endif
                q6s%qlm( -6:6, kk, i ) = qlmtemp( -6:6, jj )
                jj = jj + 1
             else
                !
                !target�����Ѥ˷׻�����qlm�ͤ���Ѥ��롣
                !
#ifdef DEBUG20050309
                write(6,*) "(FOUND)+", kk,i
#endif
                q6s%qlm( -6:6, kk, i ) = qlmtemp2( -6:6 )
             endif
          enddo
       endif
    enddo
                
    !
    !nei3�򹹿����롣
    !��������ɬ�פΤʤ���Τޤǹ������Ƥ��롣
    !
    do i=1, q6d%nei6modified
       j = q6d%newnei6idx(i)
       q6s%nei3(j) = q6d%newnei3(i)
    enddo
    do i=1, q6d%nei6modified
       j = q6d%newnei6idx(i)
       q6s%nei5(j) = q6d%newnei5(i)
    enddo

    !
    !��ư��κ�ɸ�ǡ�target��ľ�뤷��ʬ�Ҥ�neib
    !nei3����������Ƥ���Ǥʤ������ݤʤΤǤ����ˤ�äƤ�����
    !
    do jj=1, q6s%nei3(target)%nnear
       j = q6s%nei3(target)%near(jj)
       call intersection( q6s%neis(target), q6s%neis(j), q6s%neib( jj, target) )
       do kk=1, q6s%nei3(j)%nnear
          k = q6s%nei3(j)%near(kk)
          call intersection( q6s%neis(j), q6s%neis(k), q6s%neib( kk, j) )
       enddo
    enddo
    !
    !��Ȥ��target�ȷ�礷�Ƥ�������ư�ˤ���礬�ڤ줿����neib���������ɬ�פ����롣
    !
    do ii=1, lost%size
       i = lost%set(ii)
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          call intersection( q6s%neis(i), q6s%neis(j), q6s%neib( jj, i) )
       end do
    end do
  end subroutine sq6system_Integrate

  subroutine sq6system_Reject( q6s, q6d, target, nmol )
    type( sq6system ) :: q6s
    type( sq6diff )   :: q6d
    integer, intent(IN) :: target, nmol

    integer :: i,j, k,m, ii,jj,kk

    !
    !qlmsum������
    !
    do jj=1, q6d%newnei3(1)%nnear
       j = q6d%newnei3(1)%near(jj)
       do kk=1, q6d%neib( jj )%size
          k = q6d%neib( jj )%set(kk)
          q6s%nbond(k) = q6s%nbond(k) - 2
          do m=-6,6
             q6s%qlmsum( m, k ) = q6s%qlmsum( m, k ) - q6d%qlm( m, jj ) * 2
          enddo
       enddo
    enddo
    do ii=1, q6s%nei3(target)%nnear
       i = q6s%nei3(target)%near(ii)
       !
       !target-i���ζ�˵��ʬ�ҤˤĤ���
       !
       do kk=1, q6s%neib( ii, target )%size
          k = q6s%neib( ii, target )%set(kk)
          q6s%nbond(k) = q6s%nbond(k) + 2
          do m=-6,6
             q6s%qlmsum(m,k) = q6s%qlmsum(m,k) + q6s%qlm( m, ii, target ) * 2
          enddo
       enddo
    enddo
    do ii=1, q6d%plus%size
       i = q6d%plus%set(ii)
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          if ( set_exist( q6d%neis, j ) ) then
             q6s%nbond(target) = q6s%nbond(target) - 2
             do m=-6,6
                q6s%qlmsum( m, target ) = q6s%qlmsum( m, target ) - q6s%qlm( m, jj, i ) * 2
             enddo
          endif
       enddo
    enddo
    do ii=1, q6d%minus%size
       i = q6d%minus%set(ii)
       do jj=1, q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          if ( set_exist( q6s%neis(target), j ) ) then
             q6s%nbond(target) = q6s%nbond(target) + 2
             do m=-6,6
                q6s%qlmsum( m, target ) = q6s%qlmsum( m, target ) + q6s%qlm( m, jj, i ) * 2
             enddo
          endif
       enddo
    enddo
#ifdef DEBUG20050309
    call sq6system_write( q6s, nmol, STDOUT )
    write(STDOUT,*) "######################################"
#endif
  end subroutine sq6system_Reject

  subroutine sq6system_EndAStep( q6s )
    type( sq6system ) :: q6s
  end subroutine sq6system_EndAStep

  subroutine sQ6System_Write( q6s, nmol, file )
    use neighbors_module
    type( sq6system ) :: q6s
    integer :: file
    integer :: nmol

    integer :: i,j,k,ii,jj,kk
    write(file,*) "nei6========================="
    do i=1, nmol
       write(file,*) i
       call sneighbor1_write( q6s%nei6(i), file )
    enddo
    write(file,*) "nei5========================="
    do i=1, nmol
       write(file,*) i
       call sneighbor1_write( q6s%nei5(i), file )
    enddo
    write(file,*) "nei3========================="
    do i=1, nmol
       write(file,*) i
       call sneighbor1_write( q6s%nei3(i), file )
    enddo
    write(file,*) "q6()=========================", q6s%sumq6
    do i=1, nmol
       write(file,*) i, q6s%q6(i)
    enddo
    write(file,*) "qlmsum()========================="
    do i=1, nmol
       write(file,*) i, q6s%qlmsum( -6, i )
    enddo
    write(file,*) "qlm()========================="
    do i=1, nmol
       do j=1,q6s%nei3(i)%nnear
          write(file,*) i, j, q6s%qlm( -6, j, i )
       enddo
    enddo
    write(file,*) "neis()========================="
    do i=1, nmol
       write(file,*) i
       call sset_write( q6s%neis(i), file )
    enddo
    write(file,*) "neib()========================="
    do i=1, nmol
       write(file,*) i
       do j=1,q6s%nei3(i)%nnear
          write(file,*) i,j
          call sset_write( q6s%neib(j,i), file )
       enddo
    enddo
    write(file,*) "nbond()========================="
    do i=1, nmol
       write(file,*) i,q6s%nbond(i)
    enddo
    !
    !����ʬ�Ҥ�q6�˴�Ϳ���Ƥ�����Υꥹ�Ȥ��ۤ�����
    !
    do i=1, nmol
       do jj=1,q6s%nei3(i)%nnear
          j = q6s%nei3(i)%near(jj)
          do kk=1, q6s%neib(jj,i)%size
             k = q6s%neib(jj,i)%set(kk)
             if ( k == 5 ) then
                write(file,*) i,j,k
                !if ( ( i==4 .and. j==34 ) .or. ( i==34 .and. j==4 ) ) then
                !   write(file,*) "4,17",distance(q6s,4,17)
                !   write(file,*) "34,17",distance(q6s,34,17)
                !   write(file,*) "4,34",distance(q6s,4,34)
                !endif
             endif
          enddo
       enddo
    enddo
  end subroutine sQ6System_Write

  function distance( q6s, i, j )
    type( sq6system ) :: q6s
    integer :: i,j,k,kk
    real(kind=8) :: distance
    distance = 0d0
    do kk=1, q6s%nei6(i)%nnear
       k = q6s%nei6(i)%near(kk)
       if ( k == j ) then
          distance = sqrt(q6s%nei6(i)%dd(kk))
       endif
    enddo
  end function distance

  subroutine sQ6System_Compare( q1, q2, nmol, file )
    use neighbors_module
    type( sq6system ) :: q1,q2
    integer :: file, nmol

    integer :: i,j,k,ii,jj,kk
    write(file,*) "nei6++++++++++++++++++++++++="
    !
    !nei6�ϡ����⤽��˾��ӤȤ����ߤ��Ƥ���Τǡ������԰��פǤ⹽���
    !����nei6����nei5��neis�����������������Ф褤��
    !
    do i=1, nmol
       write(file,*) i
       call sneighbor1_compare( q1%nei6(i), q2%nei6(i), file )
    enddo
    write(file,*) "nei5++++++++++++++++++++++++="
    do i=1, nmol
       write(file,*) i
       call sneighbor1_compare( q1%nei5(i), q2%nei5(i), file )
    enddo
    write(file,*) "nei3++++++++++++++++++++++++="
    do i=1, nmol
       write(file,*) i
       call sneighbor1_compare( q1%nei3(i), q2%nei3(i), file )
    enddo
    write(file,*) "q6()++++++++++++++++++++++++="
    if ( 1d-12 < real_compare( q1%sumq6, q2%sumq6 ) ) then
       write(file,*) real_compare( q1%sumq6, q2%sumq6 ), "differs(sumq6)."
    endif
    do i=1, nmol
       if ( 1d-12 < real_compare( q1%q6(i), q2%q6(i) ) ) then
          write(file,*) i, real_compare( q1%q6(i), q2%q6(i) ), "differs(q6(i))."
       endif
    enddo
    write(file,*) "qlmsum()++++++++++++++++++++++++="
    do i=1, nmol
       if ( 1d-12 < complex_compare( q1%qlmsum( -6, i), q2%qlmsum( -6, i) ) ) then
          write(file,*) i, complex_compare( q1%qlmsum( -6, i), q2%qlmsum( -6, i) ), "differs(qlmsum( -6, i))."
       endif
    enddo
    write(file,*) "qlm()++++++++++++++++++++++++="
    do i=1, nmol
       k = min(q1%nei3(i)%nnear,q2%nei3(i)%nnear)
       do jj=1,k
          if ( 1d-12 < complex_compare( q1%qlm( -6, jj, i), q2%qlm( -6, jj, i) ) ) then
             write(file,*) i, jj, q1%nei3(i)%near(jj), q2%nei3(i)%near(jj), complex_compare( q1%qlm( -6, jj, i), q2%qlm( -6, jj, i) ), "differs(qlm( -6,jj,i)).",q1%qlm( -6, jj, i), q2%qlm( -6, jj, i)
          endif
       enddo
    enddo
    !      write(file,*) i, j, q6s%qlm( -6, j, i )
    !   enddo
    !enddo
    write(file,*) "neis()++++++++++++++++++++++++="
    do i=1, nmol
       write(file,*) i
       call sset_compare( q1%neis(i), q2%neis(i), file )
    enddo
    write(file,*) "neib()++++++++++++++++++++++++="
    do i=1, nmol
       write(file,*) i
       write(file,*) q1%nei3(i)%nnear, q2%nei3(i)%nnear, "(nei3(i)%nnear)."
       if ( q1%nei3(i)%nnear /= q2%nei3(i)%nnear ) then
          write(file,*) "differs."
       endif
       k = min(q1%nei3(i)%nnear,q2%nei3(i)%nnear)
       do j=1,k
          write(file,*) i,j,q1%nei3(i)%near(j),q2%nei3(i)%near(j)
          call sset_compare( q1%neib(j,i), q2%neib(j,i), file )
       enddo
    enddo
    write(file,*) "nbond()++++++++++++++++++++++++="
    do i=1, nmol
       if ( q1%nbond(i) /= q2%nbond(i) ) then
          write(file,*) q1%nbond(i), q2%nbond(i), " differs(nbond(i))."
       endif
    enddo
  end subroutine sQ6System_Compare

  function real_compare( x,y )
    real(kind=8) :: x,y,real_compare
    real_compare = abs((x-y)/y)
  end function real_compare

  function complex_compare( x,y )
    complex(kind=8) :: x,y
    real(kind=8) :: complex_compare
    complex_compare = abs((x-y)/y)
  end function complex_compare

end module q6local2
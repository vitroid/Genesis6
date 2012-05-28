! -*- f90 -*-
!���γ�(O-O-O)��MC��Ԥ�����λ�����
!
module triangle
  use set_module
  use neighbors_module
  use distrib_module
  use umbrella
  implicit none

  type sTriangle
#ifdef DEBUG
     type(sHistogram), pointer :: lasthist
#ifdef DEBUG2
     type(sHistogram), pointer :: tmphist
     type(sNeighbors) :: nei2copy(MAXMOL)
#endif
#endif
     type(sHistogram), pointer :: hist, newhist
     type(sSet) :: set
     type(sNeighbors) :: nei(MAXMOL)  ! ������ϺǶ���4ʬ�Ҥ�����
     type(sNeighbors) :: nei2(MAXMOL) ! �������5A����ζ�˵�ˤ���ʬ�ҡ�
     !
     !temporary neighbor list
     !
     type(sNeighbors) :: newnei(99)   !newnei2��Ʊ���礭����ɬ�ס�
     integer          :: newnei2idx(99)
     type(sNeighbors) :: newnei2(99)
     integer          :: nmodified
     logical          :: logarithmic
  end type sTriangle


  type sTrianglePtr
     type(sTriangle), pointer :: p
  end type sTrianglePtr

  interface new
     module procedure sTriangle_new
  end interface

contains

  subroutine sTriangle_new( tri, islog )
    type( sTriangle ), pointer :: tri
    logical, intent(IN) :: islog
    allocate( tri )
    tri%logarithmic = islog
    nullify( tri%hist )
    nullify( tri%newhist )
  end subroutine sTriangle_new

  subroutine sTriangle_prepareall( tri, nmol, ww, com, refhist, op )
    use interaction_module
    use vector_module
    type( sTriangle ) :: tri
    integer, intent(IN) :: nmol
    type( sInteraction ), intent(IN) :: ww
    type(vector3), intent(IN) :: com( * )
    type( sHistogram ), pointer :: refhist
    real(kind=8), intent(OUT) :: op

    integer :: i,j,k
    real(kind=8) :: costh
    !
    !����ʬ�Ҥζ�˵4ʬ�Ҥ��������٤�ʬ�ۤ���롣
    !1��ʬ�Τߡ�
    !
    call new( tri%hist, 21, 0.1d0, -1d0 )
    !
    !5A����ˤ�������ʬ�Ҥ�ɽ��1 MC step�δְݻ����롣
    !��˵ʬ�Ҥ����ޤ�󤯤ˤ��äƤ��ޤ�ʤ����Ȥ�����
    !
#ifdef DEBUG2
    if ( loop .ne. 1 ) then
       do i=1, nmol
          tri%nei2copy(i) = tri%nei2(i)
       enddo
    endif
#endif
    call neighbors_all2( nmol, tri%nei2, ww, com, 0, ORDER_BY_DISTANCE, 5d0 )
    !
    !nei2�򸵤ˡ���ʬ�Ҥ˺Ƕ��ܤ�4ʬ�Ҥ�ɽ���������롣ɬ�פ˱������������롣
    !
    do i=1, nmol
       call neighbor1_compress2( tri%nei(i), tri%nei2(i), 4, ORDER_BY_DISTANCE )
    enddo
    do i=1, nmol
       do j=1, tri%nei(i)%nnear
          do k=j+1, tri%nei(i)%nnear
             !
             ! O-O-O�Ѥ�;����׻����롣
             !
             costh = tri%nei(i)%dx(j)*tri%nei(i)%dx(k) + tri%nei(i)%dy(j)*tri%nei(i)%dy(k) + tri%nei(i)%dz(j)*tri%nei(i)%dz(k)
             costh = costh / sqrt ( tri%nei(i)%dd(j) * tri%nei(i)%dd(k) )
             call histogram_accumulate( tri%hist, costh, 1d0 )
          enddo
       enddo
    enddo
    op = kullback_div( tri%hist, refhist )
    !call histogram_show( tri%hist, STDOUT )
    !write(6,*) op, log(op)
    if ( tri%logarithmic ) op = log(op)
    !stop
#ifdef DEBUG
    !
    !������ѻ��ҥ��ȥ����ȡ������ǺƷ׻������ҥ��ȥ���ब���פ��뤫�ɤ���
    !tmphist�Ǵ����˰��פ��Ƥ���(accept����histogram�ι���������ʤ�)�Τˡ�
    !�������԰��פ�������Ȥ���С�nei2���礭���Ѳ��������Ƥ����ǽ�������롣
    !
    if ( associated( tri%lasthist ) ) then
       diff = kullback_div( tri%hist, tri%lasthist )
       if ( diff /= 0d0 ) then
          write(STDERR,*) diff, "DIV"
          call histogram_show( tri%hist, STDERR )
          call histogram_show( tri%lasthist, STDERR )
          call die( 0, "Histogram differs." )
       endif
       call done( tri%lasthist )
    endif
#endif /*DEBUG*/
  end subroutine sTriangle_prepareall
    

  subroutine sTriangle_Differentiate( tri, target, deltax, deltay, deltaz, refhist, op )
    type( sTriangle ) :: tri
    integer, intent(IN) :: target
    real(kind=8), intent(IN) :: deltax,deltay,deltaz
    type( sHistogram ), pointer :: refhist
    real(kind=8), intent(OUT) :: op
    
    integer :: i, j, k, ii
    logical :: verbose
    real(kind=8) :: costh
    !
    !�ޤ���5A����ζ�˵ʬ�ҥꥹ�Ȥ򹹿����롣tri%newnei2�ˤϡ���ư����ʬ�ξ���������֤���롣
    !
    call neighbors_displacement( tri%nei2, target, deltax, deltay, deltaz, tri%nmodified, tri%newnei2, tri%newnei2idx )
    !
    !target���դ�ʬ�Ҥ�����ɽ����ʤ�����
    !(�ºݤˤ�̵�̤ʷ׻���¿�������������nei2�ε�Υ���ͤ�Ŭ�ڤ����٤�С��׻��̤��Ǿ����Ǥ��롣)
    !
    do i=1, tri%nmodified
       call neighbor1_compress2( tri%newnei(i), tri%newnei2(i), 4, ORDER_BY_DISTANCE )
    enddo
    !
    !tri%nmodified���ͤϡ�target��ư�����Ȥˤ�äơ��Ѳ�������˵�ꥹ�Ȥ����ǿ������㤤�Ǥʤ���С�target�����ܿ�+1(target���Ȥ�ޤޤ�뤫��)���������ʤ�Ȼפ���
    !
#ifdef DEBUG
    if ( tri%nei2(target)%nnear+1 /= tri%nmodified ) then
       write(STDERR,*) tri%nei2(target)%nnear+1, tri%nmodified, "should be the same."
    endif
#endif
    call new( tri%set, MAXMOL )
    !
    !��ư���ζ�˵��
    !
    call join( tri%set, target )
    do i=1,tri%nei(target)%nnear
       call join( tri%set, tri%nei(target)%near(i) )
    enddo
    do ii=1, tri%nei2(target)%nnear
       !
       !tgt�ζ�˵5A�ˤ���ʬ��i������
       !
       i = tri%nei2(target)%near(ii)
       !
       !i�ζ�˵��tgt�Ϥ���Ϥ���
       !
       k = 0
       do j=1, tri%nei(i)%nnear
          if ( tri%nei(i)%near(j) == target ) then
             call join( tri%set, i )
          endif
       enddo
    enddo
    !
    !��ư��ζ�˵��
    !
    do i=1,tri%newnei(1)%nnear
       call join( tri%set, tri%newnei(1)%near(i) )
    enddo
    do i=2, tri%nmodified
       k = 0
       do j=1, tri%newnei(i)%nnear
          if ( tri%newnei(i)%near(j) == target ) then
             call join( tri%set, tri%newnei2idx(i) )
          endif
       enddo
    enddo
    !write(STDERR,*) "SET", (set%set(i), i=1, set%nmark )
    !
    !����ʬ�Ҥ���ư����ȡ��ɤ�O-O-O�Ѥ��Ѳ����뤫
    !���ˤ�äƤϡ���ư���뤳�ȤǺǶ���4ʬ�Ҥ������ǽ�������롣
    !
    !histogram��ʣ������
    !
    call histogram_dup( tri%newhist, tri%hist )
#ifdef DEBUG2
    !call histogram_dup( tmphist, hist )
#endif
    !
    !�ҥ��ȥ���फ���ö��������
    !
    !verbose = loop == 1 .and. trial == 2
    !verbose = .true.
    verbose = .false.
    if ( verbose ) call histogram_show( tri%hist, STDERR )
    call partial_accumulate( tri%newhist, tri%set, tri%nei, tri%nei2, -1d0, verbose )
    !
    !�ޤ����濴ʬ�ҤΡ��Ƕ���4ʬ�Ҥ�ɽ�򹹿����롣tri%newnei2����������������Ƥ���Τǡ�������Ȥˤ���Ф褤��
    !
    !call neighbor1_compress( tri%newnei(1), tri%newnei2(1), 4 )
    if ( verbose ) then
       write(STDERR,*) target, ">", (tri%newnei(1)%near(i), i=1, tri%newnei(1)%nnear)
    endif
    !
    !�ҥ��ȥ����˿����˲û����롣
    !
    call partial_accumulate_indexed( tri%newhist, tri%set, tri%newnei, tri%newnei2, +1d0, tri%nmodified, tri%newnei2idx, verbose )
    
#ifdef DEBUG
    if ( verbose ) then
       do i=1, nmol
          do j=1, tri%nei(i)%nnear
             do k=j+1, tri%nei(i)%nnear
                !
                ! O-O-O�Ѥ�;����׻����롣
                !
                costh = tri%nei(i)%dx(j)*tri%nei(i)%dx(k) + tri%nei(i)%dy(j)*tri%nei(i)%dy(k) + tri%nei(i)%dz(j)*tri%nei(i)%dz(k)
                costh = costh / sqrt ( tri%nei(i)%dd(j) * tri%nei(i)%dd(k) )
                !if ( 0.55d0 < costh .and. costh < 0.65d0 ) then
                !   write(STDERR,*) "Large:", nei(i)%near(j), i, nei(i)%near(k), nint( ( costh - hist%minvalue ) / hist%binwidth ) + 1
                !endif
                !if ( i .eq. target .or. nei(i)%near(j) .eq. target .or. nei(i)%near(k) .eq.target ) then
                !   write(STDERR,*) "Before:", nei(i)%near(j), i, nei(i)%near(k), nint( ( costh - hist%minvalue ) / hist%binwidth ) + 1
                !endif
             enddo
          enddo
       enddo
    endif
#endif
    
    !
    !�����ޤ��ɤ���¤�ǤϥХ��餷����ΤϤʤ���
    !
    op = kullback_div( tri%newhist, refhist )
    if ( tri%logarithmic ) op = log(op)
  end subroutine sTriangle_Differentiate
  
  subroutine sTriangle_Integrate( tri )
    type( sTriangle ) :: tri
    
    integer :: i,j
    !
    !����ɽ(nei2)�򹹿����롣
    !
    do i=1, tri%nmodified
       j = tri%newnei2idx(i)
       tri%nei2(j) = tri%newnei2(i)
    enddo
    !
    !����ɽ(nei)�򹹿����롣
    !
    do i=1, tri%nmodified
       j = tri%newnei2idx(i)
       tri%nei(j) = tri%newnei(i)
    enddo
    !
    !�ҥ��ȥ����򹹿����롣
    !tri%newhist�Υݥ��󥿤�hist�˰ܤ�����hist�����Ƥ��Ѵ��������Τ�������������ˡ���褯�狼��ʤ��Τǡ��ޤ�꤯�ɤ������򤹤롣
    !
    call done( tri%hist )
    call histogram_dup( tri%hist, tri%newhist )
    call done( tri%newhist )
    
#ifdef DEBUG2
    !
    !1trial����ѻ��ҥ��ȥ����򡢺Ʒ׻����ư��פ��뤫�ɤ������ǧ
    !
    if ( .true. ) then
       call new( tri%tmphist, 21, 0.1d0, -1d0 )
       do i=1, sys%mol(1)%nmol
          do j=1, tri%nei(i)%nnear
             do k=j+1, tri%nei(i)%nnear
                !
                ! O-O-O�Ѥ�;����׻����롣
                !
                costh = tri%nei(i)%dx(j)*tri%nei(i)%dx(k) + tri%nei(i)%dy(j)*tri%nei(i)%dy(k) + tri%nei(i)%dz(j)*tri%nei(i)%dz(k)
                costh = costh / sqrt ( tri%nei(i)%dd(j) * tri%nei(i)%dd(k) )
                call histogram_accumulate( tri%tmphist, costh, 1d0 )
             enddo
          enddo
       enddo
       if ( associated( tri%tmphist ) ) then
          diff  = kullback_div( tri%hist, tri%tmphist )
          if ( diff /= 0d0 ) then
             call histogram_show( tri%hist, STDERR )
             call histogram_show( tri%tmphist, STDERR )
             write(STDERR,*) diff, "DIV", loop, trial
             call die( 0, "TMP Histogram differs." )
          endif
       endif
       call done( tri%tmphist )
    endif
#endif /*DEBUG2*/
  end subroutine sTriangle_Integrate

  subroutine sTriangle_Reject( tri )
    type( sTriangle ) :: tri
    call done( tri%newhist )
  end subroutine sTriangle_Reject

  subroutine sTriangle_EndAStep( tri )
    type( sTriangle ) :: tri
#ifdef DEBUG
    call histogram_dup( tri%lasthist, tri%hist )
#endif
    call done( tri%hist )
  end subroutine sTriangle_EndAStep

end module triangle
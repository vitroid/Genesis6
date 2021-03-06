! -*- f90 -*-
! apply various external fields.

module external_module
  use error_module
  use vector_module
  implicit none
  !
  !$BJ,;R$N=E?4$rFCDj$N:BI8$K%P%M$G$D$J$0!#(B
  !$B0l1~!"%P%MDj?t$rA4It<+M3$KA*$Y$k$h$&$K$7$F$*$3$&$+!#(B
  !
  !$B%P%M$N%(%M%k%.!<$O(Ba * x**b
  !x$B$NC10L$O(BAngstrom, $B%(%M%k%.!<$NC10L$OFbItC10L(B(dJ/mol)
  !$B5wN%(B100A$B$"$?$j$G(B10kJ/mol$B$0$i$$$NJI$r:n$k$J$i!"(B
  !a * 100**b == 10kJ/mol = 1000dJ/mol$B$+$i(Ba$B$r7h$a$l$P$h$$!#(B
  !
  integer, parameter :: TIE_NONE=0, TIE_ACTIVE=1, BIND_NONE=0, BIND_ACTIVE=1, MAXBIND = 100, JOINTHB_ACTIVE=1
  type sTie
     integer :: mode
     integer :: n,b
     real(kind=8),  dimension(:), pointer :: a
     type(vector3), dimension(:), pointer :: r
  end type sTie

  type sBind
     integer :: mode
     integer :: b, npair, balance
     logical :: press, stretch
     integer :: pair_i(MAXBIND), pair_j(MAXBIND)
     real(kind=8) :: a
  end type sBind
#ifdef SOLVATIONTEST
  !
  !Bind$B$H;w$F$$$k$,!"$3$A$i$O?e$N?eAG;@AG$rL@<(;XDj$7$F7k9g$rGSB>E*$K7A@.$5$;$k!#HFMQ@-$O$J$$!#(B
  !
  type sJointHB
     integer :: mode
     !
     !$B%/%i%9%?$KB0$9$k?eJ,;R(B
     !
     integer,pointer :: members(:)
     !
     !$B%/%i%9%?Fb?eAG7k9g!#$3$3$K4^$^$l$k7k9g0J30$OGS=|$5$l$k!#(B
     !
     integer,pointer :: bonds(:,:)
     !
     !$B7k9gAj<j!#(B0$B$J$iAj<j$O;XDj$5$l$F$$$J$$!#(B
     !
     integer,pointer :: partner(:,:)
     !
     !$B7k9g$NDj?t$J$I(B
     !
     integer :: b, nmol, balance
     real(kind=8) :: a
  end type sJointHB
#endif

  interface new
     module procedure tie_constructor
  end interface

  interface save_binary
     module procedure tie_writebinarytie4
  end interface


contains

  subroutine Tie_Constructor(t)
    implicit none
    type(sTie),intent(inout) :: t
    t%mode = TIE_NONE
    t%n    = 0
    t%b    = 2
  end subroutine Tie_Constructor
  
  subroutine Tie_Initialize(t,n,x,y,z,a,b)
    implicit none
    type(sTie),intent(inout) :: t
    integer,intent(in) :: n,b
    real(kind=8),dimension(*),intent(in) :: x,y,z,a
    integer :: i
    t%n = n
    if ( associated( t%r ) ) deallocate( t%r )
    if ( associated( t%a ) ) deallocate( t%a )
    if ( n > 0 ) then
       t%mode = TIE_ACTIVE
       allocate(t%r(n))
       allocate(t%a(n))
       t%r(1:n)%vec(1) = x(1:n)
       t%r(1:n)%vec(2) = y(1:n)
       t%r(1:n)%vec(3) = z(1:n)
       t%a(:)          = a(1:n)
       t%b             = b
    else
       t%mode = TIE_NONE
    endif
  end subroutine Tie_Initialize
  
  function Tie_ReadTIE3(t,file)
    implicit none
    type(sTie),intent(inout) :: t
    integer,intent(in) :: file
    integer :: i,j
    logical :: Tie_ReadTIE3
    Tie_ReadTIE3=.false.
    t%mode=TIE_NONE
    read(file,*) t%n
    t%b = 2
    if(t%n.lt.0)return
    t%mode=TIE_ACTIVE
    allocate(t%r(t%n))
    do i=1,t%n
       read(file,*) (t%r(i)%vec(j),j=1,3)
    enddo
    Tie_ReadTIE3=.true.
  end function Tie_ReadTIE3
  
  function Tie_ReadTIE4(t,file)
    implicit none
    type(sTie),intent(inout) :: t
    integer,intent(in) :: file
    integer :: i,j
    logical :: Tie_ReadTIE4
    Tie_ReadTIE4=.false.
    t%mode=TIE_NONE
    read(file,*) t%n, t%b
    if(t%n.lt.0)return
    t%mode=TIE_ACTIVE
    allocate(t%r(t%n))
    allocate(t%a(t%n))
    do i=1,t%n
       read(file,*) (t%r(i)%vec(j),j=1,3),t%a(i)
    enddo
    Tie_ReadTIE4=.true.
  end function Tie_ReadTIE4
  
  function Tie_ReadBinaryTIE3(t,file)
    implicit none
    type(sTie),intent(inout) :: t
    integer,intent(in) :: file
    integer :: i,j
    logical :: Tie_ReadBinaryTIE3
    Tie_ReadBinaryTIE3=.false.
    t%mode=TIE_NONE
    read(file) t%n
    t%b = 2
    if(t%n.lt.0)return
    t%mode=TIE_ACTIVE
    allocate(t%r(t%n))
    do i=1,t%n
       read(file) (t%r(i)%vec(j),j=1,3)
    enddo
    Tie_ReadBinaryTIE3=.true.
  end function Tie_ReadBinaryTIE3
  
  function Tie_ReadBinaryTIE4(t,file)
    implicit none
    type(sTie),intent(inout) :: t
    integer,intent(in) :: file
    integer :: i,j
    logical :: Tie_ReadBinaryTIE4
    Tie_ReadBinaryTIE4=.false.
    t%mode=TIE_NONE
    read(file) t%n,t%b
    if(t%n.le.0)return
    t%mode=TIE_ACTIVE
    allocate(t%r(t%n))
    allocate(t%a(t%n))
    do i=1,t%n
       read(file) (t%r(i)%vec(j),j=1,3),t%a(i)
    enddo
    Tie_ReadBinaryTIE4=.true.
  end function Tie_ReadBinaryTIE4
  
  subroutine Interaction_Force_Tie(t,b,n,r,f,ep)
    use box_module
    use vector_module
    type(sTie),intent(in) :: t
    type(sBox),intent(in) :: b
    integer,intent(in) :: n
    type(vector3),dimension(*),intent(in)    :: r
    type(vector3),dimension(*),intent(inout) :: f
    real(kind=8),intent(out) :: ep
    real(kind=8) :: dd
    type(vector3) :: delta,cell
    integer :: i
    if(n.ne.t%n) call die( error_different_num_of_molecules, "External 1" )
    ep=0d0
    do i=1,n
       delta%vec(:) = r(i)%vec(:) - t%r(i)%vec(:)
       call Box_Renormalize(b,delta,cell)
       dd=inner_product( delta, delta )
       ep=ep + t%a(i) * dd**(t%b/2)
       f(i)%vec(:)=f(i)%vec(:) - t%a*t%b*delta%vec(:)*dd**(t%b/2-1)
    enddo
  end subroutine Interaction_Force_Tie
  
  !
  !$BJ,;RBP$r7k$S$D$1$k!#(B
  !
  subroutine Interaction_Force_Bind(bind,b,r,f,ep)
    use box_module
    use vector_module
    type(sBind),intent(in) :: bind
    type(sBox),intent(in) :: b
    type(vector3),dimension(*),intent(in)    :: r
    type(vector3),dimension(*),intent(inout) :: f
    real(kind=8),intent(out) :: ep
    real(kind=8) :: dd, d
    type(vector3) :: delta,cell
    integer :: i,j,k
    ep=0d0
    do k=1, bind%npair
       i = bind%pair_i(k)
       j = bind%pair_j(k)
       delta%vec(:) = r(i)%vec(:) - r(j)%vec(:)
       call Box_Renormalize(b,delta,cell)
       dd = inner_product( delta, delta )
       if ( ( dd < bind%balance**2 .and. bind%stretch ) .or. ( bind%balance**2 < dd .and. bind%press ) ) then
          d  = sqrt(dd)
          ep = ep + bind%a * (d - bind%balance)**bind%b
          f(i)%vec(:) = f(i)%vec(:) - bind%a * bind%b * delta%vec(:) * ( d - bind%balance )**( bind%b - 1 ) / d
          f(j)%vec(:) = f(j)%vec(:) + bind%a * bind%b * delta%vec(:) * ( d - bind%balance )**( bind%b - 1 ) / d
       endif
    enddo
  end subroutine Interaction_Force_Bind
  
!  subroutine Tie_WriteBinaryTIEP(t,f)
!    implicit none
!    integer,intent(in) :: f
!    type(sTie),intent(in) :: t
!    write(f)'@TIEP'
!    write(f) t%a
!  end subroutine Tie_WriteBinaryTIEP
  
  subroutine Tie_WriteBinaryTIE3(t,f)
    implicit none
    integer,intent(in) :: f
    type(sTie),intent(in) :: t
    integer :: i,j
    write(f) '@TIE3'
    write(f) t%n
    do i=1,t%n
       write(f) (t%r(i)%vec(j),j=1,3)
    enddo
  end subroutine Tie_WriteBinaryTIE3

  subroutine Tie_WriteBinaryTIE4(t,f)
    implicit none
    integer,intent(in) :: f
    type(sTie),intent(in) :: t
    integer :: i,j
    write(f) '@TIE4'
    if ( t%mode == TIE_NONE ) then
       write(f) 0,2
    else
       write(f) t%n,t%b
       do i=1,t%n
          write(f) (t%r(i)%vec(j),j=1,3),t%a(i)
       enddo
    endif
  end subroutine Tie_WriteBinaryTIE4
end module external_module


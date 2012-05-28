! -*- f90 -*-
module book_module
  use vector_module
  use error_module
  implicit none
  integer, parameter :: BOOK_NONE=0, BOOK_FIX=1
  integer, parameter :: BOOK_AUTOINTERVAL=2, BOOK_AUTOMERGIN=3
  type sBook
     !middle is usually set as rc+margin/2
     real(kind=8) :: outer,margin
     !Previous position to monitor displacement during book-rebuild
     type(vector3), dimension(:), pointer :: prev
     !outer booking radius and middle of booking margin radius
     !outer is usually set as rc+margin
     !real(kind=8) :: middle
     !update interval
     !
     !帳簿を作りなおす間隔
     !
     integer :: interval,countdown
     !
     !帳簿法のモード
     !
     integer :: mode
  end type sBook

  interface new
     module procedure Book_Constructor
  end interface
  interface load
     module procedure book_loader
  end interface
  interface save
     module procedure book_save
  end interface
  interface load_binary
     module procedure book_binaryloader
  end interface
  interface save_binary
     module procedure book_savebinary
  end interface
  interface active
     module procedure book_isactive
  end interface

contains

  function Book_ReadBOOK(b,file)
    integer,intent(IN) ::file
    type(sBook),intent(OUT) :: b
    logical :: Book_ReadBOOK
    read(file,*) b%interval,b%margin
#ifdef VERBOSE
    write(STDERR,*) "@BOOK",b%interval,b%margin
#endif
    if (b%interval.gt.0) then
       Book_ReadBOOK=.true.
       b%mode = BOOK_FIX
    else
       Book_ReadBOOK=.false.
       b%mode = BOOK_NONE
    endif
    return
  end function Book_ReadBOOK
  
  function Book_ReadBinaryBOOK(b,file)
    integer,intent(IN) ::file
    type(sBook),intent(OUT) :: b
    logical :: Book_ReadBinaryBOOK
    read(file) b%interval,b%margin
#ifdef VERBOSE
    write(STDERR,*) "@BOOK",b%interval,b%margin
#endif
    if (b%interval.gt.0) then
       Book_ReadBinaryBOOK=.true.
       b%mode = BOOK_FIX
    else
       Book_ReadBinaryBOOK=.false.
       b%mode = BOOK_NONE
    endif
    return
  end function Book_ReadBinaryBOOK
  
  subroutine Book_WriteBinaryBOOK(b,file)
    integer,intent(IN) ::file
    type(sBook),intent(IN) :: b
    write(file) "@BOOK"
    write(file) b%interval,b%margin
    return
  end subroutine Book_WriteBinaryBOOK
  
  subroutine Book_WriteBOOK(b,file)
    integer,intent(IN) ::file
    type(sBook),intent(IN) :: b
    write(file,1)
    write(file,*) b%interval,b%margin
  1 format( "@BOOK" )
    return
  end subroutine Book_WriteBOOK
  
  subroutine Book_Constructor(bk)
    type(sBook),intent(OUT) :: bk
    bk%mode=BOOK_NONE
    bk%countdown=0
    bk%interval=1
    bk%margin=0d0
    bk%outer=-1
  end subroutine Book_Constructor
  
  function Book_IsActive(bk)
    logical                :: Book_IsActive
    Type(sBook),intent(in) :: bk
    Book_IsActive = (bk%mode.ne.BOOK_NONE)
  end function Book_IsActive

  subroutine Book_Initialize(bk, mode, margin, interval)
    type(sBook),   intent(inout) :: bk
    integer,       intent(in)    :: mode, interval
    real(kind=8),  intent(in)    :: margin
    bk%mode=mode
    bk%countdown=0
    bk%interval=interval
    bk%margin=margin
    bk%outer=-1
  end subroutine Book_Initialize

  subroutine Book_ResetInterval(bk,outer,n,com)
    integer,intent(IN) ::n
    real(kind=8),intent(in) :: outer
    type(vector3), dimension(*), intent(IN) :: com
    type(sBook),intent(INOUT) :: bk
    integer :: i
    bk%outer=outer
    bk%countdown=0
    if(associated(bk%prev))then
       deallocate(bk%prev)
    endif
    allocate(bk%prev(n))
    bk%prev(1:n)=com(1:n)
    return
  end subroutine Book_ResetInterval
  
  subroutine Book_AssumeInterval(bk,n,com,recomm)
    real(kind=8),intent(IN) :: recomm
    integer,intent(IN) ::n
    type(vector3),dimension(*),intent(IN) :: com
    type(sBook),intent(INOUT) :: bk
    integer i
    if(recomm.lt.bk%margin)then
       bk%interval = bk%interval +2!* 1.05
    else
       bk%interval = 0.9d0*bk%interval + 0.5d0
    endif
#ifdef VERBOSE
    write(STDERR,*) bk%interval,' Book remaking interval'
#endif
     bk%prev(1:n)=com(1:n)
  end subroutine Book_AssumeInterval
  
  function Book_BestMargin(bk,n,com)
    type(sBook),intent(in) :: bk
    type(vector3), intent(in), dimension(*) :: com
    integer,intent(in) :: n
    real(kind=8) :: Book_BestMargin
    !local
    real(kind=8) :: dx,dy,dz,dd,rrmax
    integer :: i
    rrmax=0d0
    if(bk%outer < 0)then
       write(STDERR,*) "BOOKING INFORMATIONS ARE NOT SET."
       call die( 0, "Book 1" )
    endif
    do i=1,n
       dx = com(i)%vec(1) - bk%prev(i)%vec(1)
       dy = com(i)%vec(2) - bk%prev(i)%vec(2)
       dz = com(i)%vec(3) - bk%prev(i)%vec(3)
       dd = dx**2 + dy**2 + dz**2
       if(dd.gt.rrmax)then
          rrmax=dd
       endif
    enddo
  !     margin should be twice as far as relative maximum square displacement
    Book_BestMargin=dsqrt(rrmax*2d0)*2d0
#ifdef VERBOSE
    write(STDERR,*) "Recommended margin width=",dsqrt(rrmax*2d0)*2d0
#endif
  end function Book_BestMargin
  
  subroutine Book_Loader(bk,file,tag)
    type(sBook),intent(inout) :: bk
    integer,intent(IN) :: file
    character(len=5),intent(IN) :: tag
    logical :: result
    if(tag == "@BOOK")then
       result = Book_ReadBOOK(bk,file)
    endif
  end subroutine Book_Loader
  
  subroutine Book_BinaryLoader(bk,file,tag)
    type(sBook),intent(inout) :: bk
    integer,intent(IN) :: file
    character(len=5),intent(IN) :: tag
    logical :: result
    if(tag == "@BOOK")then
       result = Book_ReadBinaryBOOK(bk,file)
    endif
  end subroutine Book_BinaryLoader
  
  subroutine Book_Save(bk,file)
    type(sBook),intent(in) :: bk
    integer,intent(IN) :: file
    if(bk%mode .ne. BOOK_NONE)then
       call Book_WriteBOOK(bk,file)
    endif
  end subroutine Book_Save
  
  subroutine Book_SaveBinary(bk,file)
    type(sBook),intent(in) :: bk
    integer,intent(IN) :: file
    if(bk%mode .ne. BOOK_NONE)then
       call Book_WriteBinaryBOOK(bk,file)
    endif
  end subroutine Book_SaveBinary
end module book_module
  

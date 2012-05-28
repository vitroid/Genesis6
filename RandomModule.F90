! -*- f90 -*-
!mtrndのwrapper。mtrndが大域的変数を使っているが、擬似的にオブジェクト化する。複数同時に利用してはいけない。
!
module random_module
  use common_module
  implicit none
  integer,parameter :: BUFSIZ = 1000
  type pRandom
     real(kind=8)      :: random_numbers(BUFSIZ)
     integer           :: next
  end type pRandom
  integer,save,private :: initialized = 0
  private BUFSIZ

  interface new
     module procedure random_initialize
  end interface

  interface done
     module procedure random_done
  end interface

contains
  !
  !Box-Muller
  !
  !generate a normal random number from two given white random numbers
  function Random_Normal(avg,var,ran1,ran2)
    use physconst_module
    real(kind=8), intent(in) :: ran1,ran2,avg,var
    real(kind=8) :: r,theta
    real(kind=8) :: Random_Normal
    r=dsqrt(-2d0*dlog(ran1))
    theta=2d0*PI*ran2
    Random_Normal=r*dsin(theta)*dsqrt(var)+avg
  end function Random_Normal
  !
  !Initialize random number with given seed
  !
  subroutine Random_Initialize( rand, seed )
    use error_module
    type(pRandom), pointer :: rand
    integer, optional :: seed
    if ( 0 < initialized ) then
       write(STDERR,*) initialized
       call die( 0, "Random must not be initialized multiply." )
    endif
    initialized = initialized + 1
    allocate( rand )
    if ( present( seed ) ) call mtrngi(seed)
    rand%next = BUFSIZ + 1
  end subroutine Random_Initialize


  subroutine random_done( rand )
    type(pRandom), pointer :: rand
    initialized = initialized - 1
    deallocate( rand )
  end subroutine random_done
  !
  !Initialize random number with internal state
  !
  subroutine Random_Restore( rand, lbuf )
    type(pRandom), pointer :: rand
    integer :: lbuf(624)
    call mtrnl( lbuf )
  end subroutine Random_Restore
  !
  !Evacuate internal state
  !
  subroutine Random_Evacuate( rand, lbuf )
    type(pRandom), pointer :: rand
    integer :: lbuf(624)
    call mtrns( lbuf )
  end subroutine Random_Evacuate
  !
  !Load last state from file
  !
  subroutine Random_LoadMTRN( rand, file )
    type(pRandom), pointer :: rand
    integer :: file

    integer :: lbuf(624), i, bufsize

    call Random_Initialize( rand )
    read(file,*) rand%next
    read(file,*) bufsize  !dummy
    read(file,*) ( rand%random_numbers(i), i=1, BUFSIZ )
    read(file,*) ( lbuf(i), i=1, 624 )
    call Random_Restore( rand, lbuf )
  end subroutine Random_LoadMTRN
  !
  !Save last state to file
  !
  subroutine Random_SaveMTRN( rand, file )
    type(pRandom), pointer :: rand
    integer :: file

    integer :: lbuf(624), i, bufsize

    call writetag( file, "@MTRN" )
    write(file,*) rand%next
    write(file,*) BUFSIZ  !dummy
    do i=1, BUFSIZ
       write(file,fmt0) rand%random_numbers(i)
    enddo
    call Random_Evacuate( rand, lbuf )
    do i=1, 624
       write(file,*) lbuf(i)
    enddo
  end subroutine Random_SaveMTRN
  !
  !fill buffer with random number
  !
  subroutine Random_Prepare( rand )
    type(pRandom), pointer :: rand
    call mtrndv( rand%random_numbers, BUFSIZ )
    rand%next = 1
  end subroutine Random_Prepare
  !
  !Get next random bumber from buffer
  !
  function Random_GetNext( rand )
    type(pRandom), pointer :: rand
    real(kind=8) :: x, Random_GetNext
    if( rand%next.gt.BUFSIZ ) call Random_Prepare( rand )
    x = rand%random_numbers( rand%next )
    rand%next = rand%next + 1
    Random_GetNext = 1d0 - x
  end function Random_GetNext
  !
  !Save internal state
  !
  subroutine Random_Save( rand, file )
    type(pRandom), pointer :: rand
    integer,         intent(in) :: file
    integer :: lbuf(624)
    integer :: i
    write(file,"('@RNDV')")
    call mtrns(lbuf)
    do i=1,624
       write(file,*) lbuf(i)
    enddo
    write(file,"('@RNDB')")
    write(file,*) BUFSIZ,rand%next
    do i=1,BUFSIZ
       write(file,*) rand%random_numbers(i)
    enddo
  end subroutine Random_Save
  !
  !Recover internal state
  !
  subroutine Random_Loader(rand, file,tag)
    type(pRandom), pointer :: rand
    character(len=5),intent(in) :: tag
    integer,         intent(in) :: file
    integer :: lbuf(624)
    integer :: b,i
    if ( tag == '@RNDV' )then
       do i=1,624
          read(file,*) lbuf(i)
       enddo
       call mtrnl(lbuf)
    else if ( tag == '@RNDB' ) then
       read(file,*) b,rand%next
       if ( b /= BUFSIZ ) stop
       do i=1,b
          read(file,*) rand%random_numbers(i)
       enddo
    endif
  end subroutine Random_Loader

  subroutine shuffle(rand,size,list)
    type(pRandom), pointer :: rand
    integer :: size
    integer :: list(size)
    integer :: mark(size)
    integer :: i,j,pick,last,firstslot
    mark(:) = 0
    last = 0
    firstslot = 0
    do i=1,size
       findnew: do
          j = Random_GetNext( rand ) * size + 1
          if ( mark(j) == 0 ) exit findnew
       enddo findnew
       if ( firstslot == 0 ) then
          firstslot = j
       endif
       mark( j ) = 1
       pick = list(j)
       list(j) = last
       last = pick
    enddo
    list( firstslot ) = last
  end subroutine shuffle
end module random_module

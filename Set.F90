! -*- f90 -*-
!
!要求されるサイズが非常に大きい場合はhashの方が便利。
!
module set_module
  use common_module
  implicit none
  type sSet
     logical :: private_mark(MAXMOL)
     integer :: size, maxsize
     integer :: set(MAXMOL)
  end type sSet

  interface exist
     module procedure set_exist
  end interface

  interface new
     module procedure set_new
     module procedure set_initial
  end interface

contains

  subroutine set_new( set, maxsize )
    type(sSet) :: set
    integer    :: maxsize
    set%size = 0
    set%private_mark(:) = .false.
    set%maxsize = maxsize
  end subroutine set_new

  subroutine set_initial( set, maxsize, n, list, dup )
    type(sSet) :: set
    integer, intent(IN) :: n
    integer, intent(IN) :: list(*)
    logical, intent(IN) :: dup  !listの要素に重複があるか否か。
    integer    :: maxsize

    integer :: i
    call new( set, maxsize )
    if( dup ) then
       do i=1,n
          call join( set, list(i) )
       enddo
    else
       set%size = n
       do i=1,n
          set%private_mark(list(i)) = .true.
          set%set(i) = list(i)
       enddo
    endif
  end subroutine set_initial

  subroutine join( set, elem )
    use error_module
    type(sSet) :: set
    integer, intent(IN) :: elem
    if ( .not. set%private_mark(elem) ) then
       set%size = set%size + 1
       set%set(set%size) = elem
       set%private_mark(elem) = .true.
    end if
  end subroutine join

  subroutine drop( set, elem )
    use error_module
    type(sSet) :: set
    integer, intent(IN) :: elem

    integer :: i
    if ( set%private_mark(elem) ) then
       loop:do i=1, set%size
          if ( set%set(i) == elem ) then
             set%set(i) = set%set(set%size)
             set%size = set%size -1
             set%private_mark(elem) = .false.
             exit loop
          endif
       end do loop
    end if
  end subroutine drop

  subroutine intersection( set1, set2, result )
    use error_module
    type(sSet) :: set1, set2, result

    integer :: i,ii
    
    if ( set1%maxsize /= set2%maxsize ) then
       call warn( 0, "Set size differs.(intersection)" )
       write(STDERR,*) set1%maxsize,set2%maxsize
    endif
    call new( result, set1%maxsize )
    do ii=1, set1%size
       i = set1%set(ii)
       if ( exist( set2, i ) ) then
          call join( result, i )
       endif
    enddo
  end subroutine intersection
      
    
  subroutine difference( set1, set2, result )
    use error_module
    type(sSet) :: set1, set2, result

    integer :: i,ii
    
    if ( set1%maxsize /= set2%maxsize ) then
       call warn( 0, "Set size differs.(difference)" )
       write(STDERR,*) set1%maxsize,set2%maxsize
    endif
    call new( result, set1%maxsize )
    do ii=1, set1%size
       i = set1%set(ii)
       if ( .not. exist( set2, i ) ) then
          call join( result, i )
       endif
    enddo
  end subroutine difference
      


  function set_exist( set, element )
    integer, intent(IN) :: element
    type(sSet)          :: set
    logical :: set_exist
    set_exist = set%private_mark( element )
  end function set_exist

  subroutine sSet_write( set, file )
    type(sSet)          :: set
    integer             :: file

    integer :: i
    write(file,*) ( set%set(i), i=1, set%size )
  end subroutine sSet_write

  subroutine sSet_compare( s1,s2, file )
    type(sSet), intent(IN) :: s1,s2
    integer, intent(IN)    :: file

    type(sSet) :: result
    call new( result, s2%maxsize )
    call difference( s1,s2,result )
    if ( 0 < result%size ) then
       write(file,*) "Only in s1:"
       call sSet_write( result, file )
    endif
    call difference( s2,s1,result )
    if ( 0 < result%size ) then
       write(file,*) "Only in s2:"
       call sSet_write( result, file )
    endif
  end subroutine sSet_compare

end module set_module

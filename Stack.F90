module stack_module
  implicit none
  type real8stack
     real(kind=8),pointer :: value(:)
     integer          :: depth,size
  end type real8stack

  type indexedreal8stack
     integer, pointer :: index(:)
     real(kind=8),pointer :: value(:)
     integer          :: depth,size
  end type indexedreal8stack

  interface push
     module procedure real8stack_push
     module procedure indexedreal8stack_push
  end interface

  interface done
     module procedure real8stack_done
     module procedure indexedreal8stack_done
  end interface

  interface new
     module procedure real8stack_new
     module procedure indexedreal8stack_new
  end interface
  
  interface empty
     module procedure real8stack_empty
     module procedure indexedreal8stack_empty
  end interface
  
contains

  subroutine real8stack_new( stack, size )
    type(real8stack), intent(INOUT) :: stack
    integer :: size
    allocate( stack%value( size ) )
    stack%size = size
    call empty( stack )
  end subroutine real8stack_new

  subroutine indexedreal8stack_new( stack, size )
    type(indexedreal8stack), intent(INOUT) :: stack
    integer :: size
    allocate( stack%value( size ) )
    allocate( stack%index( size ) )
    stack%size = size
    call empty( stack )
  end subroutine indexedreal8stack_new

  subroutine real8stack_done( stack )
    type(real8stack), intent(INOUT) :: stack
    deallocate( stack%value )
    stack%size = 0
    call empty( stack )
  end subroutine real8stack_done

  subroutine indexedreal8stack_done( stack )
    type(indexedreal8stack), intent(INOUT) :: stack
    deallocate( stack%value )
    deallocate( stack%index )
    stack%size = 0
    call empty( stack )
  end subroutine indexedreal8stack_done

  subroutine real8stack_empty( stack )
    type(real8stack), intent(INOUT) :: stack
    stack%depth = 0
  end subroutine real8stack_empty

  subroutine indexedreal8stack_empty( stack )
    type(indexedreal8stack), intent(INOUT) :: stack
    stack%depth = 0
  end subroutine indexedreal8stack_empty

  subroutine real8stack_push( stack, value )
    type(real8stack), intent(INOUT) :: stack
    real(kind=8),     intent(IN)    :: value
    stack%depth = stack%depth + 1
    stack%value( stack%depth ) = value
  end subroutine real8stack_push
  
  subroutine indexedreal8stack_push( stack, index, value )
    type(indexedreal8stack), intent(INOUT) :: stack
    real(kind=8),     intent(IN)    :: value
    integer,          intent(IN)    :: index
    stack%depth = stack%depth + 1
    stack%value( stack%depth ) = value
    stack%index( stack%depth ) = index
  end subroutine indexedreal8stack_push
  
end module stack_module

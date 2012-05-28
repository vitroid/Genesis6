module error_module
  integer, parameter :: error_none=0
  integer, parameter :: error_id_not_found=1
  integer, parameter :: error_different_num_of_molecules=2
  integer, parameter :: error_too_many_components=3
  integer, parameter :: error_anonymous_molecule=4
  integer, parameter :: error_max=4


contains

  subroutine die( err, additional )
    integer, intent(IN) :: err
    character(LEN=*), intent(IN), optional :: additional
    if ( present( additional ) ) then
       call warn( err, additional )
    else
       call warn( err )
    endif
    stop
  end subroutine die

  subroutine warn( err, additional )
    integer, intent(IN) :: err
    character(LEN=*), intent(IN), optional :: additional
    character(len=50) :: errorMessages(100)
    if ( present( additional ) ) then
       write( STDERR,* ) additional
    endif
       
    errorMessages(error_id_not_found) = "Molecule ID not found."
    errorMessages(error_different_num_of_molecules) = "Different number of molecules for morph."
    errorMessages(error_too_many_components) = "Number of componets exceeded the limit."
    errorMessages(error_anonymous_molecule) = "Molecule name is not specified."

    if ( 0 < err .and. err .le. error_max ) then
       write( STDERR, * ) errorMessages( err )
    else
       write( STDERR, * ) "Error code ", err
    endif
  end subroutine warn


end module error_module

! -*- f90 -*-
module Umbrella
  implicit none
  integer, parameter :: umbrella_size = 8*2
  type umbrella_parabolic
     real(kind=8) :: center, coeffi
  end type umbrella_parabolic
  
contains
  
  function umbrella_potential( um, op )
    type( umbrella_parabolic ), intent(IN) :: um
    real(kind=8), intent(IN) :: op
    real(kind=8) :: umbrella_potential
    umbrella_potential = um%coeffi * ( op - um%center )**2
  end function umbrella_potential

  subroutine umbrella_ReadParam( um, file )
    type( umbrella_parabolic ), intent(INOUT) :: um
    integer, intent(IN) :: file
    read(file,*) um%center, um%coeffi
    write(6,*) um%center, um%coeffi
  end subroutine umbrella_ReadParam

end module Umbrella

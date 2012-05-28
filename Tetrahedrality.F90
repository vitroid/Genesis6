module tetrahedrality_module
  implicit none
contains
  function tetrahedrality( x, y, z, dd )
    real(kind=8) :: x(4), y(4), z(4), dd(4)
    real(kind=8) :: tetrahedrality
    real(kind=8) :: prod, sum
    real(kind=8) :: r(4)
    integer      :: i,j
    sum = 0d0
    do i=1,4
       r(i) = sqrt( dd(i) )
    enddo
    do i=1,3
       do j=i+1,4
          prod = x(i)*x(j) + y(i)*y(j) + z(i)*z(j)
          prod = prod / ( r(i)*r(j) )
          sum = sum + (prod + 1d0/3d0) ** 2
       enddo
    enddo
    tetrahedrality = 1d0 - 3d0/8d0 * sum
    !write(6,*) "TET", sum, (r(i),i=1,4)
  end function tetrahedrality
end module tetrahedrality_module

! -*- f90 -*-
module vector_module
  implicit none
  type,public :: vector3
     sequence
     real(kind=8) :: vec(3)
  end type vector3
  type vector4
     sequence
     real(kind=8) :: vec(4)
  end type vector4

  interface scale1
     module procedure vector_scale1
  end interface

  interface scalev
     module procedure vector_scalev
  end interface

  interface scale
     module procedure vector_scale1
     module procedure vector_scalev
  end interface

  interface volume
     module procedure vector_volume
  end interface

contains

  real(kind=8) function vector_volume(this)
    type(vector3),intent(in) :: this
    vector_volume = this%vec(1) * this%vec(2) * this%vec(3)
  end function vector_volume

  subroutine invert(this)
    type(vector3),intent(inout) :: this
    this%vec(:) = 1d0 / this%vec(:)
  end subroutine invert

  subroutine vector_scale1(this,ratio)
    type(vector3),intent(inout) :: this
    real(kind=8),intent(in)     :: ratio
    this%vec(:) = this%vec(:) * ratio
  end subroutine vector_scale1

  subroutine vector_scalev(this,ratio)
    type(vector3),intent(inout) :: this
    type(vector3),intent(in)     :: ratio
    this%vec(:) = this%vec(:) * ratio%vec(:)
  end subroutine vector_scalev

  function inner_product(this,ratio)
    real(kind=8)             :: inner_product
    type(vector3),intent(in) :: this
    type(vector3),intent(in) :: ratio
    real(kind=8)             :: sum
    integer                  :: i
    sum = 0
    do i=1,3
       sum = sum + this%vec(i) * ratio%vec(i)
    enddo
    inner_product = sum
  end function inner_product

  subroutine outer_product(a1,a2,result)
    type(vector3),intent(in) :: a1,a2
    type(vector3),intent(out) :: result

    result%vec(1) = a1%vec(2)*a2%vec(3) - a1%vec(3)*a2%vec(2)
    result%vec(2) = a1%vec(3)*a2%vec(1) - a1%vec(1)*a2%vec(3)
    result%vec(3) = a1%vec(1)*a2%vec(2) - a1%vec(2)*a2%vec(1)
  end subroutine outer_product

  subroutine vector_normalize( x, z )
    type(vector3),intent(in) :: x
    type(vector3),intent(out) :: z
    real(kind=8) :: y
    y = 1d0 / dsqrt(x%vec(1)**2 + x%vec(2)**2 + x%vec(3)**2)
    z = x
    call vector_scale1( z, y )
  end subroutine vector_normalize
    
  !
  !$BM?$($i$l$?(B2$B$D$N%Y%/%H%k$H$N3QEY$,$=$l$>$l(Btheta1,theta2$B$J$k$h$&$JBh(B3
  !$B$N%Y%/%H%k$r5a$a$k!#(B
  !
  subroutine thirdvector(first,second,theta1,theta2, third)
    use matrix_module
    implicit none
    type(vector3), intent(in)  :: first, second
    real(kind=8),intent(in)    :: theta1,theta2
    type(vector3), intent(out) :: third
    type(vector3)              :: first1, second1, third1, third2
    real(kind=8) :: r11,r12,r13,r21,r22,r23,r31,r32,r33,c1,c2,c3
    integer :: i
    call vector_normalize(first, first1)
    call vector_normalize(second, second1)
    !$B2r@OE*$K5a$a$i$l$=$&$J5$$,$9$k$,!"?tCM7W;;$G5a$a$F$_$k!#(B
    !x1x3 + y1y3 + z1z3 = cos theta1
    !x2x3 + y2y3 + z2z3 = cos theta2
    !x3x3 + y3y3 + z3z3 = 1d0
    !$B$r$_$?$9(B(x3,y3,z3)$B$r5a$a$l$P$h$$!#>e$N<0$r9TNs$H%Y%/%H%k$G(B
    !R v3 = c
    !$B$H=q$1$P!"(Bv3 = R^-1 c$B$H$J$k!#(Bv3$B$NCM$r(BR$B$K:F5"$5$;$l$P<}B+$9$k$@$m$&!#(B
    !$B=i4|CM$O!"(Bfirst$B$H(Bsecond$B$KD>8r$9$k%Y%/%H%k$H$9$k(B
    call outer_product( first1, second1, third1 )
    c1 = dcos(theta1)
    c2 = dcos(theta2)
    c3 = 1d0
    r11 = first1%vec(1)
    r12 = first1%vec(2)
    r13 = first1%vec(3)
    r21 = second1%vec(1)
    r22 = second1%vec(2)
    r23 = second1%vec(3)
    do i=1,50
       r31 = third1%vec(1)
       r32 = third1%vec(2)
       r33 = third1%vec(3)
       call inv33(r11,r12,r13,r21,r22,r23,r31,r32,r33)
       third%vec(1) = r11*c1 + r12*c2 + r13*c3
       third%vec(2) = r21*c1 + r22*c2 + r23*c3
       third%vec(3) = r31*c1 + r32*c2 + r33*c3
       call vector_normalize( third, third2 )
       third1%vec(:) = third1%vec(:)*0.9 + third2%vec(:)*0.1
       write(6,*) third1%vec(1),third1%vec(2),third1%vec(3)
    enddo
  end subroutine thirdvector
  
  subroutine thirdvector1(first,second,theta1,theta2, third)
    use matrix_module
    implicit none
    type(vector3), intent(in)  :: first, second
    real(kind=8),intent(in)    :: theta1,theta2
    type(vector3), intent(out) :: third
    type(vector3)              :: first1, second1, third1
    integer :: i,count
    real(kind=8)               :: a,b,c,d,e,f, costh1,costh2, dx3,dy3,dz3, len
    call vector_normalize(first, first1)
    call vector_normalize(second, second1)
    !
    !thirdvector0$B$N:F5"K!$G$O<}B+$7$J$+$C$?$N$GJL$NJ}K!$r9M$($k!#(B
    !3$BHVL\$N<0$r(B1,2$B$KBeF~$7!"O"N)(BNewton$BK!$G2r$/!#(B
    !f(y3,z3) = x1x3 + y1y3 + z1z3 - costh1 = 0
    !g(y3,z3) = x2x3 + y2y3 + z2z3 - costh2 = 0
    !$B$?$@$7(Bx3 = sqrt(1 - y2^2 - z2^2)
    !
    !df/dy3   = y1 - x1y3/x3
    !df/dz3   = z1 - x1z3/x3
    !dg/dy3   = y2 - x2y3/x3
    !dg/dz3   = z2 - x2z3/x3
    !
    !F( a, b ) + dF(a,b)/dX *$B&D(BX + dF(a,b)/dY *$B&D(BY = 0 
    !G( a, b ) + dG(a,b)/dX *$B&D(BX + dG(a,b)/dY *$B&D(BY = 0 
    !
    !$B=i4|CM$O!"(Bfirst$B$H(Bsecond$B$KD>8r$9$k%Y%/%H%k$H$9$k(B
    call outer_product( first1, second1, third1 )
    !third%vec(1) = 1d0
    !third%vec(2) = 1d0
    third1%vec(1) = third1%vec(1)
    third1%vec(2) = third1%vec(2)
    third1%vec(3) = third1%vec(3)
    !call vector_normalize(third, third1)
    costh1 = dcos(theta1)
    costh2 = dcos(theta2)
    if ( 0.1d0 < dabs( third1%vec(1) ) ) then
       do i=1,100
          !f()
          a = inner_product( first1, third1 ) - costh1
          !g()
          d = inner_product( second1, third1 ) - costh2
          !df/dy3
          b = first1%vec(2) - first1%vec(1)*third1%vec(2) / third1%vec(1)
          !df/dz3
          c = first1%vec(3) - first1%vec(1)*third1%vec(3) / third1%vec(1)
          !dg/dy3
          e = second1%vec(2) - second1%vec(1)*third1%vec(2) / third1%vec(1)
          !dg/dz3
          f = second1%vec(3) - second1%vec(1)*third1%vec(3) / third1%vec(1)
          !$B0J2<$NO"N)J}Dx<0$r2r$/!#(B
          !a + b dy3 + c dz3 = 0
          !d + e dy3 + f dz3 = 0
          dy3 = ( c*d - f*a ) / ( f*b - c*e )
          dz3 = ( a*e - b*d ) / ( f*b - c*e )
          count = 0
          do
             len = ( third1%vec(2) + dy3 )**2 + ( third1%vec(3) + dz3 )**2
             if ( len < 1d0 ) exit
             dy3 = dy3 * 0.5
             dz3 = dz3 * 0.5
             count = count + 1
          enddo
          write(STDERR,*) count
          if ( dabs( dy3 * 1e10 ) < dabs( third1%vec(2) ) .and. dabs( dz3 * 1e10 ) < dabs( third1%vec(3) ) )exit
          third1%vec(2) = third1%vec(2) + dy3
          third1%vec(3) = third1%vec(3) + dz3
          third1%vec(1) = sqrt(1d0 - third1%vec(2)**2 - third1%vec(3)**2 )
          !write(6,*) a,b,c,d,e,f
          write(STDERR,*) i, third1%vec(1),third1%vec(2),third1%vec(3)
       enddo
    else if ( 0.1d0 < dabs( third1%vec(2) ) ) then
       do i=1,100
          !f()
          a = inner_product( first1, third1 ) - costh1
          !g()
          d = inner_product( second1, third1 ) - costh2
          !df/dx3
          b = first1%vec(1) - first1%vec(2)*third1%vec(1) / third1%vec(2)
          !df/dz3
          c = first1%vec(3) - first1%vec(2)*third1%vec(3) / third1%vec(2)
          !dg/dx3
          e = second1%vec(1) - second1%vec(2)*third1%vec(1) / third1%vec(2)
          !dg/dz3
          f = second1%vec(3) - second1%vec(2)*third1%vec(3) / third1%vec(2)
          !$B0J2<$NO"N)J}Dx<0$r2r$/!#(B
          !a + b dx3 + c dz3 = 0
          !d + e dx3 + f dz3 = 0
          dx3 = ( c*d - f*a ) / ( f*b - c*e )
          dz3 = ( a*e - b*d ) / ( f*b - c*e )
          count = 0
          do
             len = ( third1%vec(1) + dx3 )**2 + ( third1%vec(3) + dz3 )**2
             if ( len < 1d0 ) exit
             dx3 = dx3 * 0.5
             dz3 = dz3 * 0.5
             count = count + 1
          enddo
          write(STDERR,*) count
          if ( dabs( dx3 * 1e10 ) < dabs( third1%vec(1) ) .and. dabs( dz3 * 1e10 ) < dabs( third1%vec(3) ) )exit
          third1%vec(1) = third1%vec(1) + dx3
          third1%vec(3) = third1%vec(3) + dz3
          third1%vec(2) = sqrt(1d0 - third1%vec(1)**2 - third1%vec(3)**2 )
          !write(6,*) a,b,c,d,e,f
          write(STDERR,*) i, third1%vec(1),third1%vec(2),third1%vec(3)
       enddo
    else
       do i=1,100
          !f()
          a = inner_product( first1, third1 ) - costh1
          !g()
          d = inner_product( second1, third1 ) - costh2
          !df/dx3
          b = first1%vec(1) - first1%vec(3)*third1%vec(1) / third1%vec(3)
          !df/dy3
          c = first1%vec(2) - first1%vec(3)*third1%vec(2) / third1%vec(3)
          !dg/dx3
          e = second1%vec(1) - second1%vec(3)*third1%vec(1) / third1%vec(3)
          !dg/dy3
          f = second1%vec(2) - second1%vec(3)*third1%vec(2) / third1%vec(3)
          !$B0J2<$NO"N)J}Dx<0$r2r$/!#(B
          !a + b dx3 + c dy3 = 0
          !d + e dx3 + f dy3 = 0
          dx3 = ( c*d - f*a ) / ( f*b - c*e )
          dy3 = ( a*e - b*d ) / ( f*b - c*e )
          count = 0
          do
             len = ( third1%vec(1) + dx3 )**2 + ( third1%vec(2) + dy3 )**2
             if ( len < 1d0 ) exit
             dx3 = dx3 * 0.5
             dy3 = dy3 * 0.5
             count = count + 1
          enddo
          write(STDERR,*) count
          if ( dabs( dx3 * 1e10 ) < dabs( third1%vec(1) ) .and. dabs( dy3 * 1e10 ) < dabs( third1%vec(2) ) )exit
          third1%vec(1) = third1%vec(1) + dx3
          third1%vec(2) = third1%vec(2) + dy3
          third1%vec(3) = sqrt(1d0 - third1%vec(1)**2 - third1%vec(2)**2 )
          !write(6,*) a,b,c,d,e,f
          write(STDERR,*) i, third1%vec(1),third1%vec(2),third1%vec(3)
       enddo
    endif
    third = third1
  end subroutine thirdvector1
end module vector_module




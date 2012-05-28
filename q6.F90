! -*- f90 -*-
module q6_module
  use physconst_module
  implicit none

contains
  function qlmsum(l,m,n,vx,vy,vz)
    integer :: l,m,n
    real(kind=8) :: vx(n),vy(n),vz(n)
    complex(kind=8) :: qlmsum
    complex(kind=8) :: sum
    !complex(kind=8) :: qlm
    integer :: i
    complex(kind=8) :: v
    sum=0d0
    !write(STDERR,*) n
    do i=1,n
       v=qlm(l,m,vx(i),vy(i),vz(i))
       sum=sum+v
    enddo
    qlmsum=sum
  end function qlmsum

  function ylm2(l,m,x,y)
    integer :: l,m
    real(kind=8) :: x
    complex(kind=8) :: y
    !complex(kind=8) :: ylm,ylm2,y
    complex(kind=8) :: ylm2
    if(m.lt.0)then
       !ylm2=(-1)**m*dconjg(ylm(l,-m,x,y))
       ylm2=(-1)**m*conjg(ylm(l,-m,x,y))
    else
       ylm2=ylm(l,m,x,y)
    endif
  end function ylm2
      
  function qlm(l,m,x,y,z)
    real(kind=8) :: x,y,z
    complex(kind=8) :: qlm
    complex(kind=8) :: w
    integer :: l,m
    !complex(kind=8) :: ylm2
    !qlm=ylm2(l,m,z,dcmplx(x,y)/sqrt(x**2+y**2))
    !
    !This causes internal error on ifc! 2004-08-06
    !
    w = cmplx(x,y) / sqrt(x**2+y**2)
    qlm=ylm2(l,m,z,w)
  end function qlm

  function ylm(l,m,costh,eiphi)
    complex(kind=8) :: ylm,v
    real(kind=8) :: costh
    complex(kind=8) :: eiphi
    integer :: l,m
    real(kind=8) :: prod(0:12)
    data prod/1d0,1d0,2d0,6d0,24d0,120d0,720d0,5040d0,40320d0,362880d0,3628800d0,39916800d0,479001600d0/
    !real(kind=8) :: plgndr
    v=sqrt((2d0*l+1d0)*prod(l-m)/(4d0*pi*prod(l+m)))*plgndr(l,m,costh)*eiphi**m
    ylm=v
  end function ylm

  FUNCTION plgndr(l,m,x)
    INTEGER :: l,m
    real(kind=8) :: plgndr,x
    INTEGER :: i,ll
    real(kind=8) :: fact,pll,pmm,pmmp1,somx2
    if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.d0) then
       write(STDERR,*) 'bad arguments in plgndr'
       stop
    endif
    pmm=1.d0
    if(m.gt.0) then
       somx2=sqrt((1.d0-x)*(1.d0+x))
       fact=1.d0
       do i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
       enddo
    endif
    if(l.eq.m) then
       plgndr=pmm
    else
       pmmp1=x*(2*m+1)*pmm
       if(l.eq.m+1) then
          plgndr=pmmp1
       else
          do ll=m+2,l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          enddo
          plgndr=pll
       endif
    endif
  end FUNCTION plgndr
!  (C) Copr. 1986-92 Numerical Recipes Software 1(-V%'2150)-3.

  !
  !d Y6m / dz, z==costh
  !
  function dy6mdz( m, z, eiphi )
    complex(kind=8) :: dy6mdz, result
    integer, intent(IN) :: m
    real(kind=8) :: z
    complex(kind=8) :: eiphi
    if ( m == -6 ) then
       result = (-3*Sqrt(3003/Pi)*z*(-1 + z**2)**2)/(32d0*eiphi**6)
    else if ( m == -5 ) then
       result = (3*Sqrt(1001/Pi)*Sqrt(1 - z**2)*(1 - 7*z**2 + 6*z**4))/(32d0*eiphi**5)
    else if ( m == -4 ) then
       result = (3*Sqrt(91/(2d0*Pi))*z*(13 - 46*z**2 + 33*z**4))/(16d0*eiphi**4)
    else if ( m == -3 ) then
       result = (-3*Sqrt(1365/Pi)*Sqrt(1 - z**2)*(1 - 15*z**2 + 22*z**4))/(32d0*eiphi**3)
    else if ( m == -2 ) then
       result = -(Sqrt(1365/Pi)*z*(19 - 102*z**2 + 99*z**4))/(32d0*eiphi**2)
    else if ( m == -1 ) then
       result =  (Sqrt(273/(2d0*Pi))*(5 - 100*z**2 + 285*z**4 - 198*z**6))/(16d0*eiphi*Sqrt(1 - z**2))
    else if ( m == 0 ) then
       result = (21*Sqrt(13/Pi)*z*(5 - 30*z**2 + 33*z**4))/16d0
    else if ( m == 1 ) then
       result = (eiphi*Sqrt(273/(2d0*Pi))*(-5 + 100*z**2 - 285*z**4 + 198*z**6))/(16d0*Sqrt(1 - z**2))
    else if ( m == 2 ) then
       result = -(eiphi**2*Sqrt(1365/Pi)*z*(19 - 102*z**2 + 99*z**4))/32d0
    else if ( m == 3 ) then
       result = (3*eiphi**3*Sqrt(1365/Pi)*Sqrt(1 - z**2)*(1 - 15*z**2 + 22*z**4))/32d0 
    else if ( m == 4 ) then
       result =  (3*eiphi**4*Sqrt(91/(2d0*Pi))*z*(13 - 46*z**2 + 33*z**4))/16d0
    else if ( m == 5 ) then
       result = (-3*eiphi**5*Sqrt(1001/Pi)*Sqrt(1 - z**2)*(1 - 7*z**2 + 6*z**4))/32d0
    else if ( m == 6 ) then
       result = (-3*eiphi**6*Sqrt(3003/Pi)*z*(-1 + z**2)**2)/32d0
    endif
    dy6mdz = result
  end function dy6mdz

  !
  !確認のため、Mathematicaの式から球面調和関数を算出する。
  !
  function y6m( m, z, eiphi )
    complex(kind=8) :: y6m, result
    integer, intent(IN) :: m
    real(kind=8) :: z
    complex(kind=8) :: eiphi
    if ( m == -6 ) then
       result = (Sqrt(3003/Pi)*(1 - z**2)**3)/(64d0*eiphi**6)
    else if ( m == -5 ) then
       result = (3*Sqrt(1001/Pi)*z*(1 - z**2)**2.5d0)/(32d0*eiphi**5)
    else if ( m == -4 ) then
       result =  (3*Sqrt(91/(2d0*Pi))*(1 - z**2)**2*(-1 + 11*z**2))/(32d0*eiphi**4)
    else if ( m == -3 ) then
       result = (Sqrt(1365/Pi)*z*(1 - z**2)**1.5d0*(-3 + 11*z**2))/(32d0*eiphi**3)
    else if ( m == -2 ) then
       result =   (Sqrt(1365/Pi)*(1 - z**2)*(1 - 18*z**2 + 33*z**4))/(64d0*eiphi**2)
    else if ( m == -1 ) then
       result = (Sqrt(273/(2d0*Pi))*z*Sqrt(1 - z**2)*(5 - 30*z**2 + 33*z**4))/(16d0*eiphi)
    else if ( m == 0 ) then
       result = (Sqrt(13/Pi)*(-5 + 105*z**2 - 315*z**4 + 231*z**6))/32d0
    else if ( m == 1 ) then
       result =  -(eiphi*Sqrt(273/(2d0*Pi))*z*Sqrt(1 - z**2)*(5 - 30*z**2 + 33*z**4))/16d0
    else if ( m == 2 ) then
       result = (eiphi**2*Sqrt(1365/Pi)*(1 - z**2)*(1 - 18*z**2 + 33*z**4))/64d0
    else if ( m == 3 ) then
       result = -(eiphi**3*Sqrt(1365/Pi)*z*(1 - z**2)**1.5d0*(-3 + 11*z**2))/32d0
    else if ( m == 4 ) then
       result = (3*eiphi**4*Sqrt(91/(2d0*Pi))*(1 - z**2)**2*(-1 + 11*z**2))/32d0
    else if ( m == 5 ) then
       result = (-3*eiphi**5*Sqrt(1001/Pi)*z*(1 - z**2)**2.5d0)/32d0
    else if ( m == 6 ) then
       result = (eiphi**6*Sqrt(3003/Pi)*(1 - z**2)**3)/64d0
    end if
    y6m = result 
  end function y6m
  !
  !d Y6m / d(exp(i phi))
  !
  function dy6mdeiphi( m, z, eiphi )
    complex(kind=8) :: dy6mdeiphi, result
    integer, intent(IN) :: m
    real(kind=8) :: z
    complex(kind=8) :: eiphi
    result = ylm2( 6, m, z, eiphi ) * m  / eiphi
    if ( m < 0 ) result = - result
    dy6mdeiphi = result
  end function dy6mdeiphi

end module q6_module

#undef DEBUG
#ifdef DEBUG
program main
  use q6_module
  implicit none
  integer :: m
  real(kind=8) :: x,y,z
  real(kind=8) :: xx,yy,zz
  complex(kind=8) :: w, ww, y6md, dy6m, y0, y1
  y = 0.3001d0
  z = 0.3001d0
  x = sqrt( 1 - z**2 - y**2 )
  w = 0.301d0  !cmplx( x, y )/sqrt(x**2+y**2)
  yy = 0.3d0
  zz = 0.3d0
  xx = sqrt( 1 - zz**2 - yy**2 )
  ww = 0.3d0 !cmplx( xx, yy )/sqrt(xx**2+yy**2)
  do m=-6,6
     !まずはwwでの微分の評価。
     y0 = ylm2( 6, m, z, w )
     y1 = y6m( m,z,w )
     y6md = ylm2( 6, m, z, ww ) - ylm2( 6, m, z, w )
     dy6m = dy6mdeiphi( m,z,w ) * (ww - w)
     write(6,*) ( dy6m - y6md ) / y6md, y1, (y1 - y0)/y0
     y6md = ylm2( 6, m, zz, w ) - ylm2( 6, m, z, w )
     dy6m = dy6mdz( m, z, w ) * (zz-z)
     write(6,*) ( dy6m - y6md ) / y6md 
     write(6,*)
  enddo
end program main
#endif
  

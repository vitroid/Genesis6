! -*- f90 -*-

!This module contains manupilations of quaternion
module quat_module
  
contains
  subroutine abc2abcd(ea,eb,ec,a,b,c,d)
    real(kind=8),intent(in) ::ea,eb,ec
    real(kind=8),intent(out) :: a,b,c,d
    real(kind=8) :: eah,ebc,ecb
    eah = ea/2d0
    ebc = (eb+ec)/2d0
    ecb = (ec-eb)/2d0
    a = dcos(eah)*dcos(ebc)
    b = dsin(eah)*dcos(ecb)
    c = dsin(eah)*dsin(ecb)
    d = dcos(eah)*dsin(ebc)
  end subroutine abc2abcd

  subroutine abcd2abc_bar(nmole, qua1, qua2, qua3, qua4, eul1, eul2, eul3)
    implicit none
    integer nmole
    real(kind=8), intent(IN)  :: qua1(*), qua2(*), qua3(*), qua4(*)
    real(kind=8), intent(OUT) :: eul1(*), eul2(*), eul3(*)
    !
    !     <<  Relation Quarternion and Euler angle  >>
    !     --- definition of quarternion in my programs ---
    !         q1 = sin(the/2)sin((psi-phi)/2)
    !         q2 = sin(the/2)cos((psi-phi)/2)
    !         q3 = cos(the/2)sin((psi+phi)/2)
    !         q4 = cos(the/2)cos((psi+phi)/2)
    !
    integer :: i
    real(kind=8) :: pi, pi2, q1o, q2o, q3o, q4o
    real(kind=8) :: q1, q2, q3, q4, coef, sthh, cthh, thh, the
    real(kind=8) :: pmh, pph, psi, phi

    pi = 4.0d0 * datan(1.0d0)
    pi2 = 2.0d0 * pi

    do i = 1, nmole
       q1o = qua1(i)
       q2o = qua2(i)
       q3o = qua3(i)
       q4o = qua4(i)
       coef = 1.0d0/dsqrt(q1o*q1o + q2o*q2o + q3o*q3o + q4o*q4o)
       q1 = coef*q1o
       q2 = coef*q2o
       q3 = coef*q3o
       q4 = coef*q4o
       sthh = dsqrt(q1*q1 + q2*q2)
       cthh = dsqrt(q3*q3 + q4*q4)
       thh = datan2(sthh, cthh)
       the = 2.0d0*thh
       if (sthh .lt. 1.0d-14) then
          psi = 0.0d0
          phi = datan2(2.0d0*q3*q4, -q3*q3+q4*q4)
       else if (cthh .lt. 1.0d-14) then
          psi = 0.0d0
          phi = datan2(-2.0d0*q1*q2, -q1*q1+q2*q2)
       else
          pmh = datan2(q1, q2)
          pph = datan2(q3, q4)
          psi = pmh + pph
          phi = pph - pmh
       endif
       if (the .lt. 0.0d0) the = 0.0d0
       if (the .gt. pi) the = pi
       psi = psi - dnint(psi/pi2)*pi2
       phi = phi - dnint(phi/pi2)*pi2
       eul1(i) = the
       eul2(i) = psi
       eul3(i) = phi
    enddo
  end subroutine abcd2abc_bar

  subroutine abcd2abc(a,b,c,d,ea,eb,ec)
    implicit none
    real(kind=8),intent(in) :: a,b,c,d
    real(kind=8),intent(out) :: ea,eb,ec
    real(kind=8) :: qa,qb,qc,qd
    real(kind=8) :: thh,p,s,costhh,sinthh,r
    real(kind=8),parameter :: pid&
         & =3.1415926535897932384626433832795028841971693993d0*2d0
    r=1d0/sqrt(a**2+b**2+c**2+d**2)
    qa=a*r
    qb=b*r
    qc=c*r
    qd=d*r
    ea=acos(2d0*(qa**2+qd**2)-1d0)
    thh=ea/2d0
    sinthh=sin(thh)
    costhh=cos(thh)
    p=acos(qa/costhh)
    s=acos(qb/sinthh)
    if(qd<0d0) p=pid-p
    if(qc>0d0)then
       ec=p+s
       eb=p-s
    else
       ec=p-s
       eb=p+s
    end if
  end subroutine abcd2abc

!qaddはインライン展開すると正常に動かないようだ。
!OCL NOEVAL
  subroutine qadd(a1,b1,c1,d1,a2,b2,c2,d2)
    real(kind=8),intent(inout) :: a1,b1,c1,d1
    real(kind=8),intent(in) :: a2,b2,c2,d2
    real(kind=8) :: a3,b3,c3,d3,dd
    a3=a1*a2-b1*b2-c1*c2-d1*d2
    b3=a1*b2+b1*a2+c1*d2-d1*c2
    c3=a1*c2+c1*a2-b1*d2+d1*b2
    d3=a1*d2+d1*a2+b1*c2-c1*b2
    a1=a3
    b1=b3
    c1=c3
    d1=d3
    !
    !右回りと左回りのうち短い方(cos(phi/2)が正な方)を選ぶ。
    !
    if ( a1 < 0 ) then
       a1 = -a1
       b1 = -b1
       c1 = -c1
       d1 = -d1
    endif
    dd=a1**2+b1**2+c1**2+d1**2
    if(abs(dd-1d0) > 0.000001d0)then
       write(STDERR,*) "qadd",dd
    endif
  end subroutine qadd

  subroutine qmul(a1,b1,c1,d1,x)
    real(kind=8),intent(inout) :: a1,b1,c1,d1
    real(kind=8),intent(in) :: x
    real(kind=8) :: phi,sine
#ifdef RIGIDDEBUG
    real(kind=8) :: dd
#endif
    if(a1 >= 1d0 .or. a1 <= -1d0)return
    phi=acos(a1)
    !ex is always positive
    sine=dsqrt(1d0-a1**2)
    if(b1.lt.0d0)then
       phi=-phi
       sine=-sine
    endif
    phi  = phi*x
    sine = dsin(phi)/sine
    a1   = dcos(phi)
    b1   = b1*sine
    c1   = c1*sine
    d1   = d1*sine
#ifdef RIGIDDEBUG
    dd   = a1**2 + b1**2 + c1**2 + d1**2
    if( dabs(dd-1d0) > 0.01d0 )then
       write( STDERR, * ) "qmul", dd
    endif
#endif
  end subroutine qmul
  
  subroutine qmul2(a1,b1,c1,d1,x)
    real(kind=8),intent(inout) :: a1,b1,c1,d1
    real(kind=8),intent(in) :: x
    real(kind=8) dd
    if(a1 >= 1d0 .or. a1 <= -1d0)return
    a1=1d0-(1d0-a1)*x**2
    b1=b1*x
    c1=c1*x
    d1=d1*x
    dd=a1**2+b1**2+c1**2+d1**2
    if(abs(dd-1d0) > 0.000001d0)then
       write(STDERR,*) "qmul2",dd
    endif
  end subroutine qmul2

  !
  !z,y,x軸の順に回転させるようなquaternionを返す。ノート参照
  !x,y,zは各軸回りの回転角の半分を与える。
  !
  subroutine qrotator(x,y,z,a,b,c,d)
    real(kind=8), intent(out) :: a,b,c,d
    real(kind=8), intent(in)  :: x,y,z
    real(kind=8) :: sx,sy,sz,cx,cy,cz
    sx = sin(x)
    sy = sin(y)
    sz = sin(z)
    cx = cos(x)
    cy = cos(y)
    cz = cos(z)
    a = cx*cy*cz - sx*sy*sz
    b = sx*cy*cz + cx*sy*sz
    c = cx*sy*cz - sx*cy*sz
    d = cx*cy*sz + sx*sy*cz
  end subroutine qrotator

  !
  !回転が微小回転とみなせる場合は、もっと単純に(1,Δx,Δy,Δz)を加える
  !だけでよい。(規格化が必要)
  !
  subroutine qinfrotator(x,y,z,a,b,c,d)
    real(kind=8), intent(out) :: a,b,c,d
    real(kind=8), intent(in)  :: x,y,z
    a = 1
    b = x
    c = y
    d = z
    call qnormalize( a, b, c, d )
  end subroutine qinfrotator
  !
  !
  !
  subroutine qnormalize(a,b,c,d)
    real(kind=8), intent(inout) :: a,b,c,d
    real(kind=8) :: r
    r = 1d0/dsqrt( a**2 + b**2 + c**2 + d**2 )
    a = a * r
    b = b * r
    c = c * r
    d = d * r
  end subroutine qnormalize

end module quat_module

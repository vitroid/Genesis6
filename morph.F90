subroutine qadd(a1,b1,c1,d1,a2,b2,c2,d2)
  implicit none
  real(kind=8),intent(inout) :: a1,b1,c1,d1
  real(kind=8),intent(in) :: a2,b2,c2,d2
  real(kind=8) :: a3,b3,c3,d3
  a3=a1*a2-b1*b2-c1*c2-d1*d2
  b3=a1*b2+b1*a2+c1*d2-d1*c2
  c3=a1*c2+c1*a2-b1*d2+d1*b2
  d3=a1*d2+d1*a2+b1*c2-c1*b2
  a1=a3
  b1=b3
  c1=c3
  d1=d3
  if ( a1 < 0 ) then
     a1 = -a1
     b1 = -b1
     c1 = -c1
     d1 = -d1
  endif
end subroutine qadd

!²óÅ¾³Ñ¤òxÇÜ¡£
subroutine qmul(a1,b1,c1,d1,x)
  implicit none
  real(kind=8),intent(inout) :: a1,b1,c1,d1
  real(kind=8),intent(in) :: x
  real(kind=8) :: phi,sine
  if(a1.ge.1d0.or.a1.le.-1d0)return
  phi=acos(a1)
  !ex is always positive
  sine=sqrt(1d0-a1**2)
  if(b1.lt.0d0)then
     phi=-phi
     sine=-sine
  endif
  phi=phi*x
  sine=sin(phi)/sine
  a1=cos(phi)
  b1=b1*sine
  c1=c1*sine
  d1=d1*sine
end subroutine qmul



subroutine die( msg )
  implicit none
  character(len=*) :: msg
  write(STDERR,*) msg
  stop
end subroutine die

program main
  implicit none
  integer, parameter :: MAXMOL=100
  character(len=5) :: tag
  character(len=8) :: lastid
  integer          :: n, lastn, i, ncompo
  integer          :: loop
  real(kind=8)     :: x(MAXMOL,2)
  real(kind=8)     :: y(MAXMOL,2)
  real(kind=8)     :: z(MAXMOL,2)
  real(kind=8)     :: a(MAXMOL,2)
  real(kind=8)     :: b(MAXMOL,2)
  real(kind=8)     :: c(MAXMOL,2)
  real(kind=8)     :: d(MAXMOL,2)
  real(kind=8)     :: t
  real(kind=8)     :: xx,yy,zz,aa,bb,cc,dd

  lastn  = 0
  ncompo = 0
  do while ( ncompo < 2 )
     read(5,*,END=90) tag
     if( tag(1:5) .eq. "@ID08" )then
        read( 5, * ) lastid
        cycle
     endif
     if ( tag(1:5) .eq. "@NX4A" ) then
        read(5,*) n
        if ( lastn .ne. 0 .and. n .ne. lastn ) then
           call die( "Size mismatch." )
        endif
        lastn = n
        ncompo = ncompo + 1
        do i=1,n
           read(5,*) x(i,ncompo),y(i,ncompo),z(i,ncompo),a(i,ncompo),b(i,ncompo),c(i,ncompo),d(i,ncompo)
        enddo
     endif
  enddo
90 continue

  do loop=0,100
     t = loop
     t = t * 0.01
     write(6,'("@ID08")')
     write(6,'(a8)') lastid
     write(6,'("@NCMP")')
     write(6,*) 1
     write(6,'("@NX4A")')
     write(6,*) n
     do i=1,n
        xx = x(i,1) + t * ( x(i,2) - x(i,1) )
        yy = y(i,1) + t * ( y(i,2) - y(i,1) )
        zz = z(i,1) + t * ( z(i,2) - z(i,1) )
        aa = a(i,2)
        bb = b(i,2)
        cc = c(i,2)
        dd = d(i,2)
        call qadd( aa,bb,cc,dd, a(i,1),-b(i,1),-c(i,1),-d(i,1) )
        call qmul( aa,bb,cc,dd, t )
        call qadd( aa,bb,cc,dd, a(i,1), b(i,1), c(i,1), d(i,1) )
        write(6,'(7(e17.10,1x))') xx,yy,zz,aa,bb,cc,dd
     enddo
  enddo
end program main

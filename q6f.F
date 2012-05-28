      function qx(l,n,vx,vy,vz)
      implicit none
      real(kind=8) :: pi
      parameter ( pi = 3.14159265358979d0 )
      real(kind=8) :: qx
      integer :: n
      real(kind=8) :: vx(n),vy(n),vz(n)
      real(kind=8) :: sum,v
      integer :: l
      complex(kind=8) :: qlmbar
      integer :: m
      sum=0d0
                                !l=6
      do m=-l,l
                                !v=cdabs(qlmbar(l,m,n,vx,vy,vz))**2
        v=abs(qlmbar(l,m,n,vx,vy,vz))**2
        sum=sum+v
      enddo
      qx=sqrt(sum*4d0*pi/(2d0*l+1d0))
      end

      function qlmbar(l,m,n,vx,vy,vz)
      implicit none
      integer :: l,m,n
      real(kind=8) :: vx(n),vy(n),vz(n)
      complex(kind=8) :: qlmbar
      complex(kind=8) :: sum
      complex(kind=8) :: qlm
      integer :: i
      complex(kind=8) :: v
      sum=0d0
      do i=1,n
        v=qlm(l,m,vx(i),vy(i),vz(i))
        sum=sum+v
      enddo
      qlmbar=sum/n
      end
      
      function ylm2(l,m,x,y)
      implicit none
      integer :: l,m
      real(kind=8) :: x
      complex(kind=8) :: y
      complex(kind=8) :: ylm2
      complex(kind=8) :: ylm
      if(m.lt.0)then
                                !ylm2=(-1)**m*dconjg(ylm(l,-m,x,y))
        ylm2=(-1)**m*conjg(ylm(l,-m,x,y))
      else
        ylm2=ylm(l,m,x,y)
      endif
      end
      
      function qlm(l,m,x,y,z)
      implicit none
      real(kind=8) :: x,y,z
      complex(kind=8) :: qlm
      complex(kind=8) :: ylm2
      complex(kind=8) :: w
      integer :: l,m
                                !This causes internal error on ifc! 2004-08-06
                                !
      w = cmplx(x,y) / sqrt(x**2+y**2)
      qlm=ylm2(l,m,z,w)
      end
      
      function ylm(l,m,costh,eiphi)
      implicit none
      complex(kind=8) :: ylm,v
      real(kind=8) :: costh
      complex(kind=8) :: eiphi
      integer :: l,m
      real(kind=8) :: prod(0:12)
      real(kind=8) :: pi
      parameter ( pi = 3.14159265358979d0 )
      data prod/1d0,1d0,2d0,6d0,24d0,120d0,720d0,5040d0,40320d0,
     :     362880d0,3628800d0,39916800d0,479001600d0/
      real(kind=8) :: plgndr
      v=sqrt((2d0*l+1d0)*prod(l-m)/(4d0*pi*prod(l+m)))*
     :     plgndr(l,m,costh)*eiphi**m
      ylm=v
      end
      
      FUNCTION plgndr(l,m,x)
      implicit none
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
      end
!  (C) Copr. 1986-92 Numerical Recipes Software 1(-V%'2150)-3.

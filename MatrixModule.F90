module matrix_module
  use error_module
  implicit none
  !
  !圧力テンソル
  !
  type sTensor
     real(kind=8) :: mat33(3,3)
  end type sTensor

contains
  subroutine InvSym33(xx,yy,zz,xy,yz,zx)
    implicit none
    real(kind=8) :: xx,yy,zz,xy,yz,zx
    real(kind=8) :: xx1,yy1,zz1,xy1,yz1,zx1
    real(kind=8) :: denom
    denom = 1d0/(xx*yy*zz-xx*yz**2-yy*zx**2-zz*xy**2+2d0*xy*yz*zx)
    xx1 = yy*zz - yz**2
    yy1 = zz*xx - zx**2
    zz1 = xx*yy - xy**2
    xy1 = zx*yz - xy*zz
    yz1 = xy*zx - yz*xx
    zx1 = yz*xy - zx*yy
    xx = xx1*denom
    yy = yy1*denom
    zz = zz1*denom
    xy = xy1*denom
    yz = yz1*denom
    zx = zx1*denom
  end subroutine InvSym33

  subroutine Inv33(a,b,c,p,q,r,x,y,z)
    implicit none
    real(kind=8) :: a,b,c,p,q,r,x,y,z
    real(kind=8) :: aa,bb,cc,pp,qq,rr,xx,yy,zz,denom
    
    denom = 1d0/(a*q*z - a*r*y - p*b*z + p*c*y + x*b*r - x*c*q)
    aa = ( q*z - r*y ) * denom
    bb = ( c*y - b*z ) * denom
    cc = ( b*r - c*q ) * denom
    pp = ( r*x - p*z ) * denom
    qq = ( a*z - c*x ) * denom
    rr = ( c*p - a*r ) * denom
    xx = ( p*y - q*x ) * denom
    yy = ( b*x - a*y ) * denom
    zz = ( a*q - b*p ) * denom
    !  | a b c |^-1                  1                  | qz-ry cy-bz br-cq |
    !  | p q r |  =    ------------------------------   | rx-pz az-cx cp-ar |
    !  | x y z |    (aqz - ary - pbz + pcy + xbr - xcq) | py-qx bx-ay aq-bp |
    a = aa
    b = bb
    c = cc
    p = pp
    q = qq
    r = rr
    x = xx
    y = yy
    z = zz
  end subroutine Inv33

  subroutine normalize(x,y,z)
    real(kind=8), intent(inout) :: x,y,z
    real(kind=8) :: r
    r = 1d0 / sqrt(x**2 + y**2 + z**2 )
    x = x * r
    y = y * r
    z = z * r
  end subroutine normalize

  logical function op(xi,yi,zi,xj,yj,zj,ax,ay,az)
    real(kind=8), intent(in) :: xi,yi,zi,xj,yj,zj
    real(kind=8), intent(out) :: ax,ay,az
    op = .false.
    !/*fprintf(stderr,"a\n");*/
    if((xi*xi+yi*yi+zi*zi<0.01).or.(xj*xj+yj*yj+zj*zj<0.01))return

    ax=yi*zj-zi*yj
    ay=zi*xj-xi*zj
    az=xi*yj-yi*xj
    
    op = ( 0.01 < ax**2 + ay**2 + az**2 )
  end function op
  !
  !回転変換行列を四元数に変換。RotMX2.cから移植。
  !
  subroutine rotmx2(ix,iy,iz,jx,jy,jz,kx,ky,kz,qa,qb,qc,qd)
    real(kind=8), intent(INOUT)  :: ix,iy,iz,jx,jy,jz,kx,ky,kz
    real(kind=8), intent(OUT) :: qa,qb,qc,qd

    real(kind=8) ::  ax,ay,az
    real(kind=8) ::  i0x,i0y,i0z
    real(kind=8) ::  x0x,x0y,x0z
    real(kind=8) ::  ox,oy,oz
    real(kind=8) ::  cosine,cosh,sinh
    real(kind=8) ::  t

    !/* i軸をx軸に移す回転の軸は、iとxの2分面となる。
    !   j軸をy軸に移す回転の軸は、jとyの2分面となる。
    !   そして、それらを同時にみたす回転の軸は、それらの交線となる。
    !   交線は、2つの面の法線のいずれとも直交する=外積である。*/

    if(.not.op(ix-1.0,iy,iz,jx,jy-1.0,jz,ax,ay,az))then
       if(.not.op(ix-1.0,iy,iz,kx,ky,kz-1.0,ax,ay,az))then
          if(.not.op(kx,ky,kz-1.0,jx,jy-1.0,jz,ax,ay,az))then
             write(STDERR,*) "outer prod error"
             call die( 0, "Matrix 1" )
          endif
       endif
    endif

    call normalize(ax,ay,az)
    
    !/*回転軸aが求まった。*/
    x0x = 1
    x0y = 0
    x0z = 0
    i0x = ix
    i0y = iy
    i0z = iz
    t = ax/(ax*ax+ay*ay+az*az)
    i0x = i0x - t*ax
    i0y = i0y - t*ay
    i0z = i0z - t*az
    x0x = x0x - t*ax
    x0y = x0y - t*ay
    x0z = x0z - t*az
    !/*check*/
    !/*fprintf(stderr,"%f %f %f\n",i0x,i0y,i0z)*/
    if(i0x*i0x+i0y*i0y+i0z*i0z<0.01)then
       x0x = 0
       x0y = 1
       x0z = 0
       i0x = jx
       i0y = jy
       i0z = jz
       t = ay/(ax*ax+ay*ay+az*az)
       i0x = i0x - t*ax
       i0y = i0y - t*ay
       i0z = i0z - t*az
       x0x = x0x - t*ax
       x0y = x0y - t*ay
       x0z = x0z - t*az
       !/*check*/
       !/*fprintf(stderr,"%f %f %f t=%f\n",i0x,i0y,i0z,t)*/
    endif

    !/*fprintf(stderr,"check: %f %f\n",i0x*ax+i0y*ay+i0z*az,x0x*ax+x0y*ay+x0z*az)*/
    call normalize(i0x,i0y,i0z)
    call normalize(x0x,x0y,x0z)
    !/*inner product to determine angle*/
    cosine=i0x*x0x+i0y*x0y+i0z*x0z
    !/*fprintf(stderr,"cos %f\n",cosine)*/
    cosh=sqrt((1.0+cosine)*0.5)
    sinh=sqrt(1.0-cosh*cosh)
    !/*outer product to determine direction*/
    ox=i0y*x0z-i0z*x0y
    oy=i0z*x0x-i0x*x0z
    oz=i0x*x0y-i0y*x0x
    if(ox*ax+oy*ay+oz*az<0)then
       sinh=-sinh
    endif
    call normalize(ax,ay,az)
    !/*printf("%f %f %f\n",ax,ay,az)*/

    qa=cosh
    qb=-sinh*ax
    qc=+sinh*ay
    qd=-sinh*az
  end subroutine rotmx2

end module matrix_module

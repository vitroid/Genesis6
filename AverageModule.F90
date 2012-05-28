! -*- f90 -*-
module rigid_average_module
  use common_module
  use vector_module

  implicit none

  type sRigidAverage
     type(vector3),dimension(MAXMOL) :: com
     type(vector4),dimension(MAXMOL) :: quat
     type(vector4),dimension(MAXMOL) :: quat0
     integer :: width,interval,count
  end type sRigidAverage

  interface new
     module procedure RigidAverage_Constructor
  end interface

contains

  subroutine RigidAverage_Constructor(av)
    type(sRigidAverage),intent(out) :: av
    av%count=0
    av%interval=0
    av%width=0
  end subroutine RigidAverage_Constructor

  subroutine Average_Rigid(av,n,com,quat)
    use quat_module
    type(sRigidAverage),intent(inout) :: av
    integer,intent(in):: n
    type(vector3),dimension(*),intent(in) :: com
    type(vector4),dimension(*),intent(in) :: quat
    real(kind=8) :: ratio,a0,b0,c0,d0
    integer :: i
    av%count=av%count+1
    ratio=dble(av%count-1)/dble(av%width)
    if(av%count==1)then
       do i=1,n
          av%quat0(i)=quat(i)
          av%quat(i)%vec(1)=1d0
          av%quat(i)%vec(2)=0d0
          av%quat(i)%vec(3)=0d0
          av%quat(i)%vec(4)=0d0
          av%com(i)=com(i)
       enddo
    else
       av%com(1:n)%vec(1)=av%com(1:n)%vec(1)+com(1:n)%vec(1)
       av%com(1:n)%vec(2)=av%com(1:n)%vec(2)+com(1:n)%vec(2)
       av%com(1:n)%vec(3)=av%com(1:n)%vec(3)+com(1:n)%vec(3)
       do i=1,n
          !,A:9$$H$k!#(B
          a0=quat(i)%vec(1)
          b0=quat(i)%vec(2)
          c0=quat(i)%vec(3)
          d0=quat(i)%vec(4)
          !Difference from the last conf.
          call qadd(a0,b0,c0,d0,av%quat(i)%vec(1),-av%quat(i)%vec(2),-av%quat(i)%vec(3),-av%quat(i)%vec(4))
          av%quat0(i)=quat(i)
          call qmul2(a0,b0,c0,d0,ratio)
          call qadd(av%quat(i)%vec(1),av%quat(i)%vec(2),av%quat(i)%vec(3),av%quat(i)%vec(4),a0,b0,c0,d0)
       enddo
    endif
    if(av%count==av%width)then
       ratio = 1d0/av%width
       do i=1,n
          av%com(i)%vec(1:3) = av%com(i)%vec(1:3) * ratio
          a0 = +av%quat(i)%vec(1)
          b0 = -av%quat(i)%vec(2)
          c0 = -av%quat(i)%vec(3)
          d0 = -av%quat(i)%vec(4)
          call qadd(a0, b0, c0, d0, quat(i)%vec(1),quat(i)%vec(2),quat(i)%vec(3)&
               & ,quat(i)%vec(4))
          av%quat(i)%vec(1) = a0
          av%quat(i)%vec(2) = b0
          av%quat(i)%vec(3) = c0
          av%quat(i)%vec(4) = d0
       enddo
       av%count=0
    endif
  end subroutine Average_Rigid

  subroutine Average_RigidWriteNX4A(av,file,n)
    type(sRigidAverage),intent(in) ::av
    integer,intent(in) :: file,n
    integer :: i
    write(file,'("@NX4A")')
    write(file,*) n
    do i=1,n
       write(file,"(7(e17.10,1x))") av%com(i)%vec(1:3),av%quat(i)%vec(1:4)
    enddo
  end subroutine Average_RigidWriteNX4A
end module rigid_average_module

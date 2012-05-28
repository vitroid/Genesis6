! -*- f90 -*-
module time_module
  implicit none
  type sTime
     real(kind=8) :: dt
  end type sTime
contains
  subroutine time_initialize(t,dt)
    implicit none
    type(sTime),intent(out) :: t
    real(kind=8),intent(in) :: dt
    t%dt=dt
    return
  end subroutine time_initialize
  
  subroutine Time_WriteBinaryDTPS(t,file)
    implicit none
    integer,intent(IN) :: file
    type(sTime),intent(IN) :: t
    write(file)"@DTPS"
    write(file)t%dt
    return
  end subroutine Time_WriteBinaryDTPS
  
  subroutine Time_WriteDTPS(t,file)
    implicit none
    integer,intent(IN) :: file
    type(sTime),intent(IN) :: t
    write(file,1)
    write(file,*)t%dt
1   format("@DTPS")
    return
  end subroutine Time_WriteDTPS
  
  subroutine Time_ReadBinaryDTPS(t,file)
    implicit none
    integer,intent(IN) :: file
    type(sTime),intent(OUT) :: t
    real(kind=8) :: dt
    read(file)dt
    call time_initialize(t,dt)
    return
  end subroutine Time_ReadBinaryDTPS
  
  subroutine Time_ReadDTPS(t,file)
    implicit none
    integer,intent(IN) :: file
    type(sTime),intent(OUT) :: t
    real(kind=8) :: dt
    read(file,*)dt
    call time_initialize(t,dt)
    return
  end subroutine Time_ReadDTPS
  
  subroutine Time_ReadBinaryDTSC(t,file)
    implicit none
    integer,intent(IN) :: file
    type(sTime),intent(OUT) :: t
    real(kind=8) :: dt
    read(file)dt
    dt=dt*1d12
    call time_initialize(t,dt)
    return
  end subroutine Time_ReadBinaryDTSC
end module time_module

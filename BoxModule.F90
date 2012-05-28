! -*- f90 -*-
!試しにモジュール化してみようか。平成14年11月1日(金)
module box_module
  use vector_module
  implicit none
  type sBox
     type(vector3) :: size,invsize
     integer :: mode
  end type sBox
  integer, parameter :: BOX_NONE=0, BOX_ORTHO=1
  integer, parameter :: BOX_BOX3=0, BOX_MDVW=1
  interface new
     module procedure box_constructor
  end interface
  interface load
     module procedure box_loader
  end interface
  interface save
     module procedure box_save
  end interface
  interface load_binary
     module procedure box_binaryloader
  end interface
  interface save_binary
     module procedure box_savebinary
  end interface
  interface renormalize
     module procedure box_renormalize
  end interface
  interface volume
     module procedure box_volume
  end interface

  interface scale1
     module procedure box_scale1
  end interface

  interface scalev
     module procedure box_scalev
  end interface

  interface scale
     module procedure box_scale1
     module procedure box_scalev
  end interface

  interface active
     module procedure box_active
  end interface

contains
  function box_active( box )
    logical :: box_active
    type(sBox) :: box
    box_active = ( box%mode .ne. box_none )
  end function box_active

  subroutine box_renormalize(b,delta,cell)
    use vector_module
    implicit none
    type(sBox),intent(in) :: b
    type(vector3),intent(inout) :: delta
    type(vector3),intent(out)   :: cell
    if ( b%mode /= BOX_NONE ) then
       cell%vec(:)  = dnint( delta%vec(:) * b%invsize%vec(:) ) * b%size%vec(:)
       delta%vec(:) = delta%vec(:) - cell%vec(:)
    endif
  end subroutine box_renormalize
  
  subroutine box_renormalize0(b,dx,dy,dz,cx,cy,cz)
    use vector_module
    implicit none
    type(sBox),intent(in) :: b
    real(kind=8),intent(inout) :: dx,dy,dz
    real(kind=8),intent(out)   :: cx,cy,cz
    cx  = dnint( dx * b%invsize%vec(1) ) * b%size%vec(1)
    cy  = dnint( dy * b%invsize%vec(2) ) * b%size%vec(2)
    cz  = dnint( dz * b%invsize%vec(3) ) * b%size%vec(3)
    dx = dx - cx
    dy = dy - cy
    dz = dz - cz
  end subroutine box_renormalize0
  
  subroutine Box_ReadBXLA(b,file)
    implicit none
    type(sBox) :: b
    integer file
    real(kind=8) :: bxl,bxli
    read(file,*) bxl
    bxli             = 1d0 / bxl
    b%size%vec(:)    = bxl
    b%invsize%vec(:) = bxli
    return
  end subroutine Box_ReadBXLA
  
  subroutine Box_ReadBinaryBXLA(b,file)
    implicit none
    type(sBox) :: b
    integer file
    real(kind=8) :: bxl,bxli
    read(file) bxl
    bxli = 1d0/bxl
    bxli             = 1d0 / bxl
    b%size%vec(:)    = bxl
    b%invsize%vec(:) = bxli
    return
  end subroutine Box_ReadBinaryBXLA
  
  subroutine Box_ReadBOX3(b,file)
    implicit none
    type(sBox) :: b
    integer file
    read(file,*) b%size%vec(:)
    b%invsize%vec(:) = b%size%vec(:)
    call invert( b%invsize )
  end subroutine Box_ReadBOX3
  
  subroutine Box_ReadBinaryBOX3(b,file)
    implicit none
    type(sBox) :: b
    integer file
    read(file) b%size%vec(:)
    b%invsize%vec(:) = b%size%vec(:)
    call invert( b%invsize )
  end subroutine Box_ReadBinaryBOX3
  
  subroutine Box_WriteBOX3(b,file)
    implicit none
    type(sBox),intent(in) :: b
    integer,intent(in)    :: file
    integer               :: i
    write(file,101)
101 format("@BOX3")
    write(file,*) (b%size%vec(i),i=1,3)
  end subroutine Box_WriteBOX3
  
  subroutine Box_WriteMDVW(b,file)
    implicit none
    type(sBox),intent(in) :: b
    integer,intent(in)    :: file
    integer               :: i
    write(file,'("-center 0.0 0.0 0.0")')
    write(file,'("-fold")')
    write(file,'("-length ''(",2(f12.5,","),f12.5,")''")') (b%size%vec(i),i=1,3)
  end subroutine Box_WriteMDVW
  
  subroutine Box_WriteBinaryBOX3(b,file)
    implicit none
    type(sBox),intent(in) :: b
    integer,intent(in)    :: file
    integer               :: i
    write(file)"@BOX3"
    write(file) (b%size%vec(i),i=1,3)
  end subroutine Box_WriteBinaryBOX3
  
  subroutine Box_Version
    write(STDERR,*) "$Id: BoxModule.F90,v 1.9 2005/02/10 06:06:45 matto Exp $"
  end subroutine Box_Version
  
  function Box_Volume(b)
    use vector_module
    implicit none
    type(sBox),intent(in) :: b
    real(kind=8)          :: Box_Volume
    Box_Volume = volume( b%size )
  end function Box_Volume
  
  subroutine Box_Scalev(b, ratio)
    implicit none
    type(sBox) :: b
    type(vector3), intent(in) :: ratio
    call scale( b%size, ratio )
    b%invsize%vec(:) = 1d0 / b%size%vec(:)
  end subroutine Box_Scalev
  
  subroutine Box_Scale1(b, ratio)
    use vector_module
    implicit none
    type(sBox) :: b
    real(kind=8), intent(in) :: ratio
    call scale( b%size, ratio )
    b%invsize%vec(:) = 1d0 / b%size%vec(:)
  end subroutine Box_Scale1
  
  subroutine Box_Constructor(b)
    implicit none
    type(sBox),intent(inout) :: b
    b%mode=BOX_NONE
  end subroutine Box_Constructor
  
  subroutine Box_Loader(b,file,tag)
    implicit none
    type(sBox) :: b
    integer file
    character(len=5) :: tag
    if(tag == "@BXLA")then
       call Box_ReadBXLA(b,file)
       b%mode=BOX_ORTHO
    endif
    if(tag == "@BOX3")then
       call Box_ReadBOX3(b,file)
       b%mode=BOX_ORTHO
    endif
  end subroutine Box_Loader
  
  subroutine Box_BinaryLoader(b,file,tag)
    implicit none
    type(sBox) :: b
    integer file
    character(len=5) :: tag
    if(tag == "@BXLA")then
       call Box_ReadBinaryBXLA(b,file)
       b%mode=BOX_ORTHO
    endif
    if(tag == "@BOX3")then
       call Box_ReadBinaryBOX3(b,file)
       b%mode=BOX_ORTHO
    endif
  end subroutine Box_BinaryLoader
  
  subroutine Box_Save(b,file,writemode)
    implicit none
    type(sBox) :: b
    integer :: file, writemode
    if(b%mode == BOX_ORTHO)then
       if ( writemode .eq. BOX_BOX3 ) then
          call Box_WriteBOX3(b,file)
       else if ( writemode .eq. BOX_MDVW ) then
          call Box_WriteMDVW(b,file)
       endif
    endif
  end subroutine Box_Save
  
  subroutine Box_SaveBinary(b,file)
    implicit none
    type(sBox) :: b
    integer file
    if(b%mode == BOX_ORTHO)then
       call Box_WriteBinaryBOX3(b,file)
    endif
  end subroutine Box_SaveBinary
end module box_module

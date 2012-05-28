! -*- f90 -*-
!
!MEMUSAGE関係の処理を集約する。
!
module monitor_module
  implicit none
  integer, parameter :: MON_NOP=0, MON_MAXIMUM=1, MON_MINIMUM=2
  integer, parameter :: MON_WARN=1, MON_HALT=2

  type sMonitor_integer
     integer      :: threshold, value
     integer      :: type
     integer      :: action
     character(len=16) :: label
  end type sMonitor_integer

  interface new
     module procedure assign_new
  end interface
  interface update
     module procedure monitor_update
  end interface

contains
  
  !
  !呼ばれる度に新しいデバッグ変数を確保し、その番号を返す。
  !
  subroutine assign_new( monitor, type, threshold, action, label, value )
    type(sMonitor_integer), intent(out) :: monitor
    integer,        intent(in)  :: type
    integer,        intent(in)  :: threshold
    integer,        intent(in)  :: action
    character(len=16), intent(in) :: label
    integer,        intent(in), optional :: value
    monitor%type      = type
    monitor%threshold = threshold
    monitor%action    = action
    monitor%label     = label
    if ( present( value ) )then
       monitor%value = value
    else
       if( type .eq. MON_MAXIMUM )then
          monitor%value  = -2147483647
       else if( type .eq. MON_MINIMUM )then
          monitor%value  =  2147483647
       endif
    endif
  end subroutine assign_new
  !
  !変数値を監査する。
  !
  subroutine monitor_update( monitor, x )
    type(sMonitor_integer), intent(inout) :: monitor
    integer,        intent(in)    :: x
    integer :: over
    over= 0
    if ( monitor%type .eq. MON_MAXIMUM )then
       if ( monitor%value .lt. x ) then
          monitor%value = x
          if ( monitor%threshold .lt. x ) then
             over = 1
          endif
       endif
    else if ( monitor%type .eq. MON_MINIMUM )then
       if ( x .lt. monitor%value ) then
          monitor%value = x
          if ( x .lt. monitor%threshold ) then
             over = 1
          endif
       endif
    endif
    if ( over .ne. 0 ) then
       write( STDERR, * ) "Overflow:",x,"exceeds",monitor%threshold
       call showstatus( monitor )
       if ( monitor%action .eq. MON_HALT ) then
          stop
       endif
    endif
  end subroutine monitor_update

  subroutine showstatus( monitor )
    type(sMonitor_integer), intent(in) :: monitor
    write( STDERR, '(a16,1x,i8,"/",i8)') monitor%label, monitor%value, monitor%threshold
  end subroutine showstatus

end module monitor_module

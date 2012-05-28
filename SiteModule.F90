! -*- f90 -*-
!
!$BAj8_:nMQE@(B($B%5%$%H(B)$B$r0l3g$7$F07$&%b%8%e!<%k!#(B
!
!$B9dBNJ,;R$O!"J,;RA4BN$r0l$D$N9dBN$H$7$F07$$!"1?F0$O9dBN$N2sE>$HJB?J$H(B
!$B$7$F5-=R$9$k!#0lJ}!"J,;R4VAj8_:nMQ$O!"9dBNFb:BI8$K8GDj$5$l$?Aj8_:nMQ(B
!$BE@(B($B%5%$%H(B)$BF1;N$NAjBP:BI8$G7W;;$5$l$k!#(BGenesis6$B$G$O!"7O$K4^$^$l$kA4J,(B
!$B;R$NA4%5%$%H$r!"C10l$NG[Ns$H$7$F07$&$3$H$G!":G$b7W;;;~4V$,$+$+$kNO$N(B
!$B7W;;$N=hM}$r$G$-$k$@$18zN(2=$G$-$k$h$&$K@_7W$7$F$$$k!#(B
!
!$B=@FpJ,;R$N>l9g$O!"Aj8_:nMQE@0lE@0lE@$,<ANL$r;}$D$b$N$H9M$($k!#$3$N$?(B
!$B$a!"9dBNJ,;R$K4^$^$l$k%5%$%H$H=@FpJ,;R$N%5%$%H$G$O!"B0@-$,0[$J$k(B($B=@Fp(B
!$BJ,;R$N<AE@$O!"0l$D0l$D$,$$$o$PJ,;R$H$7$F$U$k$^$&$N$G!"<ANL$dB.EY!"2C(B
!$BB.EY$H$$$C$?B0@-$rJ];}$7$J$1$l$P$J$i$J$$!#9dBNJ,;R$N>l9g$O9dBNA4BN$N(B
!$B2sE>B.EY$d2CB.EY$H$$$C$?B0@-$O$"$k$,%5%$%H$K$O$=$N$h$&$JB0@-$OITMW(B)$B!#(B
!$B=@FpJ,;RFCM-$NB0@-$O(BFlexModule.F90$B$GIUM?$7$F$$$k!#9dBNJ,;R$NB0@-$O(B
!RigidModule.F90$B$G07$&!#(B
!
!$BAj8_:nMQ$r%+%C%H%*%U$9$k>l9g!"0lHLE*$K$OJ,;R$N=E?4$+$i$N5wN%$r4p=`$K(B
!$B$7$F!"Aj8_:nMQ$r%9%/%j!<%K%s%0$9$k!#$3$N>l9g!"=E?4$KBP$7$F$bNO$,F/$/(B
!$B$3$H$K$J$k!#(B($B&U(B(r)=s(r)$B&U(B0(r)$B$r(Br$B$GHyJ,$7$FNO$r5a$a$k(B)$B$=$3$G!"(BTIP4P$B$J(B
!$B$I$N9dBNJ,;R$r07$&>l9g$K$O!"=E?4%5%$%H$r2>A[E*$J%5%$%H$H$7$FDI2C$7!"(B
!$B=E?4$KF/$/NO$rB>$N%5%$%HF1MM$K;;=P$7!"9dBN1?F0$K4T85$5$;$F$$$k!#(B
!
module site_module
  use common_module
  implicit none

  type sSite
     sequence
     !$BA4%5%$%H$N8D?t(B
     integer :: nsite
     !$B8eH>$O9dBNJ,;R$N(Bsite$B$,F~$k!#(B
     integer :: nflexsite
     real(kind=8),pointer :: x(:), y(:), z(:)
     real(kind=8),pointer :: fx(:), fy(:), fz(:)
  end type sSite
  interface new
     module procedure Site_Constructor
  end interface
contains
  subroutine Site_constructor(si, n)
    type(sSite)        :: si
    integer,intent(IN) :: n
    if ( associated( si%x ) ) then
       deallocate( si%x )
       deallocate( si%y )
       deallocate( si%z )
       deallocate( si%fx )
       deallocate( si%fy )
       deallocate( si%fz )
    endif
    !write(STDERR,*) "SI ALLOC",n
    allocate(si%x(n))
    allocate(si%y(n))
    allocate(si%z(n))
    allocate(si%fx(n))
    allocate(si%fy(n))
    allocate(si%fz(n))
    call Site_initialize(si)
  end subroutine Site_constructor
  
  subroutine Site_initialize(si)
    type(sSite) :: si
    si%nsite=0
    si%nflexsite=0
  end subroutine Site_initialize
  
  subroutine Site_resetforce(si)
    type(sSite) :: si
    !write(STDERR,*) "SI",si%nsite
    si%fx(1:si%Nsite)=0d0
    si%fy(1:si%Nsite)=0d0
    si%fz(1:si%Nsite)=0d0
  end subroutine Site_resetforce
end module site_module
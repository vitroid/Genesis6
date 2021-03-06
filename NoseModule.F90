! -*- f90 -*-
! Nose's constant temperature method.
! Originally from Genesis2

!$BJ?@.#1#2G/#47n#6F|(B($BLZ(B)
!$B3HD%(BHamiltonian$B$O!"$A$c$s$H(Bs0$B$G7W;;$7$J$$$H@:EY$,=P$J$$!#$=$N$?$a!"8=(B
!$B:_$NBg;(GD$J7W;;J}K!$G$O8m:9$,Bg$-$$$h$&$K8+$($k$,!"3HD%(BHamiltonian$B$N(B
!$B@:EY$N0-$5$O7W;;@:EY<+BN$K$O(Bfeedback$B$5$l$J$$$N$G!"L5;k$7$F9=$o$J$$!#(B
!$B@:EY$,5$$K$J$k$J$i!"(BNOSEDBG$B$NItJ,$r;2>H$7$F!"(BNkTlog(s0)$B$G3HD%7O$N%]%F(B
!$B%s%7%c%k$r7W;;$9$k$Y$-!#(B

module nose_module
  implicit none
  type sNose
     logical :: active
     real(kind=8) :: q,dump,temp
     real(kind=8) :: zeta0,zeta1,zeta2,zeta3,zeta4,zetasum
#ifdef NOSEDBG
     real(kind=8) :: s0,s1,s2,s3,s4
#endif
  end type sNose
  interface new
     module procedure Nose_Constructor
  end interface
contains
  subroutine Nose_Constructor(nose)
    type(sNose),intent(inout) :: nose
    nose%active=.false.
    nose%zetasum=0d0
#ifdef NOSEDBG
    nose%s0=1d0
    nose%s1=0d0
    nose%s2=0d0
    nose%s3=0d0
    nose%s4=0d0
#endif
    return
  end subroutine Nose_Constructor
  
#ifdef NOSEDBG
  subroutine Nose_Correct(nose,temp,t)
#else
  subroutine Nose_Correct(nose,temp)
#endif
    use common_module
    use physconst_module
    real(kind=8) temp
    type(sNose),intent(inout) :: nose
    real(kind=8) dz2,zetaf
#ifdef NOSEDBG
    type(sTime),intent(in) ::t
    real(kind=8) :: sf,ds1
#endif
    zetaf = (temp-nose%temp)/nose%q-nose%dump*nose%zeta0
    dz2 = zetaf - nose%zeta1
    nose%zeta0 = nose%zeta0 + Gear4_0*dz2
    nose%zeta1 = zetaf
    nose%zeta2 = nose%zeta2 + Gear4_2*dz2
    nose%zeta3 = nose%zeta3 + Gear4_3*dz2
    nose%zeta4 = nose%zeta4 + Gear4_4*dz2
    nose%zetasum=nose%zetasum+nose%zeta0
#ifdef NOSEDBG
    sf = nose%s0*nose%zeta0*t%dt
    ds1 = sf-nose%s1
    nose%s0 = nose%s0 + Gear4_0*ds1
    nose%s1 = sf
    nose%s2 = nose%s2 + Gear4_2*ds1
    nose%s3 = nose%s3 + Gear4_3*ds1
    nose%s4 = nose%s4 + Gear4_4*ds1
#endif
    return
  end subroutine Nose_Correct

  subroutine Nose_Predict(n)
    type(sNose),intent(inout) :: n
    n%zeta0  = n%zeta0 + n%zeta1 + n%zeta2 + n%zeta3 + n%zeta4
    n%zeta1  = n%zeta1 + 2d0*n%zeta2 + 3d0*n%zeta3 + 4d0*n%zeta4
    n%zeta2  = n%zeta2 + 3d0*n%zeta3 + 6d0*n%zeta4
    n%zeta3  = n%zeta3 + 4d0*n%zeta4
#ifdef NOSEDBG
    n%s0  = n%s0 + n%s1 + n%s2 + n%s3 + n%s4
    n%s1  = n%s1 + 2d0*n%s2 + 3d0*n%s3 + 4d0*n%s4
    n%s2  = n%s2 + 3d0*n%s3 + 6d0*n%s4
    n%s3  = n%s3 + 4d0*n%s4
#endif
    return
  end subroutine Nose_Predict
  
!$BJ?@.#1#2G/#47n#6F|(B($BLZ(B)NOSE$B$O(BNPT2$B$H$N8_497A<0!#:#8e$O(BNOS2$B7A<0$r;HMQ$9$k!#(B
  function Nose_ReadNOSE(nose,file)
    type(sNose),intent(inout) :: nose
    integer file
    logical :: Nose_ReadNOSE
    read(file,*) nose%temp,nose%q,nose%dump,nose%zeta0,nose%zeta1&
         & ,nose%zeta2,nose%zeta3,nose%zeta4
#ifdef VERBOSE
    write(STDERR,*) "@NOSE",nose%temp,nose%q,nose%dump,nose%zeta0&
         & ,nose%zeta1,nose%zeta2,nose%zeta3,nose%zeta4
#endif
    Nose_ReadNOSE=(nose%temp >= 0d0)
    nose%q=nose%q * 0021.8852721616497249235d0**2
    return
  end function Nose_ReadNOSE
  
  function Nose_ReadBinaryNOSE(nose,file)
    type(sNose),intent(inout) :: nose
    integer file
    logical :: Nose_ReadBinaryNOSE
    read(file) nose%temp,nose%q,nose%dump,nose%zeta0,nose%zeta1,nose&
         & %zeta2,nose%zeta3,nose%zeta4
#ifdef VERBOSE
    write(STDERR,*) "@NOSE",nose%temp,nose%q,nose%dump,nose%zeta0&
         & ,nose%zeta1,nose%zeta2,nose%zeta3,nose%zeta4
#endif
    Nose_ReadBinaryNOSE=(nose%temp >= 0d0)
    nose%q=nose%q * 0021.8852721616497249235d0**2
    return
  end function Nose_ReadBinaryNOSE
  
  function Nose_ReadNOS2(nose,file)
    type(sNose),intent(inout) :: nose
    integer file
    logical :: Nose_ReadNOS2
    read(file,*) nose%temp,nose%q,nose%dump,nose%zeta0,nose%zeta1&
         & ,nose%zeta2,nose%zeta3,nose%zeta4,nose%zetasum
    Nose_ReadNOS2=(nose%temp >= 0d0)
#ifdef VERBOSE
    write(STDERR,*) "@NOS2",nose%temp,nose%q,nose%dump,nose%zeta0&
         & ,nose%zeta1,nose%zeta2,nose%zeta3,nose%zeta4,nose%zetasum
#endif
    return
  end function Nose_ReadNOS2
  
  function Nose_ReadBinaryNOS2(nose,file)
    type(sNose),intent(inout) :: nose
    integer file
    logical :: Nose_ReadBinaryNOS2
    read(file) nose%temp,nose%q,nose%dump,nose%zeta0,nose%zeta1,nose&
         & %zeta2,nose%zeta3,nose%zeta4,nose%zetasum
    Nose_ReadBinaryNOS2=(nose%temp >= 0d0)
#ifdef VERBOSE
    write(STDERR,*) "@NOS2",nose%temp,nose%q,nose%dump,nose%zeta0&
         & ,nose%zeta1,nose%zeta2,nose%zeta3,nose%zeta4,nose%zetasum
#endif
    return
  end function Nose_ReadBinaryNOS2
  
  subroutine Nose_WriteBinaryNOS2(nose,file)
    type(sNose),intent(in) :: nose
    integer file
    write(file) "@NOS2"
    write(file) nose%temp,nose%q,nose%dump,nose%zeta0,nose%zeta1,nose&
         & %zeta2,nose%zeta3,nose%zeta4,nose%zetasum
    return
  end subroutine Nose_WriteBinaryNOS2
  
  subroutine Nose_SaveBinary(nose,file)
    type(sNose),intent(in) :: nose
    integer file
    if(nose%active)then
       call Nose_WriteBinaryNOS2(nose,file)
    endif
  end subroutine Nose_SaveBinary
  
  subroutine Nose_WriteNOS2(nose,file)
    type(sNose),intent(in) :: nose
    integer file
    write(file,1)
    write(file,*) nose%temp,nose%q,nose%dump,nose%zeta0,nose%zeta1&
         & ,nose%zeta2,nose%zeta3,nose%zeta4,nose%zetasum
1   format( "@NOS2" )
    return
  end subroutine Nose_WriteNOS2
  
  subroutine Nose_Loader(nose,file,tag)
    type(sNose),intent(inout) :: nose
    integer,intent(IN) :: file
    character(len=5),intent(IN) :: tag
    if(tag == "@NOSE")then
       nose%active=Nose_ReadNOSE(nose,file)
    endif
    if(tag == "@NOS2")then
       nose%active=Nose_ReadNOS2(nose,file)
    endif
  end subroutine Nose_Loader
  
  subroutine Nose_BinaryLoader(nose,file,tag)
    type(sNose),intent(inout) :: nose
    integer,intent(IN) :: file
    character(len=5),intent(IN) :: tag
    if(tag == "@NOSE")then
       nose%active=Nose_ReadBinaryNOSE(nose,file)
    endif
    if(tag == "@NOS2")then
       nose%active=Nose_ReadBinaryNOS2(nose,file)
    endif
  end subroutine Nose_BinaryLoader
  
  subroutine Nose_Version
    write(STDERR,*) "$Id: NoseModule.F90,v 1.2 2002/12/09 08:55:52 matto Exp $"
    return
  end subroutine Nose_Version
end module nose_module

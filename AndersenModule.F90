! -*- f90 -*-
! Andersen's constant-pressure method.

!$B;~4V9o$_$rJQ99$7$?>l9g$H!"%7%9%F%`%5%$%:$rJQ99$7$?>l9g$K!"%T%9%H%s$N?6(B
!$BF0?t$,JQ2=$7$J$$$h$&$JC10L7O$K$7$J$1$l$P$J$i$J$$!#(B
!$B8e<T$O(B64$B8D$N7O$r(B8$B8D$/$C$D$1$?>l9g$GHf3S$9$l$P$h$$!#(B
!$BJ]B8NL$,@5$7$/J]B8$5$l$F$$$k$+$b3NG'$9$k$3$H!#(B

!$BJ?@.#1#2G/#67n#7F|(B($B?e(B)$B;~4V9o$_$K4X$7$F$O3NG'$7$?!#(B
!system size$B$K4X$7$F$b40A4$G$O$J$$$K$;$h0lCW$9$k$h$&$K$J$C$?!#(B
!NPT2$B$G$OJ,;R$"$?$j$NBN@Q$G7W;;$7$?$,!"$3$3$G$OA4BN@Q$G7W;;$9$k!#(B
!$B$"$$$+$o$i$:!"J]B8NL$OJ]B8$7$J$$!#7W;;J}K!$,8m$C$F$$$k$h$&$@$,!"2?$,(B
!$BLdBj$+$h$/$o$+$i$J$$!#(B

!dump$B$K4X$7$F$O$^$@C10L$rD4@0$7$F$$$J$$!#(B
!$B=PNO%G!<%?$N7QB3@-$N%A%'%C%/!'J?@.#1#2G/#67n#8F|(B($BLZ(B)$B05NO$N7QB3@-$OLd(B
!$BBj$J$$$,!"1?F0NL$,0lHV2<$N(B2$B7e$,HyL/$K68$&!#$^$"!"$h$7$H$9$k$+!#(B

!$B05NO$NC10L$O(BPCN2/PCN3$B7A<0$G$O(Batm$B$G;XDj$9$k$,!"FbIt$N7W;;$O(Bpascal$B$KE}0l$9$Y$-$@$H;W$&!#J?@.(B15$BG/(B11$B7n(B26$BF|(B($B?e(B)
!
!PCN4$B$G$O!"(Bmass$B$NC10L$rJQ99$7$?(B(PCN2$B$NC10L$O$*$+$7$$(B)$B!#J?@.(B16$BG/(B8$B7n(B30$BF|(B($B7n(B)
!PCN2$B$N(Bmass$B$r(BPCN4$B$N(Bmass$B$KFI$_$+$($k$K$O!"(Bdof**2$B$G3d$k!#(B
!

module andersen_module
  use common_module
  use box_module

  implicit none

  type sAndersen
     real(kind=8) :: pex,mass,v0,v1,v2,v3,v4,v5  !,v00
     type(sBox)   :: b
     integer      :: mode
     logical      :: compat_pcn23
  end type sAndersen
  integer,parameter :: noandersen=0,orthorhombic=1,orz=2

  interface new
     module procedure andersen_constructor
  end interface

contains
  subroutine Andersen_Constructor(a)
    type(sAndersen),intent(inout) :: a
    A%v0=0d0
    A%v1=0d0
    A%v2=0d0
    A%v3=0d0
    A%v4=0d0
    A%v5=0d0
    a%mode=noandersen
  end subroutine Andersen_Constructor
  
  subroutine Andersen_Initialize(a,b)
    use box_module
    type(sAndersen),intent(inout) :: a
    type(sBox),intent(inout)      :: b
    !Copy initial volume info
    a%b   = b
    !a%v00 = box_volume( b )
    !a%v0  = a%v00
    a%v0 = box_volume( b )
  end subroutine Andersen_Initialize
  
  subroutine Andersen_Rescale(a,box)
    use vector_module
    use box_module
    type(sAndersen),intent(in) :: a
    type(sBox),intent(inout)   :: box
    type(vector3)              :: ratio
    real(kind=8)               :: r2
    if(a%mode == orthorhombic)then
       r2 = (a%v0/box_volume(box))**(1d0/3d0)
       call scale1( box, r2 )
    endif
    if(a%mode == orz)then
       ratio%vec(1) = 1d0
       ratio%vec(2) = 1d0
       ratio%vec(3) = a%v0/box_volume(box)
       call scalev( box, ratio )
    endif
  end subroutine Andersen_Rescale
  
  subroutine Andersen_Predict(a,b)
    type(sAndersen),intent(inout) :: a
    type(sBox),intent(inout) :: b
    a%v0  = a%v0 + a%v1 + a%v2 + a%v3 + a%v4 + a%v5
    a%v1  = a%v1 + 2d0*a%v2 + 3d0*a%v3 + 4d0*a%v4 + 5d0*a%v5
    a%v2  = a%v2 + 3d0*a%v3 + 6d0*a%v4 + 10d0*a%v5
    a%v3  = a%v3 + 4d0*a%v4 + 10d0*a%v5
    a%v4  = a%v4 + 5d0*a%v5
    call Andersen_Rescale(a,b)
  end subroutine Andersen_Predict
  
  subroutine Andersen_Correct(a,b,pin)
    use physconst_module
    type(sAndersen),intent(inout) :: a
    type(sBox),intent(inout) :: b
    !$BJ?@.(B15$BG/(B11$B7n(B26$BF|(B($B?e(B)$B8=:_!"05NO$O(Batm.$BC10L$NL7=b$O$9$Y$F(Bmass$B$K=8Ls$5$l$F$$$k!#(B
    real(kind=8),intent(in) :: pin
    real(kind=8) :: cvol,volp2
    volp2 = (pin-a%pex)/a%mass*0.5d0
    ! volp2: A^3 ; a%mass: atm/A^3 
    !write(STDERR,*) pin, a%pex, a%mass
    !write(STDERR,*) volp2, a%v2
    !stop
    cvol  = volp2 - a%v2
    a%v0    = a%v0  + Gear5_0*cvol
    a%v1    = a%v1  + Gear5_1*cvol
    a%v2    = volp2
    a%v3    = a%v3  + Gear5_3*cvol
    a%v4    = a%v4  + Gear5_4*cvol
    a%v5    = a%v5  + Gear5_5*cvol
    call Andersen_Rescale(a,b)
  end subroutine Andersen_Correct
  
  function Andersen_ReadPCN2(a,file)
    use error_module
    type(sAndersen),intent(inout) :: a
    integer file
    logical :: Andersen_ReadPCN2
    real(kind=8) :: dump
    call die( 0, "@PCN2 is not supported any more. Please use @PCN4 instead." )
    read(file,*) a%mode,a%pex,a%mass,dump
    Andersen_ReadPCN2=(a%mode.ne.noandersen)
#ifdef VERBOSE
    write(STDERR,*) "@PCN2",a%pex,a%mass,dump,a%v0,a%v1,a%v2,a%v3,a&
         & %v4,a%v5
#endif
    return
  end function Andersen_ReadPCN2
  
  function Andersen_ReadPCN4(a,file)
    use error_module
    type(sAndersen),intent(inout) :: a
    integer file
    logical :: Andersen_ReadPCN4
    read(file,*) a%mode,a%pex,a%mass
    A%v0=0d0
    A%v1=0d0
    A%v2=0d0
    A%v3=0d0
    A%v4=0d0
    A%v5=0d0
    a%compat_pcn23=.false.
    Andersen_ReadPCN4=(a%mode.ne.noandersen)
  end function Andersen_ReadPCN4
  
  function Andersen_ReadBinaryPCN3(a,file)
    type(sAndersen),intent(inout) :: a
    integer :: file
    logical :: Andersen_ReadBInaryPCN3
    real(kind=8) :: dump
    read(file) a%mode,a%pex,a%mass,dump, a%v0,a%v1,a%v2,a%v3,a%v4,a&
         & %v5
    Andersen_ReadBinaryPCN3=(a%mode.ne.noandersen)
    a%compat_pcn23=.true.
#ifdef VERBOSE
    write(STDERR,*) "@PCN3",a%pex,a%mass,dump,a%v0,a%v1,a%v2,a%v3,a&
         & %v4,a%v5
#endif
    return
  end function Andersen_ReadBinaryPCN3
  
  function Andersen_ReadBinaryPCN5(a,file)
    type(sAndersen),intent(inout) :: a
    integer :: file
    logical :: Andersen_ReadBInaryPCN5
    read(file) a%mode,a%pex,a%mass,a%v0,a%v1,a%v2,a%v3,a%v4,a&
         & %v5
    Andersen_ReadBinaryPCN5=(a%mode.ne.noandersen)
    a%compat_pcn23=.false.
#ifdef VERBOSE
    write(STDERR,*) "@PCN5",a%pex,a%mass,a%v0,a%v1,a%v2,a%v3,a&
         & %v4,a%v5
#endif
    return
  end function Andersen_ReadBinaryPCN5
  
  subroutine Andersen_WriteBinaryPCN5(a,file)
    type(sAndersen),intent(in) :: a
    integer file
    write(file) "@PCN5"
    write(file) a%mode,a%pex,a%mass,a%v0,a%v1,a%v2,a%v3,a%v4,a%v5
    return
  end subroutine Andersen_WriteBinaryPCN5
  
  subroutine Andersen_Loader(a,file,tag)
    type(sAndersen),intent(inout) :: a
    integer,intent(in) :: file
    character(len=5),intent(IN) :: tag
    logical :: result
    if(tag == "@PCN2")then
       result = Andersen_ReadPCN2(a,file)
    elseif(tag == "@PCN4")then
       result = Andersen_ReadPCN4(a,file)
    endif
  end subroutine Andersen_Loader
  
  subroutine Andersen_BinaryLoader(a,file,tag)
    type(sAndersen),intent(inout) :: a
    integer,intent(in) :: file
    character(len=5),intent(IN) :: tag
    logical :: result
    if(tag == "@PCN3")then
       result = Andersen_ReadBinaryPCN3(a,file)
    elseif(tag == "@PCN5")then
       result = Andersen_ReadBinaryPCN5(a,file)
    endif
  end subroutine Andersen_BinaryLoader
end module andersen_module

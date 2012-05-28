! -*- f90 -*-
! Andersen's constant-pressure method.

!時間刻みを変更した場合と、システムサイズを変更した場合に、ピストンの振
!動数が変化しないような単位系にしなければならない。
!後者は64個の系を8個くっつけた場合で比較すればよい。
!保存量が正しく保存されているかも確認すること。

!平成１２年６月７日(水)時間刻みに関しては確認した。
!system sizeに関しても完全ではないにせよ一致するようになった。
!NPT2では分子あたりの体積で計算したが、ここでは全体積で計算する。
!あいかわらず、保存量は保存しない。計算方法が誤っているようだが、何が
!問題かよくわからない。

!dumpに関してはまだ単位を調整していない。
!出力データの継続性のチェック：平成１２年６月８日(木)圧力の継続性は問
!題ないが、運動量が一番下の2桁が微妙に狂う。まあ、よしとするか。

!圧力の単位はPCN2/PCN3形式ではatmで指定するが、内部の計算はpascalに統一すべきだと思う。平成15年11月26日(水)
!
!PCN4では、massの単位を変更した(PCN2の単位はおかしい)。平成16年8月30日(月)
!PCN2のmassをPCN4のmassに読みかえるには、dof**2で割る。
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
    !平成15年11月26日(水)現在、圧力はatm.単位の矛盾はすべてmassに集約されている。
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

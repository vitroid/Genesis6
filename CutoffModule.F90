! -*- f90 -*-
module cutoff_module
  implicit none
  type sCutoff
     real(kind=8) :: in,out,c1,c2
     integer :: mode
  end type sCutoff
  integer, parameter :: CUTOFF_NONE=0,CUTOFF_POLY=1
  interface new
     module procedure cutoff_constructor
  end interface
  interface load
     module procedure cutoff_loader
  end interface
  interface save
     module procedure cutoff_save
  end interface
  interface load_binary
     module procedure cutoff_binaryloader
  end interface
  interface save_binary
     module procedure cutoff_savebinary
  end interface
contains
  subroutine cutoff_initialize(c,in,out)
    type(sCutoff),intent(out) :: c
    real(kind=8),intent(in),optional :: in,out
    if ( present( in ) ) then
       c%in = in
       c%out = out
       c%c1 = (1d0/(In-Out)**5)
       c%c2 = (c%C1*30d0)
    else
       c%mode = CUTOFF_NONE
    endif
  end subroutine cutoff_initialize
  
  subroutine cutoff_func(c,dd,hh,dh)
    type(sCutoff),intent(in) :: c
    real(kind=8),intent(in) :: dd
    real(kind=8),intent(out) :: hh,dh
    real(kind=8) :: u,l,uu,ll,dr
    hh = 0d0
    dh = 0d0
    if(dd.lt.c%Out**2)then
       if(dd.gt.c%In**2)then
          dr=dsqrt(dd)
          u=dr-c%Out
          l=dr-c%In
          uu=u**2
          ll=l**2
          hh=c%C1*uu*u*(10.0d0*ll-5.0d0*u*l+uu)
          dh =-c%C2*uu*ll/dr
       else
          hh=1d0
          dh=0d0
       endif
    endif
  end subroutine cutoff_func
  
  function Cutoff_ReadBinaryRCOA(co,file)
    type(sCutoff),intent(OUT) :: co
    integer,intent(in) :: file
    real(kind=8) :: rc,rcmargin
    logical :: Cutoff_ReadBinaryRCOA
    READ(file) rc,rcmargin
    if(rc.gt.0)then
       call cutoff_initialize(co,rc-rcmargin,rc)
    endif
    Cutoff_ReadBinaryRCOA=(rc.gt.0)
  end function Cutoff_ReadBinaryRCOA
  
  function Cutoff_ReadRCOA(co,file)
    type(sCutoff),intent(OUT) :: co
    integer,intent(in) :: file
    real(kind=8) :: rc,rcmargin
    logical :: Cutoff_ReadRCOA
    READ(file,*) rc,rcmargin
    if(rc.gt.0)then
       call cutoff_initialize(co,rc-rcmargin,rc)
    endif
    Cutoff_ReadRCOA=(rc.gt.0)
  end function Cutoff_ReadRCOA
  
  subroutine Cutoff_WriteBinaryRCOA(co,file)
    type(sCutoff),intent(in) :: co
    integer,intent(in) :: file
    write(file) "@RCOA"
    write(file) co%out,co%out-co%in
  end subroutine Cutoff_WriteBinaryRCOA
  
  subroutine Cutoff_WriteRCOA(co,file)
    type(sCutoff),intent(in) :: co
    integer,intent(in) :: file
    write(file,1)
    write(file,*) co%out,co%out-co%in
1   format( "@RCOA")
  end subroutine Cutoff_WriteRCOA
  
  subroutine CutOff_Loader(co,file,tag)
    type(sCutoff),intent(inout) :: co
    integer,intent(in) :: file
    character(len=5), intent(in) :: tag
    logical ret
    if(tag == "@RCOA")then
       co%mode = CUTOFF_NONE
       ret=Cutoff_ReadRCOA(co,file)
       if(ret)co%mode=CUTOFF_POLY
    endif
  end subroutine CutOff_Loader
  
  subroutine CutOff_BinaryLoader(co,file,tag)
    type(sCutoff),intent(inout) :: co
    integer,intent(in) :: file
    character(len=5), intent(in) :: tag
    if(tag == "@RCOA")then
       co%mode = CUTOFF_NONE
       if(Cutoff_ReadBinaryRCOA(co,file)) co%mode=CUTOFF_POLY
    endif
  end subroutine CutOff_BinaryLoader
  
  subroutine CutOff_Constructor(co)
    type(sCutoff),intent(INOUT) :: co
    co%mode=CUTOFF_NONE
  end subroutine CutOff_Constructor
  
  subroutine CutOff_Save(co,file)
    type(sCutoff),intent(in) :: co
    integer,intent(in) :: file
    if(co%mode == CUTOFF_POLY)then
       call CutOff_WriteRCOA(co,file)
    endif
  end subroutine CutOff_Save
  
  subroutine CutOff_SaveBinary(co,file)
    type(sCutoff),intent(in) :: co
    integer,intent(in) :: file
    if(co%mode == CUTOFF_POLY)then
       call CutOff_WriteBinaryRCOA(co,file)
    endif
  end subroutine CutOff_SaveBinary
  
  function CutOff_OuterLimit(co)
    type(sCutoff),intent(in) :: co
    real(kind=8) :: Cutoff_OuterLimit
    if(co%mode == CUTOFF_POLY)then
       CutOff_OuterLimit = co%out
    else
       CutOff_OuterLimit = -1d0
    endif
  end function CutOff_OuterLimit

  !
  !Make interaction dimmer array after compression
  !
  subroutine CutOff_Dimmer_obsolete(c, iv, site, mi, mj )
    use interaction_module
    use site_module
    use mol_module
    implicit none
    type(sInteraction),intent(inout) :: iv
    type(sCutOff)     ,intent(in)    :: c
    type(sSite)       ,intent(in)    :: site
    type(sMol)        ,intent(in)    :: mi,mj

    integer                          :: i,j,k,ii,jj
    real(kind=8)                     :: dx00,dy00,dz00
    real(kind=8)                     :: hh,dh
#ifdef VPOPTIMIZE
    real(kind=8)                     :: dr,dd
    real(kind=8)                     :: u,l,uu,ll
#endif

    if ( c%mode .eq. CUTOFF_POLY ) then
       do k=1,iv%npair
          i = iv%pair_i(k)
          j = iv%pair_j(k)
          !
          !分子を構成するサイトのうち最後のサイト
          !
          ii = i * mi%nsite+mi%offset
          jj = j * mj%nsite+mj%offset
          dx00 = site%x(ii) - site%x(jj) - iv%ox(k)
          dy00 = site%y(ii) - site%y(jj) - iv%oy(k)
          dz00 = site%z(ii) - site%z(jj) - iv%oz(k)
#ifdef VPOPTIMIZE
          !インライン展開されているはずなのに、実際には遅くなるので展開
          !して記述しておく。なぜだろう。
          hh = 0d0
          dh = 0d0
          dd=dx00**2+dy00**2+dz00**2
          if( dd .lt. c%Out**2 )then
             if( dd .gt. c%In**2 )then
                dr    = dsqrt( dd )
                u  = dr-c%Out
                l  = dr-c%In
                uu = u**2
                ll = l**2
                hh    = c%C1 * uu * u * (10.0d0 * ll - 5.0d0*u*l + uu)
                dh    = -c%C2 * uu * ll / dr
             else
                hh    = 1d0
                dh    = 0d0
             endif
          endif
#else
          call cutoff_func(c,dx00**2+dy00**2+dz00**2,hh,dh)
#endif
          iv%eratio( k ) = hh
          iv%fratio( k ) = dh
       enddo
    endif
  end subroutine CutOff_Dimmer_obsolete

  !
  !Make interaction dimmer array after compression
  !作業内容はCutOff_Dimmer()と同じだが、site表ではなく、重心位置を直接
  !引数として受けとる。
  !
  subroutine CutOff_Dimmer2(c, iv, com0, com1 )
    use interaction_module
    use mol_module
    use vector_module
    implicit none
    type(sInteraction),intent(inout) :: iv
    type(sCutOff)     ,intent(in)    :: c
    type(vector3),     intent(in)    :: com0(*), com1(*)

    integer                          :: i,j,k
    real(kind=8)                     :: ox,oy,oz
    real(kind=8)                     :: dx00,dy00,dz00
    real(kind=8)                     :: hh,dh
#ifdef VPOPTIMIZE
    real(kind=8)                     :: dr,dd
    real(kind=8)                     :: u,l,uu,ll
#endif
    if ( c%mode .eq. CUTOFF_POLY ) then
       do k=1,iv%npair
          i = iv%pair_i(k)
          j = iv%pair_j(k)
          dx00 = com0(i)%vec(1) - com1(j)%vec(1)
          dy00 = com0(i)%vec(2) - com1(j)%vec(2)
          dz00 = com0(i)%vec(3) - com1(j)%vec(3)
          ox   = iv%ox(k)
          oy   = iv%oy(k)
          oz   = iv%oz(k)
          dx00 = dx00-ox
          dy00 = dy00-oy
          dz00 = dz00-oz
#ifdef VPOPTIMIZE
          !インライン展開されているはずなのに、実際には遅くなるので展開
          !して記述しておく。なぜだろう。
          hh = 0d0
          dh = 0d0
          dd=dx00**2+dy00**2+dz00**2
          if( dd .lt. c%Out**2 )then
             if( dd .gt. c%In**2 )then
                dr    = dsqrt( dd )
                u  = dr-c%Out
                l  = dr-c%In
                uu = u**2
                ll = l**2
                hh    = c%C1*uu*u*(10.0d0*ll-5.0d0*u*l&
                     & +uu)
                dh    = -c%C2*uu*ll/dr
             else
                hh    = 1d0
                dh    = 0d0
             endif
          endif
#else
          call cutoff_func(c,dx00**2+dy00**2+dz00**2,hh,dh)
#endif
          iv%eratio( k ) = hh
          iv%fratio( k ) = dh
       enddo
    endif
  end subroutine CutOff_Dimmer2

  !
  !カットオフによる力は、カットオフ基準点に対して加わる。分子の属性と
  !して、この力をforceの計算時に求めて保管しておく。
  !大幅な変更になる可能性があるので、20031003時点のスナップショットを
  !cvsに追加した。
  

end module cutoff_module

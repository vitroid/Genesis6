! -*- f90 -*-
module mol_module
  !
  !剛体分子、柔軟分子共通の情報
  !
  implicit none

  integer, parameter :: RIGID_MODE=1, FLEX_MODE=0
  type sMol
     sequence
     !
     !1分子の自由度
     !
     integer :: dof
     !
     !グループに属する分子の数
     !
     integer :: nmol
     !
     !1分子あたりのサイト数
     !
     integer :: nsite
     !
     !offsetは最初の分子群は0、
     !次の分子群は最初の群のサイト数、という風に番号を与える。
     !
     integer :: offset
#ifdef VPOPTIMIZE
     integer,dimension(MAXSITE*MAXMOL) :: lvsite,lvmol
#endif
     !
     !識別子
     !
     character(len=256) id
     !
     !IsRigid
     !
     logical :: isRigid
     !
     !Fixed (immobile)
     !
     logical :: isFixed
     !
     !サイトごとの名前
     !
     character(len=8),pointer :: name(:) => null()
  end type sMol
  logical, parameter :: FIXED = .true.

  interface new
     module procedure mol_constructor
  end interface
contains
  subroutine Mol_Constructor(mol,site,dof,isRigid, isFixed)
    implicit none
    type(sMol),intent(INOUT) :: mol
    integer,intent(IN)       :: site,dof
    integer, intent(IN)      :: isRigid
    logical, intent(IN),optional :: isFixed
    mol%dof     = dof
    mol%nsite   = site
    mol%nmol    = 0
    mol%offset  = -1
    mol%isRigid = isRigid .eq. RIGID_MODE
    if ( present( isFixed ) ) then
       mol%isFixed = isFixed
    else
       mol%isFixed = .false.
    endif
    !allocate(mol%name(site))
  end subroutine Mol_Constructor
  
  !     サイト数を増加させるのみ。
  subroutine Mol_allocate(m,si,nmol,isRigid)
    use site_module
    implicit none
    type(sSite) :: si
    type(sMol) :: m
    integer,intent(IN) :: nmol
    integer,intent(IN) :: isRigid
    integer :: i
#ifdef VPOPTIMIZE
    integer :: k,s
#endif
    m%nmol=nmol
    if(m%nmol.gt.MAXmol)then
       write(STDERR,*) "ERROR: TOO LARGE SYSTEM"
       stop
    endif
    if(m%offset.ge.0)then
       write(STDERR,*) "Memory for the molecules is already&
            & allocated."
       stop
    endif
    if( RIGID_MODE .ne. IsRigid .and. (si%nsite.ne.si%nflexsite) )then
       write(STDERR,*) "Rigid molecules must be allocated later."
       stop
    endif
    m%offset = si%nsite
    si%nsite = si%nsite + (m%nsite*m%nmol)
    if( FLEX_MODE .eq. isRigid )si%nflexsite = si%nsite
#ifdef VPOPTIMIZE
    k=0
    do i=1,m%nmol
       do s=1,m%nsite
          k=k+1
          m%lvsite(k)=s
          m%lvmol(k)=i
       enddo
    enddo
#endif
    return
  end subroutine Mol_allocate
  
  function Mol_DoF(f)
    implicit none
    type(sMol),intent(IN) :: f
    integer :: Mol_DoF
    Mol_DoF=f%nmol*f%dof
    return
  end function Mol_DoF
  
  subroutine Mol_SetAtomName( m, name )
    type(sMol),intent(inout)     :: m
    character(len=8),intent(in)   :: name(*)
    !
    !local
    !
    integer :: site
    if ( associated( m%name ) ) then
       deallocate( m%name )
    endif
    allocate(m%name(1:m%nsite))
    do site=1, m%nsite
       m%name(site) = name(site)
    enddo
  end subroutine Mol_SetAtomName

  function Mol_CountSiteMDVW(mol)
    type(sMol),intent(IN) :: mol
    integer :: Mol_CountSiteMDVW
    integer :: n,site
    n=0
    !
    !原子名が空白から始まるサイトは表示しない。
    !
    do site=1,mol%nsite-1
       if ( mol%name(site)(1:1) .ne. " " ) n=n+1
    enddo
    Mol_CountSiteMDVW = n * mol%nmol
  end function Mol_CountSiteMDVW
    
  subroutine Mol_Version
    write(STDERR,*) "$Id: Mol.F90,v 1.2 2000/05/02 03:02:15 matto Exp&
         & $"
    return
  end subroutine Mol_Version
end module mol_module

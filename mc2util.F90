! -*- f90 -*-
#define OP
#ifdef OP
#define USE_FAKEMPI
#endif

module mc2util
  use common_module
  use distrib_module
  use mol_module
  use rigid_module
  use site_module
  use book_module
  use cutoff_module
  use grid_module
  use box_module
  use oplsmeoh_module
  use tip4p_module
  use tip5p_module
  use nvde_module
  use standard_interaction_module
  use umbrella
  use triangle
  use graph2
  use random_module
#ifdef OP
  use q6local2
#endif
  implicit none

  integer, parameter :: MAXUMBRELLA=10
  integer, parameter :: OP_TRIANGLE=1, OP_Q6T=2, OP_EP=3

  type sSystem
     !integer :: seed
     type(sHistogram), pointer :: refhist
     !
     !interval between saving snapshots
     !
     integer          :: snpo,nlog
     integer          :: nloop
     type(sMol) :: mol( MAXCOMPO )
     type(sSite) :: site
     type(sRigid) :: rigid( MAXCOMPO )
     type(sStdInt) :: stdint(MAXCOMPO)
     integer :: dof
     type(sBook) :: book
     type(sGrid) :: grid
     type(sBox)    :: box
     type(sCutoff) :: cutoff
     !
     !@NCMPで明示的に指定される成分数
     !
     integer :: incompo
     !
     !成分の数(とりあえず上限は2)
     !
     integer      :: ncompo
     !
     !beta in 1/K
     !
     real(kind=8) :: beta, temp
     !
     !node毎の設定ファイル
     !
     character(len=1024) :: infile, outfile
     integer :: nprocs, myrank
#ifdef OP
     !
     !several order parameters
     !
     integer :: numbrella
     integer :: optype(MAXUMBRELLA)
     type(sTrianglePtr) :: tri(MAXUMBRELLA)
     type(sQ6SystemPtr) :: q6s(MAXUMBRELLA)
     type(sQ6DiffPtr) :: q6d(MAXUMBRELLA)
     type(umbrella_parabolic) :: umbrella(MAXUMBRELLA)
     character(len=8) :: id08(MAXUMBRELLA)
     type(pGraph), pointer    :: rxgraph
     type(pPair), pointer     :: rxpair
     type(pRandom), pointer   :: rand
     integer                  :: nexchange
     integer                  :: ninnerloop
     real(kind=8), pointer    :: delta(:), rtratio(:)
#endif
  end type sSystem

  interface new
     module procedure sSystem_new
  end interface

contains

  subroutine sSystem_new( sys )
    type(sSystem), intent(INOUT) :: sys
#ifdef OP
    sys%numbrella = 0
    nullify( sys%rxgraph )
    nullify( sys%rand )
#endif
    call new(sys%site, MAXsite_total)
    !call new(an)
    !call new(av)
    call new(sys%book)
    call new(sys%grid)
    call new(sys%box)
    call new(sys%cutoff)
    
    sys%beta = 1d0/400d0
    sys%dof=0d0
    sys%ncompo=0
    sys%incompo=0
    sys%snpo = 0
    !sys%seed = 5678
    sys%nexchange = 30
    sys%ninnerloop = 100
  end subroutine sSystem_new

  subroutine ReadParams( sys, file )
    integer, intent(IN) :: file
    type(sSystem), intent(INOUT) :: sys


    character(len=5) :: tag
    integer :: i,j,nmol
    character(len=8) :: lastid

    integer :: nprocs, ierr, seed
    lastid="        "
    do
       read(file,'(a5)',END=999) tag
       write(STDERR,*) tag
       if(tag.eq."@DELT")then
          read(file,*) nmol
          allocate( sys%delta(nmol) )
          allocate( sys%rtratio(nmol) )
          do i=1,nmol
             read( file,* ) sys%delta(i), sys%rtratio(i)
          enddo
          cycle
       endif
       
       if(tag.eq."@SEED")then
          read(file,*) seed
          call new( sys%rand, seed )
          cycle
       endif
       if(tag.eq."@MTRN")then
          call random_loadMTRN( sys%rand, file )
          cycle
       endif
       if(tag.eq."@INNR")then
          read(file,*) sys%ninnerloop
          cycle
       endif
       if(tag.eq."@HIST")then
          call histogram_load( sys%refhist, file )
          cycle
       endif
       !
       !Triangle OP
       !
       if(tag.eq."@REPQ")then
          sys%numbrella = sys%numbrella + 1
          sys%id08( sys%numbrella ) = "TRIANGLE"
          call new( sys%tri( sys%numbrella )%p, .false. )
          sys%optype( sys%numbrella ) = OP_TRIANGLE
          call umbrella_ReadParam( sys%umbrella( sys%numbrella ), file )
          cycle
       endif
       !
       !General Umbrella
       !
       if(tag.eq."@GEUM")then
          sys%numbrella = sys%numbrella + 1
          read( file, * ) sys%id08( sys%numbrella )
          if ( sys%id08( sys%numbrella ) == "TRIANGLE" ) then
             call new( sys%tri( sys%numbrella )%p, .false. )
             sys%optype( sys%numbrella ) = OP_TRIANGLE
          else if ( sys%id08( sys%numbrella ) == "LOGTRIAN" ) then
             call new( sys%tri( sys%numbrella )%p, .true. )
             sys%optype( sys%numbrella ) = OP_TRIANGLE
          else if ( sys%id08( sys%numbrella ) == "Q6LOCAL_" ) then
             call new( sys%q6s( sys%numbrella )%p )
             allocate( sys%q6d( sys%numbrella )%p )
             sys%optype( sys%numbrella ) = OP_Q6T
          else if ( sys%id08( sys%numbrella ) == "EPOP____" ) then
             sys%optype( sys%numbrella ) = OP_EP
          endif
          call umbrella_ReadParam( sys%umbrella( sys%numbrella ), file )
          cycle
       endif
       !
       !Replica exchangeable pairs 
       !
       if( tag == "@RXPR" ) then
          call new( sys%rxgraph )
          call graph_read( sys%rxgraph, file )
          cycle
       endif
       if( tag == "@NRXG" ) then
          read(file,*) sys%nexchange
          cycle
       endif
       if(tag.eq."@SNPO")then
          read(file,*) sys%snpo
          cycle
       endif
       if(tag.eq."@NLOG")then
          read(file,*) sys%nlog
          cycle
       endif
       if(tag.eq."@MCLP")then
          read(file,*) sys%nloop
          cycle
       endif
       if(tag.eq."@TEMP")then
          read(file,*) sys%temp
          sys%beta = 1d0 / sys%temp
          cycle
       endif
       if(tag.eq."@MPIN")then
          read(file,*) nprocs
#ifdef USE_FAKEMPI
          call fakempif_new( nprocs, ierr )
          call fakempif_comm_rank( sys%myrank, ierr )
          call fakempif_comm_size( sys%nprocs, ierr )
#endif
          write(STDERR,*) "NODE", sys%myrank, sys%nprocs, nprocs
          do i=0, sys%myrank
             read(file,*) sys%infile
          enddo
          !
          !末尾の空白を落し、拡張子を付加する。
          !
          write(STDERR,*)"CONFIG FILE 2:", sys%infile
          sys%outfile = trim( sys%infile ) // ".out"
          cycle
       endif
       if(tag.eq."@ID08")then
          read(file,*) lastid
          cycle
       endif
       if(tag.eq."@NCMP")then
          read(file,*) sys%incompo
          if(sys%incompo.gt.MAXCOMPO)then
             call die( error_too_many_components, "mctest 1" )
          endif
          cycle
       endif
       if(tag.eq.'@NX4A')then
          if(lastid.eq."        ")then
             write(STDERR,*) "Error: Anonymous data"
             stop
          endif
          sys%ncompo = sys%ncompo + 1
          if(sys%ncompo.gt.MAXCOMPO)then
             write(STDERR,*) "Number of componets exceeded the&
                  & limit."
             stop
          endif
          sys%mol(sys%ncompo)%id=lastid
          if(lastid.eq."TIP4P   ")then
             call new(sys%mol(sys%ncompo),5,6, RIGID_MODE, .not. FIXED )
             call Rigid_TIP4P_Constructor(sys%rigid(sys%ncompo))
             call TIP4P_SetInteraction(sys%stdint(sys%ncompo))
          else if(lastid.eq."TIP5P   ")then
             call new(sys%mol(sys%ncompo),6,6, RIGID_MODE, .not. FIXED )
             call Rigid_TIP5P_Constructor(sys%rigid(sys%ncompo))
             call TIP5P_SetInteraction(sys%stdint(sys%ncompo))
          else if(lastid.eq."NVDE____")then
             call new(sys%mol(sys%ncompo),7,6, RIGID_MODE, .not. FIXED )
             call Rigid_NvdE_Constructor(sys%rigid(sys%ncompo))
             call NvdE_SetInteraction(sys%stdint(sys%ncompo))
          else if (lastid.eq."OPLSMEOH")then
             call new(sys%mol(sys%ncompo),4,6, RIGID_MODE, .not. FIXED )
             call Rigid_OPLSMeOH_Constructor(sys%rigid(sys%ncompo))
             call OPLSMeOH_SetInteraction(sys%stdint(sys%ncompo))
          else
             write(STDERR,*) "Unknown rigid body: ",lastid
             stop
          endif
          call Rigid_ReadNX4A(sys%rigid(sys%ncompo),sys%mol(sys%ncompo),sys%site, file)
          sys%dof = sys%dof + Mol_DoF(sys%mol(sys%ncompo))
          !
          !指定された成分数を読みこんだら入力をうちきる。
          !
          if(sys%ncompo==sys%incompo)then
             exit
          endif
          lastid=""
          cycle
       endif
       call Box_Loader(sys%box,file,tag)
       call Cutoff_Loader(sys%cutoff,file,tag)
       call Grid_Loader(sys%grid,file,tag)
       call Book_Loader(sys%book,file,tag)
    enddo
999 continue
  end subroutine ReadParams

end module mc2util

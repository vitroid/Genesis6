module MDreplicaexchange
  use random_module
  use graph2
  implicit none
  type sMDReplica
     !beta in 1/K
     !
     real(kind=8) :: beta, temp
     !
     !node毎の設定ファイル
     !
     character(len=1024) :: infile, outfile, outtraj, intraj
     integer :: nprocs, myrank
     type(pRandom), pointer :: rand
     type(pGraph), pointer  :: rxgraph
     type(pPair), pointer     :: rxpair
     integer :: ninnerloop, nexchange
  end type sMDReplica

  interface new
     module procedure MDReplica_new
  end interface

contains

  subroutine MDReplica_new( rep )
    type(sMDReplica),intent(INOUT)  :: rep
    rep%beta = 1d0/400d0
    rep%temp = 400d0
  end subroutine MDReplica_new
  
  subroutine MDReplica_Loader( rep, file, tag )
    integer, intent(IN) :: file
    type(sMDReplica), intent(INOUT) :: rep
    !local
    character(len=5), intent(IN) :: tag
    integer :: ierr, seed
    integer :: nprocs, i
    
    write(STDERR,*) tag
    if(tag.eq."@SEED")then
       read(file,*) seed
       call new( rep%rand, seed )
    endif
    if(tag.eq."@MTRN")then
       call random_loadMTRN( rep%rand, file )
    endif
    if(tag.eq."@INNR")then
       read(file,*) rep%ninnerloop
    endif
    !
    !Replica exchangeable pairs 
    !
    if( tag == "@RXPR" ) then
       call new( rep%rxgraph )
       call graph_read( rep%rxgraph, file )
    endif
    if( tag == "@NRXG" ) then
       read(file,*) rep%nexchange
    endif
    if(tag.eq."@TEMP")then
       read(file,*) rep%temp
       rep%beta = 1d0 / rep%temp
    endif
    if(tag.eq."@MPIN")then
       read(file,*) nprocs
       call fakempif_new( nprocs, ierr )
       call fakempif_comm_rank( rep%myrank, ierr )
       call fakempif_comm_size( rep%nprocs, ierr )
       write(STDERR,*) "NODE", rep%myrank, rep%nprocs, nprocs
       do i=0, rep%myrank
          read(file,*) rep%infile
       enddo
       !
       !末尾の空白を落し、拡張子を付加する。
       !
       write(STDERR,*)"CONFIG FILE 2:", rep%infile
       rep%outfile = trim( rep%infile ) // ".log"
       rep%outtraj = trim( rep%infile ) // ".traj"
       rep%intraj = trim( rep%infile ) // ".intraj"
    endif
  end subroutine MDReplica_Loader

end module MDreplicaexchange

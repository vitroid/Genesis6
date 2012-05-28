! -*- f90 -*-
!
!ToDo
!  1. Flexible molecule�ؤ��б�(������Argon�η׻����Ǥ��ʤ�)
!
program main
  use common_module
  use rigid_module
  use flex_module
  use mol_module
  use site_module
  use box_module
  use grid_module
  use cutoff_module
  use rigidquench_module
  use quench_module
  use standard_interaction_module
  use rigid_setup_module

  implicit none

  type(sRigid), target  :: rigid(MAXCOMPO)
  type(sFlex)  ,target  :: flex( MAXCOMPO )
  type(sMol)  , target  :: mol(MAXCOMPO)
  type(sCutoff),target  :: co
  type(sSite)  ,target  :: site
  type(sBox)   ,target  :: box
  type(sGrid)  ,target  :: g
  type(sStdInt),target  :: si(MAXCOMPO)
  !real(kind=8) :: rc
  character(len=5) :: tag
  !,fupdated
  !real(kind=8) :: range
  real(kind=8) :: error, emin
  integer :: err
  integer :: niter
  integer :: i,j,k
  integer :: offset
  integer :: dof
  logical :: isFixed
  !real(kind=8) :: ea,eb,ec,aaa,bbb,ccc,ddd
  real(kind=8), allocatable :: p(:)
  character(len=8) :: id
  integer :: incompo, ncompo
  integer, parameter :: mode = RQ_QUAT
  !
  ! counter for energy and force
  !
  integer :: ecount, fcount
  !
  !����ʬ�Ҥμ����μ���
  !
  real(kind=8) :: molmass(MAXSITE)
  !
  !Fast calculation
  !
  logical :: ffast
  !     initialization
  !     ���ޤ����������������Ȥꤢ��������������٤���¤�Τ򤹤٤Ƥ��Τޤ�
  !     �錄�����Ȥˤ��롣
  !
  !ʬ�Ҵ���ߺ��Ѥ����㷸����ݤ���
  !
  real(kind=8) :: ivscale( MAXCOMPO, MAXCOMPO ), value
  !
  !��ߺ��Ѥ����㷸����ݤ��롣����ͤ�1
  !
  ivscale(:,:) = 1d0
  call new(site, MAXsite_TOTAL)
  call new(g)
  call new(box)
  call new(co)

  isFixed = .false.
  id = "        "
  !
  !outer loop
  !
  !
  !Iteration ( overridden by @NITE )
  !
  niter = 1000
  !
  !Number of components specified by @NCMP
  !
  incompo = 0
  ncompo  = 0
  ffast   = .false.
  do
     !
     !Read data from stdin
     !
     do
        read(STDIN,'(a5)',END=999) tag
        write(STDERR,*) tag
        if(tag.eq."@ID08")then
           read(STDIN,*) id
           cycle
        endif
        if(tag.eq."@FAST")then
           read(STDIN,*) i
           ffast = ( i .ne. 0 )
           cycle
        endif
        !
        !i���ܤ���ʬ��j���ܤ���ʬ����ߺ��Ѥ����������ݤ��롣
        !
        if(tag.eq."@IVSC")then
           read(STDIN,*) i,j,value
           ivscale(i,j) = value
           ivscale(j,i) = value
           cycle
        endif
        !
        !1�ʤ鼡����ʬ�κ�ɸ�����
        !
        if(tag.eq."@FIXC")then
           read(STDIN,*) isFixed
           cycle
        endif
        if(tag.eq."@NITE")then
           read(STDIN,*) niter
           cycle
        endif
        !
        !incompo is effective if specified once.
        !
        if(tag.eq."@NCMP")then
           read(STDIN,*) incompo
           if ( ncompo .gt. 0 ) then
              !
              !���Ǥ��ɤߤ�����ǡ���������ʤ�quench���Ƥ��ޤ���
              !
              exit
           endif
        endif
        call load(box,STDIN,tag)
        call Cutoff_Loader(co,STDIN,tag)
        call Grid_Loader(g,STDIN,tag)
        if(tag.eq.'@NX4A'.or.tag.eq.'@NX3A'.or.tag.eq.'@WTG3')then
           !
           !Increment number of components
           !
           ncompo   = ncompo + 1
           if ( ncompo .gt. MAXCOMPO ) then
              write(STDERR,*) "Too many cmponents. Stop."
              stop
           endif
           mol(ncompo)%id = id
           !
           !id���Ȥ˰㤦��ˡ�ǽ������
           !
           err = rigid_setup( rigid(ncompo), id, isFixed, mol(ncompo), si(ncompo) )
           if ( err  .ne. error_none ) then
              call die( err, "MixtureQuench 1" )
           endif
#ifdef DEBUG
           if(id.eq."TIP4P   ")then
              call new(mol(ncompo),5,6,.true., isFixed )
              call Rigid_TIP4P_Constructor(rigid(ncompo))
              call Mol_SetAtomName(mol(ncompo),tip4pName)
              call TIP4P_SetInteraction(si(ncompo))
           else if(id.eq."TIP5P   ")then
              call new(mol(ncompo),6,6,.true., isFixed )
              call Rigid_TIP5P_Constructor(rigid(ncompo))
              call Mol_SetAtomName(mol(ncompo),tip5pName)
              call TIP5P_SetInteraction(si(ncompo))
           else if (id.eq.NvdE_ID08)then
              call new(mol(ncompo),NVDESITE,6,.true., isFixed )
              call Rigid_NvdE_Constructor(rigid(ncompo))
              call Mol_SetAtomName(mol(ncompo),NvdEName)
              call NvdE_SetInteraction(si(ncompo))
           else  if(id.eq."OPLSMEOH")then
              call new(mol(ncompo),4,6,.true., isFixed )
              call Rigid_OPLSMeOH_Constructor(rigid(ncompo))
              call Mol_SetAtomName(mol(ncompo),meohName)
              call OPLSMeOH_SetInteraction(si(ncompo))
           else if (id.eq."SPC_E   ")then
              call new(mol(ncompo),4,6,.true., isFixed )
              call Rigid_SPCE_Constructor(rigid(ncompo))
              call Mol_SetAtomName(mol(ncompo),waterName)
              call SPCE_SetInteraction(si(ncompo))
           else if (id.eq."ST2_____")then
              call new(mol(ncompo),6,6,.true., isFixed )
              call Rigid_ST2_Constructor(rigid(ncompo))
              call Mol_SetAtomName(mol(ncompo),st2Name)
              call ST2_SetInteraction(si(ncompo))
           else if (id.eq."EPM2CO2 ")then
              call new(mol(ncompo),4,5,.true., isFixed )
              call Rigid_EPM2CO2_Constructor(rigid(ncompo))
              call Mol_SetAtomName(mol(ncompo),co2Name)
              call EPM2CO2_SetInteraction(si(ncompo))
           else if (id.eq."FLATFTHF")then
              call new(mol(ncompo),14,6,.true., isFixed )
              call Rigid_FLATFTHF_Constructor(rigid(ncompo))
              call Mol_SetAtomName(mol(ncompo),thfName)
              call FLATFTHF_SetInteraction(si(ncompo))
           else if (id.eq."CJTHF___")then
              call new(mol(ncompo),UATHFSITE,6,.true., isFixed )
              call Rigid_CJTHF_Constructor(rigid(ncompo))
              call Mol_SetAtomName(mol(ncompo),uathfName)
              call CJTHF_SetInteraction(si(ncompo))
           else if ( id .eq. "SSPKMET_" )then
              call new(mol(ncompo),6,6,.true., isFixed )
              call Rigid_SSPKMethane_Constructor(rigid(ncompo))
              call Mol_SetAtomName(mol(ncompo),methaneName)
              call SSPKMethane_SetInteraction(si(ncompo))
           else if ( id .eq. "LJAR    " )then
              call new(mol(ncompo),2,3,.false., isFixed )
              call Flex_SetMass(flex(ncompo),mol(ncompo),molmass)
              call Mol_SetAtomName(mol(ncompo),argonName)
              call Argon_SetInteraction(si(ncompo))
           else if ( id .eq. "LJME____" )then
              call new(mol(ncompo),2,3,.false., isFixed )
              call UAMEthane_GetMass( molmass )
              call Mol_SetAtomName(mol(ncompo),uamethaneName)
              call UAMethane_SetInteraction(si(ncompo))
           else
              write(STDERR,*) "Error: Unknown molecule type",id
              stop
           endif
#endif
           !
           !�ǡ������ɤߤ���
           !
           if(tag.eq.'@NX4A'.or.tag.eq.'@WTG3')then
              call Rigid_ReadNX4A(rigid(ncompo),mol(ncompo),site,STDIN)
           else if ( tag.eq.'@NX3A' ) then
              call Rigid_ReadNX3A(rigid(ncompo),mol(ncompo),site,STDIN)
           endif
           isFixed = .false.
           if ( ncompo .eq. incompo ) then
              !
              !incompo�Ĥ���ʬ���ɤߤ������quench����
              !
              exit
           endif
        endif
     enddo
999  if (ncompo.eq.0) then
        stop
     endif
     !
     !Quench�ν����������å�ʬ��ʤɤ������Ǽ�ưŪ�˹Ԥ���
     !
     !call rigidquench_initialize(ncompo,rigid,mol,site,g,co,box,mode)
     if ( ffast ) then
        call rigidquench_initialize(ncompo,rigid,mol,site,g,co,box,mode)
     else
        call rigidquench_initialize(ncompo,rigid,mol,site,g,co,box,mode,si,ivscale)
     endif
     !
     !����������ݡ�������������ͳ��ʬ
     !
     dof = 0
     do i=1,ncompo
        dof = dof + mol(i)%nmol * 7
     enddo
     allocate( p(dof) )
     !
     !��ɸ�ǡ����򥷥ꥢ�饤��(ľ��)
     !
     offset = 0
     do i=1,ncompo
        if( mode .eq. RQ_QUAT ) then
           offset = rigid_serialize_positionq(rigid(i), mol(i), p, offset)
        else
           offset = rigid_serialize_position(rigid(i), mol(i), p, offset)
        endif
     enddo
     !
     !Quench�¹�
     !
     call frprmn(p, offset, error, niter, emin)
     !
     !���ͥ륮����ɽ��
     !
     write(STDERR,*) niter,emin,error
     write(STDOUT,'("@ETOT")')
     !
     !�����ͥ륮����K�ǽ���
     !
     write(STDOUT,*) emin * I2J * J2K
     !
     !ľ�󲽤��줿�ǡ������ɸ���᤹��
     !
     offset = 0
     do i=1,ncompo
        if ( mode .eq. RQ_QUAT ) then
           offset = rigid_unserialize_positionq(rigid(i), mol(i), p, offset)
        else
           offset = rigid_unserialize_position(rigid(i), mol(i), p, offset)
        endif
     enddo
     !
     !��ɸ�ν���
     !
     write(STDOUT,'("@NCMP")')
     write(STDOUT,*) ncompo
     do i=1,ncompo
        do j=i, ncompo
           if ( ivscale(i,j) .ne. 1d0 ) then
              write(STDOUT,'("@IVSC")')
              write(STDOUT,*) i,j,ivscale(i,j)
           endif
        enddo
     enddo
     call save( box, STDOUT, BOX_BOX3 )
     do i=1,ncompo
        call save( rigid(i), mol(i), STDOUT, RIGID_NX4A )
     enddo
     call getcount( ecount, fcount )
     write(STDERR,*) ecount, fcount
     !
     !�������γ���
     !
     deallocate(p)
     !
     !�����Ȥ򤹤٤��˴����롣(�����ȿ���0�Ȥ���)
     !
     call new( site, MAXsite_total )
     ncompo = 0
  enddo
end program main

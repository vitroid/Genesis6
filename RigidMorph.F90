!
!���빽¤�����̤ι�¤��Ϣ³�ѷ������롣
!
module rigidmorph_module
  use rigid_module
  implicit none
  integer, parameter :: rigidmorph_none = 0, rigidmorph_active = 1
  type srigidmorph
     integer      :: mode
     integer      :: nmol
     type(sRigid) :: initial, final
     !�ƥ����Ȥ��Ѱ̥٥��ȥ�
     real(kind=8) :: tx(1000), ty(1000), tz(1000)
  end type srigidmorph

  interface save_binary
     module procedure rigidmorph_writebinarynx4a
  end interface
contains

  subroutine RigidMorph_Interpolate( rigidmorph, ratio, rigid )
    use quat_module
    type( sRigidMorph ), intent( IN )  :: rigidmorph
    real(kind=8),        intent( IN )  :: ratio
    type( sRigid ),      intent( OUT ) :: rigid
    integer        :: i,n
    real( kind=8 ) :: a0,b0,c0,d0
    real( kind=8 ) :: a1,b1,c1,d1
    do i=1,rigidmorph%nmol
       rigid%mol(i)%com%vec(1:3) = rigidmorph%initial%mol(i)%com%vec(1:3) + ratio * ( rigidmorph%final%mol(i)%com%vec(1:3) - rigidmorph%initial%mol(i)%com%vec(1:3) )
       a0 = rigidmorph%initial%mol(i)%quat(1)%vec(1)
       b0 = rigidmorph%initial%mol(i)%quat(1)%vec(2)
       c0 = rigidmorph%initial%mol(i)%quat(1)%vec(3)
       d0 = rigidmorph%initial%mol(i)%quat(1)%vec(4)
       a1 = rigidmorph%final%mol(i)%quat(1)%vec(1)
       b1 = rigidmorph%final%mol(i)%quat(1)%vec(2)
       c1 = rigidmorph%final%mol(i)%quat(1)%vec(3)
       d1 = rigidmorph%final%mol(i)%quat(1)%vec(4)
       call qadd( a1,b1,c1,d1, a0,-b0,-c0,-d0 )
       !if ( i.le.3 ) then
       !   write(10,'(a1,6(1x,e17.10))') 'l', 0,0,0,b1,c1,d1
       !   write(10,'(a1,3(1x,e17.10),1x,a1,i1)') 't', b1,c1,d1, "Q", i
       !endif
       call qmul( a1,b1,c1,d1, ratio )
       call qadd( a1,b1,c1,d1, a0,b0,c0,d0 )
       !if ( i.le.3 ) then
       !   write(10,'(a1,6(1x,e17.10))') 'l', 0,0,0,b1,c1,d1
       !   write(10,'(a1,3(1x,e17.10),1x,a1,i1)') 't', b1,c1,d1, "R", i
       !endif
       rigid%mol(i)%quat(1)%vec(1) = a1
       rigid%mol(i)%quat(1)%vec(2) = b1
       rigid%mol(i)%quat(1)%vec(3) = c1
       rigid%mol(i)%quat(1)%vec(4) = d1
    enddo
  end subroutine RigidMorph_Interpolate

  !
  !mix�������Ѳ����������γ���ߺ��������Ѱ̤�Ͽ���Ƥ�����
  !��ͳ���ͥ륮���׻���®���Ǥ��뤫�⤷��ʤ���ʿ��16ǯ8��19��(��)
  !�����Ϥ����ä���ʬ����Ф�ä�®���Ǥ��뤬���ޤ��ϼ¾ڼ¸��Ȥ������Ȥ�
  !��ʬ�ǽ񤯡�(���ʤ���˽)
  !
  !ʿ��16ǯ8��19��(��)�ʤ�Ȥ������� 
  !
  subroutine RigidMorph_SetTangentVectors0( rigidmorph, ratio, rigid, mol, site )
    use quat_module
    use mol_module
    use site_module
    type( sRigidMorph ), intent( INOUT )  :: rigidmorph
    real(kind=8),        intent( IN )  :: ratio
    type( sRigid ),      intent( INOUT ) :: rigid
    type( sMol ),        intent( INOUT ) :: mol
    type( sSite ),       intent( INOUT ) :: site
    type( sRigid ) :: middle
    integer        :: i,j
    call Rigid_SetSitePosition( rigid, mol, site )
    !do i=1,14
    !   write(STDERR,*) i,site%x( i + mol%offset )
    !enddo
    do i=1, mol%nmol * mol%nsite
       rigidmorph%tx(i) = site%x( i + mol%offset )
       rigidmorph%ty(i) = site%y( i + mol%offset )
       rigidmorph%tz(i) = site%z( i + mol%offset )
    enddo
    call RigidMorph_Interpolate( rigidmorph, ratio - 0.01, middle )
    middle%molx = rigid%molx
    middle%moly = rigid%moly
    middle%molz = rigid%molz
    call Rigid_SetSitePosition( middle, mol, site )
    !do i=1,14
    !   write(STDERR,*) i,site%x( i + mol%offset )
    !enddo
    do i=1, mol%nmol * mol%nsite
       rigidmorph%tx(i) = ( rigidmorph%tx(i) - site%x( i + mol%offset ) ) / 0.01d0
       rigidmorph%ty(i) = ( rigidmorph%ty(i) - site%y( i + mol%offset ) ) / 0.01d0
       rigidmorph%tz(i) = ( rigidmorph%tz(i) - site%z( i + mol%offset ) ) / 0.01d0
       if ( i .le. 3 ) then
          !write(11,'(a1,6(1x,e17.10))') 'l', 0,0,0,n1,n2,n3
          !write(11,'(a1,6(1x,e17.10))') 'l', 0,0,0,rigid%mol(i)%intrax(j),rigid%mol(i)%intray(j),rigid%mol(i)%intraz(j)
          !write(10,'(a1,6(1x,e17.10))') 'l', 0,0,0,rigidmorph%tx(i),rigidmorph%ty(i),rigidmorph%tz(i)
          !write(11,'(a1,3(1x,e17.10),1x,a1)') 't', n1,n2,n3, "n"
          !write(11,'(a1,3(1x,e17.10),1x,a1)') 't', rigid%mol(i)%intrax(j),rigid%mol(i)%intray(j),rigid%mol(i)%intraz(j), "i"
          !write(10,'(a1,3(1x,e17.10),1x,a1,i1)') 't', rigidmorph%tx(i),rigidmorph%ty(i),rigidmorph%tz(i), "T", i
       endif
    enddo
  end subroutine RigidMorph_SetTangentVectors0
  
  !
  !��ꥹ�ޡ��Ȥʼ��Ƿ׻����롣
  !ʿ��16ǯ8��20��(��)�ɤ����Ƥ⤦�ޤ��׻��Ǥ��ʤ����͸��Ҥ��鼴���������ˡ�򡢲������㤤���Ƥ���Τ��������������֤���ä����ʤ��Τǡ����̾�μ�ˡ��Ȥ����Ȥˤ��롣
  !
  subroutine RigidMorph_SetTangentVectors( rigidmorph, rigid, mol )
    use quat_module
    use mol_module
    type( sRigidMorph ), intent( INOUT )  :: rigidmorph
    type( sRigid ),      intent( IN ) :: rigid
    type( sMol ),        intent( IN ) :: mol
    integer        :: i,j,k
    real( kind=8 ) :: a0,b0,c0,d0
    real( kind=8 ) :: a1,b1,c1,d1
    real( kind=8 ) :: phi, n1,n2,n3, d
    k = 0
    do i=1,rigidmorph%nmol
       a0 = rigidmorph%initial%mol(i)%quat(1)%vec(1)
       b0 = rigidmorph%initial%mol(i)%quat(1)%vec(2)
       c0 = rigidmorph%initial%mol(i)%quat(1)%vec(3)
       d0 = rigidmorph%initial%mol(i)%quat(1)%vec(4)
       a1 = rigidmorph%final%mol(i)%quat(1)%vec(1)
       b1 = rigidmorph%final%mol(i)%quat(1)%vec(2)
       c1 = rigidmorph%final%mol(i)%quat(1)%vec(3)
       d1 = rigidmorph%final%mol(i)%quat(1)%vec(4)
       call qadd( a1,b1,c1,d1, a0,-b0,-c0,-d0 )
       ! 
       !��ž���򡢼������ȳ��٤�ʬ�䤹�롣
       !
       d = sqrt( b1**2 + c1**2 + d1**2 )
       if ( d .le. 0d0 ) then
          !write(STDERR,*) "D:", d
          n1 = 0d0
          n2 = 0d0
          n3 = 0d0
       else
          phi = 2d0 * acos( a1 )
          n1 =-c1 / d * phi
          n2 =-d1 / d * phi
          n3 =-b1 / d * phi
          !write(STDERR,*) "phi:", phi, n1,n2,n3,b1,c1,d1
          !
          !ʿ��16ǯ8��20��(��)
          !���󤶤��ϫ������̡���Τ褦���б�����ǧ���줿��
          !
       endif
       do j=1, mol%nsite
          k = k + 1
          rigidmorph%tx(k) = n2 * rigid%mol(i)%intraz(j) - n3 * rigid%mol(i)%intray(j)
          rigidmorph%ty(k) = n3 * rigid%mol(i)%intrax(j) - n1 * rigid%mol(i)%intraz(j)
          rigidmorph%tz(k) = n1 * rigid%mol(i)%intray(j) - n2 * rigid%mol(i)%intrax(j)
          !if ( k .le. 3 ) then
          !   write(10,'(a1,6(1x,e17.10))') 'l', 0,0,0,n1,n2,n3
          !   write(10,'(a1,6(1x,e17.10))') 'l', 0,0,0,rigid%mol(i)%intrax(j),rigid%mol(i)%intray(j),rigid%mol(i)%intraz(j)
          !   write(10,'(a1,6(1x,e17.10))') 'l', 0,0,0,rigidmorph%tx(k),rigidmorph%ty(k),rigidmorph%tz(k)
          !   write(10,'(a1,3(1x,e17.10),1x,a1)') 't', n1,n2,n3, "n"
          !   write(10,'(a1,3(1x,e17.10),1x,a1,i1)') 't', rigid%mol(i)%intrax(j),rigid%mol(i)%intray(j),rigid%mol(i)%intraz(j), "i", k
          !   write(10,'(a1,3(1x,e17.10),1x,a1,i1)') 't', rigidmorph%tx(k),rigidmorph%ty(k),rigidmorph%tz(k), "t", k
          !endif
          rigidmorph%tx(k) = rigidmorph%tx(k) + rigidmorph%final%mol(i)%com%vec(1) - rigidmorph%initial%mol(i)%com%vec(1)
          rigidmorph%ty(k) = rigidmorph%ty(k) + rigidmorph%final%mol(i)%com%vec(2) - rigidmorph%initial%mol(i)%com%vec(2)
          rigidmorph%tz(k) = rigidmorph%tz(k) + rigidmorph%final%mol(i)%com%vec(3) - rigidmorph%initial%mol(i)%com%vec(3)
          !
          !�¿ʤΤۤ�������ʤ���
          !
          !rigidmorph%tx(k) = rigidmorph%final%mol(i)%com%vec(1) - rigidmorph%initial%mol(i)%com%vec(1)
          !rigidmorph%ty(k) = rigidmorph%final%mol(i)%com%vec(2) - rigidmorph%initial%mol(i)%com%vec(2)
          !rigidmorph%tz(k) = rigidmorph%final%mol(i)%com%vec(3) - rigidmorph%initial%mol(i)%com%vec(3)
       enddo
    enddo
  end subroutine RigidMorph_SetTangentVectors
  
  subroutine RigidMorph_WriteBinaryNX4A( rigidmorph, mol, output )
    type( sRigidMorph ), intent( IN )  :: rigidmorph
    type( sMol ), intent(IN)           :: mol
    integer, intent(IN)                :: output
    character(len=8) :: id
    if ( rigidmorph%mode .eq. rigidmorph_active ) then
       id = mol%id(1:8)
       write( output ) "@ID08"
       write( output ) id
       write( output ) '@RMMC'
       call Rigid_WriteBinaryNX4A( rigidmorph%initial, rigidmorph%nmol, output )
       call Rigid_WriteBinaryNX4A( rigidmorph%final,   rigidmorph%nmol, output )
    endif
  end subroutine RigidMorph_WriteBinaryNX4A

  subroutine RigidMorph_ReadBinaryNX4A( rigidmorph, rigid, mol, site, input )
    use mol_module
    use error_module
    use site_module
    type( sRigidMorph ), intent( INOUT )  :: rigidmorph
    type( sRigid ), intent( INOUT )       :: rigid
    type( sMol ), intent( INOUT )         :: mol
    type( sSite ), intent( INOUT )        :: site
    integer, intent(IN)                   :: input
    character(len=5) :: tag
    integer          :: nmol
    ! Ignore leading tags ( @NX4A )
    tag = '     '
    do while ( tag .ne. '@NX4A' )
       read(input) tag
    enddo
    call Rigid_ReadBinaryNX4A( rigid, mol, site, input )
    rigidmorph%initial = rigid
    rigidmorph%nmol    = mol%nmol
    tag = '     '
    do while ( tag .ne. '@NX4A' )
       read(input) tag
    enddo
    !read(input) tag
    read(input) nmol
    if ( nmol .ne. mol%nmol ) then
       call die ( error_different_num_of_molecules, "RigidMorph 1" )
    endif
    call Rigid_ReadOverBinaryNX4A( rigidmorph%final, nmol, input )
    rigidmorph%mode = rigidmorph_active
  end subroutine RigidMorph_ReadBinaryNX4A

  subroutine RigidMorph_ReadNX4A( rigidmorph, rigid, mol, site, input )
    use mol_module
    use error_module
    use site_module
    type( sRigidMorph ), intent( INOUT )  :: rigidmorph
    type( sRigid ), intent( INOUT )       :: rigid
    type( sMol ), intent( INOUT )         :: mol
    type( sSite ), intent( INOUT )        :: site
    integer, intent(IN)                   :: input
    character(len=5) :: tag
    integer          :: nmol
    ! Ignore leading tags ( @NX4A )
    tag = '     '
    do while ( tag .ne. '@NX4A' )
       read(input,*) tag
    enddo
    !read(input,*) tag
    call Rigid_ReadNX4A( rigid, mol, site, input )
    rigidmorph%initial = rigid
    rigidmorph%nmol    = mol%nmol
    tag = '     '
    do while ( tag .ne. '@NX4A' )
       read(input,*) tag
    enddo
    read(input,*) nmol
    if ( nmol .ne. mol%nmol ) then
       call die ( error_different_num_of_molecules, "RigidMorph 1" )
    endif
    call Rigid_ReadOverNX4A( rigidmorph%final, nmol, input )
    rigidmorph%mode = rigidmorph_active
  end subroutine RigidMorph_ReadNX4A

end module rigidmorph_module
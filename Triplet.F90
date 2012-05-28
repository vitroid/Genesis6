!3-body interaction handler

module triplet_module
  implicit none
  integer, parameter :: MAXNEI = 40
  type sTriplet
     integer :: num, size
     ! members of the triplet
     integer, pointer :: center(:), left(:), right(:)
     ! corresponding pair #
     integer, pointer :: centerLeft(:), centerRight(:)
  end type sTriplet
contains

  !create triplet list from pair list
  !(Single component only)
  subroutine ComposeTriplets( iv, rc, nmol, triplet )
    use error_module
    use interaction_module
    type(sInteraction),intent(in), target :: iv
    real(kind=8), intent(IN)              :: rc
    integer, intent(IN)                   :: nmol
    type( sTriplet ), intent(INOUT)       :: triplet
    ! neighbor list
    integer :: nnei( nmol ), nei( MAXNEI, nmol ), pairID( MAXNEI, nmol )
    integer :: i,j,k, jj, kk, maxnnei, num, size, jid, kid
    nnei(:) = 0
    ! compose neighbor list
    do k=1,iv%Npair
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       nnei( i ) = nnei( i ) + 1
       if ( MAXNEI < nnei(i) ) call die (0, "triplet 1")
       nei( nnei( i ), i ) = j
       pairID( nnei( i ), i ) = k
       nnei( j ) = nnei( j ) + 1
       if ( MAXNEI < nnei(j) ) call die (0, "triplet 2")
       nei( nnei( j ), j ) = i
       pairID( nnei( j ), j ) = -k
    enddo
    maxnnei = 0
    do i=1, nmol
       if ( maxnnei < nnei( i ) ) then
          maxnnei = nnei( i )
       endif
    enddo
    ! Estimated size of triplet list
    size = nmol * maxnnei * ( maxnnei - 1 ) / 2
    ! locate arrays only when size grows
    !write(STDERR,*) "LASTSIZE", triplet%size, size
    if ( triplet%size < size ) then
       if ( triplet%size .ne.0 ) then
          deallocate( triplet%center )
          deallocate( triplet%left )
          deallocate( triplet%right )
          deallocate( triplet%centerleft )
          deallocate( triplet%centerright )
       endif
       allocate( triplet%center( size ) )
       allocate( triplet%left( size ) )
       allocate( triplet%right( size ) )
       allocate( triplet%centerleft( size ) )
       allocate( triplet%centerright( size ) )
       triplet%size = size
    endif
    ! Pickup triplets
    num = 0
    do i = 1, nmol
       do j = 1, nnei( i )
          jj = nei( j, i )
          jid = pairID( j, i )
          do k = j+1, nnei( i )
             kk = nei( k, i )
             kid = pairID( k, i )
             num = num + 1
             triplet%center( num ) = i
             triplet%left( num )   = jj
             triplet%right( num )  = kk
             triplet%centerleft( num )   = jid
             triplet%centerright( num )  = kid
             !write(STDERR,'(6(i12,1x))') num, i, jj,kk,jid,kid
          enddo
       enddo
    enddo
    triplet%num = num
    !write(STDERR,*) "num", num
    !do i=1,num
    !   write(STDERR,'(6(i12,1x))') i, triplet%center(i),triplet%left(i), triplet%right(i),triplet%centerleft(i), triplet%centerright(i)
    !enddo
  end subroutine ComposeTriplets
  
end module triplet_module

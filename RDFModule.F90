module rdf_module
  use distrib_module
  implicit none

contains

  subroutine RDF_Initialize( rdf, nisite, iname, njsite, jname, nbin, binwidth )
    type(sRDF), intent(OUT) :: rdf
    integer, intent(IN) :: nisite, njsite
    real(kind=8), intent(in) :: binwidth
    integer, intent(in) :: nbin
    character(len=8), intent(IN) :: iname(*), jname(*)
    integer :: site,typ, i,j
    rdf%nitype = 0
    rdf%njtype = 0
    allocate(rdf%itype(nisite))
    allocate(rdf%jtype(njsite))
    allocate(rdf%iname(nisite))
    allocate(rdf%jname(njsite))
    allocate(rdf%imult(nisite))
    allocate(rdf%jmult(njsite))
    rdf%imult(:) = 0
    rdf%jmult(:) = 0
    do site=1,nisite
       do typ=1, rdf%nitype
          if ( rdf%iname(typ) .eq. iname(site) ) then
             goto 1
          endif
       enddo
       rdf%nitype = rdf%nitype + 1
       rdf%iname(rdf%nitype) = iname(site)
       typ = rdf%nitype
1      continue
       rdf%itype(site) = typ
       rdf%imult(typ) = rdf%imult(typ) + 1
    enddo
    do site=1,njsite
       do typ=1,rdf%njtype
          if ( rdf%jname(typ) .eq. jname(site) ) then
             goto 2
          endif
       enddo
       rdf%njtype = rdf%njtype + 1
       rdf%jname(rdf%njtype) = jname(site)
       typ = rdf%njtype
2      continue
       rdf%jtype(site) = typ
       rdf%jmult(typ) = rdf%jmult(typ) + 1
    enddo
    allocate(rdf%hist(rdf%nitype, rdf%njtype))
    do i=1,rdf%nitype
       do j=1,rdf%njtype
          call new( rdf%hist(i,j)%p, nbin, binwidth, 0d0 )
       enddo
    enddo
  end subroutine RDF_Initialize

  subroutine rdf_normalize( rdf, nloop, box, mi, mj )
    use mol_module
    use box_module
    type(sRDF), intent(INOUT) :: rdf
    integer, intent(IN) :: nloop
    type(sBox), intent(IN) :: box
    type(sMol), intent(IN) :: mi
    type(sMol), intent(IN), optional :: mj
    integer :: accum, i, j
    real(kind=8) :: nmol
    real(kind=8) :: density
    if ( present( mj ) ) then
       nmol = mj%nmol
    else
       nmol = mi%nmol - 1
       nmol = nmol * 0.5d0
    end if
    do i=1, rdf%nitype
       do j=1, rdf%njtype
          density = rdf%jmult(j) * nmol / volume(box)
          accum = nloop * mi%nmol * rdf%imult(i)
          write(6,'("#",a8,a1,a8)') rdf%iname(i), "-", rdf%jname(j) 
          call comrdf_normalize( rdf%hist(i,j)%p, density, accum )
          write(6,'()')
       enddo
    enddo
  end subroutine rdf_normalize

  !subroutine rdf_show( rdf )
  !  type(sRDF), intent(IN) :: rdf
  !  integer :: accum, i, j
  !  real(kind=8) :: density
  !  do i=1, rdf%nitype
  !     do j=1, rdf%njtype
  !        call histogram_show( rdf%hist(i,j) )
  !     enddo
  !  enddo
  !end subroutine rdf_show

  !
  !動径分布関数を計算する。
  ! moved from StdInteraction.F90
  !
  subroutine radialdistrib(iv,mi,mj,site,rdf)
    use standard_interaction_module
    use physconst_module
    use interaction_module
    use mol_module
    use site_module
    use distrib_module
    implicit none
    type(sInteraction),intent(in),   target :: iv
    type(sMol)        ,intent(in)    :: mi,mj
    type(sSite)       ,intent(inout),target :: site
    type(sRDF)        ,intent(inout) :: rdf
    real(kind=8) :: dx,dy,dz
    real(kind=8) :: radius
    integer :: i,j,k
    integer :: ii,jj,isite,jsite,itype,jtype
    real(kind=8),dimension(:),allocatable :: ep
    do isite = 1, mi%nsite
       do jsite = 1, mj%nsite
          itype = rdf%itype(isite)
          jtype = rdf%jtype(jsite)
          do k=1,iv%Npair
             i = iv%pair_i(k)
             j = iv%pair_j(k)
             !     write(STDERR,*) i,j
             ii = (i-1)*mi%nsite+mi%offset + isite
             jj = (j-1)*mj%nsite+mj%offset + jsite
             dx = site%x(ii) - site%x(jj) - iv%ox(k)
             dy = site%y(ii) - site%y(jj) - iv%oy(k)
             dz = site%z(ii) - site%z(jj) - iv%oz(k)
             radius = dsqrt( dx**2 + dy**2 + dz**2 )
             call Histogram_Accumulate(rdf%hist(itype, jtype)%p, radius, 1d0)
          enddo
       enddo
    enddo
  end subroutine radialdistrib

end module rdf_module

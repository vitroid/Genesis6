! -*- f90 -*-
!Formerly known as Montecarlo_module
module pairinfo_module
  use common_module
  type sPairInfo
     real(kind=8)  :: ep(MAXMOL,MAXMOL)
     !real(kind=8)  :: vr(MAXMOL,MAXMOL)
     real(kind=8)  :: epsum !, vrsum
  end type sPairInfo

contains

  subroutine pairenergy(iv,mi,mj,site,pairinfo,custom, si,sj)
    use interaction2_module
    use mol_module
    use site_module
    use interaction_module
    use standard_interaction_module
    type(sInteraction),intent(in),   target :: iv
    type(sMol)        ,intent(in)    :: mi,mj
    type(sStdInt)     ,intent(in),optional :: si,sj
    type(sSite)       ,intent(inout),target :: site
    type(sPairInfo)   ,intent(out)   :: pairinfo
    logical, intent(IN)                   :: custom
    real(kind=8) :: epsum
    real(kind=8) :: e0
    integer :: i,j,k
    real(kind=8),allocatable :: ep(:)
    allocate( ep(iv%npair) )
    if ( present( si ) .and. present( sj ) ) then
       call PairPotential(iv,mi,mj,site,ep, custom, si, sj )
    else
       call PairPotential(iv,mi,mj,site,ep, custom )
    endif
    epsum = 0d0
    do k=1,iv%Npair
       e0 = ep(k)
       i = iv%pair_i(k)
       j = iv%pair_j(k)
       pairinfo%ep(i,j) = e0
       pairinfo%ep(j,i) = e0
       epsum = epsum + e0
    enddo
    pairinfo%epsum = epsum
    deallocate(ep)
  end subroutine pairenergy

end module pairinfo_module

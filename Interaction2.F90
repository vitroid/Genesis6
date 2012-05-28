!Integrate force calculation

module interaction2_module
  use oplsmeoh_tip4p_module
  use st2_module
  use nvde_module
  use tip4p_module
  use oplsmeoh_module
  use sjbenzene_module
  use swsilicon_module

  use physconst_module
  use interaction_module
  use mol_module
  use site_module
  use standard_interaction_module
  implicit none
contains
  subroutine force( iv, mi, mj, site, ep, vrmat, custom, si, sj )
    type(sInteraction),intent(in), target :: iv
    type(sMol)        ,intent(in)         :: mi,mj
    type(sStdInt),intent(in),optional     :: si,sj
    type(sSite)  ,intent(inout),target    :: site
    real(kind=8),intent(out)              :: vrmat(3,3)
    real(kind=8),intent(inout)            :: ep(*)
    logical, intent(IN)                   :: custom

    if ( .not. ( present( si ) .and. present( sj ) ) .or. custom ) then
       !
       !Custom == Non standard == Individually defined
       !
       if ( mi%id(1:8) .eq. "ST2_____" .and. mj%id(1:8) .eq. "ST2_____" )then
          call force_st2(iv,mi,mj,site,ep,vrmat)
          return
       else if ( mi%id(1:8) .eq. "NVDE____" .and. mj%id(1:8) .eq. "NVDE____" )then
          call force_nvde(iv,mi,mj,site,ep,vrmat)
          return
       else if ( mi%id(1:8) .eq. "SJBENZEN" .and. mj%id(1:8) .eq. "SJBENZEN" )then
          call force_sjbenzene(iv,mi,mj,site,ep,vrmat)
          !write(STDERR,*) "SJBENZEN-FORCE"
          return
       else if ( mi%id(1:8) .eq. "SWSILICO" .and. mj%id(1:8) .eq. "SWSILICO" )then
          call force_swsilicon(iv,mi,mj,site,ep,vrmat)
          return
       else
          if ( mi%id(1:8) .eq. "TIP4P   " )then
             if ( mj%id(1:8) .eq. "TIP4P   " )then
                call interaction_force_tip4p( iv,mi,mj,site,ep,vrmat)
                return
             else if ( mj%id(1:8) .eq. "OPLSMEOH" )then
                call interaction_force_tip4p_oplsmeoh( iv,mi,mj,site,ep,vrmat)
                return
             endif
          else if ( mi%id(1:8) .eq. "OPLSMEOH" )then
             if ( mj%id(1:8) .eq. "TIP4P   " )then
                call interaction_force_tip4p_oplsmeoh( iv,mj,mi,site,ep,vrmat)
                return
             else if ( mj%id(1:8) .eq. "OPLSMEOH" )then
                call interaction_force_oplsmeoh( iv,mj,mi,site,ep,vrmat)
                return
             endif
          endif
       endif
       write(STDERR,*) "UNKNOWN COMBINATION OF ", mi%id(1:8), " - ", mj%id(1:8)
       write(STDERR,*) "FALL BACK TO GENERAL PROCEDURE."
    endif
    call StdForce( iv, mi, mj, si, sj, site, ep, vrmat )
  end subroutine force
       

  subroutine PairPotential( iv, mi, mj, site, ep, custom, si, sj )
    type(sInteraction),intent(in), target :: iv
    type(sMol)        ,intent(in)         :: mi,mj
    type(sStdInt),intent(in),optional     :: si,sj
    type(sSite)  ,intent(inout),target    :: site
    real(kind=8),intent(inout)            :: ep(*)
    logical, intent(IN)                   :: custom
    !logical, intent(IN),optional                   :: verbose
    real(kind=8)                          :: vrmat(3,3) ! dummy

    !write(STDERR,*) mi%id(1:8), "-", mj%id(1:8), custom
    if ( .not. ( present( si ) .and. present( sj ) ) .or. custom ) then
       !
       !Custom == Non standard == Individually defined
       !
       if ( mi%id(1:8) .eq. "ST2_____" .and. mj%id(1:8) .eq. "ST2_____" )then
          call PairPotential_st2(iv,mi,mj,site,ep )
          return
       else if ( mi%id(1:8) .eq. "NVDE____" .and. mj%id(1:8) .eq. "NVDE____" )then
          call force_nvde(iv,mi,mj,site,ep, vrmat )
          return
       else if ( mi%id(1:8) .eq. "SJBENZEN" .and. mj%id(1:8) .eq. "SJBENZEN" )then
          call force_sjbenzene(iv,mi,mj,site,ep, vrmat )
          !write(STDERR,*) "SJBENZEN-POTENTIAL"
          return
       else if ( mi%id(1:8) .eq. "SWSILICO" .and. mj%id(1:8) .eq. "SWSILICO" )then
          call force_swsilicon(iv,mi,mj,site,ep,vrmat)
          !write(STDERR,*) "SJBENZEN-FORCE"
          return
       else
          if ( mi%id(1:8) .eq. "TIP4P   " )then
             if ( mj%id(1:8) .eq. "TIP4P   " )then
                call interaction_force_tip4p( iv,mi,mj,site,ep, vrmat )
                return
             else if ( mj%id(1:8) .eq. "OPLSMEOH" )then
                call interaction_force_tip4p_oplsmeoh( iv,mi,mj,site,ep, vrmat )
                return
             endif
          else if ( mi%id(1:8) .eq. "OPLSMEOH" )then
             if ( mj%id(1:8) .eq. "TIP4P   " )then
                call interaction_force_tip4p_oplsmeoh( iv,mj,mi,site,ep, vrmat )
                return
             endif
          endif
       endif
       write(STDERR,*) "UNKNOWN COMBINATION OF ", mi%id(1:8), " - ", mj%id(1:8)
       write(STDERR,*) "FALL BACK TO GENERAL PROCEDURE."
    endif
    call StdPairPotential( iv, mi, mj, si, sj, site, ep )
    !write(6,*) "i2",ep(17)
  end subroutine PairPotential
       
end module interaction2_module

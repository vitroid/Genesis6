! -*- f90 -*-
module rigid_setup_module
contains
  function rigid_setup( rigid, id, isFixed, mol, si )
    use st2_module
    use oplsmeoh_module
    use tip4p_module
    use tip4p_constant_module
    use tip5p_module
    use nvde_module
    use epm2co2_module
    use flatfthf_module
    use sspkmethane_module
    use cjthf_module
    use spce_module
    use sjbenzene_module
    use dumbbell_module
    use piwater_module

    !    use stack_module
    use error_module
    use rigid_module
    use standard_interaction_module
    implicit none
    integer :: rigid_setup
    type(sRigid), intent(INOUT)  :: rigid
    character(len=8), intent(in) :: id
    logical, intent( IN )        :: isFixed
    type(sMol), intent(INOUT)    :: mol
    type(sStdInt), intent(INOUT) :: si
    
    rigid_setup = 0
    if(id.eq."TIP4P   ")then
       call new(mol,5,6,RIGID_MODE , isFixed )
       call Rigid_TIP4P_Constructor(rigid)
       call Mol_SetAtomName(mol,tip4pName)
       call TIP4P_SetInteraction(si)
    else if (id.eq."TIP5P   ")then
       call new(mol,6,6,RIGID_MODE, isFixed )
       call Rigid_TIP5P_Constructor(rigid)
       call Mol_SetAtomName(mol,tip5pName)
       call TIP5P_SetInteraction(si)
    else if ( id .eq. NvdE_ID08 )then
       call new(mol,NVDESITE,6,RIGID_MODE, isFixed )
       call Rigid_NvdE_Constructor(rigid)
       call Mol_SetAtomName(mol,NvdEName)
       call NvdE_SetInteraction(si)
    else if (id.eq."OPLSMEOH")then
       call new(mol,4,6,RIGID_MODE, isFixed )
       call Rigid_OPLSMeOH_Constructor(rigid)
       call Mol_SetAtomName(mol,meohName)
       call OPLSMeOH_SetInteraction(si)
    else if ( id .eq. SPCE_ID08 )then
       call new(mol,4,6,RIGID_MODE, isFixed )
       call Rigid_SPCE_Constructor(rigid)
       call Mol_SetAtomName(mol,waterName)
       call SPCE_SetInteraction(si)
    else if (id.eq."ST2_____")then
       call new(mol,ST2SITE,6,RIGID_MODE, isFixed )
       call Rigid_ST2_Constructor(rigid)
       call Mol_SetAtomName(mol,st2Name)
       call ST2_SetInteraction(si)
    else if (id.eq."EPM2CO2 ")then
       call new(mol,4,5,RIGID_MODE, isFixed )
       call Rigid_EPM2CO2_Constructor(rigid)
       call Mol_SetAtomName(mol,co2Name)
       call EPM2CO2_SetInteraction(si)
    else if (id.eq.piwater_id08)then
       call new(mol,4,5,RIGID_MODE, isFixed )
       call Rigid_PIWATER_Constructor(rigid)
       call Mol_SetAtomName(mol,waterName)
       call PIWATER_SetInteraction(si)
    else if (id.eq."FLATFTHF")then
       call new(mol,14,6,RIGID_MODE, isFixed )
       call Rigid_FLATFTHF_Constructor(rigid)
       call Mol_SetAtomName(mol,thfName)
       call FLATFTHF_SetInteraction(si)
    else if (id.eq."CJTHF___")then
       call new(mol,UATHFSITE,6,RIGID_MODE, isFixed )
       call Rigid_CJTHF_Constructor(rigid)
       call Mol_SetAtomName(mol,uathfName)
       call CJTHF_SetInteraction(si)
    else if (id.eq."SSPKMET_")then
       call new(mol,6,6,RIGID_MODE, isFixed )
       call Rigid_SSPKMethane_Constructor(rigid)
       call Mol_SetAtomName(mol,methaneName)
       call SSPKMethane_SetInteraction(si)
    else if (id.eq. SJBENZENE_ID08 )then
       call new(mol,13,6,RIGID_MODE, isFixed )
       call Rigid_SJBenzene_Constructor(rigid)
       call Mol_SetAtomName(mol,SjbenzeneName)
       call SJBenzene_SetInteraction(si)
    else if (id.eq. DUMBBELL_ID08 )then
       call new(mol,DUMBBELLSITE,5,RIGID_MODE, isFixed )
       call Rigid_Dumbbell_Constructor(rigid)
       call Mol_SetAtomName(mol,DumbbellName)
       call Dumbbell_SetInteraction(si)
    else if (id.eq. DUMBWATER_ID08 )then
       call new(mol,DUMBWATERSITE,6,RIGID_MODE, isFixed )
       call Rigid_SPCE_Constructor(rigid)
       call Mol_SetAtomName(mol,WaterName)
       call DumbWater_SetInteraction(si)
    else
       rigid_setup = error_id_not_found
    endif
  end function rigid_setup
end module rigid_setup_module

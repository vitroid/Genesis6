module flex_setup_module
contains
  function flex_setup( id, isFixed, mol, si, molmass )
    use error_module
    use swsilicon_module
    use argon_module
    use standard_interaction_module
    use mol_module
    implicit none
    integer :: flex_setup
    character(len=8), intent(in) :: id
    logical, intent( IN )        :: isFixed
    type(sMol), intent(INOUT)    :: mol
    type(sStdInt), intent(INOUT) :: si
    real(kind=8), intent(OUT)    :: molmass(*)
    
    flex_setup = 0
    if (id.eq."LJAR    ")then
       call new(mol,2,3,FLEX_MODE , isFixed )
       call Argon_GetMass( molmass )
       call Mol_SetAtomName(mol,argonName)
       call Argon_SetInteraction(si)
    else if (id.eq."LJME____")then
       call new(mol,2,3,FLEX_MODE , isFixed )
       call UAMethane_GetMass( molmass )
       call Mol_SetAtomName(mol,uaMethaneName)
       call UAMethane_SetInteraction(si)
    else if ( id .eq. REDUCEAR_ID08 )then
       call new(mol,2,3,FLEX_MODE , isFixed )
       call IdealGas_GetMass( molmass )
       call Mol_SetAtomName(mol,ArgonName)
       call ReduceAr_SetInteraction(si)
    else if ( id .eq. IDEALGAS_ID08 )then
       call new(mol,2,3,FLEX_MODE , isFixed )
       call IdealGas_GetMass( molmass )
       call Mol_SetAtomName(mol,IdealGasName)
       call IdealGas_SetInteraction(si)
    else if ( id .eq. "SWSILICO" ) then
       call new(mol,2,3,FLEX_MODE , isFixed )
       call SWSilicon_GetMass( molmass )
       call Mol_SetAtomName(mol,SWSiliconName)
       call SWSilicon_SetInteraction( si )
    else
       flex_setup = error_id_not_found
    endif
  end function flex_setup
end module flex_setup_module

! Auto-generated package activation for: chapman_emis_dep
nullify(chem_chapman_emis_dep_inActive)
call mpas_pool_get_package(packages, 'chem_chapman_emis_dep_inActive', &
                           chem_chapman_emis_dep_inActive)
if (associated(chem_chapman_emis_dep_inActive)) then
   chem_chapman_emis_dep_inActive = .false.
end if

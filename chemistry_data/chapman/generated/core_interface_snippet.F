! Auto-generated package activation for: chapman
nullify(chem_chapman_inActive)
call mpas_pool_get_package(packages, 'chem_chapman_inActive', &
                           chem_chapman_inActive)
if (associated(chem_chapman_inActive)) then
   chem_chapman_inActive = .false.
end if

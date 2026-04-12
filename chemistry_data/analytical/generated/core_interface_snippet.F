! Auto-generated package activation for: analytical
nullify(chem_analytical_inActive)
call mpas_pool_get_package(packages, 'chem_analytical_inActive', &
                           chem_analytical_inActive)
if (associated(chem_analytical_inActive)) then
   chem_analytical_inActive = .false.
end if

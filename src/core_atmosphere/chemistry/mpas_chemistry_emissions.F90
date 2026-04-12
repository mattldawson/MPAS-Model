! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Surface emission rate computation for MICM EMISSION reactions.
! Discovers emitting species from MICM species properties at init-time;
! computes volumetric emission rates [mol/m³/s] for the surface level.
!
! Configuration is per-mechanism via species property:
!   "__mpas_surface_emission_flux [molec cm-2 s-1]"
!
module mpas_chemistry_emissions

   use mpas_kind_types,  only : RKIND
   use mpas_log,         only : mpas_log_write
   use iso_fortran_env,  only : real64
   use mpas_chemistry_utils, only : AVOGADRO

   implicit none

   private
   public :: emissions_init, emissions_set_rates, emissions_cleanup

   !> Number of species with surface emissions
   integer, save :: n_emitted = 0

   !> MICM rate_parameters index for each emission reaction (1-based)
   integer, allocatable, save :: emitted_rp_idx(:)

   !> Surface emission flux [molec cm⁻² s⁻¹] for each emission species
   real (kind=real64), allocatable, save :: emitted_flux(:)

contains

   !> Discover emission species from MICM and build rate param mappings.
   !!
   !! For each MICM species with the "__mpas_surface_emission_flux [molec cm-2 s-1]"
   !! property, looks up the corresponding EMIS.<species_name> rate parameter index.
   subroutine emissions_init(micm_solver, micm_state, errmsg, errcode)

      use musica_micm,  only : micm_t
      use musica_state, only : state_t
      use musica_util,  only : error_t

      type(micm_t),  pointer, intent(in) :: micm_solver
      type(state_t), pointer, intent(in) :: micm_state
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errcode

      type(error_t)     :: error
      integer           :: n_species, i, n_emis
      character(len=128) :: species_name
      real (kind=real64) :: flux_val

      ! Temporaries (allocatable, sized to n_species)
      integer, allocatable :: tmp_rp_idx(:)
      real (kind=real64), allocatable :: tmp_flux(:)

      errmsg  = ''
      errcode = 0
      n_emis  = 0

      n_species = micm_state%species_ordering%size()
      allocate(tmp_rp_idx(n_species), tmp_flux(n_species))

      do i = 1, n_species
         species_name = micm_state%species_ordering%name(i)

         ! Check for surface emission property
         flux_val = micm_solver%get_species_property_double( &
            trim(species_name), &
            '__mpas_surface_emission_flux [molec cm-2 s-1]', error)

         if (.not. error%is_success()) cycle  ! No emission for this species

         n_emis = n_emis + 1
         tmp_flux(n_emis) = flux_val

         ! Look up MICM rate parameter: EMIS.<species_name>
         tmp_rp_idx(n_emis) = micm_state%rate_parameters_ordering%index( &
            'EMIS.' // trim(species_name), error)
         if (.not. error%is_success()) then
            errmsg = '[CheMPAS] Cannot find rate param EMIS.' &
                     // trim(species_name) // ': ' // error%message()
            errcode = 1; return
         end if

         call mpas_log_write('[CheMPAS]   Emission: ' // trim(species_name) &
                             // ' (surface flux)')
      end do

      ! Copy to module arrays
      n_emitted = n_emis
      if (n_emis > 0) then
         allocate(emitted_rp_idx(n_emis))
         allocate(emitted_flux(n_emis))
         emitted_rp_idx(1:n_emis) = tmp_rp_idx(1:n_emis)
         emitted_flux(1:n_emis) = tmp_flux(1:n_emis)
      end if

      deallocate(tmp_rp_idx, tmp_flux)

      call mpas_log_write('[CheMPAS] Emissions: $i species', &
                          intArgs=(/n_emitted/))

   end subroutine emissions_init


   !> Set emission rate parameters in MICM state for all grid cells.
   !!
   !! Surface emissions inject into the bottom model level only.
   !! Rate [mol/m³/s] = flux [molec/cm²/s] × 1e4 [cm²/m²] / (N_A [molec/mol] × Δz [m])
   subroutine emissions_set_rates(nCellsSolve, nVertLevels, zgrid, &
                                   rate_params, rp_gc_stride, rp_var_stride)

      integer,           intent(in)    :: nCellsSolve, nVertLevels
      real (kind=RKIND), intent(in)    :: zgrid(:,:)
      real (kind=real64), intent(inout) :: rate_params(:)
      integer,           intent(in)    :: rp_gc_stride, rp_var_stride

      integer :: iCell, k, i_cell, flat_idx, s
      real (kind=real64) :: dz, rate

      if (n_emitted == 0) return

      i_cell = 0
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            i_cell = i_cell + 1
            do s = 1, n_emitted
               flat_idx = (i_cell - 1) * rp_gc_stride &
                        + (emitted_rp_idx(s) - 1) * rp_var_stride + 1
               if (k == nVertLevels) then
                  ! Bottom (surface) level — compute emission rate
                  dz = abs(real(zgrid(nVertLevels, iCell) &
                          - zgrid(nVertLevels + 1, iCell), real64))
                  rate = emitted_flux(s) * 1.0e4_real64 / (AVOGADRO * dz)
                  rate_params(flat_idx) = rate
               else
                  ! All other levels — zero emission
                  rate_params(flat_idx) = 0.0_real64
               end if
            end do
         end do
      end do

   end subroutine emissions_set_rates


   !> Deallocate emission arrays.
   subroutine emissions_cleanup()

      if (allocated(emitted_rp_idx)) deallocate(emitted_rp_idx)
      if (allocated(emitted_flux))   deallocate(emitted_flux)
      n_emitted = 0

   end subroutine emissions_cleanup

end module mpas_chemistry_emissions

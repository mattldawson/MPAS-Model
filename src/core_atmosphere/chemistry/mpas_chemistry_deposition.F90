! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Surface dry deposition rate computation for MICM FIRST_ORDER_LOSS reactions.
! Discovers deposited species from MICM species properties at init-time;
! computes first-order loss rates [s⁻¹] for the surface level.
!
! Configuration is per-mechanism via species property:
!   "__mpas_surface_deposition_velocity [cm s-1]"
!
module mpas_chemistry_deposition

   use mpas_kind_types,  only : RKIND
   use mpas_log,         only : mpas_log_write
   use iso_fortran_env,  only : real64

   implicit none

   private
   public :: deposition_init, deposition_set_rates, deposition_cleanup

   !> Number of species with surface deposition
   integer, save :: n_deposited = 0

   !> MICM rate_parameters index for each deposition reaction (1-based)
   integer, allocatable, save :: deposited_rp_idx(:)

   !> Surface deposition velocity [cm s⁻¹] for each deposited species
   real (kind=real64), allocatable, save :: deposited_vel(:)

contains

   !> Discover deposited species from MICM and build rate param mappings.
   !!
   !! For each MICM species with the "__mpas_surface_deposition_velocity [cm s-1]"
   !! property, looks up the corresponding LOSS.<species_name> rate parameter index.
   subroutine deposition_init(micm_solver, micm_state, errmsg, errcode)

      use musica_micm,  only : micm_t
      use musica_state, only : state_t
      use musica_util,  only : error_t

      type(micm_t),  pointer, intent(in) :: micm_solver
      type(state_t), pointer, intent(in) :: micm_state
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errcode

      type(error_t)     :: error
      integer           :: n_species, i, n_dep
      character(len=128) :: species_name
      real (kind=real64) :: vel_val

      ! Temporaries (allocatable, sized to n_species)
      integer, allocatable :: tmp_rp_idx(:)
      real (kind=real64), allocatable :: tmp_vel(:)

      errmsg  = ''
      errcode = 0
      n_dep   = 0

      n_species = micm_state%species_ordering%size()
      allocate(tmp_rp_idx(n_species), tmp_vel(n_species))

      do i = 1, n_species
         species_name = micm_state%species_ordering%name(i)

         ! Check for deposition velocity property
         vel_val = micm_solver%get_species_property_double( &
            trim(species_name), &
            '__mpas_surface_deposition_velocity [cm s-1]', error)

         if (.not. error%is_success()) cycle  ! No deposition for this species

         n_dep = n_dep + 1
         tmp_vel(n_dep) = vel_val

         ! Look up MICM rate parameter: LOSS.<species_name>
         tmp_rp_idx(n_dep) = micm_state%rate_parameters_ordering%index( &
            'LOSS.' // trim(species_name), error)
         if (.not. error%is_success()) then
            errmsg = '[CheMPAS] Cannot find rate param LOSS.' &
                     // trim(species_name) // ': ' // error%message()
            errcode = 1; return
         end if

         call mpas_log_write('[CheMPAS]   Deposition: ' // trim(species_name) &
                             // ' (surface velocity)')
      end do

      ! Copy to module arrays
      n_deposited = n_dep
      if (n_dep > 0) then
         allocate(deposited_rp_idx(n_dep))
         allocate(deposited_vel(n_dep))
         deposited_rp_idx(1:n_dep) = tmp_rp_idx(1:n_dep)
         deposited_vel(1:n_dep) = tmp_vel(1:n_dep)
      end if

      deallocate(tmp_rp_idx, tmp_vel)

      call mpas_log_write('[CheMPAS] Deposition: $i species', &
                          intArgs=(/n_deposited/))

   end subroutine deposition_init


   !> Set deposition rate parameters in MICM state for all grid cells.
   !!
   !! Surface deposition applies to the bottom model level only.
   !! Rate [s⁻¹] = v_d [cm/s] × 0.01 [m/cm] / Δz [m]
   subroutine deposition_set_rates(nCellsSolve, nVertLevels, zgrid, &
                                    rate_params, rp_gc_stride, rp_var_stride)

      integer,           intent(in)    :: nCellsSolve, nVertLevels
      real (kind=RKIND), intent(in)    :: zgrid(:,:)
      real (kind=real64), intent(inout) :: rate_params(:)
      integer,           intent(in)    :: rp_gc_stride, rp_var_stride

      integer :: iCell, k, i_cell, flat_idx, s
      real (kind=real64) :: dz, rate

      if (n_deposited == 0) return

      i_cell = 0
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            i_cell = i_cell + 1
            do s = 1, n_deposited
               flat_idx = (i_cell - 1) * rp_gc_stride &
                        + (deposited_rp_idx(s) - 1) * rp_var_stride + 1
               if (k == nVertLevels) then
                  ! Bottom (surface) level — compute deposition rate
                  dz = abs(real(zgrid(nVertLevels, iCell) &
                          - zgrid(nVertLevels + 1, iCell), real64))
                  rate = deposited_vel(s) * 0.01_real64 / dz
                  rate_params(flat_idx) = rate
               else
                  ! All other levels — zero deposition
                  rate_params(flat_idx) = 0.0_real64
               end if
            end do
         end do
      end do

   end subroutine deposition_set_rates


   !> Deallocate deposition arrays.
   subroutine deposition_cleanup()

      if (allocated(deposited_rp_idx)) deallocate(deposited_rp_idx)
      if (allocated(deposited_vel))    deallocate(deposited_vel)
      n_deposited = 0

   end subroutine deposition_cleanup

end module mpas_chemistry_deposition

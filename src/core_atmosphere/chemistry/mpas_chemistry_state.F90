! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! State copy utilities: MPAS ↔ MICM data transfer with unit conversion.
! Mechanism-agnostic: uses species arrays from mpas_chemistry_species.
!
module mpas_chemistry_state

   use mpas_kind_types, only : RKIND
   use iso_fortran_env, only : real64
   use mpas_chemistry_species, only : n_advected, advected_mpas_idx, &
                                       advected_micm_idx, advected_molar_mass, &
                                       n_constant, constant_micm_idx, constant_vmr, &
                                       MW_AIR

   implicit none

   private
   public :: update_micm_from_mpas, update_mpas_from_micm

contains

   !> Copy MPAS state into MICM state arrays for a batch of grid cells.
   !!
   !! Iterates generic species arrays for advected and constant species.
   !! Advected species: mmr [kg/kg] * rho [kg/m3] / Mw [kg/mol] → mol/m3
   !! Constant species: VMR * air_conc [mol/m3] → mol/m3
   subroutine update_micm_from_mpas(nCellsSolve, nVertLevels, scalars,       &
                                     temperature, pressure, rho_dry,          &
                                     photo_rates,                             &
                                     n_photo_rxns, photo_mapping,            &
                                     conditions, concentrations, rate_params, &
                                     sp_gc_stride, sp_var_stride,            &
                                     rp_gc_stride, rp_var_stride,            &
                                     offset, batch_size)

      use musica_state, only : conditions_t

      integer,                    intent(in)    :: nCellsSolve, nVertLevels
      real (kind=RKIND),          intent(in)    :: scalars(:,:,:)
      real (kind=RKIND),          intent(in)    :: temperature(:,:)
      real (kind=RKIND),          intent(in)    :: pressure(:,:)
      real (kind=RKIND),          intent(in)    :: rho_dry(:,:)
      real (kind=real64),         intent(in)    :: photo_rates(:,:,:)
      integer,                    intent(in)    :: n_photo_rxns
      integer,                    intent(in)    :: photo_mapping(:)
      type(conditions_t),         intent(inout) :: conditions(:)
      real (kind=real64),         intent(inout) :: concentrations(:)
      real (kind=real64),         intent(inout) :: rate_params(:)
      integer,                    intent(in)    :: sp_gc_stride, sp_var_stride
      integer,                    intent(in)    :: rp_gc_stride, rp_var_stride
      integer,                    intent(in)    :: offset, batch_size

      integer :: iCell, k, i_cell, i_local, flat_idx, r, s
      real (kind=real64) :: rho_d, air_conc

      i_cell = 0
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            i_cell = i_cell + 1
            if (i_cell <= offset .or. i_cell > offset + batch_size) cycle

            i_local = i_cell - offset  ! 1-based local index within batch

            rho_d = real(rho_dry(k, iCell), real64)
            air_conc = rho_d / MW_AIR  ! mol/m3

            ! Conditions
            conditions(i_local)%temperature = real(temperature(k, iCell), real64)
            conditions(i_local)%pressure    = real(pressure(k, iCell), real64)
            conditions(i_local)%air_density = air_conc

            ! Advected species: mmr → mol/m3
            do s = 1, n_advected
               flat_idx = (i_local - 1) * sp_gc_stride &
                        + (advected_micm_idx(s) - 1) * sp_var_stride + 1
               concentrations(flat_idx) = &
                  real(scalars(advected_mpas_idx(s), k, iCell), real64) &
                  * rho_d / advected_molar_mass(s)
            end do

            ! Constant species: VMR * air_conc → mol/m3
            do s = 1, n_constant
               flat_idx = (i_local - 1) * sp_gc_stride &
                        + (constant_micm_idx(s) - 1) * sp_var_stride + 1
               concentrations(flat_idx) = constant_vmr(s) * air_conc
            end do

            ! Rate parameters (photolysis rates)
            do r = 1, n_photo_rxns
               flat_idx = (i_local - 1) * rp_gc_stride &
                        + (photo_mapping(r) - 1) * rp_var_stride + 1
               rate_params(flat_idx) = photo_rates(k, iCell, r)
            end do
         end do
      end do

   end subroutine update_micm_from_mpas


   !> Copy MICM state back to MPAS scalars after chemistry solve.
   !!
   !! Only advected species are copied back (constant species are read-only).
   !! mol/m3 * Mw [kg/mol] / rho [kg/m3] → mmr [kg/kg]
   subroutine update_mpas_from_micm(nCellsSolve, nVertLevels, scalars,       &
                                     rho_dry,                                 &
                                     concentrations,                          &
                                     sp_gc_stride, sp_var_stride,            &
                                     offset, batch_size)

      integer,                    intent(in)    :: nCellsSolve, nVertLevels
      real (kind=RKIND),          intent(inout) :: scalars(:,:,:)
      real (kind=RKIND),          intent(in)    :: rho_dry(:,:)
      real (kind=real64),         intent(in)    :: concentrations(:)
      integer,                    intent(in)    :: sp_gc_stride, sp_var_stride
      integer,                    intent(in)    :: offset, batch_size

      integer :: iCell, k, i_cell, i_local, flat_idx, s
      real (kind=real64) :: rho_d, mmr_val

      i_cell = 0
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            i_cell = i_cell + 1
            if (i_cell <= offset .or. i_cell > offset + batch_size) cycle

            i_local = i_cell - offset

            rho_d = real(rho_dry(k, iCell), real64)

            ! Advected species: mol/m3 → mmr
            do s = 1, n_advected
               flat_idx = (i_local - 1) * sp_gc_stride &
                        + (advected_micm_idx(s) - 1) * sp_var_stride + 1
               mmr_val = concentrations(flat_idx) * advected_molar_mass(s) / rho_d
               scalars(advected_mpas_idx(s), k, iCell) = max(real(mmr_val, RKIND), 0.0_RKIND)
            end do
         end do
      end do

   end subroutine update_mpas_from_micm

end module mpas_chemistry_state

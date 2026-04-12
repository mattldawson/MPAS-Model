! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! State copy utilities: MPAS ↔ MICM data transfer with unit conversion
!
module mpas_chemistry_state

   use mpas_kind_types, only : RKIND
   use iso_fortran_env, only : real64

   implicit none

   private
   public :: update_micm_from_mpas, update_mpas_from_micm

   ! Number of Chapman reactive species
   integer, parameter, public :: N_CHAPMAN_SPECIES = 3

   ! Molar masses [kg/mol]
   real (kind=real64), parameter, public :: MW_O3  = 0.048_real64
   real (kind=real64), parameter, public :: MW_O   = 0.016_real64
   real (kind=real64), parameter, public :: MW_O1D = 0.016_real64
   real (kind=real64), parameter, public :: MW_O2  = 0.032_real64
   real (kind=real64), parameter, public :: MW_N2  = 0.028_real64
   real (kind=real64), parameter, public :: MW_AIR = 0.029_real64

   ! Volume mixing ratios of constant species
   real (kind=real64), parameter, public :: VMR_O2 = 0.2095_real64
   real (kind=real64), parameter, public :: VMR_N2 = 0.7808_real64

contains

   !> Copy MPAS state into MICM state arrays for a batch of grid cells
   !!
   !! @param[in]     nCellsSolve    Number of owned cells
   !! @param[in]     nVertLevels    Number of vertical levels
   !! @param[in]     scalars        MPAS scalars(nScalars, nVertLevels, nCells)
   !! @param[in]     temperature    Temperature(nVertLevels, nCells) [K]
   !! @param[in]     pressure       Pressure(nVertLevels, nCells) [Pa]
   !! @param[in]     rho_dry        Dry air density(nVertLevels, nCells) [kg/m3]
   !! @param[in]     photo_rates    Photolysis rates(nVertLevels, nCells, n_photo) [s-1]
   !! @param[in]     idx_o3, idx_o, idx_o1d  MPAS scalar indices for Chapman species
   !! @param[in]     micm_idx_o3, micm_idx_o, micm_idx_o1d, micm_idx_o2, micm_idx_n2
   !!                               MICM species ordering indices (1-based)
   !! @param[in]     n_photo_rxns   Number of photolysis rate parameters
   !! @param[in]     photo_mapping  Mapping from TUV-x ordering to MICM ordering
   !! @param[in]     conditions     MICM conditions array
   !! @param[inout]  concentrations MICM concentrations flat array
   !! @param[inout]  rate_params    MICM rate_parameters flat array
   !! @param[in]     sp_strides     Species strides (grid_cell, variable)
   !! @param[in]     rp_strides     Rate parameter strides (grid_cell, variable)
   !! @param[in]     offset         Offset into the batch (0-based)
   !! @param[in]     batch_size     Number of grid cells in this batch
   subroutine update_micm_from_mpas(nCellsSolve, nVertLevels, scalars,       &
                                     temperature, pressure, rho_dry,          &
                                     photo_rates,                             &
                                     idx_o3, idx_o, idx_o1d,                  &
                                     micm_idx_o3, micm_idx_o, micm_idx_o1d,  &
                                     micm_idx_o2, micm_idx_n2,               &
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
      integer,                    intent(in)    :: idx_o3, idx_o, idx_o1d
      integer,                    intent(in)    :: micm_idx_o3, micm_idx_o
      integer,                    intent(in)    :: micm_idx_o1d
      integer,                    intent(in)    :: micm_idx_o2, micm_idx_n2
      integer,                    intent(in)    :: n_photo_rxns
      integer,                    intent(in)    :: photo_mapping(:)
      type(conditions_t),         intent(inout) :: conditions(:)
      real (kind=real64),         intent(inout) :: concentrations(:)
      real (kind=real64),         intent(inout) :: rate_params(:)
      integer,                    intent(in)    :: sp_gc_stride, sp_var_stride
      integer,                    intent(in)    :: rp_gc_stride, rp_var_stride
      integer,                    intent(in)    :: offset, batch_size

      integer :: iCell, k, i_cell, i_local, flat_idx, r
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

            ! Species concentrations: mmr [kg/kg] * rho [kg/m3] / Mw [kg/mol] = mol/m3
            flat_idx = (i_local - 1) * sp_gc_stride + (micm_idx_o3 - 1) * sp_var_stride + 1
            concentrations(flat_idx) = real(scalars(idx_o3, k, iCell), real64) * rho_d / MW_O3

            flat_idx = (i_local - 1) * sp_gc_stride + (micm_idx_o - 1) * sp_var_stride + 1
            concentrations(flat_idx) = real(scalars(idx_o, k, iCell), real64) * rho_d / MW_O

            flat_idx = (i_local - 1) * sp_gc_stride + (micm_idx_o1d - 1) * sp_var_stride + 1
            concentrations(flat_idx) = real(scalars(idx_o1d, k, iCell), real64) * rho_d / MW_O1D

            ! Constant species: O2, N2 diagnosed from air density
            flat_idx = (i_local - 1) * sp_gc_stride + (micm_idx_o2 - 1) * sp_var_stride + 1
            concentrations(flat_idx) = VMR_O2 * air_conc

            flat_idx = (i_local - 1) * sp_gc_stride + (micm_idx_n2 - 1) * sp_var_stride + 1
            concentrations(flat_idx) = VMR_N2 * air_conc

            ! Rate parameters (photolysis rates)
            do r = 1, n_photo_rxns
               flat_idx = (i_local - 1) * rp_gc_stride &
                        + (photo_mapping(r) - 1) * rp_var_stride + 1
               rate_params(flat_idx) = photo_rates(k, iCell, r)
            end do
         end do
      end do

   end subroutine update_micm_from_mpas


   !> Copy MICM state back to MPAS scalars after chemistry solve
   !!
   !! @param[in]     nCellsSolve    Number of owned cells
   !! @param[in]     nVertLevels    Number of vertical levels
   !! @param[inout]  scalars        MPAS scalars(nScalars, nVertLevels, nCells)
   !! @param[in]     rho_dry        Dry air density(nVertLevels, nCells) [kg/m3]
   !! @param[in]     idx_o3, idx_o, idx_o1d  MPAS scalar indices
   !! @param[in]     micm_idx_o3, micm_idx_o, micm_idx_o1d  MICM species indices
   !! @param[in]     concentrations MICM concentrations flat array
   !! @param[in]     sp_gc_stride, sp_var_stride  Species strides
   !! @param[in]     offset         Offset into the batch (0-based)
   !! @param[in]     batch_size     Number of grid cells in this batch
   subroutine update_mpas_from_micm(nCellsSolve, nVertLevels, scalars,       &
                                     rho_dry,                                 &
                                     idx_o3, idx_o, idx_o1d,                  &
                                     micm_idx_o3, micm_idx_o, micm_idx_o1d,  &
                                     concentrations,                          &
                                     sp_gc_stride, sp_var_stride,            &
                                     offset, batch_size)

      integer,                    intent(in)    :: nCellsSolve, nVertLevels
      real (kind=RKIND),          intent(inout) :: scalars(:,:,:)
      real (kind=RKIND),          intent(in)    :: rho_dry(:,:)
      integer,                    intent(in)    :: idx_o3, idx_o, idx_o1d
      integer,                    intent(in)    :: micm_idx_o3, micm_idx_o
      integer,                    intent(in)    :: micm_idx_o1d
      real (kind=real64),         intent(in)    :: concentrations(:)
      integer,                    intent(in)    :: sp_gc_stride, sp_var_stride
      integer,                    intent(in)    :: offset, batch_size

      integer :: iCell, k, i_cell, i_local, flat_idx
      real (kind=real64) :: rho_d, mmr_val

      i_cell = 0
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            i_cell = i_cell + 1
            if (i_cell <= offset .or. i_cell > offset + batch_size) cycle

            i_local = i_cell - offset

            rho_d = real(rho_dry(k, iCell), real64)

            ! mol/m3 * Mw [kg/mol] / rho [kg/m3] = mmr [kg/kg]
            flat_idx = (i_local - 1) * sp_gc_stride + (micm_idx_o3 - 1) * sp_var_stride + 1
            mmr_val = concentrations(flat_idx) * MW_O3 / rho_d
            scalars(idx_o3, k, iCell) = max(real(mmr_val, RKIND), 0.0_RKIND)

            flat_idx = (i_local - 1) * sp_gc_stride + (micm_idx_o - 1) * sp_var_stride + 1
            mmr_val = concentrations(flat_idx) * MW_O / rho_d
            scalars(idx_o, k, iCell) = max(real(mmr_val, RKIND), 0.0_RKIND)

            flat_idx = (i_local - 1) * sp_gc_stride + (micm_idx_o1d - 1) * sp_var_stride + 1
            mmr_val = concentrations(flat_idx) * MW_O1D / rho_d
            scalars(idx_o1d, k, iCell) = max(real(mmr_val, RKIND), 0.0_RKIND)
         end do
      end do

   end subroutine update_mpas_from_micm

end module mpas_chemistry_state

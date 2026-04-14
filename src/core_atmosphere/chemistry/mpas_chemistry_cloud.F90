! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Cloud chemistry support module.
! Reads cloud water species, prescribed cloud profile, and default
! concentrations for non-advected aqueous equilibrium species.
! Sets MICM concentrations each timestep for cloud-water and
! equilibrium initial-guess species.
!
module mpas_chemistry_cloud

   use mpas_kind_types, only : RKIND
   use mpas_log,        only : mpas_log_write
   use iso_fortran_env, only : real64

   implicit none

   private
   public :: cloud_init, cloud_set_state, cloud_cleanup

   !> Cloud water species — MICM index (1-based); -1 if not used
   integer, save :: cloud_water_micm_idx = -1

   !> Molar mass of cloud water species [kg/mol]
   real (kind=real64), save :: cloud_water_mw = 0.01802_real64

   !> Prescribed cloud: pressure bounds [Pa] and LWC [kg/kg]
   logical, save :: use_prescribed_cloud = .false.
   real (kind=real64), save :: prescribed_p_top = 0.0_real64
   real (kind=real64), save :: prescribed_p_bot = 0.0_real64
   real (kind=real64), save :: prescribed_lwc   = 0.0_real64

   !> Number of default-concentration species (aqueous equilibrium)
   integer, save :: n_default = 0

   !> MICM species indices for default-concentration species
   integer, allocatable, save :: default_micm_idx(:)

   !> Default concentration values [mol/m3]
   real (kind=real64), allocatable, save :: default_conc(:)

contains

   !> Initialize cloud chemistry from support files.
   subroutine cloud_init(config_dir, micm_state, errmsg, errcode)

      use musica_state, only : state_t
      use musica_util,  only : error_t

      character(len=*),  intent(in)  :: config_dir
      type(state_t), pointer, intent(in) :: micm_state
      character(len=*),  intent(out) :: errmsg
      integer,           intent(out) :: errcode

      character(len=512) :: fpath, line
      character(len=128) :: species_name
      real (kind=real64) :: conc_val
      type(error_t) :: error
      integer :: iunit, ios, idx, n_loaded
      logical :: file_exists

      ! Temporaries for default concentrations
      integer, parameter :: MAX_DEFAULT = 50
      integer :: tmp_idx(MAX_DEFAULT)
      real (kind=real64) :: tmp_conc(MAX_DEFAULT)

      errmsg  = ''
      errcode = 0

      ! --- Read cloud water species name ---
      fpath = trim(config_dir) // '/mpas_cloud_water.txt'
      inquire(file=trim(fpath), exist=file_exists)
      if (.not. file_exists) then
         call mpas_log_write('[CheMPAS-Cloud] No mpas_cloud_water.txt — cloud chemistry disabled')
         return
      end if

      open(newunit=iunit, file=trim(fpath), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         errmsg = '[CheMPAS-Cloud] Cannot open ' // trim(fpath)
         errcode = 1; return
      end if
      read(iunit, '(A)', iostat=ios) line
      close(iunit)
      line = adjustl(line)

      if (len_trim(line) == 0) then
         call mpas_log_write('[CheMPAS-Cloud] mpas_cloud_water.txt is empty — cloud chemistry disabled')
         return
      end if

      ! Look up cloud water species in MICM ordering
      species_name = trim(line)
      idx = micm_state%species_ordering%index(trim(species_name), error)
      if (.not. error%is_success()) then
         errmsg = '[CheMPAS-Cloud] Cloud water species "' // trim(species_name) &
                  // '" not found in MICM species ordering'
         errcode = 1; return
      end if
      cloud_water_micm_idx = idx
      call mpas_log_write('[CheMPAS-Cloud] Cloud water species: ' &
                          // trim(species_name) // ' → MICM index $i', &
                          intArgs=(/idx/))

      ! --- Read prescribed cloud profile ---
      fpath = trim(config_dir) // '/prescribed_cloud.txt'
      inquire(file=trim(fpath), exist=file_exists)
      if (file_exists) then
         open(newunit=iunit, file=trim(fpath), status='old', action='read', iostat=ios)
         if (ios == 0) then
            do
               read(iunit, '(A)', iostat=ios) line
               if (ios /= 0) exit
               line = adjustl(line)
               if (len_trim(line) == 0) cycle
               if (line(1:1) == '#') cycle
               read(line, *, iostat=ios) prescribed_p_top, prescribed_p_bot, prescribed_lwc
               if (ios == 0) then
                  use_prescribed_cloud = .true.
                  call mpas_log_write('[CheMPAS-Cloud] Prescribed cloud: p_top=$r Pa, p_bot=$r Pa, LWC=$r kg/kg', &
                                      realArgs=(/real(prescribed_p_top, RKIND), &
                                                 real(prescribed_p_bot, RKIND), &
                                                 real(prescribed_lwc, RKIND)/))
               end if
               exit
            end do
            close(iunit)
         end if
      end if

      ! --- Read default concentrations for equilibrium species ---
      fpath = trim(config_dir) // '/default_concentrations.txt'
      inquire(file=trim(fpath), exist=file_exists)
      if (.not. file_exists) then
         call mpas_log_write('[CheMPAS-Cloud] No default_concentrations.txt')
         return
      end if

      n_loaded = 0
      open(newunit=iunit, file=trim(fpath), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call mpas_log_write('[CheMPAS-Cloud] Cannot open default_concentrations.txt')
         return
      end if

      do
         read(iunit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle

         read(line, *, iostat=ios) species_name, conc_val
         if (ios /= 0) cycle

         idx = micm_state%species_ordering%index(trim(species_name), error)
         if (.not. error%is_success()) then
            call mpas_log_write('[CheMPAS-Cloud] Default conc species not in MICM: ' &
                                // trim(species_name))
            cycle
         end if

         if (n_loaded >= MAX_DEFAULT) exit
         n_loaded = n_loaded + 1
         tmp_idx(n_loaded) = idx
         tmp_conc(n_loaded) = conc_val
         call mpas_log_write('[CheMPAS-Cloud]   Default: ' // trim(species_name) &
                             // ' → MICM index $i', intArgs=(/idx/))
      end do
      close(iunit)

      n_default = n_loaded
      if (n_loaded > 0) then
         allocate(default_micm_idx(n_loaded))
         allocate(default_conc(n_loaded))
         default_micm_idx(1:n_loaded) = tmp_idx(1:n_loaded)
         default_conc(1:n_loaded) = tmp_conc(1:n_loaded)
      end if

      call mpas_log_write('[CheMPAS-Cloud] Initialized: $i default conc species', &
                          intArgs=(/n_default/))

   end subroutine cloud_init


   !> Set cloud water and default concentrations in MICM state.
   !!
   !! For prescribed cloud: if cell pressure is within cloud bounds,
   !! set cloud water = LWC * rho_dry / MW_H2O; else set to 0.
   !!
   !! Default concentration species are always set to their configured
   !! values as initial guesses for the equilibrium solver.
   subroutine cloud_set_state(nCellsSolve, nVertLevels, pressure, rho_dry, &
                               concentrations, sp_gc_stride, sp_var_stride)

      integer,            intent(in)    :: nCellsSolve, nVertLevels
      real (kind=RKIND),  intent(in)    :: pressure(:,:)
      real (kind=RKIND),  intent(in)    :: rho_dry(:,:)
      real (kind=real64), intent(inout) :: concentrations(:)
      integer,            intent(in)    :: sp_gc_stride, sp_var_stride

      integer :: iCell, k, i_cell, flat_idx, s
      real (kind=real64) :: rho_d, cloud_conc

      if (cloud_water_micm_idx < 1 .and. n_default == 0) return

      i_cell = 0
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            i_cell = i_cell + 1

            ! Set cloud water concentration
            if (cloud_water_micm_idx >= 1) then
               flat_idx = (i_cell - 1) * sp_gc_stride &
                        + (cloud_water_micm_idx - 1) * sp_var_stride + 1
               rho_d = real(rho_dry(k, iCell), real64)

               if (use_prescribed_cloud) then
                  if (real(pressure(k, iCell), real64) >= prescribed_p_top .and. &
                      real(pressure(k, iCell), real64) <= prescribed_p_bot) then
                     ! LWC [kg/kg] * rho [kg/m3] / MW [kg/mol] → mol/m3
                     cloud_conc = prescribed_lwc * rho_d / cloud_water_mw
                  else
                     cloud_conc = 0.0_real64
                  end if
               else
                  ! No prescribed cloud — would need MPAS qc (future work)
                  cloud_conc = 0.0_real64
               end if
               concentrations(flat_idx) = cloud_conc
            end if

            ! Set default concentrations for equilibrium species
            do s = 1, n_default
               flat_idx = (i_cell - 1) * sp_gc_stride &
                        + (default_micm_idx(s) - 1) * sp_var_stride + 1
               concentrations(flat_idx) = default_conc(s)
            end do
         end do
      end do

   end subroutine cloud_set_state


   !> Deallocate module arrays.
   subroutine cloud_cleanup()
      if (allocated(default_micm_idx)) deallocate(default_micm_idx)
      if (allocated(default_conc))     deallocate(default_conc)
      cloud_water_micm_idx = -1
      n_default = 0
      use_prescribed_cloud = .false.
   end subroutine cloud_cleanup

end module mpas_chemistry_cloud

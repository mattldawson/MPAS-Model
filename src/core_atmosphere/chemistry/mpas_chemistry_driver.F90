! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Mechanism-agnostic chemistry driver for MPAS-A.
! Discovers species from MICM at init-time; TUV-x is optional.
!
module mpas_chemistry_driver

   use mpas_kind_types,     only : RKIND, StrKIND
   use mpas_derived_types,  only : domain_type, block_type, mpas_pool_type, &
                                   MPAS_NOW, MPAS_LOG_ERR, MPAS_LOG_CRIT
   use mpas_pool_routines,  only : mpas_pool_get_config, mpas_pool_get_subpool, &
                                   mpas_pool_get_dimension, mpas_pool_get_array
   use mpas_log,            only : mpas_log_write
   use mpas_timekeeping,    only : mpas_get_clock_time, mpas_get_time, &
                                   MPAS_Time_type
   use iso_fortran_env,     only : real64

   use mpas_chemistry_micm,    only : micm_setup, micm_solve, micm_cleanup, &
                                      micm_state, micm_solver_ptr
   use mpas_chemistry_utils,   only : compute_solar_zenith_angle, &
                                      compute_earth_sun_distance
   use mpas_chemistry_state,   only : update_micm_from_mpas, update_mpas_from_micm
   use mpas_chemistry_species, only : chem_species_init, chem_species_cleanup, &
                                      n_advected, advected_mpas_idx, &
                                      advected_micm_idx, advected_molar_mass, &
                                      n_constant, constant_micm_idx, constant_vmr, &
                                      n_tuvx_profiles, tuvx_profiles
   use mpas_chemistry_emissions,  only : emissions_init, emissions_set_rates, &
                                         emissions_cleanup
   use mpas_chemistry_deposition, only : deposition_init, deposition_set_rates, &
                                         deposition_cleanup

   implicit none

   private
   public :: chemistry_init, chemistry_timestep, chemistry_finalize

   ! Module-level state
   logical, save :: chemistry_enabled = .false.
   logical, save :: tuvx_enabled = .false.
   real (kind=real64), save :: chem_dt = 0.0_real64

   ! Photolysis mapping: photo_mapping(r) = MICM rate_parameters index
   ! for TUV-x photolysis reaction r
   integer, allocatable, save :: photo_mapping(:)
   integer, save :: n_photo_rxns_local = 0

   ! R_v / R_d for temperature from moist potential temperature
   real (kind=RKIND), parameter :: rvord = 461.51_RKIND / 287.04_RKIND

contains

   !> Initialize chemistry — called once during model startup.
   subroutine chemistry_init(domain)

      use musica_util, only : error_t

      type(domain_type), intent(inout) :: domain

      type(mpas_pool_type), pointer :: mesh, state
      logical, pointer :: config_chemistry_enabled
      character(len=StrKIND), pointer :: config_chemistry_config_path
      real (kind=RKIND), pointer :: config_chemistry_dt
      real (kind=RKIND), pointer :: config_chemistry_surface_albedo
      integer, pointer :: nCellsSolve, nVertLevels

      character(len=512) :: micm_config_path, tuvx_config_path
      character(len=512) :: errmsg
      logical :: file_exists
      integer :: errcode, r, n_grid_cells
      type(error_t) :: error

      call mpas_pool_get_config(domain % blocklist % configs, &
                                'config_chemistry_enabled', config_chemistry_enabled)
      if (.not. config_chemistry_enabled) then
         call mpas_log_write('[CheMPAS] Chemistry disabled')
         return
      end if

      chemistry_enabled = .true.
      call mpas_log_write('[CheMPAS] Initializing chemistry...')

      call mpas_pool_get_config(domain % blocklist % configs, &
                                'config_chemistry_dt', config_chemistry_dt)
      chem_dt = real(config_chemistry_dt, real64)

      call mpas_pool_get_config(domain % blocklist % configs, &
                                'config_chemistry_config_path', config_chemistry_config_path)
      call mpas_pool_get_config(domain % blocklist % configs, &
                                'config_chemistry_surface_albedo', &
                                config_chemistry_surface_albedo)

      ! Derive sub-paths from config directory
      ! MICM: try v1 (top-level config.json), fall back to v0 (micm/config.json)
      micm_config_path = trim(config_chemistry_config_path) // '/config.json'
      inquire(file=trim(micm_config_path), exist=file_exists)
      if (.not. file_exists) then
         micm_config_path = trim(config_chemistry_config_path) // '/micm/config.json'
      end if

      ! TUV-x: optional, check if tuvx/config.json exists
      tuvx_config_path = trim(config_chemistry_config_path) // '/tuvx/config.json'

      ! Get mesh dimensions from first block
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

      n_grid_cells = nCellsSolve * nVertLevels

      ! --- Setup MICM ---
      call micm_setup(trim(micm_config_path), n_grid_cells, errmsg, errcode)
      if (errcode /= 0) then
         call mpas_log_write(trim(errmsg), messageType=MPAS_LOG_CRIT)
         return
      end if

      call mpas_log_write('[CheMPAS] MICM: $i species, $i rate params, $i grid cells', &
                          intArgs=(/micm_state%number_of_species, &
                                    micm_state%number_of_rate_parameters, &
                                    micm_state%number_of_grid_cells/))

      ! --- Discover species from MICM ---
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call chem_species_init(state, micm_solver_ptr, micm_state, errmsg, errcode)
      if (errcode /= 0) then
         call mpas_log_write(trim(errmsg), messageType=MPAS_LOG_CRIT)
         return
      end if

      ! --- Discover emission/deposition species from MICM ---
      call emissions_init(micm_solver_ptr, micm_state, errmsg, errcode)
      if (errcode /= 0) then
         call mpas_log_write(trim(errmsg), messageType=MPAS_LOG_CRIT)
         return
      end if

      call deposition_init(micm_solver_ptr, micm_state, errmsg, errcode)
      if (errcode /= 0) then
         call mpas_log_write(trim(errmsg), messageType=MPAS_LOG_CRIT)
         return
      end if

      ! --- Setup TUV-x (optional) ---
      inquire(file=trim(tuvx_config_path), exist=tuvx_enabled)
      if (tuvx_enabled) then
         call tuvx_init(tuvx_config_path, nVertLevels, &
                         real(config_chemistry_surface_albedo, real64), &
                         errmsg, errcode)
         if (errcode /= 0) then
            call mpas_log_write(trim(errmsg), messageType=MPAS_LOG_CRIT)
            return
         end if

         ! Build photolysis mapping: TUV-x reaction name → MICM rate param index
         call build_photo_mapping(errmsg, errcode)
         if (errcode /= 0) then
            call mpas_log_write(trim(errmsg), messageType=MPAS_LOG_CRIT)
            return
         end if
      else
         call mpas_log_write('[CheMPAS] TUV-x disabled (no config path)')
         n_photo_rxns_local = 0
         allocate(photo_mapping(0))
      end if

      if (chem_dt <= 0.0_real64) then
         call mpas_log_write('[CheMPAS] config_chemistry_dt <= 0; using dynamics dt')
      end if

      call mpas_log_write('[CheMPAS] Chemistry initialization complete')

   end subroutine chemistry_init


   !> Initialize TUV-x subsystem.
   subroutine tuvx_init(config_path, nVertLevels, surface_albedo, errmsg, errcode)

      use mpas_chemistry_tuvx, only : tuvx_setup, n_photo_rxns, photo_ordering

      character(len=*), intent(in) :: config_path
      integer, pointer, intent(in) :: nVertLevels
      real (kind=real64), intent(in) :: surface_albedo
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errcode

      character(len=64), allocatable :: profile_names(:)
      integer :: i

      ! Build profile name array from tuvx_profiles descriptors
      allocate(profile_names(n_tuvx_profiles))
      do i = 1, n_tuvx_profiles
         profile_names(i) = tuvx_profiles(i)%profile_name
      end do

      call tuvx_setup(trim(config_path), nVertLevels, n_tuvx_profiles, &
                       profile_names, surface_albedo, errmsg, errcode)
      deallocate(profile_names)
      if (errcode /= 0) return

      n_photo_rxns_local = n_photo_rxns
      call mpas_log_write('[CheMPAS] TUV-x: $i photolysis reactions', &
                          intArgs=(/n_photo_rxns_local/))

   end subroutine tuvx_init


   !> Build photolysis rate mapping: TUV-x ordering → MICM rate_parameters.
   subroutine build_photo_mapping(errmsg, errcode)

      use musica_util, only : error_t
      use mpas_chemistry_tuvx, only : n_photo_rxns, photo_ordering

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errcode

      type(error_t) :: error
      integer :: r

      errmsg  = ''
      errcode = 0

      allocate(photo_mapping(n_photo_rxns))
      do r = 1, n_photo_rxns
         photo_mapping(r) = micm_state%rate_parameters_ordering%index( &
            'PHOTO.' // trim(photo_ordering%name(r)), error)
         if (.not. error%is_success()) then
            errmsg = '[CheMPAS] Cannot map photo rxn ' // trim(photo_ordering%name(r)) &
                     // ': ' // error%message()
            errcode = 1; return
         end if
      end do

   end subroutine build_photo_mapping


   !> Run one chemistry timestep — called every dynamics timestep.
   subroutine chemistry_timestep(domain, dt, itimestep)

      type(domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in)    :: dt
      integer,           intent(in)    :: itimestep

      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: mesh, state, diag
      type(MPAS_Time_type) :: currTime

      ! Mesh & state arrays (pointers into MPAS pools)
      real (kind=RKIND), dimension(:,:), pointer :: zgrid
      real (kind=RKIND), dimension(:,:), pointer :: zz
      real (kind=RKIND), dimension(:,:), pointer :: rho_zz
      real (kind=RKIND), dimension(:,:), pointer :: theta_m
      real (kind=RKIND), dimension(:,:), pointer :: exner
      real (kind=RKIND), dimension(:,:), pointer :: pressure_base
      real (kind=RKIND), dimension(:,:), pointer :: pressure_p
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      integer, pointer :: nCellsSolve, nVertLevels, index_qv

      ! Pre-computed 2D working arrays  (allocatable to avoid large stack)
      real (kind=RKIND), allocatable :: temp_2d(:,:)    ! (nVL, nCells)
      real (kind=RKIND), allocatable :: pres_2d(:,:)    ! (nVL, nCells)
      real (kind=RKIND), allocatable :: rho_2d(:,:)     ! (nVL, nCells)

      ! All-cell photolysis rates
      real (kind=real64), allocatable :: photo_all(:,:,:)  ! (nVL, nCells, n_photo)

      real (kind=real64) :: sza_r64, earth_sun_r64, solve_dt
      real (kind=RKIND)  :: sza, earth_sun, qv_k

      character(len=512) :: errmsg
      integer :: errcode, iCell, k, ierr
      integer :: year, julday, hour, minute, second
      real (kind=RKIND) :: ut_hours

      if (.not. chemistry_enabled) return

      ! Chemistry solve time step
      if (chem_dt > 0.0_real64) then
         solve_dt = chem_dt
      else
         solve_dt = real(dt, real64)
      end if

      ! Current simulation time → SZA
      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_get_time(currTime, YYYY=year, DoY=julday, H=hour, M=minute, S=second)
      ut_hours = real(hour, RKIND) + real(minute, RKIND) / 60.0_RKIND &
               + real(second, RKIND) / 3600.0_RKIND

      ! Earth-Sun distance (same for all columns)
      earth_sun = compute_earth_sun_distance(julday)
      earth_sun_r64 = real(earth_sun, real64)

      ! Process each block
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
         call mpas_pool_get_subpool(block % structs, 'state', state)
         call mpas_pool_get_subpool(block % structs, 'diag', diag)

         call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

         call mpas_pool_get_array(mesh, 'zgrid', zgrid)
         call mpas_pool_get_array(mesh, 'zz', zz)
         call mpas_pool_get_array(mesh, 'latCell', latCell)
         call mpas_pool_get_array(mesh, 'lonCell', lonCell)

         call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)
         call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
         call mpas_pool_get_array(state, 'scalars', scalars, 1)

         call mpas_pool_get_array(diag, 'exner', exner)
         call mpas_pool_get_array(diag, 'pressure_base', pressure_base)
         call mpas_pool_get_array(diag, 'pressure_p', pressure_p)

         call mpas_pool_get_dimension(state, 'index_qv', index_qv)

         ! Pre-compute temperature, pressure, dry density for all cells
         allocate(temp_2d(nVertLevels, nCellsSolve))
         allocate(pres_2d(nVertLevels, nCellsSolve))
         allocate(rho_2d(nVertLevels, nCellsSolve))
         allocate(photo_all(nVertLevels, nCellsSolve, max(n_photo_rxns_local, 1)))

         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               rho_2d(k, iCell) = rho_zz(k, iCell) * zz(k, iCell)

               if (associated(index_qv) .and. index_qv > 0) then
                  qv_k = scalars(index_qv, k, iCell)
               else
                  qv_k = 0.0_RKIND
               end if

               temp_2d(k, iCell) = theta_m(k, iCell) * exner(k, iCell) &
                                 / (1.0_RKIND + rvord * qv_k)
               pres_2d(k, iCell) = pressure_base(k, iCell) + pressure_p(k, iCell)
            end do
         end do

         ! === Phase 1: TUV-x photolysis (if enabled) ===
         photo_all(:,:,:) = 0.0_real64

         if (tuvx_enabled .and. n_photo_rxns_local > 0) then
            call run_tuvx_photolysis(nCellsSolve, nVertLevels, scalars, &
                                     temp_2d, rho_2d, zgrid, latCell, lonCell, &
                                     julday, ut_hours, earth_sun_r64, &
                                     photo_all)
         end if

         ! === Phase 2: Fill MICM state from MPAS ===
         call update_micm_from_mpas(nCellsSolve, nVertLevels, scalars,         &
                                     temp_2d, pres_2d, rho_2d,                 &
                                     photo_all,                                &
                                     n_photo_rxns_local, photo_mapping,        &
                                     micm_state%conditions,                    &
                                     micm_state%concentrations,                &
                                     micm_state%rate_parameters,               &
                                     micm_state%species_strides%grid_cell,     &
                                     micm_state%species_strides%variable,      &
                                     micm_state%rate_parameters_strides%grid_cell, &
                                     micm_state%rate_parameters_strides%variable,  &
                                     0, nCellsSolve * nVertLevels)

         ! === Phase 2b: Set emission rate parameters ===
         call emissions_set_rates(nCellsSolve, nVertLevels, zgrid, &
                                  micm_state%rate_parameters, &
                                  micm_state%rate_parameters_strides%grid_cell, &
                                  micm_state%rate_parameters_strides%variable)

         ! === Phase 2c: Set deposition rate parameters ===
         call deposition_set_rates(nCellsSolve, nVertLevels, zgrid, &
                                   micm_state%rate_parameters, &
                                   micm_state%rate_parameters_strides%grid_cell, &
                                   micm_state%rate_parameters_strides%variable)

         ! === Phase 3: Solve chemistry ===
         call micm_solve(solve_dt, errmsg, errcode)
         if (errcode /= 0) then
            call mpas_log_write('[CheMPAS] MICM solve error: ' // trim(errmsg), &
                                messageType=MPAS_LOG_ERR)
         end if

         ! === Phase 4: Copy results back to MPAS ===
         call update_mpas_from_micm(nCellsSolve, nVertLevels, scalars,         &
                                     rho_2d,                                    &
                                     micm_state%concentrations,                 &
                                     micm_state%species_strides%grid_cell,      &
                                     micm_state%species_strides%variable,       &
                                     0, nCellsSolve * nVertLevels)

         deallocate(temp_2d, pres_2d, rho_2d, photo_all)

         block => block % next
      end do

   end subroutine chemistry_timestep


   !> Run TUV-x photolysis for all columns.
   subroutine run_tuvx_photolysis(nCellsSolve, nVertLevels, scalars, &
                                   temp_2d, rho_2d, zgrid, latCell, lonCell, &
                                   julday, ut_hours, earth_sun_r64, &
                                   photo_all)

      use mpas_chemistry_tuvx, only : tuvx_run_column

      integer, intent(in) :: nCellsSolve, nVertLevels
      real (kind=RKIND), intent(in) :: scalars(:,:,:)
      real (kind=RKIND), intent(in) :: temp_2d(:,:)
      real (kind=RKIND), intent(in) :: rho_2d(:,:)
      real (kind=RKIND), intent(in) :: zgrid(:,:)
      real (kind=RKIND), intent(in) :: latCell(:), lonCell(:)
      integer, intent(in) :: julday
      real (kind=RKIND), intent(in) :: ut_hours
      real (kind=real64), intent(in) :: earth_sun_r64
      real (kind=real64), intent(inout) :: photo_all(:,:,:)

      real (kind=RKIND)  :: zgrid_col(nVertLevels+1), rho_col(nVertLevels)
      real (kind=RKIND)  :: temp_col(nVertLevels)
      real (kind=real64) :: photo_col(nVertLevels, n_photo_rxns_local)
      real (kind=real64) :: sza_r64
      real (kind=RKIND)  :: sza
      character(len=512) :: errmsg
      integer :: errcode, iCell, k

      do iCell = 1, nCellsSolve
         sza = compute_solar_zenith_angle(latCell(iCell), lonCell(iCell), &
                                           julday, ut_hours)
         sza_r64 = real(sza, real64)

         rho_col(1:nVertLevels)     = rho_2d(1:nVertLevels, iCell)
         temp_col(1:nVertLevels)    = temp_2d(1:nVertLevels, iCell)
         zgrid_col(1:nVertLevels+1) = zgrid(1:nVertLevels+1, iCell)

         call tuvx_run_column(nVertLevels, zgrid_col, temp_col, &
                               rho_col, scalars(:, :, iCell), &
                               sza_r64, earth_sun_r64, &
                               photo_col, errmsg, errcode)
         if (errcode /= 0) then
            call mpas_log_write('[CheMPAS] TUV-x error cell $i: ' // &
                                trim(errmsg), intArgs=(/iCell/))
         else
            photo_all(1:nVertLevels, iCell, 1:n_photo_rxns_local) = &
               photo_col(1:nVertLevels, 1:n_photo_rxns_local)
         end if
      end do

   end subroutine run_tuvx_photolysis


   !> Finalize chemistry — called at model shutdown.
   subroutine chemistry_finalize()

      use mpas_chemistry_tuvx, only : tuvx_cleanup

      if (.not. chemistry_enabled) return

      call mpas_log_write('[CheMPAS] Finalizing chemistry...')
      call micm_cleanup()
      if (tuvx_enabled) call tuvx_cleanup()
      call chem_species_cleanup()
      call emissions_cleanup()
      call deposition_cleanup()
      if (allocated(photo_mapping)) deallocate(photo_mapping)
      chemistry_enabled = .false.
      call mpas_log_write('[CheMPAS] Chemistry finalized')

   end subroutine chemistry_finalize

end module mpas_chemistry_driver

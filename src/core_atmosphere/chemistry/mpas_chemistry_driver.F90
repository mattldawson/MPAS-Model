! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Chemistry driver for MPAS-A: initializes MICM + TUV-x and orchestrates
! per-timestep photolysis and chemical kinetics for Chapman mechanism.
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

   use mpas_chemistry_micm, only : micm_setup, micm_solve, micm_cleanup, &
                                   micm_state
   use mpas_chemistry_tuvx, only : tuvx_setup, tuvx_run_column, tuvx_cleanup, &
                                   n_photo_rxns, photo_ordering
   use mpas_chemistry_utils, only : compute_solar_zenith_angle, &
                                    compute_earth_sun_distance
   use mpas_chemistry_state, only : update_micm_from_mpas, update_mpas_from_micm

   implicit none

   private
   public :: chemistry_init, chemistry_timestep, chemistry_finalize

   ! Module-level state
   logical, save :: chemistry_enabled = .false.
   real (kind=real64), save :: chem_dt = 0.0_real64

   ! MPAS scalar indices for Chapman species
   integer, save :: idx_o3  = 0
   integer, save :: idx_o   = 0
   integer, save :: idx_o1d = 0

   ! MICM species indices (1-based Fortran)
   integer, save :: micm_idx_o3  = 0
   integer, save :: micm_idx_o   = 0
   integer, save :: micm_idx_o1d = 0
   integer, save :: micm_idx_o2  = 0
   integer, save :: micm_idx_n2  = 0

   ! Photolysis mapping: photo_mapping(r) = MICM rate_parameters index
   ! for TUV-x photolysis reaction r
   integer, allocatable, save :: photo_mapping(:)

   ! R_v / R_d for temperature from moist potential temperature
   real (kind=RKIND), parameter :: rvord = 461.51_RKIND / 287.04_RKIND

contains

   !> Initialize chemistry — called once during model startup.
   subroutine chemistry_init(domain)

      use musica_util, only : error_t

      type(domain_type), intent(inout) :: domain

      type(mpas_pool_type), pointer :: mesh, state
      logical, pointer :: config_chemistry_enabled
      character(len=StrKIND), pointer :: config_micm_config_path
      character(len=StrKIND), pointer :: config_tuvx_config_path
      character(len=StrKIND), pointer :: config_tuvx_micm_mapping_path
      real (kind=RKIND), pointer :: config_chemistry_dt
      integer, pointer :: nCellsSolve, nVertLevels

      character(len=512) :: errmsg
      integer :: errcode, r, n_grid_cells
      integer, pointer :: idx_ptr
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
                                'config_micm_config_path', config_micm_config_path)
      call mpas_pool_get_config(domain % blocklist % configs, &
                                'config_tuvx_config_path', config_tuvx_config_path)
      call mpas_pool_get_config(domain % blocklist % configs, &
                                'config_tuvx_micm_mapping_path', &
                                config_tuvx_micm_mapping_path)

      ! Get mesh dimensions from first block
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

      n_grid_cells = nCellsSolve * nVertLevels

      ! Get MPAS scalar indices for Chapman species
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
      call mpas_pool_get_dimension(state, 'index_o3', idx_ptr)
      idx_o3 = idx_ptr
      call mpas_pool_get_dimension(state, 'index_o', idx_ptr)
      idx_o = idx_ptr
      call mpas_pool_get_dimension(state, 'index_o1d', idx_ptr)
      idx_o1d = idx_ptr

      call mpas_log_write('[CheMPAS] MPAS scalar indices: O3=$i, O=$i, O1D=$i', &
                          intArgs=(/idx_o3, idx_o, idx_o1d/))

      ! --- Setup MICM ---
      call micm_setup(trim(config_micm_config_path), n_grid_cells, errmsg, errcode)
      if (errcode /= 0) then
         call mpas_log_write(trim(errmsg), messageType=MPAS_LOG_CRIT)
         return
      end if

      ! Get MICM species indices
      micm_idx_o3  = micm_state%species_ordering%index('O3', error)
      micm_idx_o   = micm_state%species_ordering%index('O', error)
      micm_idx_o1d = micm_state%species_ordering%index('O1D', error)
      micm_idx_o2  = micm_state%species_ordering%index('O2', error)
      micm_idx_n2  = micm_state%species_ordering%index('N2', error)

      call mpas_log_write('[CheMPAS] MICM species: O3=$i O=$i O1D=$i O2=$i N2=$i', &
                          intArgs=(/micm_idx_o3, micm_idx_o, micm_idx_o1d, &
                                    micm_idx_o2, micm_idx_n2/))
      call mpas_log_write('[CheMPAS] MICM: $i species, $i rate params, $i grid cells', &
                          intArgs=(/micm_state%number_of_species, &
                                    micm_state%number_of_rate_parameters, &
                                    micm_state%number_of_grid_cells/))

      ! --- Setup TUV-x ---
      call tuvx_setup(trim(config_tuvx_config_path), nVertLevels, errmsg, errcode)
      if (errcode /= 0) then
         call mpas_log_write(trim(errmsg), messageType=MPAS_LOG_CRIT)
         return
      end if

      call mpas_log_write('[CheMPAS] TUV-x: $i photolysis reactions', &
                          intArgs=(/n_photo_rxns/))

      ! --- Build photolysis mapping ---
      ! Convention: MICM rate parameters named "PHOTO.<tuvx_reaction_name>"
      allocate(photo_mapping(n_photo_rxns))
      do r = 1, n_photo_rxns
         photo_mapping(r) = micm_state%rate_parameters_ordering%index( &
            'PHOTO.' // trim(photo_ordering%name(r)), error)
         if ((.not. error%is_success())) then
            call mpas_log_write('[CheMPAS] Cannot map photo rxn $i: ' // &
                                error%message(), intArgs=(/r/), &
                                messageType=MPAS_LOG_CRIT)
            return
         end if
      end do

      if (chem_dt <= 0.0_real64) then
         call mpas_log_write('[CheMPAS] config_chemistry_dt <= 0; using dynamics dt')
      end if

      call mpas_log_write('[CheMPAS] Chemistry initialization complete')

   end subroutine chemistry_init


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

      ! Per-column arrays (stack — max supported nVertLevels = 200)
      real (kind=RKIND)  :: o3_col(200), zgrid_col(201)
      real (kind=RKIND)  :: rho_col(200)
      real (kind=real64) :: photo_col(200, 10)

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
         allocate(photo_all(nVertLevels, nCellsSolve, n_photo_rxns))

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

         ! === Phase 1: TUV-x photolysis for each column ===
         photo_all(:,:,:) = 0.0_real64

         do iCell = 1, nCellsSolve
            sza = compute_solar_zenith_angle(latCell(iCell), lonCell(iCell), &
                                              julday, ut_hours)
            sza_r64 = real(sza, real64)

            o3_col(1:nVertLevels)      = scalars(idx_o3, 1:nVertLevels, iCell)
            rho_col(1:nVertLevels)     = rho_2d(1:nVertLevels, iCell)
            zgrid_col(1:nVertLevels+1) = zgrid(1:nVertLevels+1, iCell)

            call tuvx_run_column(nVertLevels, zgrid_col, temp_2d(:,iCell), &
                                  rho_col, o3_col, sza_r64, earth_sun_r64, &
                                  photo_col, errmsg, errcode)
            if (errcode /= 0) then
               call mpas_log_write('[CheMPAS] TUV-x error cell $i: ' // &
                                   trim(errmsg), intArgs=(/iCell/))
            else
               photo_all(1:nVertLevels, iCell, 1:n_photo_rxns) = &
                  photo_col(1:nVertLevels, 1:n_photo_rxns)
            end if
         end do

         ! === Phase 2: Fill MICM state from MPAS ===
         call update_micm_from_mpas(nCellsSolve, nVertLevels, scalars,         &
                                     temp_2d, pres_2d, rho_2d,                 &
                                     photo_all,                                &
                                     idx_o3, idx_o, idx_o1d,                   &
                                     micm_idx_o3, micm_idx_o, micm_idx_o1d,   &
                                     micm_idx_o2, micm_idx_n2,                &
                                     n_photo_rxns, photo_mapping,             &
                                     micm_state%conditions,                    &
                                     micm_state%concentrations,                &
                                     micm_state%rate_parameters,               &
                                     micm_state%species_strides%grid_cell,     &
                                     micm_state%species_strides%variable,      &
                                     micm_state%rate_parameters_strides%grid_cell, &
                                     micm_state%rate_parameters_strides%variable,  &
                                     0, nCellsSolve * nVertLevels)

         ! === Phase 3: Solve chemistry ===
         call micm_solve(solve_dt, errmsg, errcode)
         if (errcode /= 0) then
            call mpas_log_write('[CheMPAS] MICM solve error: ' // trim(errmsg), &
                                messageType=MPAS_LOG_ERR)
         end if

         ! === Phase 4: Copy results back to MPAS ===
         call update_mpas_from_micm(nCellsSolve, nVertLevels, scalars,         &
                                     rho_2d,                                    &
                                     idx_o3, idx_o, idx_o1d,                    &
                                     micm_idx_o3, micm_idx_o, micm_idx_o1d,    &
                                     micm_state%concentrations,                 &
                                     micm_state%species_strides%grid_cell,      &
                                     micm_state%species_strides%variable,       &
                                     0, nCellsSolve * nVertLevels)

         deallocate(temp_2d, pres_2d, rho_2d, photo_all)

         block => block % next
      end do

   end subroutine chemistry_timestep


   !> Finalize chemistry — called at model shutdown.
   subroutine chemistry_finalize()

      if (.not. chemistry_enabled) return

      call mpas_log_write('[CheMPAS] Finalizing chemistry...')
      call micm_cleanup()
      call tuvx_cleanup()
      if (allocated(photo_mapping)) deallocate(photo_mapping)
      chemistry_enabled = .false.
      call mpas_log_write('[CheMPAS] Chemistry finalized')

   end subroutine chemistry_finalize

end module mpas_chemistry_driver

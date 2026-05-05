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

   !> Default concentration values [mol/L water]; converted to mol/m^3 cell
   !! at apply time using the local LWC (see cloud_set_state).
   real (kind=real64), allocatable, save :: default_conc(:)

   !> Liquid water density [kg/m^3] at ~277 K (cloud-water reference).
   !! Used to convert default_conc from mol/L water → mol/m^3 cell.
   real (kind=real64), parameter :: RHO_H2O_KG_PER_M3 = 997.0_real64

   !> Number of ALL aqueous species (for concentration floor)
   integer, save :: n_aqueous = 0

   !> MICM species indices for all aqueous species
   integer, allocatable, save :: aqueous_micm_idx(:)

   !> Aqueous-species name prefix, derived at init time from the
   !! cloud-water species name in mpas_cloud_water.txt (everything up to and
   !! including the final '.'). Empty string disables prefix-based discovery.
   character(len=128), save :: aqueous_prefix = ''

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

      ! Derive the aqueous-species name prefix from the cloud-water species
      ! name: take everything up to and including the final '.' (e.g. for
      ! "CLOUD.AQUEOUS.H2O" the prefix is "CLOUD.AQUEOUS."). All other
      ! aqueous species in the mechanism share this prefix by convention; it
      ! is the discovery key for the concentration floor. No mechanism-
      ! specific string lives in the source code — only in the config file.
      block
         integer :: dot_pos
         dot_pos = index(trim(species_name), '.', back=.true.)
         if (dot_pos > 0) then
            aqueous_prefix = species_name(1:dot_pos)
            call mpas_log_write('[CheMPAS-Cloud] Aqueous species prefix: ' // trim(aqueous_prefix))
         else
            aqueous_prefix = ''
            call mpas_log_write('[CheMPAS-Cloud] Cloud water species name has no "."; ' &
                                // 'aqueous-species discovery disabled')
         end if
      end block

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

      ! --- Discover ALL aqueous species (sharing aqueous_prefix) for floor ---
      call discover_aqueous_species(micm_state)

   end subroutine cloud_init


   !> Find all species whose names start with the aqueous-species prefix
   !! (derived from mpas_cloud_water.txt) and store their MICM indices for
   !! the concentration floor.
   subroutine discover_aqueous_species(micm_state)

      use musica_state, only : state_t

      type(state_t), pointer, intent(in) :: micm_state

      integer, parameter :: MAX_AQ = 200
      integer :: tmp_aq(MAX_AQ)
      integer :: ns, i, aq_count, plen
      character(len=256) :: sname

      plen = len_trim(aqueous_prefix)
      if (plen == 0) then
         n_aqueous = 0
         call mpas_log_write('[CheMPAS-Cloud] No aqueous prefix configured; floor disabled')
         return
      end if

      ns = micm_state%species_ordering%size()
      aq_count = 0
      do i = 1, ns
         sname = micm_state%species_ordering%name(i)
         if (len_trim(sname) >= plen) then
            if (sname(1:plen) == aqueous_prefix(1:plen)) then
               if (aq_count < MAX_AQ) then
                  aq_count = aq_count + 1
                  tmp_aq(aq_count) = micm_state%species_ordering%index(i)
               end if
            end if
         end if
      end do
      n_aqueous = aq_count
      if (aq_count > 0) then
         allocate(aqueous_micm_idx(aq_count))
         aqueous_micm_idx(1:aq_count) = tmp_aq(1:aq_count)
      end if
      call mpas_log_write('[CheMPAS-Cloud] Found $i aqueous species for floor', &
                          intArgs=(/n_aqueous/))

   end subroutine discover_aqueous_species


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
      logical, save :: diag_printed = .false.
      logical, save :: first_call = .true.
      integer :: n_cloud_cells
      real (kind=real64) :: max_cloud_conc
      logical :: in_cloud

      if (cloud_water_micm_idx < 1 .and. n_default == 0) return

      ! After the first chemistry timestep, cloud water and all aqueous
      ! condensed-phase species are advected by MPAS (see advected_species.txt).
      ! We only seed them once at the very start of the run; thereafter the
      ! values flowing in from MPAS scalar transport are authoritative and
      ! must not be reset (resetting them every step throws away dissolved
      ! S(IV) accumulated during chemistry, suppressing SO4 production).
      ! All we still do per-step is floor any aqueous species to a tiny
      ! positive value to avoid division-by-zero in dissolved-reaction rates.
      if (.not. first_call) then
         i_cell = 0
         do iCell = 1, nCellsSolve
            do k = 1, nVertLevels
               i_cell = i_cell + 1
               do s = 1, n_aqueous
                  flat_idx = (i_cell - 1) * sp_gc_stride &
                           + (aqueous_micm_idx(s) - 1) * sp_var_stride + 1
                  concentrations(flat_idx) = max(concentrations(flat_idx), 1.0e-30_real64)
               end do
            end do
         end do
         return
      end if

      n_cloud_cells = 0
      max_cloud_conc = 0.0_real64

      i_cell = 0
      do iCell = 1, nCellsSolve
         do k = 1, nVertLevels
            i_cell = i_cell + 1
            in_cloud = .false.

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
                     n_cloud_cells = n_cloud_cells + 1
                     if (cloud_conc > max_cloud_conc) max_cloud_conc = cloud_conc
                     in_cloud = .true.
                  else
                     cloud_conc = 0.0_real64
                  end if
               else
                  ! No prescribed cloud — would need MPAS qc (future work)
                  cloud_conc = 0.0_real64
               end if
               ! Floor to tiny value to prevent division-by-zero in
               ! dissolved reactions (rate / solvent^n)
               concentrations(flat_idx) = max(cloud_conc, 1.0e-30_real64)
            end if

            ! Set default concentrations ONLY in cloud cells.
            ! In non-cloud cells (H2O≈1e-30), defaults like Hp=1e-4
            ! are wildly inconsistent with dissolved equilibria and
            ! prevent constraint initialization from converging.
            !
            ! default_conc is in mol/L water (literature units). Convert
            ! to mol/m^3 cell using the local LWC volume fraction:
            !   conc_cell = conc_per_L * cloud_conc * MW_H2O / rho_H2O * 1000
            ! This makes seed values physically meaningful (e.g. Hp=1e-4
            ! mol/L → pH 4) regardless of LWC, instead of producing
            ! pH<1 in thin clouds and triggering constraint-solver
            ! divergence (Phase 8 ts1_cloud blowup root cause).
            if (in_cloud) then
               ! mol_per_L_water → mol_per_m3_cell scale factor
               ! = LWC_volume_fraction * 1000 L/m^3
               ! = (cloud_conc * MW_H2O / rho_H2O) * 1000
               do s = 1, n_default
                  flat_idx = (i_cell - 1) * sp_gc_stride &
                           + (default_micm_idx(s) - 1) * sp_var_stride + 1
                  concentrations(flat_idx) = default_conc(s) &
                       * cloud_conc * cloud_water_mw &
                       / RHO_H2O_KG_PER_M3 * 1000.0_real64
               end do
            end if

            ! Floor ALL aqueous species to avoid division-by-zero
            ! in dissolved reactions and equilibrium constraints
            do s = 1, n_aqueous
               flat_idx = (i_cell - 1) * sp_gc_stride &
                        + (aqueous_micm_idx(s) - 1) * sp_var_stride + 1
               concentrations(flat_idx) = max(concentrations(flat_idx), 1.0e-30_real64)
            end do
         end do
      end do

      if (.not. diag_printed) then
         diag_printed = .true.
         call mpas_log_write('[CheMPAS-Cloud] DIAG: cloud_water_micm_idx=$i, n_default=$i', &
                             intArgs=(/cloud_water_micm_idx, n_default/))
         call mpas_log_write('[CheMPAS-Cloud] DIAG: use_prescribed_cloud=$l, LWC=$r', &
                             logicArgs=(/use_prescribed_cloud/), &
                             realArgs=(/real(prescribed_lwc, RKIND)/))
         call mpas_log_write('[CheMPAS-Cloud] DIAG: p_top=$r Pa, p_bot=$r Pa', &
                             realArgs=(/real(prescribed_p_top, RKIND), &
                                        real(prescribed_p_bot, RKIND)/))
         call mpas_log_write('[CheMPAS-Cloud] DIAG: n_cloud_cells=$i / $i total, max_cloud_conc=$r mol/m3', &
                             intArgs=(/n_cloud_cells, nCellsSolve * nVertLevels/), &
                             realArgs=(/real(max_cloud_conc, RKIND)/))
         call mpas_log_write('[CheMPAS-Cloud] DIAG: sp_gc_stride=$i, sp_var_stride=$i', &
                             intArgs=(/sp_gc_stride, sp_var_stride/))
      end if

      ! Mark first-call seeding as done — subsequent calls only floor.
      first_call = .false.

   end subroutine cloud_set_state


   !> Deallocate module arrays.
   subroutine cloud_cleanup()
      if (allocated(default_micm_idx)) deallocate(default_micm_idx)
      if (allocated(default_conc))     deallocate(default_conc)
      if (allocated(aqueous_micm_idx)) deallocate(aqueous_micm_idx)
      cloud_water_micm_idx = -1
      n_default = 0
      n_aqueous = 0
      use_prescribed_cloud = .false.
   end subroutine cloud_cleanup

end module mpas_chemistry_cloud

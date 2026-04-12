! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! TUV-x interface for MPAS-A chemistry coupling (Chapman mechanism).
! Creates TUV-x instance with height + wavelength grids and profiles
! from MPAS host code, runs per-column photolysis rate calculations.
! Follows the CAM-SIMA / atmospheric_physics integration pattern:
!   grids and profiles are created by the host, not in the JSON config.
!
module mpas_chemistry_tuvx

   use mpas_kind_types,  only : RKIND
   use iso_fortran_env,  only : real64
   use musica_tuvx,      only : tuvx_t, grid_t, profile_t, &
                                grid_map_t, profile_map_t, radiator_map_t
   use musica_util,      only : error_t, mappings_t

   implicit none

   private
   public :: tuvx_setup, tuvx_run_column, tuvx_cleanup
   public :: n_photo_rxns, photo_ordering

   ! Module-level TUV-x objects (persistent across timesteps)
   type(tuvx_t),     pointer :: tuvx_solver               => null()
   type(grid_t),     pointer :: height_grid                => null()
   type(grid_t),     pointer :: wavelength_grid            => null()
   type(profile_t),  pointer :: temperature_profile        => null()
   type(profile_t),  pointer :: dry_air_profile            => null()
   type(profile_t),  pointer :: o2_profile                 => null()
   type(profile_t),  pointer :: o3_profile                 => null()
   type(profile_t),  pointer :: surface_albedo_profile     => null()
   type(profile_t),  pointer :: et_flux_profile            => null()
   type(mappings_t), pointer :: photo_ordering             => null()

   integer :: n_photo_rxns = 0
   integer :: n_tuvx_layers = 0    ! = nVertLevels + 1

   ! Wavelength grid: 102 bins, 103 edges (120-750 nm)
   ! Same grid as CAM-SIMA / atmospheric_physics for TUV-x compatibility
   integer, parameter :: N_WAVELENGTH_BINS = 102

   ! Conversion constants
   real (kind=real64), parameter :: km_to_cm  = 1.0e5_real64
   real (kind=real64), parameter :: m3_to_cm3 = 1.0e6_real64
   real (kind=real64), parameter :: avogadro  = 6.02214076e23_real64
   real (kind=real64), parameter :: pi64      = 3.14159265358979323846_real64

   ! Gas species scale heights for exo-layer [km]
   real (kind=real64), parameter :: SCALE_HEIGHT_AIR = 8.01_real64
   real (kind=real64), parameter :: SCALE_HEIGHT_O2  = 8.01_real64
   real (kind=real64), parameter :: SCALE_HEIGHT_O3  = 4.5_real64

   ! Max SZA for photolysis calculations [degrees]
   real (kind=real64), parameter :: MAX_SZA_DEG = 110.0_real64

   ! Default surface albedo for JW baroclinic wave test (no land model)
   real (kind=real64), parameter :: DEFAULT_SURFACE_ALBEDO = 0.1_real64

contains

   !> Initialize TUV-x for Chapman photolysis with MPAS grid.
   !! Creates height + wavelength grids and all profiles from host code,
   !! following the CAM-SIMA / atmospheric_physics integration pattern.
   subroutine tuvx_setup(config_path, nVertLevels, errmsg, errcode)

      character(len=*), intent(in)  :: config_path
      integer,          intent(in)  :: nVertLevels
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errcode

      type(grid_map_t),      pointer :: grids     => null()
      type(profile_map_t),   pointer :: profiles  => null()
      type(radiator_map_t),  pointer :: radiators => null()
      type(error_t) :: error

      ! Wavelength grid: 103 edges (102 bins), 120-750 nm
      ! Same grid as CAM-SIMA (from atmospheric_physics photolysis lookup tables)
      real (kind=real64), target :: wl_edges(N_WAVELENGTH_BINS + 1)
      real (kind=real64), target :: wl_midpoints(N_WAVELENGTH_BINS)

      ! Extraterrestrial solar flux [photon cm-2 s-1 nm-1] on the 102 bins
      real (kind=real64) :: et_flux_per_nm(N_WAVELENGTH_BINS)
      real (kind=real64), target :: et_flux_per_bin(N_WAVELENGTH_BINS)

      ! Surface albedo (uniform across wavelength bins)
      real (kind=real64), target :: albedo_edges(N_WAVELENGTH_BINS + 1)

      integer :: i

      errmsg  = ''
      errcode = 0
      n_tuvx_layers = nVertLevels + 1  ! half-layer at surface + exo-layer at top

      ! ===== Create height grid from host =====
      height_grid => grid_t("height", "km", n_tuvx_layers, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to create height grid: ' // error%message()
         errcode = 1; return
      end if

      ! ===== Create wavelength grid from host =====
      call get_wavelength_edges(wl_edges)
      do i = 1, N_WAVELENGTH_BINS
         wl_midpoints(i) = 0.5_real64 * (wl_edges(i) + wl_edges(i+1))
      end do

      wavelength_grid => grid_t("wavelength", "nm", N_WAVELENGTH_BINS, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to create wavelength grid: ' // error%message()
         errcode = 1; return
      end if
      call wavelength_grid%set_edges(wl_edges, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to set wavelength edges: ' // error%message()
         errcode = 1; return
      end if
      call wavelength_grid%set_midpoints(wl_midpoints, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to set wavelength midpoints: ' // error%message()
         errcode = 1; return
      end if

      ! ===== Populate grid map =====
      grids => grid_map_t(error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to create grid map: ' // error%message()
         errcode = 1; return
      end if
      call grids%add(height_grid, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to add height grid: ' // error%message()
         errcode = 1; return
      end if
      call grids%add(wavelength_grid, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to add wavelength grid: ' // error%message()
         errcode = 1; return
      end if

      ! ===== Create profiles from host =====
      profiles => profile_map_t(error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to create profile map: ' // error%message()
         errcode = 1; return
      end if

      ! --- Temperature profile (on height grid) ---
      temperature_profile => profile_t("temperature", "K", height_grid, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to create temperature profile: ' // error%message()
         errcode = 1; return
      end if
      call profiles%add(temperature_profile, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to add temperature profile: ' // error%message()
         errcode = 1; return
      end if

      ! --- Dry air profile (on height grid) ---
      dry_air_profile => profile_t("air", "molecule cm-3", height_grid, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to create dry air profile: ' // error%message()
         errcode = 1; return
      end if
      call profiles%add(dry_air_profile, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to add dry air profile: ' // error%message()
         errcode = 1; return
      end if

      ! --- O2 profile (on height grid) ---
      o2_profile => profile_t("O2", "molecule cm-3", height_grid, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to create O2 profile: ' // error%message()
         errcode = 1; return
      end if
      call profiles%add(o2_profile, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to add O2 profile: ' // error%message()
         errcode = 1; return
      end if

      ! --- O3 profile (on height grid) ---
      o3_profile => profile_t("O3", "molecule cm-3", height_grid, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to create O3 profile: ' // error%message()
         errcode = 1; return
      end if
      call profiles%add(o3_profile, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to add O3 profile: ' // error%message()
         errcode = 1; return
      end if

      ! --- Surface albedo profile (on wavelength grid) ---
      surface_albedo_profile => profile_t("surface albedo", "none", &
                                          wavelength_grid, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to create surface albedo profile: ' // error%message()
         errcode = 1; return
      end if
      call profiles%add(surface_albedo_profile, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to add surface albedo profile: ' // error%message()
         errcode = 1; return
      end if

      ! --- Extraterrestrial flux profile (on wavelength grid) ---
      et_flux_profile => profile_t("extraterrestrial flux", "photon cm-2 s-1", &
                                   wavelength_grid, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to create ET flux profile: ' // error%message()
         errcode = 1; return
      end if
      call profiles%add(et_flux_profile, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to add ET flux profile: ' // error%message()
         errcode = 1; return
      end if

      ! ===== Radiator map (empty -- no clouds/aerosols for JW case) =====
      radiators => radiator_map_t(error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to create radiator map: ' // error%message()
         errcode = 1; return
      end if

      ! ===== Construct TUV-x solver =====
      tuvx_solver => tuvx_t(trim(config_path), grids, profiles, radiators, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to create TUV-x solver: ' // error%message()
         errcode = 1; return
      end if

      ! TUV-x owns the data now -- release construction-time maps
      deallocate(grids);    nullify(grids)
      deallocate(profiles); nullify(profiles)
      deallocate(radiators); nullify(radiators)

      ! ===== Re-obtain handles from TUV-x =====
      grids => tuvx_solver%get_grids(error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to get TUV-x grids: ' // error%message()
         errcode = 1; return
      end if
      height_grid => grids%get("height", "km", error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to get height grid: ' // error%message()
         errcode = 1; return
      end if
      wavelength_grid => grids%get("wavelength", "nm", error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to get wavelength grid: ' // error%message()
         errcode = 1; return
      end if
      deallocate(grids); nullify(grids)

      profiles => tuvx_solver%get_profiles(error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to get TUV-x profiles: ' // error%message()
         errcode = 1; return
      end if
      temperature_profile => profiles%get("temperature", "K", error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to get temperature profile: ' // error%message()
         errcode = 1; return
      end if
      dry_air_profile => profiles%get("air", "molecule cm-3", error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to get air profile: ' // error%message()
         errcode = 1; return
      end if
      o2_profile => profiles%get("O2", "molecule cm-3", error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to get O2 profile: ' // error%message()
         errcode = 1; return
      end if
      o3_profile => profiles%get("O3", "molecule cm-3", error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to get O3 profile: ' // error%message()
         errcode = 1; return
      end if
      surface_albedo_profile => profiles%get("surface albedo", "none", error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to get surface albedo profile: ' // error%message()
         errcode = 1; return
      end if
      et_flux_profile => profiles%get("extraterrestrial flux", &
                                      "photon cm-2 s-1", error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to get ET flux profile: ' // error%message()
         errcode = 1; return
      end if
      deallocate(profiles); nullify(profiles)

      ! ===== Set surface albedo (constant for JW test) =====
      albedo_edges(:) = DEFAULT_SURFACE_ALBEDO
      call surface_albedo_profile%set_edge_values(albedo_edges, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to set surface albedo: ' // error%message()
         errcode = 1; return
      end if

      ! ===== Set extraterrestrial flux =====
      ! Flux per nm -> flux per bin (multiply by bin width)
      call get_extraterrestrial_flux(et_flux_per_nm)
      do i = 1, N_WAVELENGTH_BINS
         et_flux_per_bin(i) = et_flux_per_nm(i) * (wl_edges(i+1) - wl_edges(i))
      end do
      call et_flux_profile%set_midpoint_values(et_flux_per_bin, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to set ET flux: ' // error%message()
         errcode = 1; return
      end if

      ! ===== Get photolysis rate ordering =====
      photo_ordering => tuvx_solver%get_photolysis_rate_constants_ordering(error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] Failed to get photolysis ordering: ' // error%message()
         errcode = 1; return
      end if
      n_photo_rxns = photo_ordering%size()

   end subroutine tuvx_setup


   !> Run TUV-x for one MPAS column -- compute photolysis rates.
   subroutine tuvx_run_column(nVertLevels, zgrid, temperature, rho_dry, &
                               o3_mmr, sza, earth_sun_dist,             &
                               photo_rates, errmsg, errcode)
      use mpas_chemistry_state, only : MW_AIR, MW_O3, MW_O2, VMR_O2

      integer,             intent(in)  :: nVertLevels
      real (kind=RKIND),   intent(in)  :: zgrid(:)           ! (nVertLevels+1) [m]
      real (kind=RKIND),   intent(in)  :: temperature(:)     ! (nVertLevels) [K]
      real (kind=RKIND),   intent(in)  :: rho_dry(:)         ! (nVertLevels) [kg/m3]
      real (kind=RKIND),   intent(in)  :: o3_mmr(:)          ! (nVertLevels) [kg/kg]
      real (kind=real64),  intent(in)  :: sza                ! radians
      real (kind=real64),  intent(in)  :: earth_sun_dist     ! AU
      real (kind=real64),  intent(out) :: photo_rates(:,:)   ! (nVertLevels, n_photo_rxns)
      character(len=*),    intent(out) :: errmsg
      integer,             intent(out) :: errcode

      ! TUV-x grid: nVertLevels+2 interfaces, nVertLevels+1 layers
      real (kind=real64), target :: tuvx_interfaces(nVertLevels + 2)
      real (kind=real64), target :: tuvx_midpoints(nVertLevels + 1)
      real (kind=real64) :: height_deltas(nVertLevels + 1)

      ! TUV-x profiles
      real (kind=real64), target :: temp_edges(nVertLevels + 2)

      ! TUV-x output
      real (kind=real64), target :: photo_out(nVertLevels + 2, n_photo_rxns)
      real (kind=real64), target :: heating_out(nVertLevels + 2, n_photo_rxns)

      real (kind=real64) :: zmid_km(nVertLevels)
      real (kind=real64) :: zint_km(nVertLevels + 1)
      real (kind=real64) :: sza_deg
      real (kind=RKIND)  :: ones(nVertLevels), vmr_o2_arr(nVertLevels)
      type(error_t) :: error
      integer :: k

      errmsg  = ''
      errcode = 0

      ! Check if dark (SZA > threshold)
      sza_deg = sza * 180.0_real64 / pi64
      if (sza_deg > MAX_SZA_DEG .or. sza_deg < 0.0_real64) then
         photo_rates(:,:) = 0.0_real64
         return
      end if

      ! Heights in km -- MPAS zgrid is surface-to-top (same as TUV-x)
      do k = 1, nVertLevels + 1
         zint_km(k) = real(zgrid(k), real64) * 1.0e-3_real64
      end do
      do k = 1, nVertLevels
         zmid_km(k) = 0.5_real64 * (zint_km(k) + zint_km(k+1))
      end do

      ! TUV-x interfaces = [surface, MPAS midpoints, model top]
      tuvx_interfaces(1) = zint_km(1)
      do k = 1, nVertLevels
         tuvx_interfaces(k+1) = zmid_km(k)
      end do
      tuvx_interfaces(nVertLevels+2) = zint_km(nVertLevels+1)

      ! TUV-x midpoints
      tuvx_midpoints(1) = 0.5_real64 * (tuvx_interfaces(1) + tuvx_interfaces(2))
      do k = 2, nVertLevels
         tuvx_midpoints(k) = zint_km(k)  ! MPAS interfaces -> TUV-x midpoints
      end do
      tuvx_midpoints(nVertLevels+1) = 0.5_real64 * &
         (tuvx_interfaces(nVertLevels+1) + tuvx_interfaces(nVertLevels+2))

      ! Height deltas [km]
      do k = 1, nVertLevels + 1
         height_deltas(k) = tuvx_interfaces(k+1) - tuvx_interfaces(k)
      end do

      ! Set height grid
      call height_grid%set_edges(tuvx_interfaces, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] set height edges: ' // error%message()
         errcode = 1; return
      end if
      call height_grid%set_midpoints(tuvx_midpoints, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] set height midpoints: ' // error%message()
         errcode = 1; return
      end if

      ! Set temperature profile (edge values)
      temp_edges(1) = real(temperature(1), real64)
      do k = 1, nVertLevels
         temp_edges(k+1) = real(temperature(k), real64)
      end do
      temp_edges(nVertLevels+2) = real(temperature(nVertLevels), real64)
      call temperature_profile%set_edge_values(temp_edges, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] set temperature: ' // error%message()
         errcode = 1; return
      end if

      ! Set gas species profiles
      ! Air: mmr=1 (total dry air)
      ones(:) = 1.0_RKIND
      call set_gas_profile(dry_air_profile, nVertLevels, rho_dry, &
                           ones, MW_AIR, height_deltas, &
                           SCALE_HEIGHT_AIR, .false., errmsg, errcode)
      if (errcode /= 0) return

      ! O2: convert VMR to MMR then use standard formula
      ! MMR = VMR * MW_species / MW_air
      vmr_o2_arr(:) = real(VMR_O2 * MW_O2 / MW_AIR, RKIND)
      call set_gas_profile(o2_profile, nVertLevels, rho_dry, &
                           vmr_o2_arr, MW_O2, height_deltas, &
                           SCALE_HEIGHT_O2, .false., errmsg, errcode)
      if (errcode /= 0) return

      ! O3: from MPAS scalars
      call set_gas_profile(o3_profile, nVertLevels, rho_dry, &
                           o3_mmr, MW_O3, height_deltas, &
                           SCALE_HEIGHT_O3, .true., errmsg, errcode)
      if (errcode /= 0) return

      ! Run TUV-x
      photo_out(:,:)   = 0.0_real64
      heating_out(:,:) = 0.0_real64

      call tuvx_solver%run(sza, earth_sun_dist, photo_out, heating_out, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] TUV-x run failed: ' // error%message()
         errcode = 1; return
      end if

      ! Map TUV-x layers to MPAS midpoints (layer k -> MPAS level k)
      ! Layer nVL+1 (exo-layer) is not used.
      do k = 1, nVertLevels
         photo_rates(k, :) = max(photo_out(k, :), 0.0_real64)
      end do

   end subroutine tuvx_run_column


   !> Set a gas species profile from MPAS data.
   !! Converts mmr -> molecule/cm3 edge values and column densities.
   subroutine set_gas_profile(prof, nVertLevels, rho_dry, mmr, molar_mass, &
                               height_deltas, scale_height, is_o3, errmsg, errcode)

      type(profile_t),    intent(inout) :: prof
      integer,            intent(in)    :: nVertLevels
      real (kind=RKIND),  intent(in)    :: rho_dry(:)       ! (nVertLevels) [kg/m3]
      real (kind=RKIND),  intent(in)    :: mmr(:)           ! (nVertLevels) [kg/kg]
      real (kind=real64), intent(in)    :: molar_mass       ! [kg/mol]
      real (kind=real64), intent(in)    :: height_deltas(:) ! (nVertLevels+1) [km]
      real (kind=real64), intent(in)    :: scale_height     ! [km]
      logical,            intent(in)    :: is_o3
      character(len=*),   intent(out)   :: errmsg
      integer,            intent(out)   :: errcode

      real (kind=real64), target :: edges(nVertLevels + 2)
      real (kind=real64), target :: densities(nVertLevels + 1)
      real (kind=real64) :: conc(nVertLevels)
      type(error_t) :: error
      integer :: k

      errmsg  = ''
      errcode = 0

      ! mmr -> molecule/cm3:
      ! n [molec/cm3] = mmr * rho [kg/m3] / M [kg/mol] * N_A / 1e6
      do k = 1, nVertLevels
         conc(k) = real(mmr(k), real64) * real(rho_dry(k), real64) &
                 / molar_mass * avogadro / m3_to_cm3
      end do

      ! Edge values (surface-to-top, no inversion)
      edges(1) = conc(1)
      do k = 1, nVertLevels
         edges(k+1) = conc(k)
      end do
      edges(nVertLevels+2) = conc(nVertLevels)

      ! Layer column densities [molecule/cm2]
      if (is_o3) then
         ! Arithmetic mean (better for rapidly varying species)
         do k = 1, nVertLevels + 1
            densities(k) = height_deltas(k) * km_to_cm * &
                           0.5_real64 * (edges(k) + edges(k+1))
         end do
      else
         ! Geometric mean (standard for well-mixed species)
         do k = 1, nVertLevels + 1
            densities(k) = height_deltas(k) * km_to_cm * &
                           sqrt(max(edges(k), 0.0_real64)) * &
                           sqrt(max(edges(k+1), 0.0_real64))
         end do
      end if

      call prof%set_edge_values(edges, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] gas profile edges: ' // error%message()
         errcode = 1; return
      end if

      call prof%set_layer_densities(densities, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] gas profile densities: ' // error%message()
         errcode = 1; return
      end if

      call prof%calculate_exo_layer_density(scale_height, error)
      if ((.not. error%is_success())) then
         errmsg = '[CheMPAS] exo-layer density: ' // error%message()
         errcode = 1; return
      end if

   end subroutine set_gas_profile


   !> Clean up TUV-x resources.
   subroutine tuvx_cleanup()

      if (associated(photo_ordering)) then
         deallocate(photo_ordering); nullify(photo_ordering)
      end if
      if (associated(height_grid)) then
         deallocate(height_grid); nullify(height_grid)
      end if
      if (associated(wavelength_grid)) then
         deallocate(wavelength_grid); nullify(wavelength_grid)
      end if
      if (associated(temperature_profile)) then
         deallocate(temperature_profile); nullify(temperature_profile)
      end if
      if (associated(dry_air_profile)) then
         deallocate(dry_air_profile); nullify(dry_air_profile)
      end if
      if (associated(o2_profile)) then
         deallocate(o2_profile); nullify(o2_profile)
      end if
      if (associated(o3_profile)) then
         deallocate(o3_profile); nullify(o3_profile)
      end if
      if (associated(surface_albedo_profile)) then
         deallocate(surface_albedo_profile); nullify(surface_albedo_profile)
      end if
      if (associated(et_flux_profile)) then
         deallocate(et_flux_profile); nullify(et_flux_profile)
      end if
      if (associated(tuvx_solver)) then
         deallocate(tuvx_solver); nullify(tuvx_solver)
      end if

   end subroutine tuvx_cleanup


   !> Wavelength grid edges [nm] -- 103 interfaces for 102 bins.
   !! Matches the CAM-SIMA / atmospheric_physics photolysis wavelength grid.
   subroutine get_wavelength_edges(edges)
      real (kind=real64), intent(out) :: edges(N_WAVELENGTH_BINS + 1)

      edges = (/ &
         120.0_real64, 121.4_real64, 121.9_real64, 123.5_real64, 124.3_real64, &
         125.5_real64, 126.3_real64, 127.1_real64, 130.1_real64, 131.1_real64, &
         135.0_real64, 140.0_real64, 145.0_real64, 150.0_real64, 155.0_real64, &
         160.0_real64, 165.0_real64, 168.0_real64, 171.0_real64, 173.0_real64, &
         174.4_real64, 175.4_real64, 177.0_real64, 178.6_real64, 180.2_real64, &
         181.8_real64, 183.5_real64, 185.2_real64, 186.9_real64, 188.7_real64, &
         190.5_real64, 192.3_real64, 194.2_real64, 196.1_real64, 198.0_real64, &
         200.0_real64, 202.0_real64, 204.1_real64, 206.2_real64, 208.0_real64, &
         211.0_real64, 214.0_real64, 217.0_real64, 220.0_real64, 223.0_real64, &
         226.0_real64, 229.0_real64, 232.0_real64, 235.0_real64, 238.0_real64, &
         241.0_real64, 244.0_real64, 247.0_real64, 250.0_real64, 253.0_real64, &
         256.0_real64, 259.0_real64, 263.0_real64, 267.0_real64, 271.0_real64, &
         275.0_real64, 279.0_real64, 283.0_real64, 287.0_real64, 291.0_real64, &
         295.0_real64, 298.5_real64, 302.5_real64, 305.5_real64, 308.5_real64, &
         311.5_real64, 314.5_real64, 317.5_real64, 322.5_real64, 327.5_real64, &
         332.5_real64, 337.5_real64, 342.5_real64, 347.5_real64, 350.0_real64, &
         355.0_real64, 360.0_real64, 365.0_real64, 370.0_real64, 375.0_real64, &
         380.0_real64, 385.0_real64, 390.0_real64, 395.0_real64, 400.0_real64, &
         405.0_real64, 410.0_real64, 415.0_real64, 420.0_real64, 430.0_real64, &
         440.0_real64, 450.0_real64, 500.0_real64, 550.0_real64, 600.0_real64, &
         650.0_real64, 700.0_real64, 750.0_real64 /)

   end subroutine get_wavelength_edges


   !> Extraterrestrial solar flux [photon cm-2 s-1 nm-1] on the 102 bins.
   !! Standard solar spectrum matching the CAM-SIMA wavelength grid.
   subroutine get_extraterrestrial_flux(flux)
      real (kind=real64), intent(out) :: flux(N_WAVELENGTH_BINS)

      flux = (/ &
         1.67555e+14_real64, 1.61720e+14_real64, 1.54341e+14_real64, &
         1.52118e+14_real64, 1.58469e+14_real64, 1.66036e+14_real64, &
         1.29769e+14_real64, 1.47316e+14_real64, 1.56274e+14_real64, &
         1.60325e+14_real64, 1.70225e+14_real64, 1.61097e+14_real64, &
         1.71742e+14_real64, 1.70514e+14_real64, 1.24876e+14_real64, &
         1.68376e+14_real64, 1.60490e+14_real64, 1.57979e+14_real64, &
         1.66505e+14_real64, 1.52378e+14_real64, 1.97811e+14_real64, &
         1.75901e+14_real64, 1.54751e+14_real64, 1.98793e+14_real64, &
         2.02567e+14_real64, 1.89537e+14_real64, 1.68532e+14_real64, &
         1.60546e+14_real64, 1.13346e+14_real64, 2.05967e+14_real64, &
         1.77864e+14_real64, 1.62888e+14_real64, 2.14804e+14_real64, &
         1.75501e+14_real64, 1.86444e+14_real64, 2.32774e+14_real64, &
         2.30828e+14_real64, 2.24982e+14_real64, 2.02219e+14_real64, &
         2.47972e+14_real64, 2.00891e+14_real64, 2.44825e+14_real64, &
         1.99913e+14_real64, 1.57765e+14_real64, 1.65729e+14_real64, &
         2.16080e+14_real64, 2.09073e+14_real64, 2.45738e+14_real64, &
         2.55917e+14_real64, 1.91251e+14_real64, 2.47308e+14_real64, &
         2.10889e+14_real64, 1.41346e+14_real64, 1.32255e+14_real64, &
         1.98997e+14_real64, 1.85347e+14_real64, 2.08594e+14_real64, &
         1.88650e+14_real64, 1.78574e+14_real64, 2.41000e+14_real64, &
         2.40832e+14_real64, 2.75942e+14_real64, 1.88910e+14_real64, &
         9.70730e+13_real64, 2.19069e+14_real64, 2.74779e+14_real64, &
         1.29952e+14_real64, 2.08327e+14_real64, 3.08963e+14_real64, &
         3.33270e+14_real64, 3.32895e+14_real64, 3.63642e+14_real64, &
         3.65966e+14_real64, 3.37217e+14_real64, 3.26648e+14_real64, &
         3.41746e+14_real64, 3.32766e+14_real64, 3.17377e+14_real64, &
         3.75741e+14_real64, 3.52127e+14_real64, 3.10827e+14_real64, &
         3.77464e+14_real64, 3.72567e+14_real64, 3.66392e+14_real64, &
         3.63313e+14_real64, 3.63562e+14_real64, 3.87293e+14_real64, &
         3.50809e+14_real64, 3.55652e+14_real64, 3.60092e+14_real64, &
         3.73018e+14_real64, 3.82393e+14_real64, 3.37355e+14_real64, &
         3.65658e+14_real64, 3.78703e+14_real64, 3.63957e+14_real64, &
         3.65457e+14_real64, 3.38550e+14_real64, 3.43225e+14_real64, &
         3.19808e+14_real64, 2.46645e+14_real64, 3.67134e+14_real64 /)

   end subroutine get_extraterrestrial_flux

end module mpas_chemistry_tuvx

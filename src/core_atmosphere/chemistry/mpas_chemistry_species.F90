! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Runtime species discovery for mechanism-agnostic chemistry driver.
! Queries MICM for species names and properties, builds MPAS↔MICM
! index mappings and molar mass tables.
!
module mpas_chemistry_species

   use mpas_kind_types,    only : RKIND, StrKIND
   use mpas_derived_types, only : mpas_pool_type, MPAS_LOG_ERR, MPAS_LOG_CRIT
   use mpas_pool_routines, only : mpas_pool_get_dimension
   use mpas_log,           only : mpas_log_write
   use iso_fortran_env,    only : real64

   implicit none

   private
   public :: chem_species_init, chem_species_cleanup
   public :: n_advected, advected_mpas_idx, advected_micm_idx, advected_molar_mass
   public :: n_constant, constant_micm_idx, constant_vmr
   public :: tuvx_o3_mpas_idx

   !> Number of advected chemistry species (get MPAS tracers)
   integer, save :: n_advected = 0

   !> MPAS scalar indices for advected species
   integer, allocatable, save :: advected_mpas_idx(:)

   !> MICM species ordering indices for advected species (1-based)
   integer, allocatable, save :: advected_micm_idx(:)

   !> Molar masses [kg/mol] for advected species
   real (kind=real64), allocatable, save :: advected_molar_mass(:)

   !> Number of constant species (diagnosed from air composition)
   integer, save :: n_constant = 0

   !> MICM species ordering indices for constant species
   integer, allocatable, save :: constant_micm_idx(:)

   !> Volume mixing ratios for constant species
   real (kind=real64), allocatable, save :: constant_vmr(:)

   !> MPAS scalar index for O3 (needed for TUV-x profile); 0 if no O3
   integer, save :: tuvx_o3_mpas_idx = 0

   ! Molar mass of dry air [kg/mol]
   real (kind=real64), parameter, public :: MW_AIR = 0.029_real64

contains

   !> Discover species from MICM ordering and build index mappings.
   !!
   !! For each species in MICM's species_ordering:
   !!   - If it has a "__mpas_constant_vmr" property → constant species
   !!   - Otherwise → advected species (needs MPAS scalar)
   !!   - Molar mass from "molecular weight [kg mol-1]" property
   !!   - MPAS index from pool dimension "index_<lowercase_name>"
   subroutine chem_species_init(state_pool, micm_solver, micm_state, errmsg, errcode)

      use musica_micm,  only : micm_t
      use musica_state, only : state_t
      use musica_util,  only : error_t, mappings_t

      type(mpas_pool_type), pointer, intent(in) :: state_pool
      type(micm_t),         pointer, intent(in) :: micm_solver
      type(state_t),        pointer, intent(in) :: micm_state
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errcode

      type(error_t) :: error
      integer :: n_species, i, mpas_idx
      integer, pointer :: idx_ptr
      character(len=128) :: species_name, mpas_name
      real (kind=real64) :: vmr_val, mw_val

      ! Temporary arrays (sized to max possible)
      integer :: tmp_adv_mpas(100), tmp_adv_micm(100)
      real (kind=real64) :: tmp_adv_mw(100)
      integer :: tmp_con_micm(100)
      real (kind=real64) :: tmp_con_vmr(100)
      integer :: n_adv, n_con

      errmsg  = ''
      errcode = 0

      n_species = micm_state%species_ordering%size()
      call mpas_log_write('[CheMPAS] Discovering $i MICM species...', &
                          intArgs=(/n_species/))

      n_adv = 0
      n_con = 0
      tuvx_o3_mpas_idx = 0

      do i = 1, n_species
         species_name = micm_state%species_ordering%name(i)

         ! Check if this is a constant species (has __mpas_constant_vmr)
         vmr_val = micm_solver%get_species_property_double( &
            trim(species_name), '__mpas_constant_vmr', error)

         if (error%is_success()) then
            ! Constant species — diagnosed from VMR
            n_con = n_con + 1
            if (n_con > 100) then
               errmsg = '[CheMPAS] Too many constant species (max 100)'
               errcode = 1; return
            end if
            tmp_con_micm(n_con) = micm_state%species_ordering%index( &
               trim(species_name), error)
            tmp_con_vmr(n_con) = vmr_val

            call mpas_log_write('[CheMPAS]   Constant: ' // trim(species_name))
            cycle
         end if

         ! Query molar mass — if absent, skip (shouldn't happen)
         mw_val = micm_solver%get_species_property_double( &
            trim(species_name), 'molecular weight [kg mol-1]', error)
         if (.not. error%is_success()) then
            call mpas_log_write('[CheMPAS]   Skipping ' // trim(species_name) &
                                // ' (no molar mass)')
            cycle
         end if

         ! Advected species — look up MPAS scalar index
         mpas_name = to_lower(trim(species_name))
         nullify(idx_ptr)
         call mpas_pool_get_dimension(state_pool, 'index_' // trim(mpas_name), idx_ptr)
         if (.not. associated(idx_ptr)) then
            errmsg = '[CheMPAS] MICM species "' // trim(species_name) &
                     // '" has no MPAS scalar (index_' // trim(mpas_name) &
                     // '). Check Registry.xml packages.'
            errcode = 1; return
         end if
         mpas_idx = idx_ptr
         if (mpas_idx <= 0) then
            errmsg = '[CheMPAS] MICM species "' // trim(species_name) &
                     // '": MPAS index is 0 (inactive package?)'
            errcode = 1; return
         end if

         n_adv = n_adv + 1
         if (n_adv > 100) then
            errmsg = '[CheMPAS] Too many advected species (max 100)'
            errcode = 1; return
         end if
         tmp_adv_mpas(n_adv) = mpas_idx
         tmp_adv_micm(n_adv) = micm_state%species_ordering%index( &
            trim(species_name), error)
         tmp_adv_mw(n_adv) = mw_val

         ! Track O3 for TUV-x
         if (trim(species_name) == 'O3') then
            tuvx_o3_mpas_idx = mpas_idx
         end if

         call mpas_log_write('[CheMPAS]   Advected: ' // trim(species_name) &
                             // ' → index_' // trim(mpas_name))
      end do

      ! Copy to module arrays
      n_advected = n_adv
      if (n_adv > 0) then
         allocate(advected_mpas_idx(n_adv))
         allocate(advected_micm_idx(n_adv))
         allocate(advected_molar_mass(n_adv))
         advected_mpas_idx(1:n_adv) = tmp_adv_mpas(1:n_adv)
         advected_micm_idx(1:n_adv) = tmp_adv_micm(1:n_adv)
         advected_molar_mass(1:n_adv) = tmp_adv_mw(1:n_adv)
      end if

      n_constant = n_con
      if (n_con > 0) then
         allocate(constant_micm_idx(n_con))
         allocate(constant_vmr(n_con))
         constant_micm_idx(1:n_con) = tmp_con_micm(1:n_con)
         constant_vmr(1:n_con) = tmp_con_vmr(1:n_con)
      end if

      call mpas_log_write('[CheMPAS] Species: $i advected, $i constant', &
                          intArgs=(/n_advected, n_constant/))

   end subroutine chem_species_init


   !> Deallocate species arrays.
   subroutine chem_species_cleanup()

      if (allocated(advected_mpas_idx))  deallocate(advected_mpas_idx)
      if (allocated(advected_micm_idx))  deallocate(advected_micm_idx)
      if (allocated(advected_molar_mass)) deallocate(advected_molar_mass)
      if (allocated(constant_micm_idx))  deallocate(constant_micm_idx)
      if (allocated(constant_vmr))       deallocate(constant_vmr)
      n_advected = 0
      n_constant = 0
      tuvx_o3_mpas_idx = 0

   end subroutine chem_species_cleanup


   !> Convert a string to lowercase.
   pure function to_lower(str) result(lower_str)
      character(len=*), intent(in) :: str
      character(len=len(str))      :: lower_str
      integer :: i, ic
      lower_str = str
      do i = 1, len(str)
         ic = ichar(str(i:i))
         if (ic >= ichar('A') .and. ic <= ichar('Z')) then
            lower_str(i:i) = char(ic + 32)
         end if
      end do
   end function to_lower

end module mpas_chemistry_species

! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Runtime chemistry scalar allocation.
! Reads a species list from an advected_species.txt file and extends
! the Registry-defined scalars (and scalars_tend) var_arrays to include
! chemistry species.  This module uses only MPAS framework APIs so it
! can be compiled by both core_atmosphere and core_init_atmosphere.
!
module mpas_chemistry_allocate

   use mpas_kind_types,    only : RKIND, StrKIND
   use mpas_derived_types, only : block_type, mpas_pool_type, field3DReal, &
                                  MPAS_LOG_ERR, MPAS_LOG_CRIT, &
                                  MPAS_POOL_SILENT, att_lists_type
   use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_subpool, &
                                  mpas_pool_get_dimension, mpas_pool_add_dimension, &
                                  mpas_pool_get_field, &
                                  mpas_pool_get_error_level, mpas_pool_set_error_level
   use mpas_attlist,       only : mpas_add_att
   use mpas_log,           only : mpas_log_write

   implicit none

   private
   public :: chemistry_allocate_scalars

   integer, parameter :: MAX_SPECIES = 500
   integer, parameter :: MAX_NAME_LEN = 128

contains

   !> Extend the scalars (and optionally scalars_tend) var_arrays with
   !! chemistry species read from an advected_species.txt file.
   !!
   !! This must be called from setup_block AFTER generate_structs but
   !! BEFORE any I/O or data initialization.
   !!
   !! @param[inout] block       The block whose pools are being extended
   !! @param[in]    extend_tend If .true., also extend the scalars_tend var_array
   !! @param[out]   errmsg      Error message (empty on success)
   !! @param[out]   errcode     0 on success
   subroutine chemistry_allocate_scalars(block, extend_tend, errmsg, errcode)

      type(block_type), pointer, intent(inout) :: block
      logical, intent(in)    :: extend_tend
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errcode

      character(len=StrKIND), pointer :: config_path
      logical, pointer :: chem_enabled
      character(len=512) :: species_file
      character(len=MAX_NAME_LEN) :: species_names(MAX_SPECIES)
      integer :: n_species, err_level
      logical :: file_exists

      errmsg  = ''
      errcode = 0

      ! Check if chemistry config path is set (non-empty means chemistry requested)
      nullify(config_path)
      err_level = mpas_pool_get_error_level()
      call mpas_pool_set_error_level(MPAS_POOL_SILENT)
      call mpas_pool_get_config(block % domain % configs, &
                                'config_chemistry_config_path', config_path)
      call mpas_pool_set_error_level(err_level)

      if (.not. associated(config_path)) return
      if (len_trim(config_path) == 0) return

      ! Build path to species list file
      species_file = trim(config_path) // '/advected_species.txt'
      inquire(file=trim(species_file), exist=file_exists)
      if (.not. file_exists) then
         ! No species file — nothing to extend (mechanism may have no advected species)
         return
      end if

      ! Read species names
      call read_species_file(trim(species_file), species_names, n_species, errmsg, errcode)
      if (errcode /= 0) return
      if (n_species == 0) return

      call mpas_log_write('[CheMPAS] Extending scalars with $i chemistry species from ' &
                          // trim(species_file), intArgs=(/n_species/))

      ! Extend the scalars var_array in the state pool
      call extend_var_array(block, 'state', 'scalars', 'num_scalars', &
                            species_names, n_species, 2, errmsg, errcode)
      if (errcode /= 0) return

      ! Optionally extend scalars_tend in the tend pool
      if (extend_tend) then
         call extend_var_array(block, 'tend', 'scalars_tend', 'num_scalars_tend', &
                               species_names, n_species, 1, errmsg, errcode)
         if (errcode /= 0) return
      end if

   end subroutine chemistry_allocate_scalars


   !> Read species names from a text file (one name per line).
   subroutine read_species_file(filepath, names, n, errmsg, errcode)

      character(len=*), intent(in) :: filepath
      character(len=MAX_NAME_LEN), intent(out) :: names(MAX_SPECIES)
      integer, intent(out) :: n
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errcode

      integer :: iunit, ios
      character(len=MAX_NAME_LEN) :: line

      errmsg  = ''
      errcode = 0
      n = 0

      open(newunit=iunit, file=filepath, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         errmsg = '[CheMPAS] Cannot open species file: ' // trim(filepath)
         errcode = 1
         return
      end if

      do
         read(iunit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle           ! skip blank lines
         if (line(1:1) == '#') cycle               ! skip comments
         n = n + 1
         if (n > MAX_SPECIES) then
            errmsg = '[CheMPAS] Too many species (max 500)'
            errcode = 1
            close(iunit)
            return
         end if
         names(n) = trim(line)
      end do

      close(iunit)

   end subroutine read_species_file


   !> Extend an existing var_array field with additional constituents.
   !!
   !! For each new species:
   !!   - A dimension 'index_<name>' is added to the pool
   !!   - The constituentNames and attLists arrays are reallocated
   !!   - The num dimension is updated in-place
   subroutine extend_var_array(block, pool_name, field_name, dim_name, &
                               species_names, n_species, n_time_levels, &
                               errmsg, errcode)

      type(block_type), pointer, intent(inout) :: block
      character(len=*), intent(in) :: pool_name
      character(len=*), intent(in) :: field_name
      character(len=*), intent(in) :: dim_name
      character(len=MAX_NAME_LEN), intent(in) :: species_names(MAX_SPECIES)
      integer, intent(in) :: n_species
      integer, intent(in) :: n_time_levels
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errcode

      type(mpas_pool_type), pointer :: sub_pool
      type(field3DReal), pointer :: fld
      integer, pointer :: num_dim_ptr
      integer :: old_num, new_num, i, t, idx
      character(len=StrKIND), dimension(:), pointer :: old_names, new_names
      type(att_lists_type), dimension(:), pointer :: old_atts, new_atts
      character(len=MAX_NAME_LEN) :: dim_key

      errmsg  = ''
      errcode = 0

      nullify(sub_pool)
      call mpas_pool_get_subpool(block % structs, pool_name, sub_pool)
      if (.not. associated(sub_pool)) then
         errmsg = '[CheMPAS] Pool not found: ' // trim(pool_name)
         errcode = 1
         return
      end if

      ! Get current dimension size
      nullify(num_dim_ptr)
      call mpas_pool_get_dimension(sub_pool, dim_name, num_dim_ptr)
      if (.not. associated(num_dim_ptr)) then
         errmsg = '[CheMPAS] Dimension not found: ' // trim(dim_name)
         errcode = 1
         return
      end if
      old_num = num_dim_ptr
      new_num = old_num + n_species

      ! Add index dimensions for each new species
      do i = 1, n_species
         dim_key = 'index_' // trim(species_names(i))
         idx = old_num + i
         call mpas_pool_add_dimension(sub_pool, trim(dim_key), idx)
      end do

      ! Update the num dimension in-place
      num_dim_ptr = new_num

      ! Extend constituentNames and attLists for each time level
      do t = 1, n_time_levels
         nullify(fld)
         call mpas_pool_get_field(sub_pool, field_name, fld, t)
         if (.not. associated(fld)) then
            errmsg = '[CheMPAS] Field not found: ' // trim(field_name) &
                     // ' time level ' // char(ichar('0') + t)
            errcode = 1
            return
         end if

         ! --- Extend constituentNames ---
         old_names => fld % constituentNames
         allocate(new_names(new_num))
         if (associated(old_names)) then
            new_names(1:old_num) = old_names(1:old_num)
            deallocate(old_names)
         end if
         do i = 1, n_species
            new_names(old_num + i) = trim(species_names(i))
         end do
         fld % constituentNames => new_names

         ! --- Extend attLists ---
         old_atts => fld % attLists
         allocate(new_atts(new_num))
         if (associated(old_atts)) then
            ! Move existing attList pointers (no deep copy needed)
            do i = 1, old_num
               new_atts(i) % attList => old_atts(i) % attList
            end do
            deallocate(old_atts)
         end if
         ! Initialize new attLists
         do i = 1, n_species
            idx = old_num + i
            allocate(new_atts(idx) % attList)
            call mpas_add_att(new_atts(idx) % attList, 'units', 'kg kg^{-1}')
            call mpas_add_att(new_atts(idx) % attList, 'long_name', &
                              'Chemistry species: ' // trim(species_names(i)))
         end do
         fld % attLists => new_atts

      end do

      call mpas_log_write('[CheMPAS] Extended ' // trim(field_name) // ': $i → $i constituents', &
                          intArgs=(/old_num, new_num/))

   end subroutine extend_var_array

end module mpas_chemistry_allocate

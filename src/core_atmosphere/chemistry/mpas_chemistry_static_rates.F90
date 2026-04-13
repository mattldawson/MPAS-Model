! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Static rate parameters for USER_DEFINED reactions not computed by TUV-x.
! Reads a CSV file at init-time and sets constant rate values every timestep.
!
! CSV format:  USER.<name>,<rate_value>
! Lines starting with '#' are comments.
!
module mpas_chemistry_static_rates

   use mpas_kind_types, only : RKIND
   use mpas_log,        only : mpas_log_write
   use iso_fortran_env, only : real64

   implicit none

   private
   public :: static_rates_init, static_rates_set, static_rates_cleanup

   !> Number of static rate parameters loaded
   integer, save :: n_static = 0

   !> MICM rate_parameters index for each static rate (1-based)
   integer, allocatable, save :: static_rp_idx(:)

   !> Constant rate value for each static rate parameter
   real (kind=real64), allocatable, save :: static_rate_val(:)

contains

   !> Read static rate parameters from CSV and look up MICM indices.
   !!
   !! Skips entries already mapped by TUV-x photolysis (they are set
   !! dynamically each timestep and would be overwritten anyway).
   subroutine static_rates_init(config_dir, micm_state, &
                                 n_photo_rxns, photo_mapping, &
                                 errmsg, errcode)

      use musica_state, only : state_t
      use musica_util,  only : error_t

      character(len=*),  intent(in)  :: config_dir
      type(state_t), pointer, intent(in) :: micm_state
      integer,           intent(in)  :: n_photo_rxns
      integer,           intent(in)  :: photo_mapping(:)
      character(len=*),  intent(out) :: errmsg
      integer,           intent(out) :: errcode

      character(len=512) :: csv_path, line
      character(len=128) :: param_name
      real (kind=real64) :: rate_val
      type(error_t) :: error
      integer :: iunit, ios, n_loaded, n_skipped, idx, comma_pos, i
      logical :: file_exists

      ! Temporaries
      integer, allocatable :: tmp_idx(:)
      real (kind=real64), allocatable :: tmp_val(:)

      errmsg  = ''
      errcode = 0
      n_loaded = 0
      n_skipped = 0

      csv_path = trim(config_dir) // '/static_rate_params.csv'
      inquire(file=trim(csv_path), exist=file_exists)
      if (.not. file_exists) then
         call mpas_log_write('[CheMPAS] No static_rate_params.csv — skipping')
         return
      end if

      ! Count non-comment, non-blank lines for allocation
      open(newunit=iunit, file=trim(csv_path), status='old', &
           action='read', iostat=ios)
      if (ios /= 0) then
         errmsg = '[CheMPAS] Cannot open ' // trim(csv_path)
         errcode = 1; return
      end if

      i = 0
      do
         read(iunit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         i = i + 1
      end do
      close(iunit)

      if (i == 0) then
         call mpas_log_write('[CheMPAS] static_rate_params.csv is empty')
         return
      end if

      allocate(tmp_idx(i), tmp_val(i))

      ! Re-read and parse
      open(newunit=iunit, file=trim(csv_path), status='old', &
           action='read', iostat=ios)

      do
         read(iunit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle

         ! Parse: USER.<name>,<value>
         comma_pos = index(line, ',')
         if (comma_pos < 2) cycle
         param_name = line(1:comma_pos-1)
         read(line(comma_pos+1:), *, iostat=ios) rate_val
         if (ios /= 0) cycle

         ! Look up in MICM rate_parameters ordering
         idx = micm_state%rate_parameters_ordering%index( &
            trim(param_name), error)
         if (.not. error%is_success()) then
            ! Not in this mechanism — skip silently
            n_skipped = n_skipped + 1
            cycle
         end if

         ! Skip if this index is already mapped by TUV-x photolysis
         if (is_tuvx_mapped(idx, n_photo_rxns, photo_mapping)) then
            n_skipped = n_skipped + 1
            cycle
         end if

         n_loaded = n_loaded + 1
         tmp_idx(n_loaded) = idx
         tmp_val(n_loaded) = rate_val
      end do
      close(iunit)

      ! Copy to module-level arrays
      n_static = n_loaded
      if (n_loaded > 0) then
         allocate(static_rp_idx(n_loaded))
         allocate(static_rate_val(n_loaded))
         static_rp_idx(1:n_loaded) = tmp_idx(1:n_loaded)
         static_rate_val(1:n_loaded) = tmp_val(1:n_loaded)
      end if

      deallocate(tmp_idx, tmp_val)

      call mpas_log_write('[CheMPAS] Static rate params: $i loaded, $i skipped', &
                          intArgs=(/n_loaded, n_skipped/))

   end subroutine static_rates_init


   !> Set static rate parameters in MICM state for all grid cells.
   !! These are constant (same value at every level and cell).
   subroutine static_rates_set(n_grid_cells, rate_params, &
                                rp_gc_stride, rp_var_stride)

      integer,            intent(in)    :: n_grid_cells
      real (kind=real64), intent(inout) :: rate_params(:)
      integer,            intent(in)    :: rp_gc_stride, rp_var_stride

      integer :: gc, s, flat_idx

      if (n_static == 0) return

      do gc = 1, n_grid_cells
         do s = 1, n_static
            flat_idx = (gc - 1) * rp_gc_stride &
                     + (static_rp_idx(s) - 1) * rp_var_stride + 1
            rate_params(flat_idx) = static_rate_val(s)
         end do
      end do

   end subroutine static_rates_set


   !> Check if a rate parameter index is mapped by TUV-x photolysis.
   logical function is_tuvx_mapped(idx, n_photo, photo_mapping)
      integer, intent(in) :: idx, n_photo
      integer, intent(in) :: photo_mapping(:)
      integer :: r

      is_tuvx_mapped = .false.
      do r = 1, n_photo
         if (photo_mapping(r) == idx) then
            is_tuvx_mapped = .true.
            return
         end if
      end do
   end function is_tuvx_mapped


   !> Deallocate module arrays.
   subroutine static_rates_cleanup()
      if (allocated(static_rp_idx))  deallocate(static_rp_idx)
      if (allocated(static_rate_val)) deallocate(static_rate_val)
      n_static = 0
   end subroutine static_rates_cleanup

end module mpas_chemistry_static_rates

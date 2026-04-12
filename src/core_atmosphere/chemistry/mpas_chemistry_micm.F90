! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! MICM interface for MPAS-A chemistry coupling (Chapman mechanism).
! Creates MICM solver instance, manages state, and solves chemistry in batches.
!
module mpas_chemistry_micm

   use mpas_kind_types,  only : RKIND
   use iso_fortran_env,  only : real64
   use musica_micm,      only : micm_t, RosenbrockStandardOrder, solver_stats_t
   use musica_state,     only : state_t, conditions_t
   use musica_util,      only : error_t, mappings_t, string_t

   implicit none

   private
   public :: micm_setup, micm_solve, micm_cleanup
   public :: micm_state, n_micm_species, n_micm_rate_params

   ! Module-level MICM objects
   type(micm_t),  pointer :: micm_solver => null()
   type(state_t), pointer :: micm_state  => null()

   ! Sizing info
   integer :: n_micm_species       = 0
   integer :: n_micm_rate_params   = 0
   integer :: max_grid_cells       = 0

contains

   !> Initialize MICM solver and allocate state.
   subroutine micm_setup(config_path, n_grid_cells, errmsg, errcode)

      character(len=*), intent(in)  :: config_path
      integer,          intent(in)  :: n_grid_cells
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errcode

      type(error_t) :: error

      errmsg  = ''
      errcode = 0

      ! Create MICM solver
      micm_solver => micm_t(trim(config_path), RosenbrockStandardOrder, error)
      if (.not. error%is_success()) then
         errmsg = '[CheMPAS] Failed to create MICM solver: ' // error%message()
         errcode = 1; return
      end if
      if (.not. associated(micm_solver)) then
         errmsg = '[CheMPAS] MICM solver is null after construction'
         errcode = 1; return
      end if

      ! Query maximum grid cells
      max_grid_cells = micm_solver%get_maximum_number_of_grid_cells()

      ! Create state for the requested number of grid cells
      ! (capped at max if needed — batching handled by driver)
      micm_state => micm_solver%get_state(min(n_grid_cells, max_grid_cells), error)
      if (.not. error%is_success()) then
         errmsg = '[CheMPAS] Failed to create MICM state: ' // error%message()
         errcode = 1; return
      end if
      if (.not. associated(micm_state)) then
         errmsg = '[CheMPAS] MICM state is null after construction'
         errcode = 1; return
      end if

      n_micm_species     = micm_state%number_of_species
      n_micm_rate_params = micm_state%number_of_rate_parameters

   end subroutine micm_setup


   !> Solve chemistry for all grid cells in the state.
   !! The caller must have already filled conditions, concentrations,
   !! and rate_parameters arrays in micm_state.
   subroutine micm_solve(dt, errmsg, errcode)

      real (kind=real64), intent(in)  :: dt   ! time step [s]
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errcode

      type(error_t)        :: error
      type(string_t)       :: solver_state
      type(solver_stats_t) :: stats

      errmsg  = ''
      errcode = 0

      call micm_solver%solve(dt, micm_state, solver_state, stats, error)
      if (.not. error%is_success()) then
         errmsg = '[CheMPAS] MICM solve failed: ' // error%message()
         errcode = 1; return
      end if

   end subroutine micm_solve


   !> Clean up MICM resources.
   subroutine micm_cleanup()

      if (associated(micm_state)) then
         deallocate(micm_state); nullify(micm_state)
      end if
      if (associated(micm_solver)) then
         deallocate(micm_solver); nullify(micm_solver)
      end if

   end subroutine micm_cleanup

end module mpas_chemistry_micm

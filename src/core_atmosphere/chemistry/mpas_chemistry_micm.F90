! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! MICM interface for MPAS-A chemistry coupling.
! Creates MICM solver instance, manages state, and solves chemistry in batches.
! Auto-detects solver type: DAE4 when cloud chemistry is present (mpas_cloud_water.txt),
! otherwise standard Rosenbrock.
!
module mpas_chemistry_micm

   use mpas_kind_types,  only : RKIND
   use mpas_log,         only : mpas_log_write
   use iso_fortran_env,  only : real64
   use musica_micm,      only : micm_t, RosenbrockStandardOrder, &
                                RosenbrockDAE4StandardOrder, solver_stats_t, &
                                rosenbrock_solver_parameters_t
   use musica_state,     only : state_t, conditions_t
   use musica_util,      only : error_t, mappings_t, string_t

   implicit none

   private
   public :: micm_setup, micm_solve, micm_cleanup
   public :: micm_state, micm_solver_ptr, n_micm_species, n_micm_rate_params
   public :: using_dae_solver

   ! Module-level MICM objects
   type(micm_t),  pointer :: micm_solver => null()
   type(micm_t),  pointer :: micm_solver_ptr => null()   ! public read-only alias
   type(state_t), pointer :: micm_state  => null()

   ! Sizing info
   integer :: n_micm_species       = 0
   integer :: n_micm_rate_params   = 0
   integer :: max_grid_cells       = 0

   ! Solver type flag
   logical, save :: using_dae_solver = .false.

contains

   !> Initialize MICM solver and allocate state.
   !! Auto-detects solver type: if mpas_cloud_water.txt exists in the config
   !! directory, uses RosenbrockDAE4 (for MIAM cloud chemistry).
   !! Otherwise uses standard Rosenbrock.
   subroutine micm_setup(config_path, n_grid_cells, errmsg, errcode)

      use mpas_log, only : mpas_log_write

      character(len=*), intent(in)  :: config_path
      integer,          intent(in)  :: n_grid_cells
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errcode

      type(error_t) :: error
      character(len=512) :: cloud_water_file
      logical :: has_cloud
      integer :: solver_type
      integer :: last_slash

      errmsg  = ''
      errcode = 0

      ! Derive config directory from config_path (strip trailing /config.json)
      last_slash = index(config_path, '/', back=.true.)
      if (last_slash > 0) then
         cloud_water_file = config_path(1:last_slash) // 'mpas_cloud_water.txt'
      else
         cloud_water_file = 'mpas_cloud_water.txt'
      end if
      inquire(file=trim(cloud_water_file), exist=has_cloud)

      if (has_cloud) then
         solver_type = RosenbrockDAE4StandardOrder
         using_dae_solver = .true.
         call mpas_log_write('[CheMPAS] Cloud chemistry detected → DAE4 solver')
      else
         solver_type = RosenbrockStandardOrder
         using_dae_solver = .false.
         call mpas_log_write('[CheMPAS] Gas-phase only → Rosenbrock solver')
      end if

      ! Create MICM solver
      micm_solver => micm_t(trim(config_path), solver_type, error)
      if (.not. error%is_success()) then
         errmsg = '[CheMPAS] Failed to create MICM solver: ' // error%message()
         errcode = 1; return
      end if
      if (.not. associated(micm_solver)) then
         errmsg = '[CheMPAS] MICM solver is null after construction'
         errcode = 1; return
      end if

      ! Expose solver pointer for species discovery
      micm_solver_ptr => micm_solver

      ! For DAE solvers, increase constraint initialization iterations
      ! (default 10 is too few when starting far from equilibrium)
      if (using_dae_solver) then
         block
            type(rosenbrock_solver_parameters_t) :: params
            params = micm_solver%get_rosenbrock_solver_parameters(error)
            if (error%is_success()) then
               params%constraint_init_max_iterations = 100
               params%constraint_init_tolerance = 1.0e-8_real64
               call micm_solver%set_rosenbrock_solver_parameters(params, error)
               if (error%is_success()) then
                  call mpas_log_write('[CheMPAS] DAE constraint init: max_iter=100, tol=1e-8')
               end if
            end if
         end block
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

      ! DIAGNOSTIC: report solver state on first call
      block
         logical, save :: first = .true.
         if (first) then
            first = .false.
            call mpas_log_write('[DIAG] solver_state = ' // trim(solver_state%value_))
            call mpas_log_write('[DIAG] stats: accepted=$i, rejected=$i, decompositions=$i', &
                                intArgs=(/int(stats%accepted()), int(stats%rejected()), &
                                          int(stats%decompositions())/))
            call mpas_log_write('[DIAG] stats: solves=$i, final_time=$r', &
                                intArgs=(/int(stats%solves())/), &
                                realArgs=(/real(stats%final_time(), RKIND)/))
         end if
      end block

      ! Log solver statistics in a machine-parseable format for performance analysis.
      ! [CHEM_STATS] lines are parsed by the verification notebooks.
      ! Each call covers nCellsSolve * nVertLevels grid cells solved together.
      block
         integer, save :: call_count = 0
         call_count = call_count + 1

         ! Full stats on first call so notebook can see solver-state string
         if (call_count == 1) then
            call mpas_log_write('[DIAG] micm_solve first call: solver_state=' // &
                                trim(solver_state%value_))
         end if

         ! Parseable stats on every call:
         ! [CHEM_STATS] call=N accepted=A rejected=R nsteps=S decomp=D solves=V final_t=F
         call mpas_log_write( &
            '[CHEM_STATS] call=$i accepted=$i rejected=$i nsteps=$i decomp=$i solves=$i final_t=$r', &
            intArgs=(/call_count, int(stats%accepted()), int(stats%rejected()), &
                      int(stats%number_of_steps()), int(stats%decompositions()), &
                      int(stats%solves())/), &
            realArgs=(/real(stats%final_time(), RKIND)/))
      end block

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

! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!-----------------------------------------------------------------------
!  mpas_quacs
!
!> \brief Handles communication between QUACS and MPAS
!> \details
!>  This module handles communication between QUACS and MPAS, including
!>  constructs that have pybind11 wrappers and implement the QUACS
!>  HostAdapter protocol.
!>
!>  Specific interactions include:
!>  - Variable registration w/ callbacks
!>  - In-place state data communication
!>  - QUACS hook execution
!
!-----------------------------------------------------------------------
module mpas_quacs

    use mpas_kind_types, only : RKIND, StrKIND, ShortStrKIND
    use mpas_derived_types, only : mpas_pool_type, block_type, MPAS_LOG_ERR
    use mpas_log, only : mpas_log_write
    use mpas_runtime_vars, only : MPAS_var, MPAS_var_array, MPAS_var_add_callback, MPAS_VAR_OK, &
                                  MPAS_VAR_REAL, MPAS_VAR_INTEGER, MPAS_VAR_LOGICAL

    implicit none

    private

    public :: quacs_declare_variables, &
              quacs_add_fields, &
              quacs_set_initial_conditions

    ! Where a variable has values
    integer, parameter :: QUACS_CELL      = 1, &
                          QUACS_INTERFACE = 2, &
                          QUACS_COLUMN    = 3, &
                          QUACS_DOMAIN    = 4

    type quacs_variable
        character(len=:), allocatable :: name
        integer :: location = QUACS_CELL
        integer :: vartype = MPAS_VAR_REAL
        character(len=:), allocatable :: units
        character(len=:), allocatable :: description
        logical :: transported = .false.
    end type quacs_variable

    type (quacs_variable), dimension(:), allocatable :: variables

contains

    
   !***********************************************************************
   !
   !  routine quacs_declare_variables
   !
   !> \brief   Atmosphere core setup routine
   !> \details 
   !>  This routine queries QUACS for all required state variables and
   !>  stores them in `variables`.
   !>
   !>  STUBBED OUT - needs finished
   !
   !-----------------------------------------------------------------------
    subroutine quacs_declare_variables(configs)

        use mpas_pool_routines, only : mpas_pool_get_config

        ! Arguments
        type (mpas_pool_type), intent(in) :: configs

        ! Local variables
        logical, pointer :: config_quacs_test_variables

        ! This is called for every block, so return if the variables we already created
        if (allocated(variables)) return

        nullify(config_quacs_test_variables)
        call mpas_pool_get_config(configs, 'config_quacs_test_variables', config_quacs_test_variables)
        if (.not. associated(config_quacs_test_variables)) return
        if (.not. config_quacs_test_variables) return

        allocate(variables(7))
        variables(1) = quacs_variable('quacs_test_tracer', QUACS_CELL, MPAS_VAR_REAL, &
                                      'kg kg^{-1}', 'Test tracer mixing ratio', .true.)
        variables(2) = quacs_variable('quacs_test_cell', QUACS_CELL, MPAS_VAR_REAL, &
                                      '-', 'Test variable with one value per grid cell', .false.)
        variables(3) = quacs_variable('quacs_test_interface', QUACS_INTERFACE, MPAS_VAR_REAL, &
                                      '-', 'Test variable with one value per interface', .false.)
        variables(4) = quacs_variable('quacs_test_column', QUACS_COLUMN, MPAS_VAR_REAL, &
                                      '-', 'Test variable with one value per column', .false.)
        variables(5) = quacs_variable('quacs_test_count', QUACS_CELL, MPAS_VAR_INTEGER, &
                                      '-', 'Test integer variable', .false.)
        variables(6) = quacs_variable('quacs_test_flag', QUACS_COLUMN, MPAS_VAR_LOGICAL, &
                                      '-', 'Test logical variable', .false.)
        variables(7) = quacs_variable('quacs_test_scalar', QUACS_DOMAIN, MPAS_VAR_REAL, &
                                      '-', 'Test variable with one value for the domain', .false.)
                                      
        call MPAS_var_add_callback('state', 'scalars', add_constituents)
        call MPAS_var_add_callback('tend', 'scalars_tend', add_constituents)

    end subroutine quacs_declare_variables


    !-----------------------------------------------------------------------
    !  routine add_constituents
    !
    !> \brief Adds each transported variable to scalars and scalars_tend
    !
    !-----------------------------------------------------------------------
    subroutine add_constituents(structPool, newVarArray, poolName, varArrayName)

        ! Arguments
        type (mpas_pool_type), intent(inout) :: structPool
        type (MPAS_var_array) :: newVarArray
        character(len=*), intent(in) :: poolName
        character(len=*), intent(in) :: varArrayName

        ! Local variables
        type (MPAS_var) :: newVar
        integer :: i

        do i = 1, size(variables)
            if (.not. variables(i) % transported) cycle

            newVar % vartype = MPAS_VAR_REAL
            newVar % dimensions = 'nVertLevels nCells'
            newVar % hasTimeDimension = .true.
            newVar % isPersistent = .true.
            if (trim(varArrayName) == 'scalars') then
                newVar % name = variables(i) % name
                newVar % units = variables(i) % units
                newVar % description = variables(i) % description
            else
                newVar % name = 'tend_'//trim(variables(i) % name)
                newVar % name_in_code = variables(i) % name
                newVar % units = 'kg m^{-3} s^{-1}'
                newVar % description = 'Tendency of '//trim(variables(i) % description)
            end if
            call newVarArray % add_constituent(newVar, 'chemistry')
            call newVar % reset()
        end do

    end subroutine add_constituents


    !-----------------------------------------------------------------------
    !  routine quacs_add_fields
    !
    !> \brief Adds each variable that is not transported to a 'quacs' pool
    !> \details
    !>  Must be called after the pools of a block are generated and before the
    !>  framework allocates the block's fields.
    !
    !-----------------------------------------------------------------------
    subroutine quacs_add_fields(block, ierr)

        use mpas_pool_routines, only : mpas_pool_create_pool, mpas_pool_add_subpool

        ! Arguments
        type (block_type), pointer :: block
        integer, intent(out) :: ierr

        ! Local variables
        type (mpas_pool_type), pointer :: quacsPool
        type (MPAS_var) :: newVar
        integer :: i

        ierr = MPAS_VAR_OK
        if (.not. allocated(variables)) return

        nullify(quacsPool)
        call mpas_pool_create_pool(quacsPool)
        call mpas_pool_add_subpool(block % structs, 'quacs', quacsPool)
        call mpas_pool_add_subpool(block % allStructs, 'quacs', quacsPool)

        do i = 1, size(variables)
            if (variables(i) % transported) cycle

            newVar % name = variables(i) % name
            newVar % vartype = variables(i) % vartype
            newVar % units = variables(i) % units
            newVar % description = variables(i) % description
            newVar % hasTimeDimension = .true.
            newVar % isPersistent = .true.
            newVar % nTimeLevels = 1
            select case(variables(i) % location)
            case (QUACS_CELL)
                newVar % dimensions = 'nVertLevels nCells'
            case (QUACS_INTERFACE)
                newVar % dimensions = 'nVertLevelsP1 nCells'
            case (QUACS_COLUMN)
                newVar % dimensions = 'nCells'
            case (QUACS_DOMAIN)
                newVar % dimensions = ''
            end select

            call newVar % add_to_pool(block, quacsPool, block % packages, ierr=ierr)
            if (ierr == MPAS_VAR_OK) then
                call newVar % add_to_pool(block, block % allFields, block % packages, ierr=ierr)
            end if
            if (ierr /= MPAS_VAR_OK) then
                call mpas_log_write('Could not create variable '//trim(variables(i) % name), &
                                     messagetype=MPAS_LOG_ERR)
                return
            end if
            call newVar % reset()
        end do

    end subroutine quacs_add_fields


    !-----------------------------------------------------------------------
    !  routine quacs_set_initial_conditions
    !
    !> \brief Sets initial values for the QUACS-owned variables
    !> \details
    !>  This calls the QUACS hook that sets initial conditions.
    !>  Values are set in every cell of the block, including halo cells.
    !>
    !>  STUBBED OUT - needs finished
    !
    !-----------------------------------------------------------------------
    subroutine quacs_set_initial_conditions(block)

        use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension

        ! Arguments
        type (block_type), pointer :: block

        ! Local variables
        type (mpas_pool_type), pointer :: state, quacsPool
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: cell_values, interface_values
        real (kind=RKIND), dimension(:), pointer :: column_values
        real (kind=RKIND), pointer :: scalar_value
        integer, dimension(:,:), pointer :: counts
        integer, dimension(:), pointer :: flags
        integer, pointer :: index_tracer, nCells, nVertLevels
        integer :: iCell, k

        if (.not. allocated(variables)) return

        call mpas_pool_get_subpool(block % structs, 'state', state)
        call mpas_pool_get_subpool(block % structs, 'quacs', quacsPool)

        call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
        call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(state, 'index_quacs_test_tracer', index_tracer)

        call mpas_pool_get_array(state, 'scalars', scalars, 1)
        call mpas_pool_get_array(quacsPool, 'quacs_test_cell', cell_values)
        call mpas_pool_get_array(quacsPool, 'quacs_test_interface', interface_values)
        call mpas_pool_get_array(quacsPool, 'quacs_test_column', column_values)
        call mpas_pool_get_array(quacsPool, 'quacs_test_count', counts)
        call mpas_pool_get_array(quacsPool, 'quacs_test_flag', flags)
        call mpas_pool_get_array(quacsPool, 'quacs_test_scalar', scalar_value)

        do iCell = 1, nCells
            do k = 1, nVertLevels
                scalars(index_tracer, k, iCell) = 1.0e-6_RKIND * real(k, kind=RKIND)
                cell_values(k, iCell) = real(k, kind=RKIND)
                counts(k, iCell) = k
            end do
            do k = 1, nVertLevels + 1
                interface_values(k, iCell) = real(k, kind=RKIND) - 0.5_RKIND
            end do
            column_values(iCell) = 3.0_RKIND
            flags(iCell) = 1
        end do
        scalar_value = 600.0_RKIND

    end subroutine quacs_set_initial_conditions

end module mpas_quacs

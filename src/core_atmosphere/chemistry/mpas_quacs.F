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
!> \brief Runs QUACS parameterizations on the CheMPAS state
!> \details
!>  This module communicates with QUACS via the `quacs_embed` library's C API.
!>  It is responsible for:
!>  - Starting QUACS, passing the CheMPAS-QUACS configuration file
!>  - Creating MPAS state variables for variables owned by QUACS
!>    parameterizations
!>  - Providing QUACS parameterizations access to MPAS state data
!>  - Transferring control to QUACS at each hook.
!>
!>  The CheMPAS-QUACS configuration must list every hook that CheMPAS runs,
!>  even if no parameterizations are configured for a hook.
!>
!>  Current hooks:
!>  - `initial_conditions`: run once on non-restart runs before the first
!>                          time step
!>  - `before_dynamics`: runs for each time step, before the dynamics
!>
!>  QUACS operates on every cell in every column, including halo cells.
!>  Only one block per MPI rank is supported.
!
!-----------------------------------------------------------------------
module mpas_quacs

    use, intrinsic :: iso_c_binding, only : c_associated, c_char, c_f_pointer, &
                                            c_int, c_int64_t, c_intptr_t,      &
                                            c_loc, c_null_char, c_null_ptr,    &
                                            c_ptr, c_size_t
    use mpas_kind_types, only : RKIND, R8KIND, StrKIND, ShortStrKIND
    use mpas_derived_types, only : mpas_pool_type, mpas_pool_field_info_type,  &
                                   block_type, MPAS_POOL_REAL,                 &
                                   MPAS_POOL_INTEGER, MPAS_LOG_OUT,            &
                                   MPAS_LOG_ERR, MPAS_LOG_CRIT
    use mpas_log, only : mpas_log_write
    use mpas_runtime_vars, only : MPAS_var, MPAS_var_array,                    &
                                  MPAS_var_add_callback, MPAS_VAR_OK,          &
                                  MPAS_VAR_REAL, MPAS_VAR_INTEGER,             &
                                  MPAS_VAR_LOGICAL

    implicit none

    private

    public :: quacs_declare_variables, &
              quacs_add_fields, &
              quacs_run_hook, &
              quacs_finalize

    ! The values of the following parameters must match `quacs_embed.h`.

    ! Character array sizes
    integer, parameter :: QUACS_NAME_LENGTH = 128, &
                          QUACS_DESCRIPTION_LENGTH = 512

    ! Where a variable has values.
    integer(kind=c_int), parameter :: QUACS_CELL           = 1, &
                                       QUACS_INTERFACE      = 2, &
                                       QUACS_SURFACE_CELL   = 3, &
                                       QUACS_MODEL_TOP_CELL = 4, &
                                       QUACS_COLUMN         = 5, &
                                       QUACS_DOMAIN         = 6

    ! QUACS variable types
    integer(kind=c_int), parameter :: QUACS_FLOAT = 1, &
                                       QUACS_INT   = 2, &
                                       QUACS_BOOL  = 3

    ! Variable permissions for parameterizations
    integer(kind=c_int), parameter :: QUACS_READ_ONLY = 1, &
                                       QUACS_UPDATE    = 2, &
                                       QUACS_OWN       = 3

    ! Host-application storage types. MPAS integers are 32-bit. Booleans are
    ! stored as integers that are 0 (False) or 1 (True).
    integer(kind=c_int), parameter :: QUACS_STORAGE_FLOAT32 = 1, &
                                       QUACS_STORAGE_FLOAT64 = 2, &
                                       QUACS_STORAGE_INT32   = 3, &
                                       QUACS_STORAGE_INT64   = 4, &
                                       QUACS_STORAGE_BOOL    = 5
    integer(kind=c_int), parameter :: REAL_STORAGE = &
        merge(QUACS_STORAGE_FLOAT64, QUACS_STORAGE_FLOAT32, RKIND == R8KIND)
    
#ifdef MPAS_USE_QUACS
    ! Interface with the `quacs_embed` C API

    type, bind(C) :: quacs_variable_c
        character(kind=c_char) :: name(QUACS_NAME_LENGTH)
        character(kind=c_char) :: units(QUACS_NAME_LENGTH)
        character(kind=c_char) :: description(QUACS_DESCRIPTION_LENGTH)
        integer(kind=c_int) :: shape
        integer(kind=c_int) :: dtype
        integer(kind=c_int) :: permissions
        integer(kind=c_int) :: transported ! o (False) or 1 (True)
    end type quacs_variable_c

    interface
        function quacs_initialize_c(config_path) &
            bind(C, name='quacs_initialize') result(status)
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: config_path(*)
            integer(kind=c_int) :: status
        end function quacs_initialize_c

        function quacs_map_name_c(name, host_name) &
            bind(C, name='quacs_map_name') result(status)
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: name(*)
            character(kind=c_char), intent(in) :: host_name(*)
            integer(kind=c_int) :: status
        end function quacs_map_name_c

        function quacs_variable_count_c(count) &
            bind(C, name='quacs_variable_count') result(status)
            import :: c_int
            integer(kind=c_int), intent(out) :: count
            integer(kind=c_int) :: status
        end function quacs_variable_count_c

        function quacs_variable_info_c(index, variable) &
            bind(C, name='quacs_variable_info') result(status)
            import :: c_int, quacs_variable_c
            integer(kind=c_int), value :: index
            type (quacs_variable_c), intent(out) :: variable
            integer(kind=c_int) :: status
        end function quacs_variable_info_c

        function quacs_set_grid_c(ncolumns, nlevels) &
            bind(C, name='quacs_set_grid') result(status)
            import :: c_int
            integer(kind=c_int), value :: ncolumns
            integer(kind=c_int), value :: nlevels
            integer(kind=c_int) :: status
        end function quacs_set_grid_c

        function quacs_bind_variable_c(name, data, storage, column_stride, &
            level_stride) bind(C, name='quacs_bind_variable') result(status)
            import :: c_char, c_int, c_int64_t, c_ptr
            character(kind=c_char), intent(in) :: name(*)
            type(c_ptr), value :: data
            integer(kind=c_int), value :: storage
            integer(kind=c_int64_t), value :: column_stride
            integer(kind=c_int64_t), value :: level_stride
            integer(kind=c_int) :: status
        end function quacs_bind_variable_c

        function quacs_run_hook_c(label) &
            bind(C, name='quacs_run_hook') result(status)
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: label(*)
            integer(kind=c_int) :: status
        end function quacs_run_hook_c

        function quacs_check_conformance_c(failures) &
            bind(C, name='quacs_check_conformance') result(status)
            import :: c_char, c_int
            integer(kind=c_int), intent(out) :: failures
            integer(kind=c_int) :: status
        end function quacs_check_conformance_c

        function quacs_report_c() bind(C, name='quacs_report') result(text)
            import :: c_ptr
            type(c_ptr) :: text
        end function quacs_report_c

        function quacs_last_error_c() bind(C, name='quacs_last_error') &
            result(text)
            import :: c_ptr
            type(c_ptr) :: text
        end function quacs_last_error_c

        subroutine quacs_finalize_c() bind(C, name='quacs_finalize')
        end subroutine quacs_finalize_c

        function strlen(string) bind(C, name='strlen') result(length)
            import :: c_ptr, c_size_t
            type(c_ptr), value :: string
            integer(kind=c_size_t) :: length
        end function strlen
    end interface
#endif

    ! The MPAS name of each MPAS variable that has a different name in QUACS.
    ! Each column holds the QUACS name, then the MPAS name.
    character(len=*), dimension(2, 1), parameter :: name_map = reshape( &
        [ character (len=ShortStrKIND) :: &
        'base_state_pressure', 'pressure_base'], [2, 1])

    ! A variable used by QUACS parameterizations
    type quacs_variable
        character(len=:), allocatable :: name
        character(len=:), allocatable :: units
        character(len=:), allocatable :: description
        integer :: shape
        integer :: dtype
        logical :: owned
        logical :: transported
        integer(kind=c_intptr_t) :: address = 0
    end type quacs_variable

    type(quacs_variable), dimension(:), allocatable :: variables

    logical :: running = .false.     ! indicates whether QUACS is running
    logical :: grid_is_set = .false. ! indicates whether the MPAS grid is set

contains


    !-----------------------------------------------------------------------
    !  routine quacs_declare_variables
    !
    !> \brief Starts QUACS and collects variables used by parameterizations
    !> \details
    !>  Must be called before pools are generated, so that the transported
    !>  variables can be added to `scalars` and `scalars_tend`.
    !
    !-----------------------------------------------------------------------
    subroutine quacs_declare_variables(configs)

        use mpas_pool_routines, only : mpas_pool_get_config

        ! Arguments
        type(mpas_pool_type), intent(in) :: configs

        ! Local variables
        character(len=StrKIND), pointer :: config_quacs_config_file
#ifdef MPAS_USE_QUACS
        type(quacs_variable_c) :: info
        integer(kind=c_int) :: count, i
#endif

        ! This routine is called for every block, but QUACS only needs
        ! to be started once.
        if (running) return

        nullify(config_quacs_config_file)
        call mpas_pool_get_config(configs, 'config_quacs_config_file', &
            config_quacs_config_file)
        if (.not. associated(config_quacs_config_file)) return
        if (len_trim(config_quacs_config_file) == 0) return

#ifdef MPAS_USE_QUACS
        call check(quacs_initialize_c(&
            trim(config_quacs_config_file)//c_null_char), 'starting')
        running = .true.

        do i = 1, size(name_map, 2)
            call check(quacs_map_name_c( &
                trim(name_map(1,i))//c_null_char, &
                trim(name_map(2,i))//c_null_char &
                ), 'mapping the name '//trim(name_map(1,i)))
        end do

        call check(quacs_variable_count_c(count), 'counting its variables')
        allocate(variables(count))
        do i = 1, count
            call check(quacs_variable_info_c(i - 1, info), &
                'collecting variable info')
            variables(i) % name = to_string(info % name)
            variables(i) % units = to_string(info % units)
            variables(i) % description = to_string(info % description)
            variables(i) % shape = info % shape
            variables(i) % dtype = info % dtype
            variables(i) % owned = info % permissions == QUACS_OWN
            variables(i) % transported = info % transported == 1
            if (len(variables(i) % name) > ShortStrKIND) then
                call mpas_log_write('The name of the QUACS variable '// &
                    variables(i) % name//' is too long for MPAS', &
                    messageType=MPAS_LOG_CRIT)
            end if
        end do

        if (any(variables(:) % transported)) then
            call MPAS_var_add_callback('state', 'scalars', add_constituents)
            call MPAS_var_add_callback('tend', 'scalars_tend', add_constituents)
        end if

        call mpas_log_write('QUACS is running with the configuration file '// &
            trim(config_quacs_config_file))
#else
        call mpas_log_write('config_quacs_config_file is set, but MPAS was '// &
            'not built with QUACS=true', messageType=MPAS_LOG_CRIT)
#endif

    end subroutine quacs_declare_variables


    !-----------------------------------------------------------------------
    !  routine quacs_add_fields
    !
    !> \brief Creates each non-transported QUACS-owned variable
    !> \details
    !>  The variables are added to a `quacs` pool.
    !
    !-----------------------------------------------------------------------
    subroutine quacs_add_fields(block, ierr)

        use mpas_pool_routines, only : mpas_pool_create_pool, &
                                       mpas_pool_add_subpool, &
                                       mpas_pool_get_field_info

        ! Arguments
        type(block_type), pointer :: block
        integer, intent(out) :: ierr

        ! Local variables
        type(mpas_pool_type), pointer :: quacsPool
        type(mpas_pool_field_info_type) :: info
        type(MPAS_var) :: newVar
        integer :: i

        ierr = MPAS_VAR_OK
        if (.not. allocated(variables)) return

        nullify(quacsPool)
        call mpas_pool_create_pool(quacsPool)
        call mpas_pool_add_subpool(block % structs, 'quacs', quacsPool)
        call mpas_pool_add_subpool(block % allStructs, 'quacs', quacsPool)

        do i = 1, size(variables)
            if (.not. variables(i) % owned .or. variables(i) % transported) cycle

            call mpas_pool_get_field_info(block % allFields, &
                variables(i) % name, info)
            if (info % fieldType /= -1) then
                call mpas_log_write('A QUACS parameterization owns the '// &
                    'variable '//variables(i) % name//&
                    ', which MPAS already has', messageType=MPAS_LOG_ERR)
                ierr = 1
                return
            end if

            newVar % name = variables(i) % name
            newVar % units = variables(i) % units
            newVar % description = variables(i) % description
            newVar % dimensions = mpas_dimensions(variables(i) % shape)
            newVar % hasTimeDimension = .true.
            newVar % isPersistent = .true.
            newVar % nTimeLevels = 1
            select case (variables(i) % dtype)
            case (QUACS_INT)
                newVar % vartype = MPAS_VAR_INTEGER
            case (QUACS_BOOL)
                newVar % vartype = MPAS_VAR_LOGICAL
            case (QUACS_FLOAT)
                newVar % vartype = MPAS_VAR_REAL
            case default
                call mpas_log_write('Internal error. Unsupported QUACS '// &
                    'variable type', messageType=MPAS_LOG_CRIT)
            end select

            call newVar % add_to_pool(block, quacsPool, block % packages, &
                ierr=ierr)
            if (ierr == MPAS_VAR_OK) then
                call newVar % add_to_pool(block, block % allFields, &
                block % packages, ierr=ierr)
            end if
            if (ierr /= MPAS_VAR_OK) then
                call mpas_log_write('Could not create the QUACS variable '// &
                    variables(i) % name, messageType=MPAS_LOG_ERR)
                return
            end if
            call newVar % reset()
        end do

    end subroutine quacs_add_fields


    !-----------------------------------------------------------------------
    !  routine quacs_run_hook
    !
    !> \brief Runs QUACS parameterizations for a labelled hook
    !> \details
    !>  Parameterizations operate on all cells, including halo cells.
    !
    !-----------------------------------------------------------------------
    subroutine quacs_run_hook(block, label)

        ! Arguments
        type(block_type), pointer :: block
        character(len=*), intent(in) :: label

        if (.not. running) return

#ifdef MPAS_USE_QUACS
        call bind_variables(block)
        call check(quacs_run_hook_c(label//c_null_char), &
            'running the hook '//label)
#endif

    end subroutine quacs_run_hook


    !-----------------------------------------------------------------------
    !  routine quacs_finalize
    !
    !> \brief Stops QUACS
    !> \detail
    !>  If config_quacs_check_conformance is set, this first writes the results
    !>  of the QUACS conformance checks in the configuration file to the log.
    !>  The results are written as errors if any check failed.
    !
    !-----------------------------------------------------------------------
    subroutine quacs_finalize(configs)

        use mpas_pool_routines, only : mpas_pool_get_config

        ! Arguments
        type(mpas_pool_type), intent(in) :: configs

        ! Local variables
        logical, pointer :: config_quacs_check_conformance
#ifdef MPAS_USE_QUACS
        integer(c_int) :: failures
#endif

        if (.not. running) return

#ifdef MPAS_USE_QUACS
        nullify(config_quacs_check_conformance)
        call mpas_pool_get_config(configs, 'config_quacs_check_conformance', &
            config_quacs_check_conformance)
        if (associated(config_quacs_check_conformance)) then
            if (config_quacs_check_conformance) then
                call check(quacs_check_conformance_c(failures), &
                    'checking conformance')
                call log_lines(to_string_from_c(quacs_report_c()), &
                    merge(MPAS_LOG_ERR, MPAS_LOG_OUT, failures > 0))
            end if
        end if

        call quacs_finalize_c()
        running = .false.
#endif

    end subroutine quacs_finalize


    ! ------------------------
    ! Private helper functions
    ! ------------------------


    !-----------------------------------------------------------------------
    !  routine add_constituents
    !
    !> \brief Callback function that adds variables to `scalar` and
    !>  `scalar_tend`
    !
    !-----------------------------------------------------------------------
    subroutine add_constituents(structPool, newVarArray, poolName, varArrayName)

        ! Arguments
        type(mpas_pool_type), intent(inout) :: structPool
        type(MPAS_var_array) :: newVarArray
        character(len=*), intent(in) :: poolName
        character(len=*), intent(in) :: varArrayName

        ! Local variables
        type(MPAS_var) :: newVar
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
            else if (trim(varArrayName) == 'scalars_tend') then
                newVar % name = 'tend_'//variables(i) % name
                newVar % name_in_code = variables(i) % name
                newVar % units = variables(i) % units//' s^{-1}'
                newVar % description = 'Tendency of '// &
                    variables(i) % description
            end if
            call newVarArray % add_constituent(newVar, 'chemistry')
            call newVar % reset()
        end do

    end subroutine add_constituents


    !-----------------------------------------------------------------------
    !  function mpas_dimensions
    !
    !> \brief Returns the MPAS dimensions for a QUACS variable
    !
    !-----------------------------------------------------------------------
    function mpas_dimensions(shape) result(dimensions)

        ! Arguments
        integer, intent(in) :: shape

        ! Return value
        character(len=:), allocatable :: dimensions

        select case(shape)
        case (QUACS_CELL)
            dimensions = 'nVertLevels nCells'
        case (QUACS_INTERFACE)
            dimensions = 'nVertLevelsP1 nCells'
        case (QUACS_SURFACE_CELL)
            dimensions = 'nCells'
        case (QUACS_MODEL_TOP_CELL)
            dimensions = 'nCells'
        case (QUACS_COLUMN)
            dimensions = 'nCells'
        case (QUACS_DOMAIN)
            dimensions = ''
        case default
            call mpas_log_write('Internal error. Unsupported QUACS dimension', &
                messageType=MPAS_LOG_CRIT)
        end select

    end function mpas_dimensions


#ifdef MPAS_USE_QUACS
    !-----------------------------------------------------------------------
    !  routine bind_variables
    !
    !> \brief Binds MPAS state variables to QUACS variables
    !> \details
    !>  Provides QUACS dimensions and memory location for MPAS state variables.
    !>  Rebinds variables whenever they change memory location.
    !
    !-----------------------------------------------------------------------
    subroutine bind_variables(block)

        use mpas_pool_routines, only : mpas_pool_get_dimension

        ! Arguments
        type(block_type), pointer :: block

        ! Local variables
        integer, pointer :: nCells, nVertLevels
        integer :: i

        call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
        call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', &
            nVertLevels)

        if (.not. grid_is_set) then
            if (associated(block % prev) .or. associated(block % next)) then
                call mpas_log_write('QUACS can only be used with one block '// &
                'per MPI rank', messageType=MPAS_LOG_CRIT)
            end if
            call check(quacs_set_grid_c(int(nCells, c_int), &
                int(nVertLevels, c_int)), 'setting its grid')
            grid_is_set = .true.
        end if

        do i = 1, size(variables)
            call bind_variable(block, variables(i), nCells, nVertLevels)
        end do

    end subroutine bind_variables


    !-----------------------------------------------------------------------
    !  routine bind_variable
    !
    !> \brief Binds a single MPAS variable to a QUACS variable
    !> \details
    !>  Binds variables only when they are unbound or have moved
    !
    !-----------------------------------------------------------------------
    subroutine bind_variable(block, variable, nCells, nVertLevels)

        use mpas_pool_routines, only : mpas_pool_get_field_info, &
                                       mpas_pool_get_array, &
                                       mpas_pool_get_subpool, &
                                       mpas_pool_get_dimension

        ! Arguments
        type(block_type), pointer :: block
        type(quacs_variable), intent(inout) :: variable
        integer, intent(in) :: nCells
        integer, intent(in) :: nVertLevels

        ! Local variables
        type(mpas_pool_field_info_type) :: info
        type(mpas_pool_type), pointer :: state
        real(kind=RKIND), pointer :: real0d
        real(kind=RKIND), dimension(:), pointer :: real1d
        real(kind=RKIND), dimension(:,:), pointer :: real2d
        real(kind=RKIND), dimension(:,:,:), pointer :: scalars
        integer, pointer :: int0d
        integer, dimension(:), pointer :: int1d
        integer, dimension(:,:), pointer :: int2d
        integer, pointer :: index_constituent
        type(c_ptr) :: location
        integer(kind=c_intptr_t) :: address
        integer(kind=c_int) :: storage
        integer(kind=c_int64_t) :: value_size, column_stride, level_stride
        integer :: ndims, nlevels
        logical :: is_real

        select case (variable % shape)
        case (QUACS_CELL)
            ndims = 2
            nlevels = nVertLevels
        case (QUACS_INTERFACE)
            ndims = 2
            nlevels = nVertLevels + 1
        case (QUACS_DOMAIN)
            ndims = 0
            nlevels = 0
        case default
            ndims = 1
            nlevels = 0
        end select
        
        is_real = variable % dtype == QUACS_FLOAT
        if (is_real) then
            storage = REAL_STORAGE
            value_size = storage_size(1.0_RKIND) / 8
        else
            storage = QUACS_STORAGE_INT32
            value_size = storage_size(1) / 8
        end if

        location = c_null_ptr
        column_stride = 0
        level_stride = 0

        call mpas_pool_get_field_info(block % allFields, variable % name, info)

        if (is_real .and. info % fieldType == MPAS_POOL_REAL &
            .and. info % nDims == ndims) then

            if (ndims == 0) then
                call mpas_pool_get_array(block % allFields, variable % name, &
                    real0d)
                if (associated(real0d)) location = c_loc(real0d)
            else if (ndims == 1) then
                call mpas_pool_get_array(block % allFields, variable % name, &
                    real1d)
                if (associated(real1d)) then
                    if (size(real1d) == nCells + 1) then
                        location = c_loc(real1d(1))
                        column_stride = value_size
                    end if
                end if
            else
                call mpas_pool_get_array(block % allFields, variable % name, &
                    real2d)
                if (associated(real2d)) then
                    if (size(real2d, 1) == nlevels &
                        .and. size(real2d, 2) == nCells + 1) then
                        location = c_loc(real2d(1,1))
                        column_stride = value_size * size(real2d, 1)
                        level_stride = value_size
                    end if
                end if
            end if

        else if (.not. is_real .and. info % fieldType == MPAS_POOL_INTEGER &
                .and. info % nDims == ndims) then

            if (ndims == 0) then
                call mpas_pool_get_array(block % allFields, variable % name, &
                    int0d)
                if (associated(int0d)) location = c_loc(int0d)
            else if (ndims == 1) then
                call mpas_pool_get_array(block % allFields, variable % name, &
                    int1d)
                if (associated(int1d)) then
                    if (size(int1d) == nCells + 1) then
                        location = c_loc(int1d(1))
                        column_stride = value_size
                    end if
                end if
            else
                call mpas_pool_get_array(block % allFields, variable % name, &
                    int2d)
                if (associated(int2d)) then
                    if (size(int2d, 1) == nlevels &
                        .and. size(int2d, 2) == nCells + 1) then
                        location = c_loc(int2d(1,1))
                        column_stride = value_size * size(int2d, 1)
                        level_stride = value_size
                    end if
                end if
            end if

        else if (is_real .and. info % fieldType == -1 &
                .and. variable % shape == QUACS_CELL) then

            ! The variable is not a field, so it may be a constituent of scalars
            call mpas_pool_get_subpool(block % structs, 'state', state)
            nullify(index_constituent)
            call mpas_pool_get_dimension(state, 'index_'//variable % name, &
                index_constituent)
            if (associated(index_constituent)) then
                if (index_constituent > 0) then
                    call mpas_pool_get_array(state, 'scalars', scalars, 1)
                    location = c_loc(scalars(index_constituent,1,1))
                    column_stride = value_size * size(scalars, 1) * &
                        size(scalars, 2)
                    level_stride = value_size * size(scalars, 1)
                end if
            end if

        end if

        if (.not. c_associated(location)) then
            if (is_real) then
                call mpas_log_write('QUACS uses the variable '// &
                    variable % name//', but MPAS has no real field or '// &
                    'constituent of scalars with that name and the '// &
                    'dimensions "'//mpas_dimensions(variable % shape)//'"', &
                    messageType=MPAS_LOG_CRIT)
            else
                call mpas_log_write('QUACS uses the variable '// &
                    variable % name//', but MPAS has no integer field with '// &
                    'that name and the dimensions "'// &
                    mpas_dimensions(variable % shape)//'"', &
                    messageType=MPAS_LOG_CRIT)
            end if
        end if

        address = transfer(location, address)
        if (address == variable % address) return

        call check(quacs_bind_variable_c(variable % name//c_null_char, &
            location, storage, column_stride, level_stride), &
            'binding the variable '//variable % name)
        variable % address = address

    end subroutine bind_variable


    !-----------------------------------------------------------------------
    !  routine check
    !
    !> \brief Stops MPAS with an error message on a QUACS failure 
    !
    !-----------------------------------------------------------------------
    subroutine check(status, action)

        ! Arguments
        integer(kind=c_int), intent(in) :: status
        character(len=*), intent(in) :: action

        if (status == 0) return

        call mpas_log_write('QUACS failed while '//action//':', &
            messageType=MPAS_LOG_ERR)
        call log_lines(to_string_from_c(quacs_last_error_c()), MPAS_LOG_ERR)
        call mpas_log_write('MPAS stopping because of QUACS failure', &
            messageType=MPAS_LOG_CRIT)

    end subroutine check


    !-----------------------------------------------------------------------
    !  routine log_lines
    !
    !> \brief Writes each line of a string to the log
    !
    !-----------------------------------------------------------------------
    subroutine log_lines(text, messageType)

        ! Arguments
        character(len=*), intent(in) :: text
        integer, intent(in) :: messageType

        ! Local variables
        integer :: first, length

        first = 1
        do while (first <= len(text))
            length = index(text(first:), new_line('a')) - 1
            if (length < 0) length = len(text) - first + 1
            if (length > 0) call mpas_log_write(text(first:first+length-1), &
                messageType=messageType)
            first = first + length + 1
        end do

    end subroutine log_lines


    !-----------------------------------------------------------------------
    !  function to_string
    !
    !> \brief Returns the text that is held in an array of C characters
    !
    !-----------------------------------------------------------------------
    function to_string(chars) result(text)

        ! Arguments
        character (kind=c_char), dimension(:), intent(in) :: chars

        ! Return value
        character (len=:), allocatable :: text

        ! Local variables
        integer :: i, length

        length = size(chars)
        do i = 1, size(chars)
            if (chars(i) == c_null_char) then
                length = i - 1
                exit
            end if
        end do

        allocate(character (len=length) :: text)
        do i = 1, length
            text(i:i) = chars(i)
        end do

    end function to_string


    !-----------------------------------------------------------------------
    !  function to_string_from_c
    !
    !> \brief Returns the text that a C string points to
    !
    !-----------------------------------------------------------------------
    function to_string_from_c(string) result(text)

        ! Arguments
        type (c_ptr), intent(in) :: string

        ! Return value
        character (len=:), allocatable :: text

        ! Local variables
        character (kind=c_char), dimension(:), pointer :: chars

        if (.not. c_associated(string)) then
            text = ''
            return
        end if

        call c_f_pointer(string, chars, [strlen(string)])
        text = to_string(chars)

    end function to_string_from_c
#endif

end module mpas_quacs
Module api
    use ND_Arrays
    use vpm_lib
    use vpm_vars
    use vpm_size
    ! use openmpth
    use, intrinsic :: iso_c_binding, only: c_float, c_int, c_bool, c_null_ptr, c_double, c_ptr

    integer, parameter :: MAX_STRING_LENGTH = 256
contains
    subroutine initialize(dx_pm, dy_pm, dz_pm, proj_type, bc_type, vol_type,                        &
                         num_coarse, num_nbi, num_nbj, num_nbk, remesh_type, tree_type,             &
                          max_level, omp_threads, grid_define, write_type, write_start, write_steps &
                          ) &
        bind(C, name='init')

        use pmgrid, only: DXpm, DYpm, DZpm, IDVPM, ncoarse
        use vpm_vars, only: interf_iproj, IPMWRITE, idefine, IPMWSTART, IPMWSTEPS, OMPTHREADS, &
                            nremesh, iyntree, ilevmax, ibctyp, NBI, NBJ, NBK
        use parvar, only: set_neq
        use MPI

        ! Declare the parameters to be passed in
        implicit none
        real(c_double), intent(in) :: dx_pm, dy_pm, dz_pm
        integer(c_int), intent(in) :: proj_type, bc_type, vol_type, num_coarse
        integer(c_int), intent(in) :: num_nbi, num_nbj, num_nbk
        integer(c_int), intent(in) :: remesh_type, tree_type, max_level
        integer(c_int), intent(in) :: omp_threads, grid_define, write_type
        integer(c_int), intent(in) :: write_start(10), write_steps(10)

        integer :: ierr, my_rank, i

        ! Get the rank of the process
        call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
        if (ierr .ne. 0) then
            print *, 'Error getting the rank of the process'
            stop
        end if

        ! PM GRID
        DXpm = dx_pm
        DYpm = dy_pm
        DZpm = dz_pm
        IDVPM = vol_type
        ncoarse = num_coarse

        ! VPM_VARS
        interf_iproj = proj_type
        idefine = grid_define
        IPMWRITE = write_type
        ! Check the IPMWRITE parameter and process accordingly
        if (IPMWRITE .GT. 0) then
            do i = 1, IPMWRITE ! maximum value 10
                IPMWSTART(i) = write_start(i)
                IPMWSTEPS(i) = write_steps(i)
                if (IPMWRITE .gt. 10) stop  ! maximum value of writes equal to 10
            end do
        end if

        ! VPM_SIZE
        nremesh = remesh_type
        iyntree = tree_type
        ilevmax = max_level
        ibctyp = bc_type
        NBI = num_nbi
        NBJ = num_nbj
        NBK = num_nbk

        ! OPENMPTH
        OMPTHREADS = omp_threads

        call set_neq(3)

    end subroutine initialize

    subroutine set_verbose_level(verbocity_in) bind(C, name='set_verbosity')
        use console_io, only: verbocity
        implicit none
        integer(c_int), intent(in) :: verbocity_in
        verbocity = verbocity_in
    end subroutine set_verbose_level

    subroutine finalize() bind(C, name='finalize')
        implicit none
    end subroutine finalize

    subroutine call_vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
                        RHS_pm_out, Vel_out, NTIME_in, NI_in,                   &
                        NVR_size_in, deform_out) bind(C, name='vpm')
        !  -> XP : particle positions (3 * NVR)
        !  -> QP : particle quantities (neqpm * NVR)
        !  -> UP : particle velocities (3 * NVR)
        !  -> GP : particle defromation (3*NVR) (wmega * \nabla)  u
        !  -> NVR : number of particles
        !  -> neqpm : number of equations
        !  -> WhatToDo : 0 - initialize, 1 - solve, 2 - convect, 3 - project, 4 - back, 5 - diffuse
        !  -> \nabla^2 u_i = RHS_i (NEQ)
        !  -> Velx, Vely, Velz : velocity field at grid points
        !  -> NTIME_IN : current time
        !  -> NI_in: viscocity -> DIFFUSION OF VORTICITY
        !  -> NVRM : NVR MAX
        use vpm_lib, only: vpm
        use ND_Arrays
        use, intrinsic :: ieee_arithmetic

        ! Fortran to C bindings
        implicit none
        ! Interface for the arguments
        integer(c_int), intent(inout)          :: NVR_in
        integer(c_int), intent(in)             :: neqpm_in, WhatToDo, NVR_size_in, NTIME_in
        real(c_double), intent(in)             :: NI_in
        real(c_double), intent(inout), target  :: XP_in(3, NVR_in), QP_in(neqpm_in + 1, NVR_in)
        real(c_double), intent(inout), target  :: UP_in(3, NVR_in), GP_in(3, NVR_in)
        type(ND_Array), intent(out)            :: RHS_pm_out, Vel_out
        type(ND_Array), intent(out), optional  :: deform_out

        real(c_double), pointer          :: RHS_pm_ptr(:, :, :, :), Vel_ptr(:, :, :, :)
        real(c_double), pointer          :: deform_ptr(:, :, :, :)
        integer :: ierr, my_rank

        call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
        if (ierr .ne. 0) then
            print *, 'Error getting the rank of the process'
            stop
        end if

        call vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
                 RHS_pm_ptr, Vel_ptr, NTIME_in, NI_in, NVR_size_in, &
                 deform_ptr)

        ! Copy the data back to the arrays
        ! Assign the pointers to the arrays
        if (associated(RHS_pm_ptr)) then
            RHS_pm_out = from_intrinsic(RHS_pm_ptr, shape(RHS_pm_ptr))
        end if

        if (associated(Vel_ptr)) then
            Vel_out = from_intrinsic(Vel_ptr, shape(Vel_ptr))
        end if

        if ((associated(deform_ptr)) .and. (present(deform_out))) then
            deform_out = from_intrinsic(deform_ptr, shape(deform_ptr))
        end if
    end subroutine call_vpm

    subroutine call_remesh_particles_3d(iflag, npar_per_cell, XP_arr, QP_arr, GP_arr, UP_arr, NVR_out, cuttof_value) &
        bind(C, name='remesh_particles_3d')
        use vpm_remesh, only: remesh_particles_3d
        use base_types, only: dp
        use ND_Arrays

        implicit none
        integer(c_int), intent(in) :: iflag, npar_per_cell
        type(ND_Array), intent(out), target :: XP_arr, QP_arr, GP_arr, UP_arr
        integer(c_int), intent(out) :: NVR_out
        real(c_double), intent(in), optional :: cuttof_value

        ! Local variables
        real(dp), allocatable, target, save :: XP_out(:, :), QP_out(:, :), GP_out(:, :), UP_out(:, :)

        if (present(cuttof_value)) then
            call remesh_particles_3d(iflag, npar_per_cell, XP_out, QP_out, GP_out, UP_out, NVR_out, cuttof_value)
        else
            call remesh_particles_3d(iflag, npar_per_cell, XP_out, QP_out, GP_out, UP_out, NVR_out)
        end if
        XP_arr = from_intrinsic(XP_out, shape(XP_out))
        QP_arr = from_intrinsic(QP_out, shape(QP_out))
        GP_arr = from_intrinsic(GP_out, shape(GP_out))
        UP_arr = from_intrinsic(UP_out, shape(UP_out))
    end subroutine call_remesh_particles_3d

    subroutine write_particle_mesh_solution(folder, filename) bind(C, name='write_particle_mesh_solution')
        use file_io, only: write_pm_solution, case_folder, mesh_output_file
        use pmgrid, only: velocity_pm, deform_pm, RHS_pm, SOL_pm
        use vpm_vars, only: NTIME_pm, neqpm
        use vpm_size, only: fine_grid
        implicit none
        character(kind=c_char), intent(in), optional :: folder(*), filename(*)
        integer :: NN(3), NN_bl(6)

        if (present(folder)) then
            call set_string_f_c(case_folder, folder)
        end if

        if (present(filename)) then
            call set_string_f_c(mesh_output_file, filename)
        end if

        if (allocated(deform_pm)) then
            call write_pm_solution(NTIME_pm,fine_grid, neqpm, RHS_pm, SOL_pm, velocity_pm, deform_pm)
        else
            call write_pm_solution(NTIME_pm, fine_grid, neqpm, RHS_pm, SOL_pm, velocity_pm)
        end if
    end subroutine write_particle_mesh_solution

    subroutine write_particles_stored(folder, filename) bind(C, name='write_particles')
        use file_io, only: write_particles, case_folder, particle_output_file
        use parvar, only: XP, QP, UP, GP, NVR, NVR_size
        use vpm_vars, only: neqpm, NTIME_pm
        implicit none
        character(kind=c_char), intent(in), optional :: folder(*), filename(*)

        if (present(folder)) then
            call set_string_f_c(case_folder, folder)
        end if
        if (present(filename)) then
            call set_string_f_c(particle_output_file, filename)
        end if

        call write_particles(NTIME_pm, XP, UP, QP, GP, neqpm, NVR, NVR_size)
    end subroutine write_particles_stored

    subroutine write_particles_stored_hdf5(folder, filename) bind(C, name='write_particles_hdf5')
        use file_io, only: write_particles_hdf5, case_folder, particle_output_file
        use parvar, only: XP, QP, UP, GP, NVR, NVR_size
        use vpm_vars, only: neqpm, NTIME_pm
        implicit none
        character(kind=c_char), intent(in), optional :: folder(*), filename(*)
        if (present(folder)) then
            call set_string_f_c(case_folder, folder)
        end if
        if (present(filename)) then
            call set_string_f_c(particle_output_file, filename)
        end if
        call write_particles_hdf5(NTIME_pm, XP, UP, QP, GP, neqpm, NVR, NVR_size)
    end subroutine write_particles_stored_hdf5

    subroutine write_particle_mesh_solution_hdf5(folder, filename) bind(C, name='write_particle_mesh_solution_hdf5')
        use file_io, only: write_pm_solution_hdf5, case_folder, mesh_output_file
        use pmgrid, only: deform_pm, RHS_pm, SOL_pm, velocity_pm
        use vpm_vars, only: NTIME_pm, neqpm
        use vpm_size, only: fine_grid
        implicit none
        character(kind=c_char), intent(in), optional :: folder(*), filename(*)
        integer :: NN(3), NN_bl(6)

        if (present(folder)) then
            call set_string_f_c(case_folder, folder)
        end if

        if (present(filename)) then
            call set_string_f_c(mesh_output_file, filename)
        end if

        if (allocated(deform_pm)) then
            call write_pm_solution_hdf5(NTIME_pm, fine_grid, neqpm, RHS_pm, SOL_pm, velocity_pm, &
                                        deform_pm)
        else
            call write_pm_solution_hdf5(NTIME_pm, fine_grid, neqpm, RHS_pm, SOL_pm, velocity_pm)
        end if
    end subroutine write_particle_mesh_solution_hdf5

    subroutine get_particle_positions(XP_out) bind(C, name='get_particle_positions')
        use ND_Arrays
        use parvar, only: XP
        implicit none
        type(ND_Array), intent(out) :: XP_out
        XP_out = from_intrinsic(XP, shape(XP))
    end subroutine get_particle_positions

    subroutine get_particle_strengths(QP_out) bind(C, name='get_particle_strengths')
        use iso_c_binding
        use ND_Arrays
        use parvar, only: QP
        implicit none
        type(ND_Array), intent(out) :: QP_out
        QP_out = from_intrinsic(QP, shape(QP))
    end subroutine get_particle_strengths

    subroutine get_particle_deformation(GP_out) bind(C, name='get_particle_deformation')
        use iso_c_binding
        use ND_Arrays
        use parvar, only: GP
        implicit none
        type(ND_Array), intent(out) :: GP_out
        GP_out = from_intrinsic(GP, shape(GP))
    end subroutine get_particle_deformation

    subroutine get_particle_velocities(UP_out) bind(C, name='get_particle_velocities')
        use iso_c_binding
        use ND_Arrays
        use parvar, only: UP
        implicit none
        type(ND_Array), intent(out) :: UP_out
        UP_out = from_intrinsic(UP, shape(UP))
    end subroutine get_particle_velocities

    subroutine get_neqpm(neqpm_out) bind(C, name='get_neqpm')
        use vpm_vars, only: neqpm
        implicit none
        integer(c_int) :: neqpm_out

        neqpm_out = neqpm
    end subroutine get_neqpm

    subroutine get_velocity_pm(velocity_ptr) bind(C, name='get_velocity_pm')
        use pmgrid, only: velocity_pm
        implicit none
        real(c_double), dimension(:, :, :, :), pointer :: velocity_ptr
        velocity_ptr => velocity_pm
    end subroutine get_velocity_pm

    subroutine get_deformation_pm(deform_ptr) bind(C, name='get_deformation_pm')
        use pmgrid, only: deform_pm
        implicit none
        real(c_double), dimension(:, :, :, :), pointer :: deform_ptr
        deform_ptr => deform_pm
    end subroutine get_deformation_pm

!! VPM SIZE

!! LIBRARY PRINTS
    subroutine print_projlib() bind(C, name='print_projlib')
        use projlib, only: print_projlib_info
        implicit none

        call print_projlib_info()

    end subroutine print_projlib

    subroutine print_pmgrid() bind(C, name='print_pmgrid')
        use pmgrid, only: print_pmgrid_info
        implicit none

        call print_pmgrid_info()
    end subroutine print_pmgrid

    subroutine print_vpm_vars() bind(C, name='print_vpm_vars')
        use vpm_vars, only: print_vpm_vars_info
        implicit none

        call print_vpm_vars_info()
    end subroutine print_vpm_vars

    subroutine print_vpm_size() bind(C, name='print_vpm_size')
        use vpm_size, only: print_vpm_size_info
        implicit none

        call print_vpm_size_info()
    end subroutine print_vpm_size

    subroutine print_parvar() bind(C, name='print_parvar')
        use parvar, only: print_parvar_info
        implicit none

        call print_parvar_info()
    end subroutine print_parvar

!! STRING MANIPULATION
    subroutine get_string_f_c(fstring, cptr)
        use iso_c_binding
        implicit none
        type(c_ptr), intent(inout) :: cptr
        character(len=MAX_STRING_LENGTH), intent(in) :: fstring

        character(kind=c_char, len=MAX_STRING_LENGTH), target :: fstring_t
        integer :: i, trimmed_length, full_length

        fstring_t = trim(fstring)//c_null_char
        trimmed_length = len_trim(fstring_t)
        full_length = len(fstring_t)
        do i = trimmed_length + 1, 256
            fstring_t(i:i) = c_null_char
        end do
        cptr = c_loc(fstring_t)
    end subroutine get_string_f_c

    subroutine set_string_f_c(fstring, new_string)
        use iso_c_binding
        implicit none
        character(kind=c_char), intent(in) :: new_string(*)
        character(len=MAX_STRING_LENGTH), intent(out) :: fstring
        integer :: i

        fstring = ""
        do i = 1, 256
            if (new_string(i) == c_null_char) exit
            fstring(i:i) = new_string(i)
        end do
    end subroutine set_string_f_c
end Module api

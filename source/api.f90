module api
    use ND_Arrays
    ! use vpm_lib
    ! use vpm_vars
    ! use vpm_size
    ! use openmpth
    use, intrinsic :: iso_c_binding, only: c_float, c_int, c_bool, c_null_ptr, c_double, c_ptr

    integer, parameter :: MAX_STRING_LENGTH = 256
contains
    subroutine initialize(dx_pm, dy_pm, dz_pm, proj_type, bc_type, vol_type,                        &
                          num_coarse, num_nbi, num_nbj, num_nbk, remesh_type, tree_type,            &
                          max_level, omp_threads, grid_define, write_type, write_start,             &    
                          write_steps                                                               &
    ) bind(C, name='init')

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

    subroutine finalize() bind(C, name='finalize')
        use pmgrid, only: SOL_pm, RHS_pm, velocity_pm, deform_pm
        implicit none
        if (allocated(SOL_pm)) deallocate(SOL_pm)
        if (allocated(RHS_pm)) deallocate(RHS_pm)
        if (allocated(velocity_pm)) deallocate(velocity_pm)
        if (allocated(deform_pm)) deallocate(deform_pm)
    end subroutine finalize

    subroutine set_verbose_level(verbocity_in) bind(C, name='set_verbosity')
        use console_io, only: verbocity
        implicit none
        integer(c_int), intent(in) :: verbocity_in
        verbocity = verbocity_in
    end subroutine set_verbose_level

    subroutine set_case_folder(folder) bind(C, name='set_case_folder')
        use file_io, only: case_folder
        implicit none
        character(kind=c_char), intent(in) :: folder(*)
        call set_string_f_c(case_folder, folder)
    end subroutine set_case_folder

!!! VPM API 
    subroutine call_vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
                        RHS_pm_out, Vel_out, NTIME_in, NI_in,                   &
                        NVR_size_in                                             &
    ) bind(C, name='vpm')
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
        implicit none
        ! Interface for the arguments
        integer(c_int), intent(inout)          :: NVR_in
        integer(c_int), intent(in)             :: neqpm_in, WhatToDo, NVR_size_in, NTIME_in
        real(c_double), intent(in)             :: NI_in
        real(c_double), intent(inout), target  :: XP_in(3, NVR_in), QP_in(neqpm_in + 1, NVR_in)
        real(c_double), intent(inout), target  :: UP_in(3, NVR_in), GP_in(3, NVR_in)
        type(ND_Array), intent(out)            :: RHS_pm_out, Vel_out

        real(c_double), pointer                :: RHS_pm_ptr(:, :, :, :), Vel_ptr(:, :, :, :)

        call vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
                 RHS_pm_ptr, Vel_ptr, NTIME_in, NI_in, NVR_size_in)

        ! Assign the pointers to the arrays
        if (associated(RHS_pm_ptr)) RHS_pm_out = from_intrinsic(RHS_pm_ptr, shape(RHS_pm_ptr))
        if (associated(Vel_ptr)) Vel_out = from_intrinsic(Vel_ptr, shape(Vel_ptr))
    end subroutine call_vpm

    subroutine call_vpm_project_solve(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in, &
                                      RHS_pm_out                                             &
    ) bind(C, name='vpm_project_solve')
        use vpm_lib, only: vpm_project_and_solve
        use ND_Arrays
        implicit none
        ! Interface for the arguments
        integer(c_int), intent(in)             :: neqpm_in, NVR_size_in, NTIME_in, NVR_in
        real(c_double), intent(inout), target  :: XP_in(3, NVR_in), QP_in(neqpm_in + 1, NVR_in)
        type(ND_Array), intent(out)            :: RHS_pm_out
        ! Local variables
        real(c_double), pointer                :: RHS_pm_ptr(:, :, :, :)
        call vpm_project_and_solve(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, &
                                   neqpm_in, RHS_pm_ptr)
        if (associated(RHS_pm_ptr)) RHS_pm_out = from_intrinsic(RHS_pm_ptr, shape(RHS_pm_ptr))
    end subroutine call_vpm_project_solve

    subroutine call_vpm_define(                                                             &
        NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in                               &
    ) bind(C, name='vpm_define')
        use vpm_lib, only: vpm_define_problem
        use ND_Arrays
        implicit none
        ! Interface for the arguments
        integer(c_int), intent(in)             :: neqpm_in, NVR_size_in, NTIME_in, NVR_in
        real(c_double), intent(inout), target  :: XP_in(3, NVR_in), QP_in(neqpm_in + 1, NVR_in)
        call vpm_define_problem(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in)
    end subroutine call_vpm_define

    subroutine call_vpm_solve_velocity(NTIME_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, &
                                       neqpm_in, RHS_pm_out, Vel_out                              &
    ) bind(C, name='vpm_solve_velocity')
        use vpm_lib, only: vpm_solve_velocity
        use ND_Arrays
        implicit none
        ! Interface for the arguments 
        integer(c_int), intent(in)             :: neqpm_in, NVR_size_in, NTIME_in, NVR_in
        real(c_double), intent(inout), target  :: XP_in(3, NVR_in), QP_in(neqpm_in + 1, NVR_in)
        real(c_double), intent(inout), target  :: UP_in(3, NVR_in), GP_in(3, NVR_in)
        type(ND_Array), intent(out)            :: RHS_pm_out, Vel_out
        ! Local variables
        real(c_double), pointer                :: RHS_pm_ptr(:, :, :, :), Vel_ptr(:, :, :, :)
        
        call vpm_solve_velocity(NTIME_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, &
                                neqpm_in, RHS_pm_ptr, vel_ptr)
        if (associated(RHS_pm_ptr)) RHS_pm_out = from_intrinsic(RHS_pm_ptr, shape(RHS_pm_ptr))
        if (associated(Vel_ptr)) Vel_out = from_intrinsic(Vel_ptr, shape(Vel_ptr))
    end subroutine call_vpm_solve_velocity

    subroutine call_vpm_solve_velocity_deformation(NTIME_in, XP_in, QP_in, UP_in, GP_in, NVR_in, &
                                                   NVR_size_in,                                  &
                                                   neqpm_in, RHS_pm_out, Vel_out, Deform_out     &
    ) bind(C, name='vpm_solve_velocity_deformation')
        use vpm_lib, only: vpm_solve_velocity_deformation
        use ND_Arrays
        implicit none
        ! Interface for the arguments
        integer(c_int), intent(in)             :: neqpm_in, NVR_size_in, NTIME_in, NVR_in
        real(c_double), intent(inout), target  :: XP_in(3, NVR_in), QP_in(neqpm_in + 1, NVR_in)
        real(c_double), intent(inout), target  :: UP_in(3, NVR_in), GP_in(3, NVR_in)
        type(ND_Array), intent(out)            :: RHS_pm_out, Vel_out, Deform_out
        ! Local variables
        real(c_double), pointer                :: RHS_pm_ptr(:, :, :, :), Vel_ptr(:, :, :, :)
        real(c_double), pointer                :: deform_ptr(:, :, :, :)
        
        call vpm_solve_velocity_deformation(NTIME_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, &
                                            neqpm_in, RHS_pm_ptr, vel_ptr, deform_ptr)
        if (associated(RHS_pm_ptr)) RHS_pm_out = from_intrinsic(RHS_pm_ptr, shape(RHS_pm_ptr))
        if (associated(Vel_ptr)) Vel_out = from_intrinsic(Vel_ptr, shape(Vel_ptr))
        if ((associated(deform_ptr)))  Deform_out = from_intrinsic(deform_ptr, shape(deform_ptr))
    end subroutine call_vpm_solve_velocity_deformation

    subroutine call_vpm_interpolate(NTIME_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in,  &
                                    neqpm_in, RHS_pm_out                                        &
    ) bind(C, name='vpm_interpolate')
        use vpm_lib, only: vpm_solve_velocity_deformation
        use ND_Arrays
        implicit none
        ! Interface for the arguments
        integer(c_int), intent(in)             :: neqpm_in, NVR_size_in, NTIME_in, NVR_in
        real(c_double), intent(inout), target  :: XP_in(3, NVR_in), QP_in(neqpm_in + 1, NVR_in)
        real(c_double), intent(inout), target  :: UP_in(3, NVR_in), GP_in(3, NVR_in)
        type(ND_Array), intent(out)            :: RHS_pm_out
        ! Local variables
        real(c_double), pointer                :: RHS_pm_ptr(:, :, :, :), Vel_ptr(:, :, :, :)
        real(c_double), pointer                :: deform_ptr(:, :, :, :)
        
        call vpm_solve_velocity_deformation(NTIME_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, &
                                            neqpm_in, RHS_pm_ptr, vel_ptr, deform_ptr)
        if (associated(RHS_pm_ptr)) RHS_pm_out = from_intrinsic(RHS_pm_ptr, shape(RHS_pm_ptr))
    end subroutine call_vpm_interpolate

    subroutine call_vpm_diffuse(NI_in , XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, &
                                neqpm_in, RHS_pm_out                                     &
    ) bind(C, name='vpm_diffuse')
        use vpm_lib, only: vpm_diffuse
        use ND_Arrays
        implicit none
        ! Interface for the arguments
        real(c_double), intent(in)             :: NI_in
        integer(c_int), intent(in)             :: NVR_size_in, NVR_in, neqpm_in
        real(c_double), intent(inout), target  :: XP_in(3, NVR_in), QP_in(neqpm_in + 1, NVR_in)
        real(c_double), intent(inout), target  :: UP_in(3, NVR_in), GP_in(3, NVR_in)
        type(ND_Array), intent(out)            :: RHS_pm_out
        ! Local variables
        real(c_double), pointer                :: RHS_pm_ptr(:, :, :, :)
        
        call vpm_diffuse(NI_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, RHS_pm_ptr)
        if (associated(RHS_pm_ptr)) RHS_pm_out = from_intrinsic(RHS_pm_ptr, shape(RHS_pm_ptr))
        end subroutine call_vpm_diffuse

    subroutine call_vpm_correct_vorticity(                                                      &
        XP_in, QP_in,  NVR_in, neqpm_in, NVR_size_in                                             &
    ) bind(C, name='vpm_correct_vorticity')
        use vpm_lib, only: vpm_correct_vorticity
        use ND_Arrays
        implicit none
        ! Interface for the arguments
        integer(c_int), intent(in)             :: NVR_size_in, NVR_in, neqpm_in
        real(c_double), intent(inout), target  :: XP_in(3, NVR_in), QP_in( neqpm_in + 1, NVR_in)
        ! Local variables
        call vpm_correct_vorticity(XP_in, QP_in, NVR_in, NVR_size_in)
    end subroutine call_vpm_correct_vorticity

    subroutine call_vpm_solve_pressure(                                                         &
        vorticity_ptr, velocity_ptr, pressure_ptr, density                                      &
    ) bind(C, name='vpm_solve_pressure')
        use vpm_lib, only: vpm_solve_pressure
        use ND_Arrays
        implicit none
        ! Interface for the arguments (ND_Arrays)
        type(ND_Array), intent(in)             :: vorticity_ptr, velocity_ptr
        type(ND_Array), intent(out)            :: pressure_ptr
        real(c_double), intent(in)             :: density
        ! Local variables
        real(c_double), pointer                :: vorticity(:,:,:,:), velocity(:,:,:,:)
        real(c_double), allocatable            :: pressure(:,:,:,:)

        if (vorticity_ptr%total_size .ne. 0 ) then
            call convert_to_4D_array(vorticity_ptr, vorticity)
        end if
        
        if (velocity_ptr%total_size .ne. 0 ) then
            call convert_to_4D_array(velocity_ptr, velocity)
        end if
        call vpm_solve_pressure(vorticity, velocity, pressure, density)
        if (allocated(pressure)) then
            pressure_ptr = from_intrinsic(pressure, shape(pressure))
        end if
    end subroutine call_vpm_solve_pressure

    subroutine call_remesh_particles_3d(iflag, npar_per_cell, XP_arr, QP_arr, GP_arr, UP_arr,   &
        NVRR, cuttof_value &
    )bind(C, name='remesh_particles_3d')
        use vpm_remesh, only: remesh_particles_3d, interpolate_and_remesh_particles
        use vpm_types, only: dp
        use vpm_size, only: fine_grid
        use vpm_vars, only: neqpm
        use pmgrid, only: RHS_pm
        use ND_Arrays
        use MPI

        implicit none
        integer(c_int), intent(in) :: iflag, npar_per_cell
        type(ND_Array), intent(inout), target :: XP_arr, QP_arr, GP_arr, UP_arr
        integer(c_int), intent(inout) :: NVRR
        real(c_double), intent(in), optional :: cuttof_value
        
        ! Local variables
        real(dp), allocatable, target, save :: XP_out(:, :), QP_out(:, :), GP_out(:, :), UP_out(:, :)
        real(dp), dimension(:,:), pointer :: XP_ptr,  QP_ptr, GP_ptr, UP_ptr
        integer :: NVR_size_in, my_rank, ierr

        call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

        if ( iflag == 1 )then
            call convert_to_2D_array(XP_arr, XP_ptr)
            call convert_to_2D_array(QP_arr, QP_ptr)
            call convert_to_2D_array(GP_arr, GP_ptr)
            call convert_to_2D_array(UP_arr, UP_ptr) 
            
            if (allocated(XP_out)) deallocate(XP_out)
            if (allocated(QP_out)) deallocate(QP_out)
            if (allocated(GP_out)) deallocate(GP_out)
            if (allocated(UP_out)) deallocate(UP_out)

            NVR_size_in = size(XP_ptr, 2)
            allocate(XP_out(3, NVR_size_in)); XP_out = XP_ptr
            allocate(UP_out(3, NVR_size_in)); UP_out = UP_ptr
            allocate(GP_out(3, NVR_size_in)); GP_out = GP_ptr
            allocate(QP_out(neqpm + 1, NVR_size_in)); QP_out = QP_ptr

            ! if (my_rank.eq.0) then
            !     print *, 'Number of particles before', NVRR
            !     print *, 'Number of particles after', NVR_size_in
            !     print *, 'XP -> min and max', minval(XP_out), maxval(XP_out) 
            !     print *, 'QP -> min and max', minval(QP_out), maxval(QP_out)
            ! end if

            if (present(cuttof_value)) then
                call interpolate_and_remesh_particles(npar_per_cell, XP_out, QP_out, UP_out, GP_out, NVRR, NVR_size_in, cuttof_value)
            else
                call interpolate_and_remesh_particles(npar_per_cell, XP_out, QP_out, UP_out, GP_out, NVRR, NVR_size_in)
            end if
        else
            if (present(cuttof_value)) then
                call remesh_particles_3d(RHS_pm, fine_grid, npar_per_cell, XP_out, QP_out, UP_out, GP_out, NVRR, cuttof_value)
            else
                call remesh_particles_3d(RHS_pm, fine_grid, npar_per_cell, XP_out, QP_out, UP_out, GP_out, NVRR)
            end if
        end if


        ! if (my_rank.eq.0) then
        !     print *, "EXITING REMESH"
        !     print *, 'Number of particles before', NVR_size_in
        !     print *, 'Number of particles after', NVRR
        !     print *, 'XP -> min and max', minval(XP_out), maxval(XP_out), "Shape", shape(XP_out)
        !     print *, 'QP -> min and max', minval(QP_out), maxval(QP_out), "Shape", shape(QP_out)
        ! end if
        XP_arr = from_intrinsic(XP_out, shape(XP_out))
        QP_arr = from_intrinsic(QP_out, shape(QP_out))
        GP_arr = from_intrinsic(GP_out, shape(GP_out))
        UP_arr = from_intrinsic(UP_out, shape(UP_out))
    end subroutine call_remesh_particles_3d

    subroutine vpm_timestep_information() bind(C, name='get_timestep_information')
        ! use vpm_vars, only: timestep_info 
        implicit none

    end subroutine vpm_timestep_information 

!! FILE IO
    subroutine write_particle_mesh_solution(folder, filename) bind(C, name='write_particle_mesh_solution')
        use file_io, only: write_pm_solution, case_folder, mesh_output_file
        use pmgrid, only: velocity_pm, deform_pm, RHS_pm, SOL_pm
        use vpm_vars, only: NTIME_pm, neqpm
        use vpm_size, only: fine_grid
        implicit none
        character(kind=c_char), intent(in), optional :: folder(*), filename(*)

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

        if (present(folder)) then
            call set_string_f_c(case_folder, folder)
        end if

        if (present(filename)) then
            call set_string_f_c(mesh_output_file, filename)
        end if

        if (allocated(deform_pm)) then
            call write_pm_solution_hdf5(NTIME_pm, fine_grid%NN, fine_grid%NN_bl, neqpm, RHS_pm, SOL_pm, velocity_pm, &
                                        deform_pm)
        else
            call write_pm_solution_hdf5(NTIME_pm, fine_grid%NN, fine_grid%NN_bl, neqpm, RHS_pm, SOL_pm, velocity_pm)
        end if
    end subroutine write_particle_mesh_solution_hdf5

    subroutine write_pressure_field_hdf5(folder, filename, pressure) bind(C, name='write_pressure_hdf5')
        use file_io, only: write_field_h5, case_folder, field_output_file
        use vpm_size, only: fine_grid
        implicit none
        type(ND_Array), intent(in) :: pressure
        character(kind=c_char), intent(in), optional :: folder(*), filename(*)
        real(dp), pointer :: pressure_ptr(:,:,:,:)

        if (present(folder)) then
            call set_string_f_c(case_folder, folder)
        end if
        if (present(filename)) then
            call set_string_f_c(field_output_file, filename)
        end if
        call convert_to_4D_array(pressure, pressure_ptr)
        call write_field_h5(2, fine_grid, pressure_ptr, "pressure")
    end subroutine write_pressure_field_hdf5 

!! GETTERS
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

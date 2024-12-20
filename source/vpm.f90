module vpm_lib
    ! SPEED
#ifdef USE_INTEL
    use mkl_service
#endif
    use MPI

    use vpm_vars
    use vpm_size
    use vpm_gcalc
    use vpm_functions
    ! Constants
    use constants, only: pi, pi2, pi4
    use vpm_types, only: dp
    ! Printing
    use console_io, only: vpm_print, red, blue, green, nocolor, yellow, dummy_string, tab_level, VERBOCITY
    use parvar, only: print_particle_info, print_particle_positions, associate_particles
    ! Setting Vars
    use pmgrid, only: print_velocity_stats, print_vortex_stretching_stats, &
                      set_pm_velocities_zero, set_pm_deformations_zero, &
                      associate_velocities, associate_deformations

    !  WhatToDo flags
    integer, parameter :: DEFINE_PROBLEM = 0,               &
                          SOLVE_VELOCITY = 1,               &
                          SOLVE_VELOCITY_DELATATION = 2,    &
                          PROJECT_AND_SOLVE = 3,            &
                          INTERPOLATE = 4,                  &
                          DIFFUSE = 5,                      &
                          CORRECT_VORTICITY = 6
contains

    subroutine set_num_threads()
#ifdef USE_INTEL
        if (SOLVER .eq. 1) then
            call mkl_set_num_threads(OMPTHREADS)
        end if
#endif
    end subroutine set_num_threads

    subroutine vpm( &
        XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
        RHS_pm_ptr, vel_ptr, NTIME_in, NI_in, NVR_size_in       &
    )
        use parvar, only: NVR
        implicit None

        ! ARGUEMENTS
        real(dp), intent(inout), target  :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
        integer, intent(in)              :: NVR_in
        integer, intent(in)              :: neqpm_in, WhatToDo, NVR_size_in, NTIME_in
        real(dp), intent(in)             :: NI_in
        real(dp), intent(out), pointer   :: RHS_pm_ptr(:, :, :, :)
        real(dp), intent(out), pointer   :: vel_ptr(:, :, :, :)

        ! LOCAL VARIABLES
        integer                           :: ierr, my_rank, np
        integer                           :: i, j, k

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        call set_num_threads()

        ! Print the input arguments
        if ((my_rank .eq. 0) .and. (VERBOCITY >= 1)) then
            if (WhatToDo .eq. 0) then
                write (dummy_string, "(A)") 'VPM: Define'
            else if (WhatToDo .eq. 1) then
                write (dummy_string, "(A)") 'VPM: Project Particles- Solve PM- Get Velocities on PM'
            else if (WhatToDo .eq. 2) then
                write (dummy_string, "(A)") 'VPM: Project Particles- Solve PM- Get Velocities and Deformation on PM'
            else if (WhatToDo .eq. 3) then
                write (dummy_string, "(A)") 'VPM: Project Particles- Solve PM '
            else if (WhatToDo .eq. 4) then
                write (dummy_string, "(A)") 'VPM: Project Particles- Interpolate from PM to particles '
            else if (WhatToDo .eq. 5) then
                write (dummy_string, "(A)") 'VPM: Project Particles- Diffuse Particle vorticity '
            else if (WhatToDo .eq. 6) then
                write (dummy_string, "(A)") 'VPM: Project Particles- Correct Vorticity '
            end if
            call vpm_print(dummy_string, red, 1)

            if (VERBOCITY < 1) then
                write (*, *) ""
                return
            end if

            write (*, *) achar(9), 'Input Arguments:', achar(27)//'[1;34m'
            write (*, *) achar(9), achar(9), 'NTIME_in = ', NTIME_in
            write (*, *) achar(9), achar(9), 'WhatToDo = ', WhatToDo
            write (*, *) achar(9), achar(9), 'NVR_in = ', NVR_in
            write (*, *) achar(9), achar(9), 'neqpm_in = ', neqpm_in
            write (*, *) achar(9), achar(9), 'NI_in = ', NI_in
            write (*, *) achar(9), achar(9), 'NVR_size_in = ', NVR_size_in, achar(27)//'[0m'
        end if
        ! We set these input for compatibility with old version
        ! In the new modular approach we do not need to set these
        ND = 3
        NTIME_pm = NTIME_in
        neqpm = neqpm_in
        NI = NI_in
        call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in, UP_in, GP_in)
        if (NVR .eq. 0) return
        if (associated(RHS_pm_ptr)) nullify (RHS_pm_ptr)
        call associate_velocities(vel_ptr)

        select case (WhatToDo)
        case (DEFINE_PROBLEM)
            call vpm_define_problem(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in)
        case (SOLVE_VELOCITY)
            call vpm_solve_velocity(NTIME_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, &
                                    neqpm_in, RHS_pm_ptr, vel_ptr)

        case (SOLVE_VELOCITY_DELATATION)
            call vpm_solve_velocity_deformation(NTIME_in, XP_in, QP_in, UP_in, GP_in, NVR_in, &
                                                NVR_size_in, neqpm_in, RHS_pm_ptr, vel_ptr )

        case (PROJECT_AND_SOLVE)
            call vpm_project_and_solve(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in, &
                                       RHS_pm_ptr)

        case (INTERPOLATE)
            call vpm_interpolate(XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, neqpm_in, &
                                 RHS_pm_ptr)

        case (DIFFUSE)
            call vpm_diffuse(NI_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, neqpm_in, &
                             RHS_pm_ptr)

        case (CORRECT_VORTICITY)
            call vpm_correct_vorticity(XP_in, QP_in, NVR_in, NVR_size_in)

        end select
        return
    end subroutine vpm

    subroutine vpm_define_problem( &
        NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in &
    )
        use parvar, only: NVR
        integer, intent(in) :: NTIME_in, NVR_in, NVR_size_in, neqpm_in
        real(dp), intent(in), dimension(:, :), target :: XP_in, QP_in
        integer :: ierr, my_rank, np
        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        tab_level = 0

        ! Set Input
        ND = 3
        NTIME_pm = NTIME_in
        neqpm = neqpm_in
        call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in)
        if (NVR .eq. 0) return

        call define_sizes
        call set_pm_velocities_zero
        call set_pm_deformations_zero
        call allocate_sol_and_rhs(my_rank + 1)
    end subroutine vpm_define_problem

    subroutine vpm_interpolate( &
        XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, neqpm_in, &
        RHS_pm_ptr &
    )
        use pmgrid, only: RHS_pm, velocity_pm, deform_pm
        use parvar, only: NVR
        use MPI
        implicit none

        real(dp), intent(inout), target           :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
        integer, intent(in)                       :: NVR_in
        integer, intent(in)                       :: neqpm_in, NVR_size_in
        real(dp), intent(out), pointer            :: RHS_pm_ptr(:, :, :, :)
        integer :: ierr, my_rank, np

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        tab_level = 0

        ! Set Input
        ND = 3
        neqpm = neqpm_in
        call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in, UP_in, GP_in)
        call allocate_sol_and_rhs(my_rank + 1)
        RHS_pm_ptr => RHS_pm

        if (NVR .eq. 0) then
            print *, 'No particles to interpolate'
            return
        end if

        ! SOL_pm AND RHS_pm ARE 0 BOTH ON THE BL AND PM
        ! PARTICLES TO -> PM MEANING FROM Qs We get RHS_pm
        call project_particles_parallel

        ! We assume that the problem is solved and the velocity field is calculated
        if (my_rank .eq. 0) then
            ! call calc_velocity_serial_3d(-1) ! ONLY DEFORMATION
            call calc_vortex_stretching_conservative(velocity_pm, deform_pm)
        end if
        call interpolate_particles_parallel(1)
    end subroutine vpm_interpolate

    subroutine vpm_diffuse( &
        NI_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, neqpm_in, &
        RHS_pm_ptr &
    )
        use MPI
        use pmgrid, only: RHS_pm
        use parvar, only: NVR
        implicit none
        real(dp), intent(inout), target           :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
        integer, intent(in)                       :: NVR_in
        integer, intent(in)                       :: neqpm_in, NVR_size_in
        real(dp), intent(in)                      :: NI_in
        real(dp), intent(out), pointer            :: RHS_pm_ptr(:, :, :, :)
        integer :: ierr, my_rank, np

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        tab_level = 0

        ! Set Input
        ND = 3
        neqpm = neqpm_in
        NI = NI_in
        call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in, UP_in, GP_in)
        call allocate_sol_and_rhs(my_rank + 1)
        RHS_pm_ptr => RHS_pm

        if (NVR .eq. 0) then
            print *, 'No particles to interpolate'
            return
        end if

        ! PARTICLES TO -> PM MEANING FROM Qs We get RHS_pm
        call project_particles_parallel

        call set_pm_velocities_zero
        call set_pm_deformations_zero
        if (my_rank .eq. 0) then
            ! diffusion stores -NI*grad^2 w * Vol in GP(1,:)
            ! RHS_pm = -VIS \nabla \cdot RHS_pm
            call diffuse_vort_3d ! DIFFUSION OF VORTICITY
        end if
        ! WHEN ITYPEB = 2 WE GET THE GP FROM THE SOL_pm (DEFORMATION) and from QP
        call interpolate_particles_parallel(2)
    end subroutine vpm_diffuse

    subroutine vpm_project_and_solve( &
        NTIME_in, &
        XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in, &
        RHS_pm_ptr &
    )
        use MPI
        use pmgrid, only: RHS_pm, SOL_pm, SOL_pm_bl, RHS_pm_bl
        use parvar, only: NVR
        implicit none
        real(dp), intent(inout), target  :: XP_in(:, :), QP_in(:, :)
        integer, intent(in)              :: NVR_in
        integer, intent(in)              :: neqpm_in, NVR_size_in
        real(dp), intent(out), pointer   :: RHS_pm_ptr(:, :, :, :)
        integer, intent(in)              :: NTIME_in

        integer :: ierr, my_rank, np
        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        ! Set Input
        tab_level = 0
        ND = 3
        NTIME_pm = NTIME_in
        neqpm = neqpm_in
        call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in)
        call allocate_sol_and_rhs(my_rank + 1)
        RHS_pm_ptr => RHS_pm

        if (NVR .eq. 0) then
            print *, 'No particles to interpolate'
            return
        end if

        ! PARTICLES TO -> PM MEANING FROM Qs We get RHS_pm
        if (my_rank .eq. 0) then
            write (*, *) ""
            tab_level = tab_level - 1
            write (dummy_string, "(A)") 'Solving Primary Problem'
            call vpm_print(dummy_string, blue, 0)
            tab_level = tab_level + 1
            st = MPI_WTIME()
        end if
        call project_particles_parallel
        call solve_problem(RHS_pm, SOL_pm, RHS_pm_bl, SOL_pm_bl)
        ! Zero the velocities and deformations
        call set_pm_velocities_zero
        call set_pm_deformations_zero
    end subroutine vpm_project_and_solve

    subroutine vpm_solve_velocity( &
        NTIME_in, XP_in, QP_in, UP_in, GP_in, &
        NVR_in, NVR_size_in, neqpm_in, RHS_pm_ptr, vel_ptr &
    )
        use MPI
        real(dp), intent(inout), target  :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
        integer, intent(in)              :: NVR_in
        integer, intent(in)              :: neqpm_in, NVR_size_in
        real(dp), intent(out), pointer   :: RHS_pm_ptr(:, :, :, :)
        real(dp), intent(out), pointer   :: vel_ptr(:, :, :, :)
        integer, intent(in)              :: NTIME_in
        integer :: ierr, my_rank, np

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        call vpm_project_and_solve(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in, RHS_pm_ptr)
        call associate_velocities(vel_ptr)

        if (my_rank .eq. 0) then
            !call convect_first_order(Xbound,Dpm,NN,NN_bl)
            ! FROM THE SOLUTION OF PM WE GET THE VELOCITIES ON THE GRID
            call calc_velocity ! VELOCITY STO PM
            ! SOL_pm is stil the solution of vorticity
            call print_velocity_stats
        end if
        call interpolate_particles_parallel(1) ! INTERPOLATION FROM PM TO PARTICLES
    end subroutine vpm_solve_velocity

    subroutine vpm_solve_velocity_deformation( &
        NTIME_in, &
        XP_in, QP_in, UP_in, GP_in, &
        NVR_in, NVR_size_in, neqpm_in, &
        RHS_pm_ptr, vel_ptr, deform_ptr &
    )
        use MPI
        use pmgrid, only: velocity_pm, deform_pm, RHS_pm
        use console_io, only: tab_level
        implicit none
        integer, intent(in)                        :: NTIME_in
        real(dp), intent(inout), target            :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
        integer, intent(in)                        :: NVR_in
        integer, intent(in)                        :: neqpm_in, NVR_size_in
        real(dp), intent(out), pointer             :: RHS_pm_ptr(:, :, :, :)
        real(dp), intent(out), pointer             :: vel_ptr(:, :, :, :)
        real(dp), intent(out), pointer, optional   :: deform_ptr(:, :, :, :)

        ! LOCAL VARIABLES
        real(dp), pointer       :: RHS2_pm_ptr(:, :, :, :)
        real(dp), pointer       :: SOL2_pm_ptr(:, :, :, :)

        integer                 :: ierr, my_rank

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

        ! Project the particles to the PM grid and solve the problem
        call vpm_project_and_solve(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in, RHS_pm_ptr)

        ! Calculate the velocity and deformation on the PM grid
        call associate_velocities(vel_ptr)
        if (present(deform_ptr)) call associate_deformations(deform_ptr)
        call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in, UP_in, GP_in)
        if (my_rank .eq. 0) then
            print *, ''
            ! FROM THE SOLUTION OF PM WE GET THE VELOCITIES ON THE GRID
            call calc_velocity ! VELOCITY AND DEFORMATION STO PM
            call calc_vortex_stretching_conservative(velocity_pm, deform_pm) ! VELOCITY AND DEFORMATION STO PM
        end if

        if (my_rank .eq. 0) then
            call print_velocity_stats
            call print_vortex_stretching_stats
            tab_level = tab_level - 1
            call print_timestep_information
        end if

        call interpolate_particles_parallel(1) ! INTERPOLATION FROM PM TO PARTICLES
    end subroutine vpm_solve_velocity_deformation

    subroutine vpm_correct_vorticity(       &
        XP_in, QP_in, NVR_in, NVR_size_in  &
    )
        use parvar, only: QP, XP, UP, GP, NVR, NVR_size
        use serial_vector_field_operators, only: divergence, curl
        use vpm_interpolate, only: interpolate_particle_Q
        use vpm_remesh, only: remesh_particles_3d
        use file_io, only: write_pm_solution_hdf5, mesh_output_file, write_vector_field, write_scalar_field
        use pmgrid, only: RHS_pm
        implicit none

        ! ARGUEMENTS
        real(dp), intent(inout), target  :: XP_in(:, :), QP_in(:, :)
        integer, intent(in)              :: NVR_in
        integer, intent(in)              :: NVR_size_in
        
        ! LOCAL PARAMETERS
        real(dp), allocatable, target    :: SOL2_pm(:, :, :, :), RHS2_pm(:, :, :, :)
        real(dp), allocatable, target    :: SOL2_pm_bl(:, :, :, :), RHS2_pm_bl(:, :, :, :)
        ! INTEMEIDATE VARIABLES
        real(dp), allocatable            :: div_wmega(:, :, :)
        real(dp), allocatable            :: grad_pi_fd_nodes(:, :, :, :)
        real(dp)                         :: DXpm, DYpm, DZpm    
        ! LOCAL VARIABLES
        integer                          :: ierr, my_rank, np, nb
        integer                          :: i, j, k 
        character(len=256)               :: file_name
        ! OlD parameters
        real(dp), allocatable            :: XP_old(:,:), QP_old(:,:)
        real(dp), pointer                :: vorticity(:,:,:,:) 
        
        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        
        tab_level = 0
        if (my_rank .eq. 0) then
            write (*, *) ""
            tab_level = tab_level - 1
            write (dummy_string, "(A)") 'Solving Poisson for the vorticity field correction'
            call vpm_print(dummy_string, blue, 0)
            tab_level = tab_level + 1
            st = MPI_WTIME()
            vorticity => RHS_pm(1:3, :, :, :) 
        end if

        call project_particles_parallel
        DXpm = fine_grid%Dpm(1)
        DYpm = fine_grid%Dpm(2)
        DZpm = fine_grid%Dpm(3)
        
        call associate_particles(NVR_in, NVR_size_in, neqpm, XP_in, QP_in)
        ! Allocate the solution and the RHS
        ND = 3
        neqpm = 1
        allocate (SOL2_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
        allocate (RHS2_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
        nb = my_rank + 1
        allocate (SOL2_pm_bl(neqpm, block_grids(nb)%NN(1), block_grids(nb)%NN(2), block_grids(nb)%NN(3)))
        allocate (RHS2_pm_bl(neqpm, block_grids(nb)%NN(1), block_grids(nb)%NN(2), block_grids(nb)%NN(3)))

        
        SOL2_pm = 0.d0
        RHS2_pm = 0.d0
        ! Get the divergence of the RHS of the equations:
        ! 1) 1st equation \nabla ^ 2 \pi = \nabla \cdot \omega that gives the pi field decomposition of
        !    the vorticity field. \omega = \nabla pi + \naabla x \qi => \nabla \cdot \omega = \nabla^2 pi
        if (my_rank .eq. 0) then
            ! FINITE DIFFERENCES VERSION
            RHS2_pm(1, :, :, :) = divergence(vorticity(1:3, :, :, :), DXpm, DYpm, DZpm)
            print *, 'The divergence of the vorticity field'
            print *, achar(9)//'Max divergence of the vorticity field', maxval(abs(RHS2_pm(1, :, :, :)))
            print *, achar(9)//'Total divergence of the vorticity field', sum(abs(RHS2_pm(1, :, :, :)))

            ! GREEN GAUSS THEOREM VERSION
            !     We need to perform the green gauss theorem to calculate div(omega) in the center of the cells
            !     wmega_vel = curl(velocity_pm, DXpm, DYpm, DZpm)
            !     call compute_gauss_divergence(fine_grid,vorticity(1:3,:,:,:), RHS2_pm(1, :, :, :))
            
            ! BIOT-SAVART VERSION
        end if

        ! Solve the problem
        call solve_problem(RHS2_pm, SOL2_pm, RHS2_pm_bl, SOL2_pm_bl)

        if (my_rank .eq. 0) then
            et = MPI_WTIME()
            write (dummy_string, "(A,I5,A,F8.2,A)") &
                'Total time for solving the Secondary Problem', &
                int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
            call vpm_print(dummy_string, blue, 1)
            ! Calculate the gradient of the solution to get the correction to the omega field
            ! Use the finite differences to calculate the gradient of \pi 
            grad_pi_fd_nodes = gradient(SOL2_pm(1, :, :, :), DXpm, DYpm, DZpm)
            vorticity(1:3, :, :, :) = vorticity(1:3, :, :, :) - grad_pi_fd_nodes
            ! Calculate the divergence of the corrected vorticity field
            div_wmega = divergence(vorticity, DXpm, DYpm, DZpm)
            print *, 'Correction to the vorticity field using the normal gradient'
            print *, achar(9)//'Max divergence of the corrected vorticity field', maxval(abs(div_wmega))
            print *, achar(9)//'Total divergence of the corrected vorticity field', sum(abs(div_wmega))

            deallocate (div_wmega)
            deallocate (grad_pi_fd_nodes)
            ! call interpolate_particle_Q(vorticity, XP, QP, NVR, 4, NVR_size)
        end if
        
        neqpm = 3
        if(my_rank .eq. 0) then
            allocate (XP_old, mold = XP)
            allocate (QP_old, mold = QP)
            XP_old = XP
            QP_old = QP
            call interpolate_particle_Q(vorticity, XP, QP, NVR, 4, NVR_size)
            ! Print the difference between the old and the new Q
            print *, 'Difference between the old and the new Q'
            print *, achar(9)//'Max percentage difference in Q', &
                                100 * maxval(abs(QP - QP_old))/maxval(abs(QP_old)), " %"
            print *, achar(9)//'Tot Percentage difference in Q', &
                                100 - sum(abs(QP))/sum(abs(QP_old)) * 100, " %"
            deallocate (XP_old, QP_old)
        end if
        deallocate (SOL2_pm, RHS2_pm, SOL2_pm_bl, RHS2_pm_bl)
    end subroutine vpm_correct_vorticity

    ! !> Solve the secondary problem:
    ! !> 1) Pressure Poisson Equation
    ! !> 2) Omega Correction Equation
    ! subroutine vpm_correct_vorticity_and_calculate_pressure(

    ! )
    !     ! 2) 2nd equation is (\nabla ^2 q / rho ) = - \nabla (u x \omega) - \nabla (\sigma ) / rho
    !     !    where sigma is the stress tensor which we ignore for now
    !     allocate (cross_product(3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
    !     call compute_cross_product(velocity, vorticity, cross_product)
    !     RHS2_pm(2, :, :, :) = divergence(cross_product, DXpm, DYpm, DZpm)
    !     deallocate (cross_product)
    ! end subroutine vpm_correct_vorticity_and_calculate_pressure

    subroutine project_calc_div
        use MPI
        use parvar, only: QP, XP, UP, GP, NVR, NVR_size
        use serial_vector_field_operators, only: divergence, curl
        use vpm_interpolate, only: interpolate_particle_Q
        use vpm_remesh, only: remesh_particles_3d
        use file_io, only: write_pm_solution_hdf5, mesh_output_file, write_vector_field, write_scalar_field
        use pmgrid, only: RHS_pm, DXpm, DYpm, DZpm
        implicit none

        ! LOCAL VARIABLES
        integer                           :: ierr, my_rank, np, i
        real(dp), allocatable             :: div_wmega(:, :, :)
        real(dp), allocatable             :: vorticity_old(:, :, :, :)
        real(dp), allocatable             :: Qp_old(:,:)
        
        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        if (my_rank .eq. 0) then
            allocate (vorticity_old, mold = RHS_pm)
            allocate (Qp_old, mold = QP)
            vorticity_old = RHS_pm
            Qp_old = QP
        end if

        do i = 1,10
            if (my_rank.eq.0) then
                call interpolate_particle_Q(RHS_pm(1:3, :,:,:), XP, QP, NVR, 4, NVR_size)
                div_wmega = divergence(RHS_pm(1:3,:,:,:), DXpm, DYpm, DZpm)
                print *, ""
                print *, "-------------------------------------------------"
                print *, "Interpolating particles from the grid"
                print *, achar(9)//'Max divergence of the vorticity field            ', maxval(abs(div_wmega))
                print *, achar(9)//'Total divergence of the vorticity field          ', sum(abs(div_wmega))
                print *, achar(9)//'Max vorticity                                    ', maxval(abs(RHS_pm(1:3,:,:,:)))
                print *, achar(9)//'Total vorticity                                  ', sum(abs(RHS_pm(1:3,:,:,:)))
                print *, 'Projecting particles to the grid'
            endif
            call project_particles_parallel

            if (my_rank.eq.0) then
                ! Print the difference between the old and the new Q
                print *, achar(9)//'Difference between the old and the new Q         ', sum(abs(QP - Qp_old))
                print *, achar(9)//'Difference between the old and the new vorticity ', sum(abs(RHS_pm - vorticity_old))

            endif   
        enddo
        if (my_rank .eq. 0) then
            deallocate (vorticity_old, Qp_old)
        end if
    end subroutine project_calc_div

    subroutine print_timestep_information()
        use pmgrid, only: RHS_pm, SOL_pm, DXpm, DYpm, DZpm, velocity_pm, deform_pm
        use file_io, only: solve_stats_file, case_folder
        implicit none
        real(dp), allocatable   :: div_wmega(:, :, :), laplace_LHS_pm(:, :, :, :), div_velocity(:, :, :)
        real(dp)                :: max_div_w, mean_div_w, mean_div_u, max_div_u, &
                                   total_momentum, total_kinetic_energy, total_vorticity, total_enstrophy, &
                                   CFL_x, CFL_y, CFL_z, CFL
        character(len=256)      :: filename
        logical                 :: file_exists
        integer                 :: ierr, my_rank, np, i

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        if (my_rank .eq. 0) then
            print *, ''
            ! We need to write the residuals of the solution
            div_wmega            = divergence(RHS_pm, DXpm, DYpm, DZpm)
            div_velocity         = divergence(velocity_pm, DXpm, DYpm, DZpm)
            laplace_LHS_pm       = laplacian(SOL_pm, DXpm, DYpm, DZpm)
            total_enstrophy      = sum(RHS_pm(1:3, :, :, :)**2)*DXpm*DYpm*DZpm
            total_vorticity      = sum(RHS_pm(1:3, :, :, :))*DXpm*DYpm*DZpm
            total_momentum       = sum(velocity_pm(1:3, :, :, :)) * DXpm * DYpm * DZpm
            total_kinetic_energy = sum(velocity_pm(1:3, :, :, :)**2 ) * DXpm * DYpm * DZpm

            CFL_x = maxval(abs(velocity_pm(1, :, :, :)))/DXpm
            CFL_y = maxval(abs(velocity_pm(2, :, :, :)))/DYpm
            CFL_z = maxval(abs(velocity_pm(3, :, :, :)))/DZpm
            CFL = CFL_x + CFL_y + CFL_z

            ! Print the mean, max and min div_u
            write (dummy_string, *) ""
            call vpm_print(dummy_string, nocolor, 1)
            write (dummy_string, "(A)") 'Divergence of the Ψ field'
            call vpm_print(dummy_string, yellow, 1)
            write (dummy_string, "(A, E10.4, A,E10.4, A, E10.4)") achar(9)//"div(ω)"//achar(9)// &
                achar(9)//" min : ", minval(div_wmega), &
                achar(9)//" max : ", maxval(div_wmega), &
                achar(9)//" mean: ", sum(div_wmega)/(fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
            call vpm_print(dummy_string, blue, 1)
            write (dummy_string, "(A)") 'Total Vorticity in the domain'
            call vpm_print(dummy_string, yellow, 1)
            write (dummy_string, "(A, E10.4)") achar(9)//'Total Vorticity : ', total_vorticity
            call vpm_print(dummy_string, blue, 1)
            write (dummy_string, "(A)") 'Total Enstrophy in the domain'
            call vpm_print(dummy_string, yellow, 1)
            write (dummy_string, "(A, E10.4)") achar(9)//'Total Enstrophy : ', total_enstrophy
            call vpm_print(dummy_string, blue, 1)
            write (dummy_string, "(A)") 'Residuals of the solution'
            call vpm_print(dummy_string, yellow, 1)
            do i = 1, neqpm
                ! For each equation write the laplacian - RHS_pm
                write (dummy_string, "(A, I3, A)") '   Equation =', i, ":   Δf = RHS"
                call vpm_print(dummy_string, blue, 1)
                write (dummy_string, "(A, E10.4, A, E10.4, A, E10.4)") achar(9)//'Forcing (RHS)'// &
                    achar(9)//'min : ', minval(RHS_pm(i, :, :, :)), &
                    achar(9)//'max : ', maxval(RHS_pm(i, :, :, :)), &
                    achar(9)//'mean: ', sum(RHS_pm(i, :, :, :))/(fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
                call vpm_print(dummy_string, nocolor, 1)
                write (dummy_string, "(A, E10.4, A, E10.4, A, E10.4)") achar(9)//"Solution"// &
                    achar(9)//'min : ', minval(SOL_pm(i, :, :, :)), &
                    achar(9)//'max : ', maxval(SOL_pm(i, :, :, :)), &
                    achar(9)//'mean: ', sum(SOL_pm(i, :, :, :))/(fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
                call vpm_print(dummy_string, nocolor, 1)
                write (dummy_string, "(A, E10.4, A, E10.4, A, E10.4)") achar(9)//'Res:=Δf-RHS'// &
                    achar(9)//'min : ', minval(abs(laplace_LHS_pm(i, :, :, :) - RHS_PM(i, :, :, :))), &
                    achar(9)//'max : ', maxval(abs(laplace_LHS_pm(i, :, :, :) - RHS_PM(i, :, :, :))), &
                    achar(9)//'mean: ', sum(abs(laplace_LHS_pm(i, :, :, :) - RHS_PM(i, :, :, :))) &
                    /(fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
                call vpm_print(dummy_string, nocolor, 1)
                ! Print the CFL number
                write (dummy_string, "(A, E10.4, A, E10.4, A, E10.4)") achar(9)//'CFL number / Dt'// &
                    achar(9)//'x : ', CFL_x, &
                    achar(9)//'y : ', CFL_y, &
                    achar(9)//'z : ', CFL_z
                call vpm_print(dummy_string, RED, 1)
                write (dummy_string, *) ""
                call vpm_print(dummy_string, nocolor, 1)
            end do

            ! Construct the filename
            write (filename, "(A, A)") trim(case_folder), trim(solve_stats_file)
            ! Determine the file status and write mode
            inquire(file=filename, exist=file_exists)

            if (NTIME_pm == 1 .or. NTIME_pm == 0 .or. .not. file_exists) then
                open (unit=10, file=filename, status='replace', action='write')
                ! Write CSV header
                write (10, '(A)') "Iteration, CFL/Dt, Total_Enstrophy, Total_Vorticity," &
                                //"Total_Momentum, Total_Kinetic_Energy," &
                                //"MEAN_DIV_W,MAX_DIV_W, MEAN_DIV_VEL, MAX_DIV_VEL," &
                                //"CFL_X/Dt, CFL_Y/Dt, CFL_Z/Dt" 
            else
                open (unit=10, file=filename, status='old', position='append', action='write')
            end if

            ! Write data for this iteration
            mean_div_u = sum(div_velocity) / (fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
            mean_div_w   = sum(div_wmega) / (fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
            max_div_u  = maxval(div_velocity)
            max_div_w    = maxval(div_wmega)

            write (10, '(I0,16(",",ES14.7))') NTIME_pm, CFL,                            &
                total_enstrophy, total_vorticity, total_momentum, total_kinetic_energy, &
                mean_div_w, max_div_w, mean_div_u, max_div_u,                           &
                CFL_x, CFL_y, CFL_z
            close (10)

            tab_level = tab_level - 1

            ! Deallocate the memory
            deallocate (div_wmega, laplace_LHS_pm)
        end if
    end subroutine print_timestep_information
end Module vpm_lib

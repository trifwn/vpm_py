Module vpm_lib
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
    use base_types, only: dp
    ! Printing
    use console_io, only: vpm_print, red, blue, green, nocolor, yellow, dummy_string, tab_level, VERBOCITY
    use parvar, only: print_particle_info, print_particle_positions, associate_particles
    ! Setting Vars
    use pmgrid, only: print_velocity_stats, print_vortex_stretching_stats, &
                      set_pm_velocities_zero, set_pm_deformations_zero, &
                      associate_velocities, associate_deformations

    !  WhatToDo flags
    integer, parameter :: DEFINE_PROBLEM = 0, &
                          SOLVE_VELOCITY = 1, &
                          SOLVE_VELOCITY_DELATATION = 2, &
                          PROJECT_AND_SOLVE = 3, &
                          INTERPOLATE = 4, &
                          DIFFUSE = 5
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
        RHS_pm_ptr, vel_ptr, NTIME_in, NI_in, NVR_size_in, &
        deform_ptr &
        )
        use parvar, only: NVR
        implicit None

        ! ARGUEMENTS
        real(dp), intent(inout), target  :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
        integer, intent(inout)           :: NVR_in
        integer, intent(in)              :: neqpm_in, WhatToDo, NVR_size_in, NTIME_in
        real(dp), intent(in)             :: NI_in
        real(dp), intent(out), pointer   :: RHS_pm_ptr(:, :, :, :)
        real(dp), intent(out), pointer   :: vel_ptr(:, :, :, :)
        real(dp), intent(out), pointer, optional   :: deform_ptr(:, :, :, :)

        ! LOCAL VARIABLES
        integer                           :: ierr, my_rank, np
        integer                           :: i, j, k

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        call set_num_threads()

        ! Print the input arguments
        if ((my_rank .eq. 0) .and. (VERBOCITY >= 1)) then
            if (WhatToDo .eq. 0) then
                write (dummy_string, "(A)") 'VPM: Initialize '
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
        if (present(deform_ptr)) call associate_deformations(deform_ptr)

        select case (WhatToDo)
        case (DEFINE_PROBLEM)
            call vpm_define_problem(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in)
        case (SOLVE_VELOCITY)
            call vpm_solve_velocity(NTIME_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, &
                                    neqpm_in, RHS_pm_ptr, vel_ptr)
        case (SOLVE_VELOCITY_DELATATION)
            call vpm_solve_velocity_deformation(NTIME_in, XP_in, QP_in, UP_in, GP_in, NVR_in, &
                                                NVR_size_in, neqpm_in, RHS_pm_ptr, vel_ptr, &
                                                deform_ptr)
        case (PROJECT_AND_SOLVE)
            call vpm_project_and_solve(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in, &
                                       RHS_pm_ptr)
        case (INTERPOLATE)
            call vpm_interpolate(XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, neqpm_in, &
                                 RHS_pm_ptr, deform_ptr)
        case (DIFFUSE)
            call vpm_diffuse(XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, neqpm_in, &
                             RHS_pm_ptr, deform_ptr, NI_in)
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
        tab_level = 1

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
        RHS_pm_ptr, deform_ptr &
        )
        use pmgrid, only: RHS_pm, velocity_pm, deform_pm
        use parvar, only: NVR
        use MPI
        implicit none

        real(dp), intent(inout), target           :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
        integer, intent(inout)                    :: NVR_in
        integer, intent(in)                       :: neqpm_in, NVR_size_in
        real(dp), intent(out), pointer            :: RHS_pm_ptr(:, :, :, :)
        real(dp), intent(out), pointer, optional   :: deform_ptr(:, :, :, :)
        integer :: ierr, my_rank, np

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        tab_level = 1

        ! Set Input
        ND = 3
        neqpm = neqpm_in
        call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in, UP_in, GP_in)
        if (present(deform_ptr)) call associate_deformations(deform_ptr)
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
        XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, neqpm_in, &
        RHS_pm_ptr, deform_ptr, &
        NI_in &
        )
        use MPI
        use pmgrid, only: RHS_pm
        use parvar, only: NVR
        implicit none
        real(dp), intent(inout), target           :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
        integer, intent(inout)                    :: NVR_in
        integer, intent(in)                       :: neqpm_in, NVR_size_in
        real(dp), intent(in)                      :: NI_in
        real(dp), intent(out), pointer            :: RHS_pm_ptr(:, :, :, :)
        real(dp), intent(out), pointer, optional  :: deform_ptr(:, :, :, :)
        integer :: ierr, my_rank, np

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        tab_level = 1

        ! Set Input
        ND = 3
        neqpm = neqpm_in
        NI = NI_in
        call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in, UP_in, GP_in)
        if (present(deform_ptr)) call associate_deformations(deform_ptr)
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
        integer, intent(inout)           :: NVR_in
        integer, intent(in)              :: neqpm_in, NVR_size_in
        real(dp), intent(out), pointer   :: RHS_pm_ptr(:, :, :, :)
        integer, intent(in)              :: NTIME_in

        integer :: ierr, my_rank, np
        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        ! Set Input
        tab_level = 1
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
            tab_level = tab_level + 1
            call vpm_print(dummy_string, blue, 0)
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
        integer, intent(inout)           :: NVR_in
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
        use pmgrid, only: velocity_pm, deform_pm
        use console_io, only: tab_level
        implicit none
        integer, intent(in)                        :: NTIME_in
        real(dp), intent(inout), target            :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
        integer, intent(inout)                     :: NVR_in
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
            !call convect_first_order(Xbound,Dpm,NN,NN_bl)
            ! FROM THE SOLUTION OF PM WE GET THE VELOCITIES ON THE GRID
            call calc_velocity ! VELOCITY AND DEFORMATION STO PM
            call calc_vortex_stretching_conservative(velocity_pm, deform_pm) ! VELOCITY AND DEFORMATION STO PM
        end if

        ! Solve the secondary problem:
        ! 1) Pressure Poisson Equation
        ! 2) Omega Correction Equation
        tab_level = 1
        call solve_secondary_problem(XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, RHS2_pm_ptr, SOL2_pm_ptr)

        if (my_rank .eq. 0) then
            call print_velocity_stats
            call print_vortex_stretching_stats
            tab_level = tab_level - 1
            call print_timestep_information
        end if

        call interpolate_particles_parallel(1) ! INTERPOLATION FROM PM TO PARTICLES
    end subroutine vpm_solve_velocity_deformation

    subroutine solve_secondary_problem( &
        XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, RHS_ptr, &
        SOL_ptr &
        )
        use parvar, only: QP, XP, UP, GP, NVR, NVR_size
        use pmgrid, only: RHS_pm, SOL_pm, DXpm, DYpm, DZpm, velocity_pm, deform_pm
        use serial_vector_field_operators, only: divergence, curl
        use file_io, only: write_pm_solution_hdf5, mesh_output_file, write_vector_field, write_scalar_field
        implicit none

        ! ARGUEMENTS
        real(dp), intent(inout), target  :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
        integer, intent(in)              :: NVR_in, NVR_size_in
        real(dp), intent(out), pointer   :: SOL_ptr(:, :, :, :), RHS_ptr(:, :, :, :)
        real(dp), allocatable, target    :: SOL2_pm(:, :, :, :), RHS2_pm(:, :, :, :)
        real(dp), allocatable, target    :: SOL2_pm_bl(:, :, :, :), RHS2_pm_bl(:, :, :, :)
        real(dp), allocatable            :: sol_biot(:, :, :)
        real(dp), allocatable            :: cross_product(:, :, :, :), div_wmega(:, :, :)
        real(dp), allocatable            :: grad_pi_fluxes(:, :, :, :), grad_pi_old(:, :, :, :), grad_pi_biot(:, :, :, :)
        real(dp), allocatable            :: new_omega_fluxes(:, :, :, :), new_omega_old(:, :, :, :), new_omega_biot(:, :, :, :)

        ! LOCAL VARIABLES
        integer                           :: bc_old
        integer                           :: ierr, my_rank, np
        integer                           :: i, j, k, iter
        character(len=256)                :: file_name
        integer                           :: iter_max = 10

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        iter_max = 1
        if (NTIME_pm .lt. 2) then
            iter_max = 20
        end if

        do iter = 1, iter_max
            ND = 3
            neqpm = 2
            call associate_particles(NVR_in, NVR_size_in, 3, XP_in, QP_in, UP_in, GP_in)
            if (NVR .eq. 0) return

            ! Allocate the solution and the RHS
            if (associated(SOL_ptr)) nullify (SOL_ptr)
            if (allocated(SOL2_pm)) deallocate (SOL2_pm)
            if (allocated(RHS2_pm)) deallocate (RHS2_pm)
            if (allocated(SOL2_pm_bl)) deallocate (SOL2_pm_bl)
            if (allocated(RHS2_pm_bl)) deallocate (RHS2_pm_bl)

            allocate (SOL2_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
            allocate (RHS2_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
            allocate (SOL2_pm_bl(neqpm, &
                                 block_grids(my_rank + 1)%NN(1), block_grids(my_rank + 1)%NN(2), block_grids(my_rank + 1)%NN(3)))
            allocate (RHS2_pm_bl(neqpm, &
                                 block_grids(my_rank + 1)%NN(1), block_grids(my_rank + 1)%NN(2), block_grids(my_rank + 1)%NN(3)))
            SOL2_pm = 0.d0
            RHS2_pm = 0.d0

            if (my_rank .eq. 0) then
                SOL_ptr => SOL2_pm
                RHS_ptr => RHS2_pm

                ! Get the divergence of the RHS of the equations:
                ! 1) 1st equation \nabla ^ 2 \pi = \nabla \cdot \omega that gives the pi field decomposition of
                !    the vorticity field. \omega = \nabla pi + \naabla x \qi => \nabla \cdot \omega = \nabla^2 pi
                ! We need to perform the green gauss theorem to calculate div(omega) in the center of the cells
                ! wmega_vel = curl(velocity_pm, DXpm, DYpm, DZpm)
                ! call compute_gauss_divergence(fine_grid,RHS_pm(1:3,:,:,:), RHS2_pm(1, :, :, :))
                RHS2_pm(1, :, :, :) = divergence(RHS_pm(1:3, :, :, :), DXpm, DYpm, DZpm)
                file_name = 'div_omega_gauss'
                call write_scalar_field(fine_grid, RHS2_pm(1, :, :, :), file_name)

                ! 2) 2nd equation is (\nabla ^2 q / rho ) = - \nabla (u x \omega) - \nabla (\sigma ) / rho
                !    where sigma is the stress tensor which we ignore for now
                allocate (cross_product(3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
                !$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(cross_product, velocity_pm, RHS_pm)
                !  Compute the cross product of the velocity and the vorticity
                do i = 1, fine_grid%NN(1)
                    do j = 1, fine_grid%NN(2)
                        do k = 1, fine_grid%NN(3)
                            cross_product(1, i, j, k) = velocity_pm(2, i, j, k)*RHS_pm(3, i, j, k) - &
                                                        velocity_pm(3, i, j, k)*RHS_pm(2, i, j, k)
                            cross_product(2, i, j, k) = velocity_pm(3, i, j, k)*RHS_pm(1, i, j, k) - &
                                                        velocity_pm(1, i, j, k)*RHS_pm(3, i, j, k)
                            cross_product(3, i, j, k) = velocity_pm(1, i, j, k)*RHS_pm(2, i, j, k) - &
                                                        velocity_pm(2, i, j, k)*RHS_pm(1, i, j, k)
                        end do
                    end do
                end do
                !$OMP END PARALLEL DO
                RHS2_pm(2, :, :, :) = divergence(cross_product, DXpm, DYpm, DZpm)
                deallocate (cross_product)
            end if

            if (my_rank .eq. 0) then
                write (*, *) ""
                tab_level = tab_level - 1
                write (dummy_string, "(A)") 'Solving Secondary Problem'
                call vpm_print(dummy_string, blue, 0)
                tab_level = tab_level + 1
                st = MPI_WTIME()
            end if

            ! Solve the problem
            bc_old = ibctyp
            ibctyp = 1
            call solve_problem(RHS2_pm, SOL2_pm, RHS2_pm_bl, SOL2_pm_bl)
            ibctyp = bc_old

            if (my_rank .eq. 0) then
                et = MPI_WTIME()
                write (dummy_string, "(A,I5,A,F8.2,A)") &
                    'Total time for solving the Secondary Problem', &
                    int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
                call vpm_print(dummy_string, blue, 1)

                file_name = mesh_output_file
                mesh_output_file = 'secondary_problem'
                call write_pm_solution_hdf5(NTIME_pm, fine_grid, neqpm, RHS2_pm, SOL2_pm)
                mesh_output_file = file_name

                ! Calculate the gradient of the solution to get the correction to the omega field

                ! Use the normal gradient
                grad_pi_old = gradient(SOL2_pm(1, :, :, :), DXpm, DYpm, DZpm)
                file_name = 'grad_pi_old'
                call write_vector_field(3, fine_grid, grad_pi_old, file_name)

                ! Reconstruct omega by taking gradients of pi
                allocate (grad_pi_fluxes(3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
                call compute_gradient_at_nodes(fine_grid, SOL2_pm(1, :, :, :), grad_pi_fluxes)
                file_name = 'grad_pi_fluxes'
                call write_vector_field(3, fine_grid, grad_pi_fluxes, file_name)

                ! SOLVE USING BIOT-SAVART
                ! allocate(sol_biot(fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
                ! div_wmega = divergence(RHS_pm(1:3, :, :, :), DXpm, DYpm, DZpm)
                ! file_name = 'div_wmega'
                ! call write_scalar_field(fine_grid, div_wmega, file_name)
                ! call solve_3D_poisson_Biot(fine_grid, div_wmega, sol_biot)
                ! deallocate(div_wmega)
                ! grad_pi_biot = gradient(sol_biot, DXpm, DYpm, DZpm)
                ! file_name = 'grad_pi_biot'
                ! call write_vector_field(3, fine_grid, grad_pi_biot, file_name)

                ! Print the difference between the two divergences
                print *, "Wmega Divergence Using Finite Differences"
                print *, achar(9)//'Max divergence of the vorticity field: ', maxval(abs(divergence(RHS_pm, DXpm, DYpm, DZpm)))
                print *, achar(9)//'Total divergence of the vorticity field: ', sum(abs(divergence(RHS_pm, DXpm, DYpm, DZpm)))
                print *, "Wmega Divergence Using Green-Gauss theorem"
                print *, achar(9)//'Max divergence of the vorticity field: ', maxval(abs(RHS2_pm(1, :, :, :)))
                print *, achar(9)//'Total divergence of the vorticity field: ', sum(abs(RHS2_pm(1, :, :, :)))
                print *, "Difference between the two divergences"
                print *, achar(9)//'Max difference between the two divergences: ', maxval(abs(RHS2_pm(1, :, :, :) - &
                                                                                divergence(RHS_pm(1:3, :, :, :), DXpm, DYpm, DZpm)))
                print *, achar(9)//'Total difference between the two divergences: ', sum(abs(RHS2_pm(1, :, :, :) - &
                                                                                divergence(RHS_pm(1:3, :, :, :), DXpm, DYpm, DZpm)))

                ! Print the difference between the two gradients
                print *, ""
                print *, "Difference between the Normal And Green-Gauss gradients of pi"
                print *, achar(9)//'Max difference between the two gradients: ', maxval(abs(grad_pi_fluxes - grad_pi_old))
                print *, achar(9)//'Total difference between the two gradients: ', sum(abs(grad_pi_fluxes - grad_pi_old))
                ! print *, "Difference between the Normal And Biot-Savart gradients of pi"
                ! print *, achar(9)//'Max difference between the two gradients: ', maxval(abs(grad_pi_biot - grad_pi_old))
                ! print *, achar(9)//'Total difference between the two gradients: ', sum(abs(grad_pi_biot - grad_pi_old))
                ! print *, "Difference between the Biot-Savart And Green-Gauss gradients of pi"
                ! print *, achar(9)//'Max difference between the two gradients: ', maxval(abs(grad_pi_biot - grad_pi_fluxes))
                ! print *, achar(9)//'Total difference between the two gradients: ', sum(abs(grad_pi_biot - grad_pi_fluxes))
                ! print *, ""

                print *, ""
                allocate (new_omega_fluxes(3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
                new_omega_fluxes = RHS_pm(1:3, :, :, :) - grad_pi_fluxes
                div_wmega = divergence(new_omega_fluxes, DXpm, DYpm, DZpm)
                print *, 'Correction to the vorticity field using the Green-Gauss theorem'
                print *, achar(9)//'Max divergence of the corrected vorticity field', maxval(abs(div_wmega))
                print *, achar(9)//'Total divergence of the corrected vorticity field', sum(abs(div_wmega))
                deallocate (div_wmega)

                allocate (new_omega_old(3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
                new_omega_old = RHS_pm(1:3, :, :, :) - grad_pi_old
                div_wmega = divergence(new_omega_old, DXpm, DYpm, DZpm)
                print *, 'Correction to the vorticity field using the normal gradient'
                print *, achar(9)//'Max divergence of the corrected vorticity field', maxval(abs(div_wmega))
                print *, achar(9)//'Total divergence of the corrected vorticity field', sum(abs(div_wmega))
                deallocate (div_wmega)

                ! allocate(new_omega_biot(3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
                ! new_omega_biot = RHS_pm(1:3, :, :, :) - grad_pi_biot
                ! div_wmega = divergence(new_omega_biot, DXpm, DYpm, DZpm)
                ! print *, 'Correction to the vorticity field using the Biot-Savart theorem'
                ! print *, achar(9)//'Max divergence of the corrected vorticity field', maxval(abs(div_wmega))
                ! print *, achar(9)//'Total divergence of the corrected vorticity field', sum(abs(div_wmega))
                ! print *, ""
                ! deallocate(div_wmega)

                RHS_pm(1:3, :, :, :) = new_omega_old

                deallocate (grad_pi_fluxes)
                deallocate (grad_pi_old)
                deallocate (new_omega_fluxes)
                deallocate (new_omega_old)
                ! deallocate(sol_biot)
                ! deallocate(grad_pi_biot)
                ! deallocate(new_omega_biot)
            end if
            tab_level = tab_level - 1
        end do
        neqpm = 3
    end subroutine solve_secondary_problem

    subroutine print_timestep_information()
        use pmgrid, only: RHS_pm, SOL_pm, DXpm, DYpm, DZpm
        use file_io, only: solve_stats_file, case_folder
        implicit none
        real(dp), allocatable   :: div_wmega(:, :, :), laplace_LHS_pm(:, :, :, :)
        real(dp)                :: total_vort
        real(dp)                :: max_div_w, mean_div_w, total_enstrophy, total_vorticity
        real(dp)                :: mean_residual(neqpm), max_residual(neqpm), mean_rhs(neqpm), max_rhs(neqpm)
        character(len=256)      :: filename
        integer                 :: ierr, my_rank, np, i

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        if (my_rank .eq. 0) then
            print *, ''
            ! We need to write the residuals of the solution
            div_wmega = divergence(RHS_pm, DXpm, DYpm, DZpm)
            laplace_LHS_pm = laplacian(SOL_pm, DXpm, DYpm, DZpm)
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
            total_vort = sum(RHS_pm(1:3, :, :, :))*DXpm*DYpm*DZpm
            write (dummy_string, "(A, E10.4)") achar(9)//'Total Vorticity : ', total_vort
            call vpm_print(dummy_string, blue, 1)
            write (dummy_string, "(A)") 'Total Enstrophy in the domain'
            call vpm_print(dummy_string, yellow, 1)
            total_vort = sum(RHS_pm(1:3, :, :, :)**2)*DXpm*DYpm*DZpm
            write (dummy_string, "(A, E10.4)") achar(9)//'Total Enstrophy : ', total_vort
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
                write (dummy_string, *) ""
                call vpm_print(dummy_string, nocolor, 1)
            end do

            ! Construct the filename
            write (filename, "(A, A)") trim(case_folder), trim(solve_stats_file)

            ! Determine the file status and write mode
            if (NTIME_pm == 1 .or. NTIME_pm == 0) then
                open (unit=10, file=filename, status='replace', action='write')
                ! Write CSV header
                write (10, '(A)') "Iteration,Total_Enstrophy, Total_Vorticity," &
                    //"MEAN_DIV_W,MAX_DIV_W," &
                    //"MEAN_RESIDUAL1,MAX_RESIDUAL1,MEAN_RHS1,MAX_RHS1," &
                    //"MEAN_RESIDUAL2,MAX_RESIDUAL2,MEAN_RHS2,MAX_RHS2," &
                    //"MEAN_RESIDUAL3,MAX_RESIDUAL3,MEAN_RHS3,MAX_RHS3"
            else
                open (unit=10, file=filename, status='old', position='append', action='write')
            end if

            ! Calculate required values
            total_enstrophy = sum(RHS_pm(1:3, :, :, :)**2)*DXpm*DYpm*DZpm
            total_vorticity = sum(RHS_pm(1:3, :, :, :))

            mean_div_w = sum(div_wmega)/(fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
            max_div_w = maxval(abs(div_wmega))

            ! Calculate residual statistics across all equations
            mean_residual = sum(sum(sum(abs(laplace_LHS_pm - RHS_PM), dim=4), dim=3), dim=2)/ &
                            (neqpm*fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
            max_residual = maxval(maxval(maxval(abs(laplace_LHS_pm - RHS_PM), dim=4), dim=3), dim=2)

            ! Calculate RHS statistics across all equations
            mean_rhs = sum(sum(sum(RHS_PM, dim=4), dim=3), dim=2)/ &
                       (neqpm*fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
            max_rhs = maxval(maxval(maxval(abs(RHS_PM), dim=4), dim=3), dim=2)

            ! Write data for this iteration
            write (10, '(I0,16(",",ES14.7))') NTIME_pm, &
                total_enstrophy, total_vorticity, &
                mean_div_w, max_div_w, &
                mean_residual(1), max_residual(1), mean_rhs(1), max_rhs(1), &
                mean_residual(2), max_residual(2), mean_rhs(2), max_rhs(2), &
                mean_residual(3), max_residual(3), mean_rhs(3), max_rhs(3)

            close (10)

            tab_level = tab_level - 1

            ! Deallocate the memory
            deallocate (div_wmega, laplace_LHS_pm)
        end if
    end subroutine print_timestep_information
end Module vpm_lib

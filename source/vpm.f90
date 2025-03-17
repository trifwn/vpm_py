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
    use console_io, only: vpm_print, red, blue, green, nocolor, yellow, dummy_string, tab_level, &
                        VERBOCITY, print_stats_rank3, print_stats_rank4 
    use parvar, only: print_particle_info, print_particle_positions, associate_particles
    ! Setting Vars
    use pmgrid, only: print_velocity_stats, print_vortex_stretching_stats, &
                      set_pm_velocities_zero, set_pm_deformations_zero, &
                      associate_velocities, associate_deformations
    use serial_vector_field_operators, only: divergence, curl, laplacian, gradient

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

    subroutine read_conf()
        use pmgrid, only: DXpm, DYpm, DZpm, IDVPM, ncoarse
        integer     :: i, ncell_rem
        logical     :: pmfile_exists
        integer     :: my_rank, ierr

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

        inquire (file='pm.input', exist=pmfile_exists)
        if (pmfile_exists) then
            open (1, file='pm.input')
            read (1, *) DXpm, DYpm, DZpm     ! CELL SIZES
            read (1, *) interf_iproj         ! INTERFACE PROJECTION FUNCTION 1 , 2, 3 , 4
            read (1, *) ibctyp               ! 1 - PERIODIC, 2 - INFLOW, 3 - OUTFLOW, 4 - WALL
            read (1, *) IDVPM                ! Variable/Constant Volume(0,1)
            read (1, *)
            read (1, *) ncoarse              ! NUMBER OF FINE CELLS PER COARSE CELL per dir
            read (1, *) NBI, NBJ, NBK        !  NBI x NBJ x NBK = NUM OF PROCESSORS (NP)
            read (1, *) ncell_rem   ! 0: NO REMESHING, 1: REMESHING, ncell_rem: PARTICLE PER CELL
            read (1, *) iyntree, ilevmax     ! 1: TREE 0: NO TREE, 3: NUMB OF SUBDIVISION (2^3)
            read (1, *) OMPTHREADS           ! 1 - OPENMP THREADS
            read (1, *) idefine              ! 0: FREE GRID, 1: FIXED GRID
            close(1)

            if (my_rank .eq. 0) then
                print *, 'Inputs read:'
                print *, achar(9), 'DXpm=', DXpm
                print *, achar(9), 'DYpm=', DYpm
                print *, achar(9), 'DZpm=', DZpm
                print *, achar(9), 'Processors=', NBI, NBJ, NBK
            end if
        endif
    end subroutine read_conf

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

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        call set_num_threads()

        call read_conf()
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
            write (*, *) achar(9), achar(9), 'NVR_size_in = ', NVR_size_in, achar(27)//'[0m'
        end if
        ! We set these input for compatibility with old version
        ! In the new modular approach we do not need to set these
        ND = 3
        NTIME_pm = NTIME_in
        neqpm = neqpm_in
        call associate_particles(NVR_in, NVR_size_in, XP_in, QP_in, UP_in, GP_in)
        if (NVR .eq. 0) return
        if (associated(RHS_pm_ptr)) nullify (RHS_pm_ptr)

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
            call vpm_diffuse(NI_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, RHS_pm_ptr)

        case (CORRECT_VORTICITY)
            call vpm_correct_vorticity(XP_in, QP_in, NVR_in, NVR_size_in)
        
        end select
        call associate_velocities(vel_ptr)
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
        call associate_particles(NVR_in, NVR_size_in, XP_in, QP_in)
        if (NVR .eq. 0) return

        call define_sizes
        if (my_rank .eq. 0) then
            call set_pm_velocities_zero
            call set_pm_deformations_zero
        endif 
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
        call associate_particles(NVR_in, NVR_size_in, XP_in, QP_in, UP_in, GP_in)
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
        NI_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, &
        RHS_pm_ptr &
    )
        use MPI
        use pmgrid, only: RHS_pm
        use parvar, only: NVR
        implicit none
        real(dp), intent(in)                      :: NI_in
        real(dp), intent(inout), target           :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
        integer, intent(in)                       :: NVR_in, NVR_size_in
        real(dp), intent(out), pointer          :: RHS_pm_ptr(:, :, :, :)
        integer :: ierr, my_rank, np

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        tab_level = 0

        ! Set Input
        ND = 3
        call associate_particles(NVR_in, NVR_size_in, XP_in, QP_in, UP_in, GP_in)
        call allocate_sol_and_rhs(my_rank + 1)
        RHS_pm_ptr => RHS_pm

        if (NVR .eq. 0) then
            if (size(XP_in, 2) .eq. 0) then
                print *, 'No particles to interpolate'
                return
            end if
            NVR = size(XP_in, 2)
        end if

        ! PARTICLES TO -> PM MEANING FROM Qs We get RHS_pm
        call project_particles_parallel
        call set_pm_deformations_zero
        if (my_rank .eq. 0) then
            ! diffusion stores -NI*grad^2 w * Vol in GP
            ! RHS_pm = -VIS \nabla \cdot RHS_pm
            write (dummy_string, "(A, F12.6)") 'Diffusing Vorticity with ν=', NI_in
            call vpm_print(dummy_string, blue, 0)
            call diffuse_vort_3d(NI_in) ! DIFFUSION OF VORTICITY
        end if
        ! WHEN ITYPEB = 2 WE GET THE GP FROM THE (DEFORMATION) 
        call interpolate_particles_parallel(2)
    end subroutine vpm_diffuse

    subroutine vpm_project_and_solve(                   &
        NTIME_in,                                       &
        XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in,    &
        RHS_pm_ptr                                      &
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
        call associate_particles(NVR_in, NVR_size_in, XP_in, QP_in)
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
            write (dummy_string, "(A)") 'Solving Poisson problem for the vorticity field'
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
        
        call associate_particles(NVR_in, NVR_size_in, XP_in, QP_in, UP_in, GP_in)
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
        use console_io, only: tab_level, print_timestep_information
        use file_io, only: write_timestep_information_dat
        implicit none
        integer, intent(in)                        :: NTIME_in
        real(dp), intent(inout), target            :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
        integer, intent(in)                        :: NVR_in
        integer, intent(in)                        :: neqpm_in, NVR_size_in
        real(dp), intent(out), pointer             :: RHS_pm_ptr(:, :, :, :)
        real(dp), intent(out), pointer             :: vel_ptr(:, :, :, :)
        real(dp), intent(out), pointer, optional   :: deform_ptr(:, :, :, :)

        ! LOCAL VARIABLES
        integer                 :: ierr, my_rank

        real(dp)                :: DXpm, DYpm, DZpm
        real(dp), allocatable   :: vorticity(:,:,:,:), vorticity2(:,:,:,:), error_vorticity(:,:,:,:)

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

        call associate_particles(NVR_in, NVR_size_in, XP_in, QP_in, UP_in, GP_in)
        ! Project the particles to the PM grid and solve the problem
        call vpm_project_and_solve(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in, RHS_pm_ptr)

        call associate_velocities(vel_ptr)
        if (present(deform_ptr)) call associate_deformations(deform_ptr)
        
        ! Calculate the velocity and deformation on the PM grid
        if (my_rank .eq. 0) then
            print *, ''
            ! FROM THE SOLUTION OF PM WE GET THE VELOCITIES ON THE GRID
            call calc_velocity ! VELOCITY AND DEFORMATION STO PM
            call calc_vortex_stretching_conservative(velocity_pm, deform_pm) ! VELOCITY AND DEFORMATION STO PM
        end if

        call interpolate_particles_parallel(1) ! INTERPOLATION FROM PM TO PARTICLES

        if (my_rank .eq. 0) then
            call print_velocity_stats
            call print_vortex_stretching_stats
            tab_level = tab_level - 1
            call get_timestep_information(timestep_info)
            call get_solve_info(solve_info)
            call print_timestep_information(timestep_info, solve_info)

                DXpm = fine_grid%Dpm(1)
                DYpm = fine_grid%Dpm(2)
                DZpm = fine_grid%Dpm(3)
                allocate (vorticity, mold = RHS_pm)
                vorticity(:,:,:,:) = - RHS_pm(:,:,:,:)
                ! Calculate the difference between the vorticity and the velocity
                vorticity2 = curl(velocity_pm, DXpm, DYpm, DZpm) 
                print *, 'The difference between the vorticity and the velocity curl'
                allocate (error_vorticity, mold = vorticity)

                error_vorticity(:,:,:,:) = vorticity2(:,:,:,:) - vorticity(:,:,:,:)
                call print_stats_rank4(fine_grid, vorticity2, 'curl(u)')
                call print_stats_rank4(fine_grid, vorticity, 'ω')
                call print_stats_rank4(fine_grid, error_vorticity, 'curl(u) - ω')

            if (VERBOCITY .ge. 1) then
                call print_timestep_information(timestep_info, solve_info)
                call write_timestep_information_dat(timestep_info, NTIME_in)
            endif
        end if

    end subroutine vpm_solve_velocity_deformation

    subroutine vpm_correct_vorticity(           &
        XP_in, QP_in, NVR_in, NVR_size_in       &
    )
        use parvar, only: QP, XP, NVR, NVR_size
        use serial_vector_field_operators, only: divergence, curl
        use vpm_interpolate, only: interpolate_particle_Q
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
        ! OlD parameters
        real(dp), allocatable            :: XP_old(:,:), QP_old(:,:)
        real(dp), pointer                :: vorticity(:,:,:,:) 
        
        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        
        call associate_particles(NVR_in, NVR_size_in, XP_in, QP_in)
        
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
        
        ND = 3
        neqpm = 1

        nb = my_rank + 1
        allocate (SOL2_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
        allocate (RHS2_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
        SOL2_pm = 0.d0
        RHS2_pm = 0.d0
        allocate (SOL2_pm_bl(neqpm, block_grids(nb)%NN(1), block_grids(nb)%NN(2), block_grids(nb)%NN(3)))
        allocate (RHS2_pm_bl(neqpm, block_grids(nb)%NN(1), block_grids(nb)%NN(2), block_grids(nb)%NN(3)))
        
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
            ! call compute_gauss_divergence(fine_grid,vorticity(1:3,:,:,:), RHS2_pm(1, :, :, :))
            
            ! BIOT-SAVART VERSION
            ! call 
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
        if (allocated(SOL2_pm)) deallocate (SOL2_pm)
        if (allocated(RHS2_pm)) deallocate (RHS2_pm)
        if (allocated(SOL2_pm_bl)) deallocate (SOL2_pm_bl)
        if (allocated(RHS2_pm_bl)) deallocate (RHS2_pm_bl)
        ! call project_calc_div
    end subroutine vpm_correct_vorticity

    !> Solve the secondary problem:
    !> 1) Pressure Poisson Equation
    !>    equation is (\nabla ^2 q / rho ) = - \nabla (u x \omega) - \nabla (\sigma ) / rho
    !>    where sigma is the stress tensor which we ignore for now
    subroutine vpm_solve_pressure(                                            &
       vorticity, velocity, pressure, density, viscocity, dphi_dt, p_reference           &
    )
        use vpm_interpolate, only: interpolate_particle_Q
        real(dp), pointer, intent(in)           :: vorticity(:,:,:,:)
        real(dp), pointer, intent(in)           :: velocity(:,:,:,:)
        real(dp), allocatable, intent(out)      :: pressure(:,:,:,:)
        real(dp), intent(in)                    :: density, viscocity
        real(dp), pointer, intent(in), optional :: dphi_dt(:,:,:,:)
        real(dp), value, optional               :: p_reference
        ! LOCAL VARIABLES
        real(dp), allocatable             :: SOL2_pm(:,:,:,:), RHS2_pm(:,:,:,:)
        real(dp), allocatable             :: SOL2_pm_bl(:,:,:,:), RHS2_pm_bl(:,:,:,:)

        ! Variation 1
        real(dp), allocatable             :: cross_product(:,:,:,:), div_cross(:,:,:)
        real(dp), allocatable             :: div_stress_tensor(:,:,:,:), div_div_stress_tensor(:,:,:)

        ! Variation 2
        real(dp), allocatable             :: curl_vorticity(:,:,:,:)
        real(dp), allocatable             :: curl_w_x_u(:,:,:)
        real(dp), allocatable             :: omega_squared(:,:,:)
        real(dp), allocatable             :: vorticity2(:,:,:,:), error_vorticity(:,:,:,:)

        ! real(dp), allocatable             :: lapl
        real(dp)                          :: DXpm, DYpm, DZpm
        real(dp)                          :: velocity_squared
        integer                           :: ierr, my_rank, np, nb
        integer                           :: i, j, k
        

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        tab_level = 0
        if (my_rank .eq. 0) then
            write (*, *) ""
            tab_level = tab_level - 1
            write (dummy_string, "(A, F8.2, A, F8.2)") 'Solving Poisson for the Pressure field. ρ = ', density, ' ν = ', viscocity
            call vpm_print(dummy_string, blue, 0)
            tab_level = tab_level + 1
            st = MPI_WTIME()
        end if

        ! Allocate the solution and the RHS
        ND = 3
        neqpm = 2

        DXpm = fine_grid%Dpm(1)
        DYpm = fine_grid%Dpm(2)
        DZpm = fine_grid%Dpm(3)
        
        allocate (SOL2_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
        allocate (RHS2_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
        nb = my_rank + 1
        allocate (SOL2_pm_bl(neqpm, block_grids(nb)%NN(1), block_grids(nb)%NN(2), block_grids(nb)%NN(3)))
        allocate (RHS2_pm_bl(neqpm, block_grids(nb)%NN(1), block_grids(nb)%NN(2), block_grids(nb)%NN(3)))
                
        SOL2_pm = 0.d0
        RHS2_pm = 0.d0

        if (my_rank .eq. 0) then
            vorticity(:,:,:,:) = - vorticity(:,:,:,:)
            ! Calculate the difference between the vorticity and the velocity
            vorticity2 = curl(velocity, DXpm, DYpm, DZpm) 

            print *, 'The difference between the vorticity and the velocity curl'
            allocate (error_vorticity, mold = vorticity)
            error_vorticity(:,:,:,:) = vorticity2(:,:,:,:) - vorticity(:,:,:,:)
            call print_stats_rank4(fine_grid, vorticity2, 'curl(u)')
            call print_stats_rank4(fine_grid, vorticity, 'ω')
            call print_stats_rank4(fine_grid, error_vorticity, 'curl(u) - ω')

            ! Variation 1
            allocate (cross_product, mold = vorticity)
            call compute_cross_product(velocity, vorticity, cross_product)
            div_cross = divergence(cross_product, DXpm, DYpm, DZpm) 

            ! Print the min and max values of the divergence of the stress tensor
            call print_stats_rank3(fine_grid, div_cross, 'div(u x ω)')

            if (viscocity .lt. 1.d-14) then
                print *, 'The viscocity is 0, the stress tensor is ignored. Calculating Eulerian pressure'
                RHS2_pm(1, :, :, :) =  (div_cross(:, :, :)) * density
            else
                print *, 'The viscocity is not 0, the stress tensor is considered. Calculating NS pressure'
                !  Compute the laplacian of the deviatoric stress tensor
                div_stress_tensor = laplacian(velocity, DXpm, DYpm, DZpm)
                div_div_stress_tensor = divergence(div_stress_tensor, DXpm, DYpm, DZpm) * viscocity
                call print_stats_rank3(fine_grid, div_div_stress_tensor, '\nabla \cdot \nabla \cdot \sigma')
                
                RHS2_pm(1, :, :, :) = (div_cross(:, :, :) - div_div_stress_tensor(:, :, :)) * density
            endif
            call print_stats_rank3(fine_grid, RHS2_pm(1, :, :, :), 'Pressure Poisson RHS = div(u x ω) + div(div(σ))') 
            
            ! Variation 2
            curl_vorticity = - curl(vorticity, DXpm, DYpm, DZpm)
            allocate (curl_w_x_u(fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
            allocate (omega_squared(fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
            ! Compute the dot product of the curl of the vorticity and the velocity
            do i = 1, fine_grid%NN(1)
                do j = 1, fine_grid%NN(2)
                    do k = 1, fine_grid%NN(3)
                        omega_squared(i, j, k) = sum(vorticity(:, i, j, k)**2)
                        curl_w_x_u(i, j, k) = dot_product(curl_vorticity(:, i, j, k), velocity(:, i, j, k))
                    end do
                end do
            end do
            RHS2_pm(2, :, :, :) = (- omega_squared(:,:,:) + curl_w_x_u(:,:,:)) * density

            call print_stats_rank3(fine_grid, - omega_squared(:,:,:), '- ω^2')
            call print_stats_rank3(fine_grid, curl_w_x_u(:,:,:), 'curl(w) * u')

            call print_stats_rank3(fine_grid, RHS2_pm(2, :, :, :), 'Pressure Poisson RHS = - ω^2 - curl(w) * u')

            print *
            print *
            print *
            print *

            deallocate (cross_product)
            deallocate (div_cross)
            ! deallocate (div_stress_tensor)
            ! deallocate (div_div_stress_tensor)
            deallocate (curl_vorticity)
            deallocate (curl_w_x_u)
            deallocate (omega_squared)
        endif 

        ! Solve the problem
        call solve_problem(RHS2_pm, SOL2_pm, RHS2_pm_bl, SOL2_pm_bl)

        if (my_rank .eq. 0) then
            ! if (associated(pressure)) deallocate (pressure)
            allocate (pressure(3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
            pressure = 0.d0
            ! THIS IS THE Q PART OF THE PRESSURE
            pressure(1,:,:,:) = SOL2_pm(1,:,:,:)
            call print_stats_rank3(fine_grid, pressure(1, :, :, :), 'Pressure Q - 1')
            pressure(1,:,:,:) = SOL2_pm(2,:,:,:)
            call print_stats_rank3(fine_grid, pressure(1, :, :, :), 'Pressure Q - 2')
            
            ! Calculate the pressure field based on the bernoulli equation
            IF(.NOT. PRESENT(p_reference)) p_reference = 0 
            !101325.d0 = 1 atm 
            do i = 1, fine_grid%NN(1)
                do j = 1, fine_grid%NN(2)
                    do k = 1, fine_grid%NN(3)
                        velocity_squared = sum(velocity(:,i,j,k)**2)
                        pressure(2,i,j,k) = (density * velocity_squared / 2.d0)
                        pressure(3,i,j,k) = p_reference - pressure(1,i,j,k) - pressure(2,i,j,k) 
                    end do
                end do
            end do

            if (present(dphi_dt)) then
                ! Calculate the pressure gradient
                pressure(3, :, :, :) = pressure(3, :, :, :) - dphi_dt(1, :, :, :) * density
            end if

            call print_stats_rank3(fine_grid, pressure(2, :, :, :), 'Pressure U')
            call print_stats_rank3(fine_grid, pressure(3, :, :, :), 'Pressure Bernoulli Solution')

            ! Get the residuals
            call calculate_momentum_residuals(velocity, vorticity, pressure(3, :, :, :), density, viscocity)
        end if 

        if (my_rank .eq. 0) then
            et = MPI_WTIME()
            write (dummy_string, "(A,I5,A,F8.2,A)") &
                'Total time for solving the Pressure Problem', &
                int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
            call vpm_print(dummy_string, blue, 1)
        end if
        
        neqpm = 3
        ! DEALLOCATE THE MEMORY
        if (allocated(SOL2_pm)) deallocate (SOL2_pm)
        if (allocated(RHS2_pm)) deallocate (RHS2_pm)
        if (allocated(SOL2_pm_bl)) deallocate (SOL2_pm_bl)
        if (allocated(RHS2_pm_bl)) deallocate (RHS2_pm_bl)
    end subroutine vpm_solve_pressure

    subroutine calculate_momentum_residuals(velocity, vorticity, pressure, density, viscocity)
        implicit none
        real(dp), target, intent(in)           :: velocity(:,:,:,:)
        real(dp), target, intent(in)           :: vorticity(:,:,:,:)
        real(dp), target, intent(in)           :: pressure(:,:,:)
        real(dp), intent(in)                   :: density, viscocity

        ! LOCAL VARIABLES
        real(dp)                :: DXpm, DYpm, DZpm
        real(dp), allocatable   :: residuals(:,:,:,:)
        real(dp), allocatable   :: cross_product(:,:,:,:)
        real(dp), allocatable   :: grad_p(:,:,:,:)
        real(dp), allocatable   :: div_stress_tensor(:,:,:,:)
        real(dp), allocatable   :: velocity_squared(:,:,:)
        real(dp), allocatable   :: grad_velocity_sq(:,:,:, :)
        integer                 :: i, j, k

        DXpm = fine_grid%Dpm(1)
        DYpm = fine_grid%Dpm(2)
        DZpm = fine_grid%Dpm(3)

        ! The momentum equation is given by:
        !  \frac {\partial u }{\partial t} + \nabla \frac{u^2}{2}- u \times \omega 
        !                   = \frac{1}{\rho}\nabla \cdot ( -pI + \sigma )
        ! For the Steady state, the equation becomes:
        !  \nabla \frac{u^2}{2} - u \times \omega = \frac{1}{\rho}\nabla \cdot ( -pI + \sigma )
        !  The Residuals are calculated as follows:
        ! R = \nabla \frac{u^2}{2} - u \times \omega - \frac{1}{\rho}\nabla \cdot ( -pI + \sigma )
        ! R = \nabla \frac{u^2}{2} - u \times \omega - \frac{1}{\rho} (- grad(P) + nabla^2 u)

        ! 0) Allocate the residuals
        allocate (residuals(3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))

        ! 1) u \times \omega
        allocate (cross_product, mold = vorticity)
        call compute_cross_product(velocity, vorticity, cross_product)
        cross_product = cross_product * (-1.d0)
        residuals(:, :, :, :) = cross_product(:,:,:,:)
        
        ! Print the min and max values of the cross product
        call print_stats_rank4(fine_grid, cross_product, 'Cross Product')
        deallocate (cross_product)

        ! 2) \grad p
        grad_p = gradient(pressure(:, :, :), DXpm, DYpm, DZpm) * (-1.d0/density)
        residuals(:, :, :, :) = residuals(:, :, :, :) - grad_p(:,:,:,:)

        ! Print the min and max values of the pressure gradient
        call print_stats_rank4(fine_grid, grad_p, 'Pressure Gradient')
        deallocate (grad_p)

        ! 3) \nabla \cdot \sigma = \nabla^2 u
        if (viscocity .gt. 1.d-14) then
            div_stress_tensor = laplacian(velocity, DXpm, DYpm, DZpm) * viscocity
            residuals(:, :, :, :) = residuals(:, :, :, :) - div_stress_tensor(:,:,:,:)

            ! Print the min and max values of the divergence of the stress tensor
            call print_stats_rank4(fine_grid, div_stress_tensor, 'Divergence of the Stress Tensor')
        endif

        ! 4) \nabla (u^2/2)
        allocate (velocity_squared(fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
        do i = 1, fine_grid%NN(1)
            do j = 1, fine_grid%NN(2)
                do k = 1, fine_grid%NN(3)
                    velocity_squared(i, j, k) = sum(velocity(:, i, j, k)**2)
                end do
            end do
        end do
        grad_velocity_sq = gradient(velocity_squared, DXpm, DYpm, DZpm)
        deallocate (velocity_squared)
        residuals(:, :, :, :) = residuals(:, :, :, :) - grad_velocity_sq(:,:,:,:)

        ! Print the min and max values of the gradient of the velocity squared
        call print_stats_rank4(fine_grid, grad_velocity_sq, 'Gradient of the Velocity Squared')
        deallocate (grad_velocity_sq)

        ! Print the min and max values of the residuals
        call print_stats_rank4(fine_grid, residuals, 'Residuals')
        deallocate (residuals)

    end subroutine calculate_momentum_residuals

    subroutine project_calc_div
        use MPI
        use parvar, only: QP, XP, NVR, NVR_size
        use serial_vector_field_operators, only: divergence, curl
        use vpm_interpolate, only: interpolate_particle_Q
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

        do i = 1,100
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

    subroutine get_timestep_information(timestep_information)
        use pmgrid, only: RHS_pm, DXpm, DYpm, DZpm, velocity_pm, DVpm
        use parvar, only: XP, QP, UP, NVR, NVR_size
        use vpm_size, only: fine_grid 
        use vpm_types, only: timestepInformation, solveInformation, dp
        use serial_vector_field_operators, only: divergence, laplacian
        implicit none
        type(timestepInformation), intent(inout) :: timestep_information
        real(dp), allocatable :: div_wmega(:,:,:), div_velocity(:,:,:), particle_q
        integer :: i
        
        timestep_information%n          = NTIME_pm

        timestep_information%NVR        = NVR
        timestep_information%NVR_size   = NVR_size

        timestep_information%NN        = fine_grid%NN
        timestep_information%Xbound    = fine_grid%Xbound
        timestep_information%Dpm       = fine_grid%Dpm
        
        timestep_information%solver    = SOLVER

        ! Calculate divergences and laplacian
        div_wmega                = divergence(RHS_pm(1:3, :, :, :), DXpm, DYpm, DZpm)
        timestep_information%mean_div_w = sum(div_wmega)/(fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
        timestep_information%max_div_w  = maxval(div_wmega)
        timestep_information%min_div_w  = minval(div_wmega)
        deallocate (div_wmega)

        div_velocity             = divergence(velocity_pm, DXpm, DYpm, DZpm) 
        timestep_information%mean_div_u = sum(div_velocity)/(fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
        timestep_information%max_div_u  = maxval(div_velocity)
        timestep_information%min_div_u  = minval(div_velocity)
        deallocate (div_velocity)

        timestep_information%total_enstrophy_pm       = sum(RHS_pm(1:3,:,:,:)**2)*DXpm*DYpm*DZpm
        timestep_information%total_vorticity_pm       = sum(RHS_pm(1:3,:,:,:))*DXpm*DYpm*DZpm
        timestep_information%total_momentum_x_pm      = sum(velocity_pm(1,:,:,:))*DXpm*DYpm*DZpm
        timestep_information%total_momentum_y_pm      = sum(velocity_pm(2,:,:,:))*DXpm*DYpm*DZpm 
        timestep_information%total_momentum_z_pm      = sum(velocity_pm(3,:,:,:))*DXpm*DYpm*DZpm
        timestep_information%total_kinetic_energy_pm  = sum(velocity_pm(1:3,:,:,:)**2)*DXpm*DYpm*DZpm

        timestep_information%total_vorticity_particles = sum(QP(1, :) * QP(4, :) + QP(2, :) * QP(4, :) + QP(3, :) * QP(4, :)) 
        timestep_information%total_momentum_x_particles = sum(UP(1,:) * QP(4, :))
        timestep_information%total_momentum_y_particles = sum(UP(2,:) * QP(4, :))
        timestep_information%total_momentum_z_particles = sum(UP(3,:) * QP(4, :))

        particle_q = 0
        do i = 1, NVR
            particle_q = particle_q + sum(UP(1:3,i)**2) * QP(4, i)
        end do
        timestep_information%total_kinetic_energy_particles = particle_q

        particle_q = 0
        do i = 1, NVR
            particle_q = particle_q + sum(QP(1:3,i)**2) * QP(4, i)
        end do
        timestep_information%total_enstrophy_particles = particle_q
        

        print *, 'Total enstrophy particles', timestep_information%total_enstrophy_particles
        print *, 'Total vorticity particles', timestep_information%total_vorticity_particles
        print *, 'Total momentum x particles', timestep_information%total_momentum_x_particles
        print *, 'Total momentum y particles', timestep_information%total_momentum_y_particles
        print *, 'Total momentum z particles', timestep_information%total_momentum_z_particles
    end subroutine get_timestep_information
    
    subroutine get_solve_info(solve_information)
        use pmgrid, only: RHS_pm, SOL_pm, DXpm, DYpm, DZpm
        use vpm_size, only: fine_grid
        use vpm_types, only: solveInformation, dp
        use serial_vector_field_operators, only: laplacian
        use vpm_vars, only: neqpm
        implicit none
        type(solveInformation), intent(inout) :: solve_information
        real(dp), allocatable                 :: laplace_LHS_pm(:,:,:,:)
        real(dp)                              :: sumV
        integer                               :: i

        laplace_LHS_pm = laplacian(SOL_pm, DXpm, DYpm, DZpm)
        sumV = (fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))

        if (allocated(solve_information%f_min)) deallocate (solve_information%f_min)
        if (allocated(solve_information%f_max)) deallocate (solve_information%f_max)
        if (allocated(solve_information%f_mean)) deallocate (solve_information%f_mean)
        
        if (allocated(solve_information%sol_min)) deallocate (solve_information%sol_min)
        if (allocated(solve_information%sol_max)) deallocate (solve_information%sol_max)
        if (allocated(solve_information%sol_mean)) deallocate (solve_information%sol_mean)
        
        if (allocated(solve_information%residual_min)) deallocate (solve_information%residual_min)
        if (allocated(solve_information%residual_max)) deallocate (solve_information%residual_max)
        if (allocated(solve_information%residual_mean)) deallocate (solve_information%residual_mean)

        allocate (solve_information%f_min(neqpm))
        allocate (solve_information%f_max(neqpm))
        allocate (solve_information%f_mean(neqpm))

        allocate (solve_information%sol_min(neqpm))
        allocate (solve_information%sol_max(neqpm))
        allocate (solve_information%sol_mean(neqpm))

        allocate (solve_information%residual_min(neqpm))
        allocate (solve_information%residual_max(neqpm))
        allocate (solve_information%residual_mean(neqpm))
        
        do i = 1,neqpm
            solve_information%f_min(i)    = minval(RHS_pm(i, :, :, :))
            solve_information%f_max(i)    = maxval(RHS_pm(i, :, :, :))
            solve_information%f_mean(i)   = sum(RHS_pm(i, :, :, :))/(sumV)

            solve_information%sol_min(i)  = minval(SOL_pm(i, :, :, :))
            solve_information%sol_max(i)  = maxval(SOL_pm(i, :, :, :))
            solve_information%sol_mean(i) = sum(SOL_pm(i, :, :, :))/(sumV)

            solve_information%residual_min(i)  = minval(abs(laplace_LHS_pm(i, :, :, :) - RHS_pm(i, :, :, :)))
            solve_information%residual_max(i)  = maxval(abs(laplace_LHS_pm(i, :, :, :) - RHS_pm(i, :, :, :)))
            solve_information%residual_mean(i) = sum(abs(laplace_LHS_pm(i, :, :, :) - RHS_pm(i, :, :, :)))/(sumV)
        enddo
        deallocate (laplace_LHS_pm)
    end subroutine get_solve_info
end module vpm_lib

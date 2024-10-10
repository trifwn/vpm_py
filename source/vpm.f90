Module vpm_lib
   ! SPEED 
#ifdef USE_INTEL
      use mkl_service
#endif
   use MPI

   use vpm_vars
   use vpm_size
   use vpm_gcalc
   ! Constants
   use constants, only: pi, pi2, pi4
   use base_types, only: dp
   ! Printing
   use console_io, only: vpm_print, red, blue, green, nocolor, yellow, dummy_string, tab_level, VERBOCITY, &
                         particle_output_file_suffix, vpm_write_folder, pm_output_file_suffix
   use parvar, only:     print_particle_info, print_particle_positions, associate_particles
   use pmgrid, only:     print_velocity_stats, set_pm_velocities_zero, &
                         set_pm_deformations_zero, associate_velocities, &
                         associate_deformations 
   use file_io, only:   write_particles_hdf5, write_pm_solution_hdf5

   !  WhatToDo flags
   integer, parameter :: DEFINE_PROBLEM = 0,                                        &
                         SOLVE_VELOCITY = 1,                                        &
                         SOLVE_VELOCITY_DELATATION = 2,                             &
                         PROJECT_AND_SOLVE = 3,                                     &
                         INTERPOLATE = 4,                                           &
                         DIFFUSE = 5

   integer, parameter :: SOLVER_SERIAL_PMESH = 0, SOLVER_YAPS = 1
   integer            :: SOLVER = SOLVER_YAPS
   
   ! vpm_remesh.f90 contains
   interface remesh_particles_3d
      module subroutine remesh_particles_3d(iflag, npar_per_cell, XP_out, QP_out, GP_OUT, UP_OUT, NVR_out, cutoff_value)
         integer, intent(in)  :: iflag, npar_per_cell
         real(dp), intent(out),allocatable, target, dimension(:,:) :: XP_out, QP_out, GP_OUT, UP_OUT
         integer, intent(out) :: NVR_out
         real(dp), intent(in), optional  :: cutoff_value
      end subroutine remesh_particles_3d
   end interface remesh_particles_3d

contains

subroutine vpm_define_problem(                                 &
   NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in       &
)
   use parvar, only: NVR
   integer, intent(in) :: NTIME_in, NVR_in, NVR_size_in, neqpm_in
   real(dp), intent(in), dimension(:,:), target :: XP_in, QP_in
   integer :: ierr, my_rank, np
   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
   tab_level = 1

   ! Set Input
   ND       = 3
   NTIME_pm = NTIME_in
   neqpm    = neqpm_in
   call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in)
   if (NVR .eq. 0) return

   call define_sizes
   call set_pm_velocities_zero
   call set_pm_deformations_zero
   call allocate_sol_and_rhs(my_rank + 1)
end subroutine vpm_define_problem

subroutine vpm_interpolate(                                    &
   XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, neqpm_in,  &
   RHS_pm_ptr, deform_ptr           &
)
   use pmgrid, only: RHS_pm
   use parvar, only: NVR
   use MPI
   implicit none

   real(dp), intent(inout), target           :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
   integer, intent(inout)                    :: NVR_in
   integer, intent(in)                       :: neqpm_in, NVR_size_in
   real(dp), intent(out), pointer            :: RHS_pm_ptr(:, :, :, :)
   real(dp), intent(out), pointer,optional   :: deform_ptr(:, :, :, :)
   integer :: ierr, my_rank, np

   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
   tab_level = 1

   ! Set Input
   ND       =  3
   neqpm    = neqpm_in
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
      call calc_velocity_serial_3d(-1) 
   end if
   call interpolate_particles_parallel(1)
end subroutine vpm_interpolate

!> @brief Diffuses the particle vorticity through the PM grid
!>
!> @details
!>     The subroutine performs the following steps:
!>     - Project the particles to the PM grid
!>     - Diffuse the vorticity on the PM using D(RHS_PM) = NI \nabla^2 RHS_PM
!>     - Interpolate the PM grid to the particles
subroutine vpm_diffuse(                                         &
   XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, neqpm_in,   &
   RHS_pm_ptr, deform_ptr,           & 
   NI_in                                                        &
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
   neqpm    = neqpm_in
   NI       = NI_in
   call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in, UP_in, GP_in)
   if(present(deform_ptr)) call associate_deformations(deform_ptr)
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
      ! SOL_pm = -VIS \nabla \cdot RHS_pm
      call diffuse_vort_3d ! DIFFUSION OF VORTICITY
   end if
   ! WHEN ITYPEB = 2 WE GET THE GP FROM THE SOL_pm (DEFORMATION) and from QP
   call interpolate_particles_parallel(2)
end subroutine vpm_diffuse

subroutine vpm_project_and_solve(                              &
   NTIME_in,                                                   &
   XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in,                &
   RHS_pm_ptr                                                  &
)
   use MPI
   use pmgrid, only: RHS_pm
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
   tab_level= 1
   ND = 3
   NTIME_pm = NTIME_in
   neqpm    = neqpm_in
   call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in)
   call allocate_sol_and_rhs(my_rank + 1)
   RHS_pm_ptr => RHS_pm

   if (NVR .eq. 0) then
      print *, 'No particles to interpolate'
      return
   end if

   ! PARTICLES TO -> PM MEANING FROM Qs We get RHS_pm
   call project_particles_parallel
   call solve_problem
   ! Zero the velocities and deformations
   call set_pm_velocities_zero
   call set_pm_deformations_zero
end subroutine vpm_project_and_solve

subroutine vpm_solve_velocity(                                 &
   NTIME_in, XP_in, QP_in, UP_in, GP_in,                       &
   NVR_in, NVR_size_in, neqpm_in, RHS_pm_ptr, vel_ptr          &
)
   use MPI
   real(dp), intent(inout), target  :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
   integer, intent(inout)           :: NVR_in
   integer, intent(in)              :: neqpm_in, NVR_size_in
   real(dp), intent(out), pointer   :: RHS_pm_ptr(:, :, :, :)
   real(dp), intent(out), pointer   :: vel_ptr(:,:,:,:)
   integer, intent(in)              :: NTIME_in
   integer :: ierr, my_rank, np

   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
   call vpm_project_and_solve(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in, RHS_pm_ptr)
   call associate_velocities(vel_ptr)

   if (my_rank .eq. 0) then
      !call convect_first_order(Xbound,Dpm,NN,NN_bl)
      ! FROM THE SOLUTION OF PM WE GET THE VELOCITIES ON THE GRID
      call calc_velocity_serial_3d(0) ! VELOCITY STO PM
      ! SOL_pm is stil the solution of vorticity
      call print_velocity_stats
   end if
   call interpolate_particles_parallel(1) ! INTERPOLATION FROM PM TO PARTICLES
end subroutine vpm_solve_velocity

subroutine vpm_solve_velocity_delatation(                      &
   NTIME_in,                                                   &
   XP_in, QP_in, UP_in, GP_in,                                 &
   NVR_in, NVR_size_in, neqpm_in,                              &
   RHS_pm_ptr, vel_ptr, deform_ptr                             &
)
   use MPI
   implicit none
   real(dp), intent(inout), target            :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
   integer, intent(inout)                     :: NVR_in
   integer, intent(in)                        :: neqpm_in, NVR_size_in
   real(dp), intent(out), pointer             :: RHS_pm_ptr(:, :, :, :)
   real(dp), intent(out), pointer             :: vel_ptr(:, :,:,:)
   real(dp), intent(out), pointer, optional   :: deform_ptr(:, :,:,:)
   integer, intent(in)              :: NTIME_in
   
   integer :: ierr, my_rank, np

   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
   call vpm_project_and_solve(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in, RHS_pm_ptr)
   call associate_velocities(vel_ptr)
   if (present(deform_ptr)) call associate_deformations(deform_ptr)
   
   call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in, UP_in, GP_in)
   if (my_rank .eq. 0) then
      !call convect_first_order(Xbound,Dpm,NN,NN_bl)
      ! FROM THE SOLUTION OF PM WE GET THE VELOCITIES ON THE GRID
      call calc_velocity_serial_3d(0) ! VELOCITY STO PM
      ! SOL_pm is stil the solution of vorticity
      call print_velocity_stats
   end if
   call interpolate_particles_parallel(1) ! INTERPOLATION FROM PM TO PARTICLES
end subroutine vpm_solve_velocity_delatation

#ifdef USE_INTEL
   subroutine set_num_threads()
      if (SOLVER .eq. 1) then
         call mkl_set_num_threads(OMPTHREADS)
      end if
   end subroutine set_num_threads
#else
   subroutine set_num_threads()
   end subroutine set_num_threads
#endif

subroutine vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
               RHS_pm_ptr, vel_ptr, NTIME_in, NI_in, NVR_size_in,&
               deform_ptr)
      use parvar, only:    QP, XP, UP, GP, NVR, NVR_size,            &
                           print_particle_info, print_particle_positions, associate_particles
      use pmgrid, only:    SOL_pm, IDVPM,                            &
                           print_velocity_stats, set_pm_velocities_zero, &
                           set_pm_deformations_zero, associate_velocities, &
                           associate_deformations 
      use file_io, only:   write_particles_hdf5, write_pm_solution_hdf5

      Implicit None

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
      integer                           :: i ,j, k 

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
      call set_num_threads()
      if (my_rank .eq. 0) call print_arguements
      ! We set these input for compatibility with old version
      ! In the new modular approach we do not need to set these
      ND       = 3
      NTIME_pm = NTIME_in
      neqpm    = neqpm_in
      NI       = NI_in
      call associate_particles(NVR_in, NVR_size_in, neqpm_in, XP_in, QP_in, UP_in, GP_in)
      if (NVR .eq. 0) return
      if (associated(RHS_pm_ptr)) nullify (RHS_pm_ptr)
      call associate_velocities(vel_ptr)
      if (present(deform_ptr)) call associate_deformations(deform_ptr)

      select case(WhatToDo)
         case (DEFINE_PROBLEM)
            call vpm_define_problem(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in)
         case (SOLVE_VELOCITY)
            call vpm_solve_velocity(NTIME_in, XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in,        &
                                 neqpm_in, RHS_pm_ptr, vel_ptr)
         case (SOLVE_VELOCITY_DELATATION)
            call vpm_solve_velocity_delatation(NTIME_in, XP_in, QP_in, UP_in, GP_in, NVR_in,          &
                                 NVR_size_in, neqpm_in, RHS_pm_ptr, vel_ptr,     &
                                 deform_ptr)
         case (PROJECT_AND_SOLVE)
            call vpm_project_and_solve(NTIME_in, XP_in, QP_in, NVR_in, NVR_size_in, neqpm_in,         &
                                 RHS_pm_ptr)
         case (INTERPOLATE)
            call vpm_interpolate(XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, neqpm_in,           &
                                 RHS_pm_ptr, deform_ptr)
         case (DIFFUSE)
            call vpm_diffuse(XP_in, QP_in, UP_in, GP_in, NVR_in, NVR_size_in, neqpm_in,               &
                                 RHS_pm_ptr, deform_ptr, NI_in)
      end select
   return
   contains
   
      subroutine print_arguements
         ! write ARGUEMENTS
         if ((my_rank .eq. 0).and.(VERBOCITY>=1)) then
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
      end subroutine print_arguements
end subroutine vpm

   subroutine allocate_sol_and_rhs(n_block)
      use pmgrid, only: SOL_pm, RHS_pm, SOL_pm_bl, RHS_pm_bl
      implicit none
      integer, intent(in) :: n_block
      integer, dimension(3) :: NN_block
      ! neqpm        : is the number of equations to be solved
      ! NN           : is the number of cells in each direction
      ! NN_tmp       : is the number of cells in each direction (block cells)
      ! NN_bl_tmp    : is the start and finish of the cells in each direction
      ! Xbound       : is the boundary of the domain
      NN_block(1:3) = block_grids(n_block)%NN(1:3)

      ! SOL_pm block
      if (allocated(SOL_pm_bl)) then
         deallocate (SOL_pm_bl)
      end if
      allocate (SOL_pm_bl(neqpm, NN_block(1), NN_block(2), NN_block(3)))
      SOL_pm_bl = 0.d0

      ! RHS_pm block
      if (allocated(RHS_pm_bl)) then
         deallocate (RHS_pm_bl)
      end if 
      allocate (RHS_pm_bl(neqpm, NN_block(1), NN_block(2), NN_block(3)))
      RHS_pm_bl = 0.d0

      ! RHS_pm
      if (allocated(RHS_pm)) then
         deallocate (RHS_pm)
      end if
      allocate (RHS_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3))) 
      RHS_pm = 0.d0
      
      ! SOL_pm
      if (allocated(SOL_pm)) then
         deallocate (SOL_pm)
      end if
      allocate (SOL_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
      SOL_pm = 0.d0
   end subroutine allocate_sol_and_rhs
   
   subroutine solve_problem
      use pmgrid, only: RHS_pm, SOL_pm, RHS_pm_bl, DXpm, DYpm, DZpm, &
                        print_RHS_pm, set_pm_velocities_zero, set_pm_deformations_zero 
      use vpm_mpi, only: rhsbcast, rhsscat
      use serial_vector_field_operators, only: divergence, laplacian
      implicit none
      integer :: my_rank, ierr, i
      real(dp), allocatable :: div_wmega(:, :, :), laplace_LHS_pm(:, :, :, :)
      type(grid) :: my_block_grid
      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

      if (my_rank .eq. 0) then
         write (*, *) ""
         write (dummy_string, "(A)") 'Solving Particle Mesh'
         call vpm_print(dummy_string, blue, 1)
         st = MPI_WTIME()
      endif
      tab_level = tab_level + 1
      if (my_rank .eq. 0) then
         if (VERBOCITY >=2) then
            call print_RHS_pm()
         end if
      end if
      ! Print the Computational Domain and the number of blocks
      if (my_rank .eq. 0) then
         write (dummy_string, "(A)") 'Domain information'
         call vpm_print(dummy_string, blue, 1)
         write (dummy_string, "(A)") achar(9)//'Bounds of the fine domain:'
         call vpm_print(dummy_string, yellow, 1)
         write (dummy_string, "(A,3F8.3)") achar(9)//'X:', fine_grid%Xbound(1), fine_grid%Xbound(4)
         call vpm_print(dummy_string, yellow, 1)
         write (dummy_string, "(A,3F8.3)") achar(9)//'Y:', fine_grid%Xbound(2), fine_grid%Xbound(5)
         call vpm_print(dummy_string, yellow, 1)
         write (dummy_string, "(A,3F8.3)") achar(9)//'Z:', fine_grid%Xbound(3), fine_grid%Xbound(6)
         call vpm_print(dummy_string, yellow, 1)
         write (dummy_string, "(A,3I5)") achar(9)//'Size of the fine domain:', fine_grid%NN(1), &
                                                                               fine_grid%NN(2), &
                                                                               fine_grid%NN(3) 
         call vpm_print(dummy_string, yellow, 1)
         write (dummy_string, "(A,3I5)") achar(9)//'Size of the coarse domain:', coarse_grid%NN(1), & 
                                                                                 coarse_grid%NN(2), & 
                                                                                 coarse_grid%NN(3)
         call vpm_print(dummy_string, yellow, 1)
         write (dummy_string, "(A,3I5)") achar(9)//'Size of the block domains:', block_grids(1)%NN(1), &
                                                                                 block_grids(1)%NN(2), &
                                                                                 block_grids(1)%NN(3)
         write (dummy_string, "(A,I5)") achar(9)//'Number of blocks:', NBlocks
         call vpm_print(dummy_string, yellow, 1)
      end if

      !call diffuse_vort_3d
      
      ! SCATTER PROBLEM TO ALL PROCESSORS
      call rhsbcast(RHS_pm, fine_grid%NN, neqpm) ! RHS PM -> TO ALL PROCESSORS -> RHS_pm_BL
      IF (SOLVER .ne. 0) then
         my_block_idx = my_rank + 1
         my_block_grid = block_grids(my_block_idx)
         call rhsscat(fine_grid, RHS_pm, my_block_grid, RHS_pm_bl, nb_i, nb_j, nb_k)
      end if

      if (my_rank .eq. 0) then
         et = MPI_WTIME()
         write (dummy_string, "(A)") 'Broadcasting/Scattering RHS'
         call vpm_print(dummy_string, blue, 1)
         write (dummy_string, "(A,I5,A,F8.2,A)") &
               achar(9)//'finished in:', int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
         call vpm_print(dummy_string, yellow, 1)
      end if

      ! SOLVE THE PROBLEM
      if (my_rank .eq. 0) st = MPI_WTIME()
      call pmesh_solve
      if (my_rank .eq. 0) then
         et = MPI_WTIME()
         write (dummy_string, "(A,I5,A,F8.2,A)") &
               'Total time for solving the Particle Mesh', &
               int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
         call vpm_print(dummy_string, blue, 1)
      end if

      if (my_rank .eq. 0) then
         ! We need to write the residuals of the solution
         div_wmega = divergence(RHS_pm, DXpm , DYpm, DZpm)
         laplace_LHS_pm = laplacian(SOL_pm, DXpm, DYpm, DZpm)
         ! Print the mean, max and min div_u
         write (*, *) ""
         write (dummy_string, "(A)") 'Divergence of the Psi field'
         call vpm_print(dummy_string, yellow, 1)
         write (dummy_string, "(A, E10.4, A,E10.4, A, E10.4)") "ΔPsi"// &
               achar(9)//" min : ", minval(div_wmega), &
               achar(9)//" max : ", maxval(div_wmega), &
               achar(9)//" mean: ", sum(div_wmega)/(fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
         call vpm_print(dummy_string, blue, 1)

         write (dummy_string, "(A)") 'Residuals of the solution'
         call vpm_print(dummy_string, yellow, 1)
         do i = 1,neqpm
            ! For each equation write the laplacian - RHS_pm
            write (dummy_string, "(A, I3, A)") '   Equation =', i, "Δf = RHS"
            call vpm_print(dummy_string, blue, 1)
            write (dummy_string, "(A, E10.4, A, E10.4, A, E10.4)") achar(9)//'Forcing (RHS)'//     &
                  achar(9)//'min : ', minval(RHS_pm(i,:,:,:)),                                     &
                  achar(9)//'max : ', maxval(RHS_pm(i,:,:,:)),                                     &
                  achar(9)//'mean: ', sum(RHS_pm(i,:,:,:))/(fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
            call vpm_print(dummy_string, nocolor, 1)                                                   
            write (dummy_string, "(A, E10.4, A, E10.4, A, E10.4)") achar(9)//"Solution"//          &
                  achar(9)//'min : ', minval(SOL_pm(i,:,:,:)),                                     &
                  achar(9)//'max : ', maxval(SOL_pm(i,:,:,:)),                                     &
                  achar(9)//'mean: ', sum(SOL_pm(i,:,:,:))/(fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3)) 
            call vpm_print(dummy_string, nocolor, 1)
            write (dummy_string, "(A, E10.4, A, E10.4, A, E10.4)") achar(9)//'Res:=Δf-RHS'//  &
                  achar(9)//'min : ', minval(laplace_LHS_pm(i,:,:,:)- RHS_PM(i,:,:,:)),            &
                  achar(9)//'max : ', maxval(laplace_LHS_pm(i,:,:,:)- RHS_PM(i,:,:,:)),            &
                  achar(9)//'mean: ', sum(laplace_LHS_pm(i,:,:,:) - RHS_PM(i,:,:,:))               &
                                             /(fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3))
            call vpm_print(dummy_string, nocolor, 1)                                                   
            write (*, *) ""
         enddo 
         tab_level = tab_level - 1

         ! Deallocate the memory
         deallocate (div_wmega, laplace_LHS_pm)
      end if
      
      ! CLEAR VELOCITIES
      call set_pm_velocities_zero
      call set_pm_deformations_zero
   end subroutine solve_problem

   subroutine pmesh_solve !
      use pmgrid, only: ncoarse, SOL_pm, RHS_pm, SOL_pm_bl, RHS_pm_bl
      use parvar, only: XP, QP, NVR
      use pmlib, only:  pmesh
      use yaps, only:   yaps3d
      use vpm_mpi, only: solget, velbcast
      implicit none
      integer    :: my_rank, ierr, np
      integer    :: eq_num
      integer    :: i
      type(grid) :: my_block_grid
      ! LOCAL ARGUEMENTS DUE TO LEGACY YAPS NOT ACCEPTING GRID TYPE
      integer    ::  NN_coarse(3), NN_bl_coarse(6), NN_block(3, NBlocks), NN_bl_block(6, NBlocks)
      real(dp)   ::  Xbound_coarse(6), Dpm_coarse(3), Xbound_block(6, NBlocks), Xbound(6), Dpm_fine(3)
                 
      !Yaps or Serial Pmesh

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

      IF ((SOLVER .eq. 0).or.(np .eq. 1)) THEN
         if (my_rank .eq. 0) then
            st = MPI_WTIME()
            write (dummy_string, "(A)") 'Solving PM with Serial Pmesh (One Processor)'
            call vpm_print(dummy_string, blue, 1)
         endif

         IF (my_rank .eq. 0) then
            SOL_pm(1:neqpm, :, :, :) = 0.0
            itree = iyntree
            iynbc = 1   !for infinite domain bc
            call pmesh(SOL_pm, RHS_pm, QP, XP, &
                       fine_grid%Xbound, fine_grid%DPm, fine_grid%NN, fine_grid%NN_bl, &
                       ND, ibctyp, 1, neqpm, iynbc, NVR, itree, ilevmax)
            ! call calc_velocity_serial_3d(1)
         end if
         !--------------------------------------------
         ! call velbcast_3d
      else
         if (my_rank .eq. 0) then 
            write (dummy_string, "(A)") 'Solving PM with YAPS'
            call vpm_print(dummy_string, blue, 1)
            st = MPI_WTIME()
         endif

         iret = 0
         
         ! Coarse Grid
         Dpm_coarse(1:3)    = coarse_grid%Dpm(1:3)
         Xbound_coarse(1:6) = coarse_grid%Xbound(1:6)
         NN_coarse(1:3)     = coarse_grid%NN(1:3)
         NN_bl_coarse(1:6)  = coarse_grid%NN_bl(1:6)
         
         ! Block Grid
         do i = 1, NBlocks 
            NN_bl_block(1:6, i)  = block_grids(i)%NN_bl(1:6)
            Xbound_block(1:6, i) = block_grids(i)%Xbound(1:6)
            NN_block(1:3, i)     = block_grids(i)%NN
         end do

         ! FINE
         Dpm_fine(1:3)    = fine_grid%Dpm(1:3)
         Xbound(1:6) = fine_grid%Xbound(1:6)

         call yaps3d(SOL_pm_bl, RHS_pm_bl, Xbound_block, Xbound_coarse, Dpm_fine, Dpm_coarse, NN_block, NN_bl_block, &
                     NN_coarse, NN_bl_coarse, ND, NBlocks, ibctyp, 1, neqpm, ncoarse, NBI, NBJ, NBK, nb_i, nb_j, nb_k, &
                     iret, iyntree, ilevmax, neqpm)

         ! GATHER SOLUTION
         my_block_idx = my_rank + 1
         my_block_grid = block_grids(my_block_idx)
         call solget(NBlocks, NBI, NBJ, NBK, my_block_grid, block_grids, fine_grid, SOL_pm, SOL_pm_bl) 

         if (my_rank .eq. 0) then
            write (dummy_string, "(A)") 'Final PM solution values'
            call vpm_print(dummy_string, blue, 2)
            do eq_num = 1, neqpm
               write (dummy_string, "(A, I3,A,F8.3,A,F8.3)") &
                  achar(9)//'Equation =', eq_num, "->"// &
                  achar(9)//'min(SOL_pm)', minval(SOL_pm(eq_num,:,:,:)), &
                  achar(9)//'max(SOL_pm)', maxval(SOL_pm(eq_num,:,:,:))   
               call vpm_print(dummy_string, blue, 2)
            enddo
         end if
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      !--------------------------------------------

   end subroutine pmesh_solve

   subroutine convect_first_order(DT_convection)
      use parvar, only: NVR, XP, UP, GP, QP
      real(dp), intent(in) :: DT_convection
      integer              :: i

      do i = 1, NVR
         XP(1:3, i) = XP(1:3, i) + UP(1:3, i)*DT_convection
         QP(1:neqpm, i) = QP(1:neqpm, i) + GP(1:neqpm, i)*DT_convection
      end do
   end subroutine convect_first_order

   subroutine project_particles_parallel
      use pmgrid, only:    RHS_pm, SOL_pm, IDVPM
      use parvar, only:    NVR, XP, QP, XP_scatt, QP_scatt, NVR_projtype_scatt, NVR_p
      use projlib, only:   projlibinit, project_particles_3D, project_vol3d
      use vpm_mpi, only:   particles_scat, proj_gath
      
      integer, allocatable              :: ieq(:)
      real(dp), allocatable             :: QINF(:)
      integer                           :: my_rank, np, ierr, i

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_Size(MPI_COMM_WORLD, np, ierr)

      ! FILLS RHS_pm FROM PARTICLES
      if (my_rank .eq. 0) then 
         write (dummy_string, "(A)") 'Projecting Particles to PM Routine' 
         call vpm_print(dummy_string, blue, 1)
      endif
      tab_level = tab_level + 1
      ! BCAST NVR
      call MPI_BCAST(NVR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      ! SPLIT PARTICLES on PROCESSORS
      NVR_p = NVR/np
      ! ROUND
      if (my_rank .eq. 0) then 
         NVR_p = NVR_p + mod(NVR, np)
      end if 
      
      ! Allocate
      allocate (XP_scatt(3, NVR_p), QP_scatt(neqpm+1, NVR_p), NVR_projtype_scatt(NVR_p))
      allocate (ieq(neqpm), QINF(neqpm))

      NVR_projtype_scatt = interf_iproj
      QINF = 0.d0
      do i = 1, neqpm
         ieq(i) = i
      end do
      
      ! SCATTER PARTICLES ON EACH PROCESSOR
      call particles_scat
      ! INITIALIZATION OF PROJECTION LIBRARY
      if (my_rank .eq. 0) then
         write (dummy_string, "(A)") "Initializing Projection Library"
         call vpm_print(dummy_string, blue, 1) 
      endif

      call projlibinit(fine_grid%Xbound, fine_grid%Dpm, fine_grid%NN, fine_grid%NN_bl, IDVPM, ND)

      ! PROJECT PARTICLES TO PM
      st = MPI_WTIME()
      if (my_rank .eq. 0) then
         write (dummy_string, "(A)") "Projecting Particles on PM (RHS_pm)"
         call vpm_print(dummy_string, blue, 1) 
      endif

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! From QP_scatt we get RHS_pm
      call project_particles_3D(RHS_pm, QP_scatt, XP_scatt, NVR_projtype_scatt, NVR_p, neqpm, ieq, neqpm, QINF, NVR_p)
      
      if (my_rank .eq. 0) then
         write (dummy_string, "(A)") "Gathering RHS from all processors"
         call vpm_print(dummy_string, blue, 1) 
      endif

      call proj_gath(fine_grid%NN) ! RHS IS NOW FILLED

      tab_level = tab_level - 1
      if (my_rank .eq. 0) then
         write (dummy_string, "(A)") achar(9)//"Normalizing RHS by volume"
         call vpm_print(dummy_string, blue, 1)

         ! RHS_pm(neqpm+1,:,:,:)=DVpm
         ! RHS_pm IS NOW NORMALIZED BY VOLUME (DENSITY)
         call project_vol3d(RHS_pm, neqpm, ieq, neqpm, IDVPM) 

         et = MPI_WTIME()
         write (dummy_string, "(A, I2, A, F8.3, A)") &
            achar(9)//'finished in:'//achar(9), int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
         call vpm_print(dummy_string, yellow, 1)
      end if

      ! DEALLOCATE
      deallocate (ieq, QINF)
      deallocate (XP_scatt, QP_scatt, NVR_projtype_scatt)
   end subroutine project_particles_parallel

   subroutine interpolate_particles_parallel(itypeb)
      use parvar, only:    NVR, XP_scatt, QP_scatt, UP_scatt, GP_scatt, NVR_p
      use pmgrid, only:    velocity_pm,  deform_pm, RHS_pm, SOL_pm
      use vpm_interpolate, only: back_to_particles_3D
      use vpm_mpi, only:   rhsbcast, particles_scat, particles_gath, velbcast, defbcast

      integer, intent(in)               :: itypeb
      integer, allocatable              :: ieq(:)
      real(dp), allocatable             :: QINF(:)
      integer                           :: ierr, my_rank, np, i

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_Size(MPI_COMM_WORLD, np, ierr)

      if (my_rank .eq. 0) st = MPI_WTIME()

      ! BROADCASTING
      call rhsbcast(RHS_pm, fine_grid%NN, neqpm)
      call rhsbcast(SOL_pm, fine_grid%NN, neqpm)
      if (itypeb .eq. 1) call velbcast
      if (itypeb .eq. 2) call defbcast
      call MPI_BCAST(NVR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      ! ALLOCATE
      NVR_p = NVR/np
      if (my_rank .eq. 0) NVR_p = NVR_p + mod(NVR, np)
      allocate (XP_scatt(3, NVR_p), QP_scatt(neqpm+1, NVR_p), &
                UP_scatt(3, NVR_p), GP_scatt(3, NVR_p))
      UP_scatt = 0.d0; GP_scatt = 0.d0

      ! SCATTERING XP AND QP
      call particles_scat

      if (my_rank .eq. 0) st = MPI_WTIME()

      ! WHEN ITYPEB = 1 WE GET THE UP AND GP From the velocity field
      ! WHEN ITYPEB = 2 WE GET THE GP FROM THE deformation
      call back_to_particles_3D(SOL_pm, XP_scatt, QP_scatt, UP_scatt, GP_scatt,  &
                                velocity_pm, deform_pm,                          &
                                NVR_p, interf_iproj, itypeb, NVR_p               )
      ! GATHERS XP, QP, UP, GP
      call particles_gath

      deallocate (XP_scatt, QP_scatt, UP_scatt, GP_scatt)
      if (my_rank .eq. 0) then
         et = MPI_WTIME()
         write (dummy_string, "(A)") 'Interpolating Particles'
         call vpm_print(dummy_string, blue, 1)
         write (dummy_string, "(A,I5,A,F8.2,A)") &
               achar(9)//'finished in:', int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
         call vpm_print(dummy_string, yellow, 1)
      end if

   end subroutine interpolate_particles_parallel
end Module vpm_lib
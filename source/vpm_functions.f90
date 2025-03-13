module vpm_functions
   use MPI

   use vpm_vars
   use vpm_size

   ! Constants
   use constants, only: pi, pi2, pi4
   use vpm_types, only: dp
   ! Printing
   use console_io, only: vpm_print, red, blue, green, nocolor, yellow, dummy_string, tab_level, VERBOCITY
   use parvar, only: print_particle_info, print_particle_positions, associate_particles
   use pmgrid, only: print_velocity_stats, print_vortex_stretching_stats
   implicit none

contains
   subroutine allocate_sol_and_rhs(n_block)
      use pmgrid, only: SOL_pm, RHS_pm, SOL_pm_bl, RHS_pm_bl
      use MPI
      implicit none
      integer, intent(in) :: n_block
      integer, dimension(3) :: NN_block
      integer :: my_rank, ierr
      ! neqpm        : is the number of equations to be solved
      ! NN           : is the number of cells in each direction
      ! NN_tmp       : is the number of cells in each direction (block cells)
      ! NN_bl_tmp    : is the start and finish of the cells in each direction
      ! Xbound       : is the boundary of the domain
      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      NN_block = block_grids(n_block)%NN

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

   subroutine solve_problem(RHS_pm, SOL_pm, RHS_pm_bl, SOL_pm_bl)
      use pmgrid, only: print_RHS_pm
      use vpm_mpi, only: rhsbcast, rhsscat
      use serial_vector_field_operators, only: divergence, laplacian
      implicit none
      real(dp), allocatable, target, intent(inout) :: RHS_pm(:, :, :, :), RHS_pm_bl(:, :, :, :)
      real(dp), allocatable, target, intent(inout) :: SOL_pm(:, :, :, :), SOL_pm_bl(:, :, :, :)
      integer                                      :: my_rank, ierr
      type(cartesian_grid)                         :: my_block_grid
      ! LOCAL TEMPORARY TBR

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

      if (my_rank .eq. 0) then
         write (*, *) ""
         write (dummy_string, "(A)") 'Solving Particle Mesh'
         call vpm_print(dummy_string, blue, 1)
         st = MPI_WTIME()
      end if
      tab_level = tab_level + 1
      if (my_rank .eq. 0) then
         if (VERBOCITY >= 2) then
            call print_RHS_pm(RHS_pm)
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
         write (dummy_string, "(A,3I5,A,I15)") achar(9)//'Size of the fine domain:', fine_grid%NN(1), &
            fine_grid%NN(2), &
            fine_grid%NN(3), &
            achar(9)//achar(9)//'Total Cells:', fine_grid%NN(1)*fine_grid%NN(2)*fine_grid%NN(3)
         call vpm_print(dummy_string, yellow, 1)
         write (dummy_string, "(A,3I5,A,I15)") achar(9)//'Size of the coarse domain:', coarse_grid%NN(1), &
            coarse_grid%NN(2), &
            coarse_grid%NN(3), &
            achar(9)//'Total Cells:', coarse_grid%NN(1)*coarse_grid%NN(2)*coarse_grid%NN(3)

         call vpm_print(dummy_string, yellow, 1)
         write (dummy_string, "(A,I5)") achar(9)//'Number of blocks:', NBlocks
         call vpm_print(dummy_string, yellow, 1)
         write (dummy_string, "(A,3I5,A,I15)") achar(9)//'Size of the block domains:', block_grids(1)%NN(1), &
            block_grids(1)%NN(2), &
            block_grids(1)%NN(3), &
            achar(9)//'Total Cells:', block_grids(1)%NN(1)*block_grids(1)%NN(2)*block_grids(1)%NN(3)
         call vpm_print(dummy_string, yellow, 1)
      end if

      ! call diffuse_vort_3d

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
      call pmesh_solve(RHS_pm, SOL_pm, RHS_pm_bl, SOL_pm_bl)
      if (my_rank .eq. 0) then
         et = MPI_WTIME()
         write (dummy_string, "(A,I5,A,F8.2,A)") &
            'Total time for solving the Particle Mesh', &
            int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
         call vpm_print(dummy_string, blue, 1)
      end if
   end subroutine solve_problem

   subroutine pmesh_solve(RHS_pm, SOL_pm, RHS_pm_bl, SOL_pm_bl)
      use pmgrid, only: ncoarse
      use parvar, only: XP, QP, NVR
      use pmlib, only: pmesh
      use yaps, only: yaps3d
      use vpm_mpi, only: solget, velbcast
      implicit none
      real(dp), allocatable, target, intent(inout) :: RHS_pm(:, :, :, :), RHS_pm_bl(:, :, :, :)
      real(dp), allocatable, target, intent(inout) :: SOL_pm(:, :, :, :), SOL_pm_bl(:, :, :, :)
      integer    :: my_rank, ierr, np
      integer    :: i
      type(cartesian_grid) :: my_block_grid

      ! LOCAL ARGUEMENTS DUE TO LEGACY YAPS NOT ACCEPTING GRID TYPE
      integer    ::  NN_coarse(3), NN_bl_coarse(6), NN_block(3, NBlocks), NN_bl_block(6, NBlocks)
      real(dp)   ::  Xbound_coarse(6), Dpm_coarse(3), Xbound_block(6, NBlocks), Xbound(6), Dpm_fine(3)

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

      !Yaps or Serial Pmesh
      IF ((SOLVER .eq. 0) .or. (np .eq. 1)) THEN
         if (my_rank .eq. 0) then
            st = MPI_WTIME()
            write (dummy_string, "(A)") 'Solving PM with Serial Pmesh (One Processor)'
            call vpm_print(dummy_string, blue, 1)
         end if

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
         end if

         iret = 0

         ! Coarse Grid
         Dpm_coarse(1:3) = coarse_grid%Dpm(1:3)
         Xbound_coarse(1:6) = coarse_grid%Xbound(1:6)
         NN_coarse(1:3) = coarse_grid%NN(1:3)
         NN_bl_coarse(1:6) = coarse_grid%NN_bl(1:6)

         ! Block Grid
         do i = 1, NBlocks
            NN_bl_block(1:6, i) = block_grids(i)%NN_bl(1:6)
            Xbound_block(1:6, i) = block_grids(i)%Xbound(1:6)
            NN_block(1:3, i) = block_grids(i)%NN
         end do

         ! FINE
         Dpm_fine(1:3) = fine_grid%Dpm(1:3)
         Xbound(1:6) = fine_grid%Xbound(1:6)

         call yaps3d(SOL_pm_bl, RHS_pm_bl, Xbound_block, Xbound_coarse, Dpm_fine, Dpm_coarse, NN_block, NN_bl_block, &
                     NN_coarse, NN_bl_coarse, ND, NBlocks, ibctyp, 1, neqpm, ncoarse, NBI, NBJ, NBK, nb_i, nb_j, nb_k, &
                     iret, iyntree, ilevmax, neqpm)

         ! GATHER SOLUTION
         my_block_idx = my_rank + 1
         my_block_grid = block_grids(my_block_idx)
         call solget(NBlocks, NBI, NBJ, NBK, my_block_grid, block_grids, fine_grid, SOL_pm, SOL_pm_bl)
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

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
      use pmgrid, only: RHS_pm, IDVPM
      use parvar, only: NVR, XP_scatt, QP_scatt, NVR_projtype_scatt, NVR_p
      use projlib, only: projlibinit, project_particles_3D, project_vol3d
      use vpm_mpi, only: particles_scat, proj_gath

      integer, allocatable              :: ieq(:)
      real(dp), allocatable             :: QINF(:)
      integer                           :: my_rank, np, ierr, i

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_Size(MPI_COMM_WORLD, np, ierr)

      ! FILLS RHS_pm FROM PARTICLES
      if (my_rank .eq. 0) then
         write (dummy_string, "(A)") 'Projecting Particles to PM Routine'
         call vpm_print(dummy_string, blue, 1)
      end if
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
      allocate (XP_scatt(3, NVR_p), QP_scatt(neqpm + 1, NVR_p), NVR_projtype_scatt(NVR_p))
      allocate (ieq(neqpm), QINF(neqpm))

      NVR_projtype_scatt = interf_iproj
      QINF = 0.d0
      do i = 1, neqpm
         ieq(i) = i
      end do

      ! SCATTER PARTICLES ON EACH PROCESSOR
      call particles_scat
      ! INITIALIZATION OF PROJECTION LIBRARY
      call projlibinit(fine_grid, IDVPM, ND)
      ! PROJECT PARTICLES TO PM
      st = MPI_WTIME()
      if (my_rank .eq. 0) then
         write (dummy_string, "(A)") "Projecting Particles on PM (RHS_pm)"
         call vpm_print(dummy_string, blue, 1)
      end if

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! From QP_scatt we get RHS_pm
      call project_particles_3D(RHS_pm, QP_scatt, XP_scatt, NVR_projtype_scatt, NVR_p, neqpm, ieq, neqpm, QINF, NVR_p)
      call proj_gath ! RHS IS NOW FILLED

      if (my_rank .eq. 0) then
         write (dummy_string, "(A)") "Normalizing RHS by volume"
         call vpm_print(dummy_string, blue, 1)
         ! RHS_pm IS NOW NORMALIZED BY VOLUME (DENSITY)
         call project_vol3d(RHS_pm, neqpm, ieq, neqpm, IDVPM)
      end if

      if (my_rank .eq. 0) then
         et = MPI_WTIME()
         write (dummy_string, "(A, I2, A, F8.3, A)") &
            achar(9)//'finished in:'//achar(9), int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
         call vpm_print(dummy_string, yellow, 1)
      end if
      tab_level = tab_level - 1
      ! DEALLOCATE
      deallocate (ieq, QINF)
      deallocate (XP_scatt, QP_scatt, NVR_projtype_scatt)
   end subroutine project_particles_parallel

   subroutine interpolate_particles_parallel(itypeb)
      use parvar, only: NVR, XP_scatt, QP_scatt, UP_scatt, GP_scatt, NVR_p, QP, NVR_size
      use pmgrid, only: velocity_pm, deform_pm, RHS_pm
      use vpm_interpolate, only: back_to_particles_3D
      use vpm_mpi, only: rhsbcast, particles_scat, particles_gath, velbcast, defbcast

      integer, intent(in)               :: itypeb
      integer                           :: ierr, my_rank, np

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_Size(MPI_COMM_WORLD, np, ierr)

      if (my_rank .eq. 0) st = MPI_WTIME()

      ! BROADCASTING
      if (itypeb .eq. 1) call velbcast
      call rhsbcast(RHS_pm, fine_grid%NN, neqpm)
      call defbcast
      call MPI_BCAST(NVR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      ! ALLOCATE
      NVR_p = NVR/np
      if (my_rank .eq. 0) NVR_p = NVR_p + mod(NVR, np)
      allocate (XP_scatt(3, NVR_p), QP_scatt(neqpm + 1, NVR_p), &
                UP_scatt(3, NVR_p), GP_scatt(3, NVR_p))
      XP_scatt = 0.d0; QP_scatt = 0.d0
      UP_scatt = 0.d0; GP_scatt = 0.d0

      ! SCATTERING XP AND QP
      call particles_scat

      ! WHEN ITYPEB = 1 WE GET THE UP AND GP From the velocity field and the deformation
      ! WHEN ITYPEB = 2 WE GET THE GP FROM THE deformation
      call back_to_particles_3D(XP_scatt, QP_scatt, UP_scatt, GP_scatt, &
                                velocity_pm, deform_pm, &
                                NVR_p, interf_iproj, itypeb, NVR_p)
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

end module vpm_functions

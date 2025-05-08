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
      IF ((SOLVER .eq. 0) ) then!.or. (np .eq. 1)) THEN
         ! ---------------------------------------------------------
         ! Serial Pmesh solver
         ! ---------------------------------------------------------

         if (my_rank .eq. 0) then
            st = MPI_WTIME()
            write (dummy_string, "(A)") 'Solving PM with Serial Pmesh (One Processor)'
            call vpm_print(dummy_string, blue, 1)
         end if

         IF (my_rank .eq. 0) then
            SOL_pm(:, :, :, :) = 0.0
            ! itree = iyntree
            iynbc = 1   !for infinite domain bc
            call pmesh(SOL_pm, RHS_pm, QP, XP, &
                       fine_grid%Xbound, fine_grid%DPm, fine_grid%NN, fine_grid%NN_bl, &
                       ND, ibctyp, 1, neqpm, iynbc, NVR, iyntree, ilevmax)
            ! call calc_velocity_serial_3d(1)
         end if
         !--------------------------------------------
         ! call velbcast_3d
      elseif (SOLVER .eq. 1) then
         ! ---------------------------------------------------------
         ! YAPS parallel solver
         ! ---------------------------------------------------------

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
      elseif (SOLVER .eq. 2) then
         ! ---------------------------------------------------------
         ! Mudpack serial solver 
         ! ---------------------------------------------------------
         if (my_rank .eq. 0) then
            write (dummy_string, "(A)") 'Solving PM with Mudpack'
            call vpm_print(dummy_string, blue, 1)
            st = MPI_WTIME()

            call solve_mudpack

            et = MPI_WTIME()
            write (dummy_string, "(A,I5,A,F8.2,A)") &
               'Total time for solving the Particle Mesh', &
               int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
            call vpm_print(dummy_string, blue, 1)
         end if
      end if

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      contains 

         subroutine solve_mudpack
            use mudpack_sp 
            use console_io, only: print_stats_rank4, print_stats_rank3
            implicit none
            !--------------------------------------------------------------------
            ! Arguments:
            !   SOL_pm(3, nx, ny, nz)  : solution array (to be updated)
            !   RHS_pm(3, nx, ny, nz)  : right-hand side array
            !--------------------------------------------------------------------
            
            ! Get domain dimensions from fine_grid
            integer :: nx, ny, nz
            ! Local arrays for solution and RHS for one system
            real(dp), allocatable   :: phi(:,:,:)
            real(dp), allocatable   :: rhs(:,:,:)
            real(dp), allocatable   :: work(:)
            integer :: llwork
            ! Parameter arrays for Mudpack
            real(dp), dimension(8)  :: fprm
            integer, dimension(22)  :: iprm
            integer, dimension(4)   :: mgopt
            ! Local variables
            integer :: ierror, system
            integer :: ii, jj, kk
            integer :: ixp, jyq, kzr, iex, jey, kez
            integer :: temp, e
            
            nx = fine_grid%NN(1)
            ny = fine_grid%NN(2)
            nz = fine_grid%NN(3)
            
            ! Compute required workspace size (based on tmud3sp.f example)
            llwork = 3 * (7 * (nx+2) * (ny+2) * (nz+2)) / 2
            
            ! Initialize mgopt to zeros
            mgopt(1) = 2
            mgopt(2) = 2
            mgopt(3) = 1
            mgopt(4) = 3
            
            ! Set domain boundaries from fine_grid%Xbound
            ! Assumed ordering: [xa, xb, yc, yd, ze, zf]
            fprm(1)   = fine_grid%Xbound(1)
            fprm(2)   = fine_grid%Xbound(4)

            fprm(3)   = fine_grid%Xbound(2)
            fprm(4)   = fine_grid%Xbound(5)
            
            fprm(5)   = fine_grid%Xbound(3)
            fprm(6)   = fine_grid%Xbound(6)

            ! Set up fprm parameters
            fprm(7)   = 1e-8     ! Tolerance
            fprm(8)   = 0.0      ! Output 
            
            ! Set up iprm parameters
            ! iprm(1) is used to control the call (0: initialization, 1: solve)
            iprm(2) = 1   ! Lower X boundary: Dirichlet
            iprm(3) = 1   ! Upper X boundary: Dirichlet
            iprm(4) = 1   ! Lower Y boundary: Dirichlet
            iprm(5) = 1   ! Upper Y boundary: Dirichlet
            iprm(6) = 1   ! Lower Z boundary: Dirichlet
            iprm(7) = 1   ! Upper Z boundary: Dirichlet
            
            ! Set grid factors for a single grid level (trivial multigrid)
               ! Decompose nx: find odd part (ixp) and exponent+1 (iex)
               temp = nx - 1
               e = 0
               do while (temp .gt. 2 .and. mod(temp, 2) == 0)
                  temp = temp / 2
                  e = e + 1
               end do
               ixp = temp
               iex = e + 1

               ! Decompose ny: find odd part (jyq) and exponent+1 (jey)
               temp = ny - 1
               e = 0
               do while (temp .gt. 2 .and. mod(temp, 2) == 0)
                  temp = temp / 2
                  e = e + 1
               end do
               jyq = temp
               jey = e + 1

               ! Decompose nz: find odd part (kzr) and exponent+1 (kez)
               temp = nz - 1
               e = 0
               do while (temp .gt. 2 .and. mod(temp, 2) == 0)
                  temp = temp / 2
                  e = e + 1
               end do
               kzr = temp
               kez = e + 1

               ! Output the results
               ! print*, "Decomposition results:"
               ! print*, "nx = ", nx, " -> ixp = ", ixp, ", iex = ", iex
               ! print*, "ny = ", ny, " -> jyq = ", jyq, ", jey = ", jey
               ! print*, "nz = ", nz, " -> kzr = ", kzr, ", kez = ", kez

               ! Verification: iprm(14:16) should match nx, ny, nz respectively
               ! print*, "Verification:"
               ! print*, "iprm(14) = ixp * 2^(iex-1) + 1 = ", ixp * 2**(iex-1) + 1
               ! print*, "iprm(15) = jyq * 2^(jey-1) + 1 = ", jyq * 2**(jey-1) + 1
               ! print*, "iprm(16) = kzr * 2^(kez-1) + 1 = ", kzr * 2**(kez-1) + 1

            iprm(8)  = ixp
            iprm(9)  = jyq
            iprm(10) = kzr
            iprm(11) = iex
            iprm(12) = jey
            iprm(13) = kez
            iprm(14) = ixp*(2**(iex-1)) + 1
            iprm(15) = jyq*(2**(jey-1)) + 1 
            iprm(16) = kzr*(2**(kez-1)) + 1
            iprm(17) = 1      ! Use the initial guess (phi initially zero)
            iprm(18) = 1000   ! Maximum number of iterations/cycles
            iprm(19) = 0      ! Relaxation method (e.g., Gauss-Seidel)
            iprm(20) = llwork ! Length of the work array
            iprm(21) = 0      ! (Output: minimum required workspace; set by mud3sp)
            iprm(22) = 0      ! (Output: iteration count; set by mud3sp)
                        
            ! Allocate workspace and local solution arrays
            allocate(phi(nx, ny, nz))
            allocate(rhs(nx, ny, nz))
            allocate(work(llwork))
            
            do system = 1, neqpm

               do kk = 1, nz
                  do ii = 1, nx
                     do jj = 1, ny
                        rhs(ii, jj, kk) = RHS_pm(system, ii, jj, kk)
                        phi(ii, jj, kk) = 0.0_dp
                     end do
                  end do
               end do

               
               !----------------------------------------------------------------
               ! First call: initialization/discretization step

               iprm(1) = 0
               call mud3sp(iprm, fprm, work, cfx, cfy, cfz, bndc, rhs, phi, mgopt, ierror)
               if (ierror /= 0) then
                  write(*,*) 'Error during mud3sp initialization for system ', system, &
                             ' ierror = ', ierror
               endif
               
               !----------------------------------------------------------------
               ! Second call: iterative solve
               iprm(1) = 1
               call mud3sp(iprm, fprm, work, cfx, cfy, cfz, bndc, rhs, phi, mgopt, ierror)
               if (ierror /= 0) then
                  write(*,*) 'Error during mud3sp solve for system ', system, ' ierror = ', ierror
               endif

               ! Store the computed solution back to SOL_pm
               do kk = 1, nz
                  do jj = 1, ny
                     do ii = 1, nx
                        SOL_pm(system, ii, jj, kk) = phi(ii, jj, kk)
                     end do
                  end do
               end do

               ! call print_stats_rank3_sp(fine_grid, phi, 'Solution for system ')
               print * ,''
            end do
            deallocate(work, phi, rhs)

            call print_stats_rank3(fine_grid, RHS_pm, 'RHS for system')
            call print_stats_rank4(fine_grid, SOL_pm, 'Solution for system')
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
         end subroutine solve_mudpack

         !------------------------------------------------------------------
         ! Subroutine cfx(x,cxx,cx,cex) which provides the
         !     known real coefficients of the x derivative terms for the pde
         !     at any grid point x.  the name chosen in the calling routine
         !     may be different where the coefficient routine must be declared
         !     external.
         subroutine cfx(x, cxx, cx, cex)
            real(dp), intent(in) :: x
            real(dp), intent(out) :: cxx, cx, cex
            ! For the laplacian operator, the coefficients are constant
            cxx = 1.0
            cx = 0.0
            cex = 0.0
         end subroutine cfx

         subroutine cfy(y, cyy, cy, cey)
            real(dp), intent(in) :: y
            real(dp), intent(out) :: cyy, cy, cey
            ! For the laplacian operator, the coefficients are constant
            cyy = 1.0
            cy = 0.0
            cey = 0.0
         end subroutine cfy

         subroutine cfz(z, czz, cz, cez)
            real(dp), intent(in) :: z
            real(dp), intent(out) :: czz, cz, cez
            ! For the laplacian operator, the coefficients are constant
            czz = 1.0
            cz = 0.0
            cez = 0.0
         end subroutine cfz

         subroutine bndc(kbdy, xory, yorz, alfa, phi_val)
            ! Subroutine bndyc(kbdy,xory,yorz,alfa,gbdy).
            !       which are used to input mixed boundary conditions to mud3sp.
            !       the boundaries are numbered one thru six and the form of
            !       conditions are described below.
            !                                          
            integer, intent(in) :: kbdy, xory, yorz, alfa
            real, intent(inout) :: phi_val
            ! For Dirichlet boundary conditions with fixed value 0, nothing is required.
            phi_val = 0.0_dp
         end subroutine bndc
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

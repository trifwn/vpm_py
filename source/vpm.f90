
Module openmpth
   integer                       ::OMPTHREADS
end Module openmpth

Module vpm_lib
   use vpm_vars
   use vpm_size
   ! Constants
   use constants, only: pi, pi2, pi4
   use base_types, only: dp
   ! SPEED 
   use openmpth
   use MPI
   ! Printing
   use io, only: vpm_print, red, blue, green, nocolor, yellow, dummy_string, tab_level, VERBOCITY, &
                 particle_output_file_suffix, vpm_write_folder, pm_output_file_suffix

#ifdef USE_INTEL
      use mkl_service
#endif

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
      module subroutine remesh_particles_3d(iflag, npar_per_cell, XP_out, QP_out, GP_OUT, UP_OUT, NVR_out)
         integer, intent(in)  :: iflag, npar_per_cell
         real(dp), intent(out),allocatable, target, dimension(:,:) :: XP_out, QP_out, GP_OUT, UP_OUT
         integer, intent(out) :: NVR_out
      end subroutine remesh_particles_3d
   end interface remesh_particles_3d

   interface cell3d_interp_euler
      module function cell3d_interp_euler(F, N, M) result(FC)
         implicit none
         real(dp), dimension(8), intent(in) :: F
         integer, intent(in) :: N, M
         real(dp), dimension(N) :: FC
      end function cell3d_interp_euler
   end interface cell3d_interp_euler

   ! vpm_gcalc.f90 contains
   interface calc_velocity_serial_3d
      module subroutine calc_velocity_serial_3d(idcalc)
         integer, intent(in) :: idcalc
      end subroutine calc_velocity_serial_3d
   end interface calc_velocity_serial_3d

   interface diffuse_vort_3d
      module subroutine diffuse_vort_3d
      end subroutine diffuse_vort_3d
   end interface diffuse_vort_3d

   interface calc_antidiffusion
      module subroutine calc_antidiffusion
      end subroutine calc_antidiffusion
   end interface calc_antidiffusion

   ! vpm_mpi.f90 contains
   interface rhsbcast
      module subroutine rhsbcast(RHS_pm, NN, neq)
         integer, intent(in)                                            :: NN(3), neq
         real(dp), dimension(neq, NN(1), NN(2), NN(3)), intent(inout)   :: RHS_pm
      end subroutine rhsbcast
   end interface rhsbcast

   interface rhsscat
      module subroutine rhsscat(BLOCKS, NN_tmp, NNbl, NNbl_bl, NN_bl, nb_i, nb_j, RHS_pm_bl)
         integer, intent(in) ::BLOCKS, NNbl(3, BLOCKS), NNbl_bl(6, BLOCKS), nb_i, nb_j, NN_bl(6), NN_tmp(3)
         real(dp), intent(out) ::RHS_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
      end subroutine rhsscat
   end interface rhsscat

   interface solget
      module subroutine solget(BLOCKS, NBI, NBJ, NN_tmp, NNbl, NNbl_bl, NN_bl, SOL_pm_bl)
         integer, intent(in) ::BLOCKS, NNbl(3, BLOCKS), NNbl_bl(6, BLOCKS), NBI, NBJ, NN_bl(6), NN_tmp(3)
         real(dp), intent(in)  :: SOL_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
      end subroutine solget
   end interface solget

   interface rhsscat_3d
      module subroutine rhsscat_3d(BLOCKS, NN_tmp, NNbl, NNbl_bl, NN_bl, nb_i, nb_j, nb_k, RHS_pm_bl)
         integer, intent(in) ::BLOCKS, NNbl(3, BLOCKS), NNbl_bl(6, BLOCKS), nb_i, nb_j, nb_k, NN_bl(6), NN_tmp(3)
         real(dp), intent(out) ::RHS_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
      end subroutine rhsscat_3d
   end interface rhsscat_3d

   interface solget_3d
      module subroutine solget_3d(BLOCKS, NBI, NBJ, NBK, NN_tmp, NNbl, NNbl_bl, NN_bl, SOL_pm_bl)
         integer, intent(in) :: BLOCKS, NNbl(3, BLOCKS), NNbl_bl(6, BLOCKS), NBI, NBJ, NBK, NN_bl(6), NN_tmp(3)
         real(dp), intent(in)  :: SOL_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3))
      end subroutine solget_3d
   end interface solget_3d

   interface velbcast_3d
      module subroutine velbcast_3d
      end subroutine velbcast_3d
   end interface velbcast_3d

   interface defbcast_3d
      module subroutine defbcast_3d
      end subroutine defbcast_3d
   end interface defbcast_3d

   interface particles_scat
      module subroutine particles_scat
      end subroutine particles_scat
   end interface particles_scat

   interface particles_gath
      module subroutine particles_gath
      end subroutine particles_gath
   end interface particles_gath

   interface proj_gath
      module subroutine proj_gath(NN)
         integer, intent(in) :: NN(3)
      end subroutine proj_gath
   end interface proj_gath

   interface proj_gath_new
      module subroutine proj_gath_new(NN)
         integer, intent(in) :: NN(3)
      end subroutine proj_gath_new
   end interface proj_gath_new

   ! vpm_interpolate.f90 contains
   interface back_to_particles_3D
      module subroutine back_to_particles_3D(SOL_pm, XP, QP, UP, GP, &
                                 velvrx_pm, velvry_pm, velvrz_pm,    &
                                 deformx_pm, deformy_pm, deformz_pm, &
                                 NVR, iproj, itype, NVRM)
         Implicit None
         integer, intent(in)     :: NVR, iproj, NVRM, itype
         real(dp), intent(in)    :: SOL_pm(neqpm, NN(1), NN(2), NN(3))
         real(dp), intent(in)    :: velvrx_pm(NN(1), NN(2), NN(3))
         real(dp), intent(in)    :: velvry_pm(NN(1), NN(2), NN(3))
         real(dp), intent(in)    :: velvrz_pm(NN(1), NN(2), NN(3))
         real(dp), intent(in)    :: deformx_pm(NN(1), NN(2), NN(3))
         real(dp), intent(in)    :: deformy_pm(NN(1), NN(2), NN(3))
         real(dp), intent(in)    :: deformz_pm(NN(1), NN(2), NN(3))
         real(dp), intent(inout) :: QP(neqpm+1, NVRM), XP(3, NVRM), UP(3, NVRM), GP(3, NVRM)
      end subroutine back_to_particles_3D
   end interface

   interface interpolate_particle_Q
      module subroutine interpolate_particle_Q(Q_pm, XP, QP, NVR, iproj)
         Implicit None
         integer, intent(in)     :: NVR, iproj
         real(dp), intent(in)    :: Q_pm(neqpm, NN(1), NN(2), NN(3))
         real(dp), intent(in)    :: XP(3, NVR)
         real(dp), intent(out) :: QP(neqpm +1, NVR) 
      end subroutine interpolate_particle_Q
   end interface

contains

   !> This subroutine performs the VPM (Vortex Particle Method) calculation.
   !! Inputs:
   !!   XP_in: Array (3,NVR) of particle positions.
   !!   QP_in: Array (neq+1,NVR) of particle charges (quantities).
   !!   UP_in: Array (3,NVR) of particle velocities.
   !!   GP_in: Array (3,NVR) of particle deformations. (wmega * \nabla)  u
   !!   NVR_in: Number of particles.
   !!   neqpm_in: Number of equations to be solved.
   !!   WhatToDo: Action to be performed.
   !!   NTIME_in: Current time.
   !!   NI_in: Viscosity for diffusion of vorticity.
   !!   NVRM_in: Maximum number of particles.
   !! Outputs:
   !!   RHS_pm_ptr: Pointer that will be associated with the right hand side of the pm equation:
   !!               \nabla^2 u_i = RHS_i (NEQ)
   !!   velx_ptr: Pointer that will be associated with the x-component of the velocity field at 
   !!               pm grid points.
   !!   vely_ptr: Pointer that will be associated with the y-component of the velocity field at 
   !!               pm grid points.
   !!   velz_ptr: Pointer that will be associated with the z-component of the velocity field at 
   !!               pm grid points.
   !! Notes:
   !!   WhatToDo is an integer that specifies the action to be performed. More specifically:
   !!     0 - Define the problem.
   !!     1 - Solve problem and get velocities on PM from PM solution.
   !!     2 - Solve problem and get velocities and deformation on PM from PM solution.
   !!     3 - Project particles and solve problem.
   !!     4 - Interpolate from PM to particles.
   !!     5 - Diffuse particle vorticity.
   !!   For actions 1 through 4, the subroutine assumes the problem has already been defined.
   !!   The problem is solved in actions 1, 2, and 3. The solution is stored in SOL_pm.
   !!   The subroutine is parallelized using MPI.
   !!   IMPORTANT: OUTPUTS are guaranteed only on rank 0.
   subroutine vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
                  RHS_pm_ptr, velx_ptr, vely_ptr, velz_ptr, NTIME_in, NI_in, NVRM_in,&
                  deformx_ptr, deformy_ptr, deformz_ptr)
      use pmeshpar, only:  ND, SOL_pm, IDVPM
      use parvar, only:    QP, XP, UP, GP, NVR, print_particle_info, print_particle_positions
      use pmgrid, only:    velvrx_pm, velvry_pm, velvrz_pm, RHS_pm,  &
                           deformx_pm, deformy_pm, deformz_pm,       &
                           NXpm_coarse, NYpm_coarse, NZpm_coarse,    &
                           Nblocks, print_RHS_pm
      use pmlib, only:     pmesh
      use projlib, only:   projlibinit, project_particles_3D, project_vol3d
      use yaps, only:      yaps3d
      use openmpth, only:  OMPTHREADS

      Implicit None

      ! ARGUEMENTS
      real(dp), intent(inout), target  :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
      integer, intent(inout)           :: NVR_in
      integer, intent(in)              :: neqpm_in, WhatToDo, NVRM_in, NTIME_in
      real(dp), intent(in)             :: NI_in
      real(dp), intent(out), pointer   :: RHS_pm_ptr(:, :, :, :)
      real(dp), intent(out), pointer   :: velx_ptr(:, :, :), vely_ptr(:, :, :), velz_ptr(:, :, :)
      real(dp), intent(out), pointer, optional   :: deformx_ptr(:, :, :), deformy_ptr(:, :, :), deformz_ptr(:, :, :)

      ! LOCAL VARIABLES
      integer                           :: n_block
      real(dp), allocatable             :: SOL_pm_bl(:, :, :, :), RHS_pm_bl(:, :, :, :)

      integer                           :: ierr, my_rank, np
      integer                           :: i ,j, k 

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
      call print_arguements
      tab_level = 1
#ifdef USE_INTEL
   if (SOLVER .eq. 1) then
      call mkl_set_num_threads(OMPTHREADS)
   end if
#endif

      ! INP
      ND = 3
      NVR_size = NVRM_in
      NTIME_pm = NTIME_in
      neqpm    = neqpm_in
      NI       = NI_in
      n_block = my_rank + 1
      
      if (associated(XP))           nullify (XP)
      if (associated(QP))           nullify (QP)
      if (associated(UP))           nullify (UP)
      if (associated(GP))           nullify (GP)
      if (associated(RHS_pm_ptr))   nullify (RHS_pm_ptr)
      if (associated(velx_ptr))     nullify (velx_ptr)
      if (associated(vely_ptr))     nullify (vely_ptr)
      if (associated(velz_ptr))     nullify (velz_ptr) 
      if (my_rank .eq. 0) then
         NVR = NVR_in
         XP => XP_in
         QP => QP_in; 
         UP => UP_in; 
         GP => GP_in
      end if
      velx_ptr => velvrx_pm
      vely_ptr => velvry_pm
      velz_ptr => velvrz_pm

      if (present(deformx_ptr)) then
         deformx_ptr => deformx_pm
         deformy_ptr => deformy_pm
         deformz_ptr => deformz_pm
      end if

      ! BCAST NVR
      call MPI_BCAST(NVR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (NVR .eq. 0) return
      
      ! -> WhatToDo : 0 - Define the problem
      ! DEFINE SIZES + ALLOCATIONS
      if (WhatToDo .eq. 0) then
         call define_sizes
         call set_pm_velocities_zero
         call set_pm_deformations_zero
         ! Change to reflect if the domain actually changed or not. If the same maybe dont zero velocities?
         call allocate_sol_and_rhs
         return
      end if

      ! Make the necessary allocations
      call allocate_sol_and_rhs

      ! SOL_pm AND RHS_pm ARE 0 BOTH ON THE BL AND PM
      ! PARTICLES TO -> PM MEANING FROM Qs We get RHS_pm
      call project_particles_parallel 

      !  -> WhatToDo: 4 - interpolate from PM to particles
      if (WhatToDo .eq. 4) then
         if (my_rank .eq. 0) then
            call calc_velocity_serial_3d(-1) 
            ! ONLY DELATATION THETA STO PM
            ! WE LOST SOL_pm
            ! SOL_pm IS NOW FILLED WITH THETA (DEFORMATION)

            ! MAYBE OUTSIDE LOOP
            if (IPMWRITE .GT. 0) then
               do i = 1, IPMWRITE
                  if (NTIME_pm .ge. IPMWSTART(i) &
                  .and. NTIME_pm .le. (IPMWSTART(i) + IPMWSTEPS(i)))  then
                     call write_pm_solution(NTIME_pm, NN, NN_bl, RHS_pm, SOL_pm, &
                                            velvrx_pm, velvry_pm, velvrz_pm, &
                                            deformx_pm, deformy_pm, deformz_pm)
                  end if
               end do
            end if
         end if

         ! WHEN ITYPEB = 1 WE GET THE UP AND GP From the velocity
         call interpolate_particles_parallel(1)
         return
      end if

      !  -> WhatToDo : 5 - diffuse
      if (WhatToDo .eq. 5) then
         call set_pm_velocities_zero
         call set_pm_deformations_zero

         if (my_rank .eq. 0) then
            !diffusion stores -NI*grad^2 w * Vol in GP(1,:)
            ! SOL_pm = -VIS \nabla \cdot RHS_pm
            call diffuse_vort_3d ! DIFFUSION OF VORTICITY
         end if
         ! WHEN ITYPEB = 2 WE GET THE GP FROM THE SOL_pm (DEFORMATION) and from QP
         call interpolate_particles_parallel(2)
         ! deallocate (SOL_pm, RHS_pm)
         ! deallocate (SOL_pm_bl, RHS_pm_bl)
         return
      end if

      !------------------------------ FOR WHATTODO = 1, 2, 3 --------------------------------
      call solve_problem
      velvrx_pm = 0.d0; velvry_pm = 0.d0; velvrz_pm = 0.d0
      deformx_pm = 0.d0; deformy_pm = 0.d0; deformz_pm = 0.d0
      
      !  -> WhatToDo :  1 - GET VELOCITIES ON PM FROM PM SOLUTION
      if (WhatToDo .eq. 1) then
         if (my_rank .eq. 0) then
            st = MPI_WTIME()
            write (*, *) ""
            write (dummy_string, "(A)") 'Calculating Velocities on PM using FD'
            call vpm_print(dummy_string, blue, 1)

            !call convect_first_order(Xbound,Dpm,NN,NN_bl)
            ! FROM THE SOLUTION OF PM WE GET THE VELOCITIES ON THE GRID
            call calc_velocity_serial_3d(0) ! VELOCITY STO PM
            ! SOL_pm is stil the solution of vorticity
            et = MPI_WTIME()
            write (dummy_string, "(A,I5,A,F8.2,A)") &
                  achar(9)//'finished in:', int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
            call vpm_print(dummy_string, yellow, 1)

            call print_velocity_stats
         end if
      end if

      !  -> WhatToDo : 2 - GET VELOCITIES AND DEFORMATION ON PM FROM PM SOLUTION
      if (WhatToDo .eq. 2) then
         if (my_rank .eq. 0) then
            st = MPI_WTIME()
            write (*, *) ""
            write (dummy_string, "(A)") 'Calculating Velocities and Deformation on PM using FD'
            call vpm_print(dummy_string, blue, 1)

            !call convect_first_order(Xbound,Dpm,NN,NN_bl)
            call calc_velocity_serial_3d(1) 
            ! VELOCITY AND THETA STO PM
            ! SOL_pm IS NOW FILLED WITH THETA (DEFORMATION)
            et = MPI_WTIME()
            
            write (dummy_string, "(A,I5,A,F8.2,A)") &
                  achar(9)//'finished in:', int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
            call vpm_print(dummy_string, yellow, 1)

            call print_velocity_stats
            ! call calc_antidiffusion
         end if
         call interpolate_particles_parallel(1) ! INTERPOLATION FROM PM TO PARTICLES
      end if
      ! SAVE THE SOLUTION
      if (my_rank.eq.0) then
         if (IPMWRITE .GT. 0) then
            do i = 1, IPMWRITE
               if ((NTIME_pm .ge. IPMWSTART(i)) .and. (NTIME_pm .le. (IPMWSTART(i) + IPMWSTEPS(i)))) then
                  st = MPI_WTIME()
                  write (dummy_string, "(A)") 'Writing particles and PM solution'
                  call vpm_print(dummy_string, blue, 1)
                  tab_level = tab_level + 1

                  call write_pm_solution(NTIME_pm, NN, NN_bl, RHS_pm, SOL_pm,        &
                                         velvrx_pm, velvry_pm, velvrz_pm,     &
                                         deformx_pm, deformy_pm, deformz_pm)
                  call write_particles(NTIME_pm, XP, UP, QP, NVR)
                  
                  et = MPI_WTIME()
                  write (dummy_string, "(A,I3,A,F8.3,A)") &
                     'finished in:', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
                  call vpm_print(dummy_string, yellow, 1)
                  tab_level = tab_level - 1
               end if
            end do
         end if
      end if
   return      

   contains
      subroutine solve_problem
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
         !call diffuse_vort_3d
         
         ! SCATTER PROBLEM TO ALL PROCESSORS
         call rhsbcast(RHS_pm, NN, neqpm) ! RHS PM -> TO ALL PROCESSORS -> RHS_pm_BL
         IF (SOLVER .ne. 0) then
            n_block = my_rank + 1
            NN_tmp(1:3) = NNbl(1:3, n_block)
            call rhsscat_3d(BLOCKS, NN_tmp, NNbl, NNbl_bl, NN_bl, nb_i, nb_j, nb_k, RHS_pm_bl) !
         END IF
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
         tab_level = tab_level - 1
         
         ! CLEAR VELOCITIES
         call set_pm_velocities_zero
         call set_pm_deformations_zero
      end subroutine solve_problem

      subroutine pmesh_solve !
         integer :: rank
         !Yaps or Serial Pmesh

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
               Nblocks = 1
               call pmesh(SOL_pm, RHS_pm, QP, XP, Xbound, DPm, NN, NN_bl, ND, Nblocks, ibctyp, 1, neqpm, &
                          iynbc, NVR, itree, ilevmax)
               ! call calc_velocity_serial_3d(1)
            END IF
            !--------------------------------------------
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            ! call velbcast_3d
         ELSE
            if (my_rank .eq. 0) then 
               write (dummy_string, "(A)") 'Solving PM with YAPS'
               call vpm_print(dummy_string, blue, 1)
               st = MPI_WTIME()
            endif

            iret = 0
            call yaps3d(SOL_pm_bl, RHS_pm_bl, Xbound_bl, Xbound_coarse, Dpm, Dpm_coarse, NNbl, NNbl_bl, &
                        NN_coarse, NN_bl_coarse, ND, BLOCKS, ibctyp, 1, neqpm, ncoarse, NBI, NBJ, NBK, nb_i, nb_j, nb_k, &
                        iret, iyntree, ilevmax, neqpm)

            ! GATHER SOLUTION
            n_block = my_rank + 1
            NN_tmp(1:3) = NNbl(1:3, n_block)
            call solget_3d(BLOCKS, NBI, NBJ, NBK, NN_tmp, NNbl, NNbl_bl, NN_bl, SOL_pm_bl) 
            !if (my_rank.eq.0) call calc_velocity_serial_3d(1)
            ! call velbcast_3d

            if (my_rank .eq. 0) then
               write (dummy_string, "(A)") 'Final PM solution values'
               call vpm_print(dummy_string, blue, 2)
            END IF
            do rank = 0, np-1
               if (my_rank .eq. rank) then
                  write (dummy_string, "(A, I3,A,F8.3,A,F8.3)") &
                     'np=', my_rank, &
                     achar(9)//'min(SOL_pm)', maxval(abs(SOL_pm(neqpm, :, :, :))), &
                     achar(9)//'max(SOL_pm)', maxval(abs(SOL_pm(neqpm, :, :, :)))   
                  call vpm_print(dummy_string, blue, 2)
               end if
               call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            end do
         END IF
         !--------------------------------------------

      end subroutine pmesh_solve

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
            write (*, *) achar(9), achar(9), 'NVRM_in = ', NVRM_in, achar(27)//'[0m'
         end if
      end subroutine print_arguements
   
      subroutine allocate_sol_and_rhs
         ! NN_tmp       : is the number of cells in each direction
         ! NN_bl_tmp    : is the start and finish of the cells in each direction
         ! Xbound   : is the boundary of the domain
         NN_tmp(1:3) = NNbl(1:3, n_block)

         ! SOL_pm block
         if (allocated(SOL_pm_bl)) then
            deallocate (SOL_pm_bl)
         end if
         allocate (SOL_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3)))
         SOL_pm_bl = 0.d0

         ! RHS_pm block
         if (allocated(RHS_pm_bl)) then
            deallocate (RHS_pm_bl)
         end if 
         allocate (RHS_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3)))
         RHS_pm_bl = 0.d0

         ! RHS_pm
         if (allocated(RHS_pm)) then
            deallocate (RHS_pm)
         end if
         allocate (RHS_pm(neqpm, NN(1), NN(2), NN(3))) 
         RHS_pm = 0.d0
         RHS_pm_ptr => RHS_pm
         
         ! SOL_pm
         if (allocated(SOL_pm)) then
            deallocate (SOL_pm)
         end if
         allocate (SOL_pm(neqpm, NN(1), NN(2), NN(3)))
         SOL_pm = 0.d0
      end subroutine allocate_sol_and_rhs
   end subroutine vpm

   subroutine set_pm_velocities_zero
      use pmgrid, only: velvrx_pm, velvry_pm, velvrz_pm, NXpm_coarse, NYpm_coarse, NZpm_coarse
      if (allocated(velvrx_pm)) then
         deallocate (velvrx_pm, velvry_pm, velvrz_pm)
         allocate (velvrx_pm(NXpm_coarse, NYpm_coarse, NZpm_coarse))
         allocate(velvry_pm(NXpm_coarse, NYpm_coarse, NZpm_coarse))
         allocate(velvrz_pm(NXpm_coarse, NYpm_coarse, NZpm_coarse))
         velvrx_pm = 0.d0; 
         velvry_pm = 0.d0; 
         velvrz_pm = 0.d0
      else
         allocate (velvrx_pm(NXpm_coarse, NYpm_coarse, NZpm_coarse))
         allocate(velvry_pm(NXpm_coarse, NYpm_coarse, NZpm_coarse))
         allocate(velvrz_pm(NXpm_coarse, NYpm_coarse, NZpm_coarse))
         velvrx_pm = 0.d0; 
         velvry_pm = 0.d0; 
         velvrz_pm = 0.d0
      end if
   end subroutine set_pm_velocities_zero

   subroutine set_pm_deformations_zero
      use pmgrid, only: deformx_pm, deformy_pm, deformz_pm, NXpm_coarse, NYpm_coarse, NZpm_coarse
      if (allocated(deformx_pm)) then
         deallocate (deformx_pm, deformy_pm, deformz_pm)
         allocate (deformx_pm(NXpm_coarse, NYpm_coarse, NZpm_coarse))
         allocate(deformy_pm(NXpm_coarse, NYpm_coarse, NZpm_coarse))
         allocate(deformz_pm(NXpm_coarse, NYpm_coarse, NZpm_coarse))
         deformx_pm = 0.d0; 
         deformy_pm = 0.d0; 
         deformz_pm = 0.d0
      else
         allocate (deformx_pm(NXpm_coarse, NYpm_coarse, NZpm_coarse))
         allocate(deformy_pm(NXpm_coarse, NYpm_coarse, NZpm_coarse))
         allocate(deformz_pm(NXpm_coarse, NYpm_coarse, NZpm_coarse))
         deformx_pm = 0.d0; 
         deformy_pm = 0.d0; 
         deformz_pm = 0.d0
      end if
   end subroutine set_pm_deformations_zero

   !> Defines Coarse and Fine GRID
   subroutine define_sizes
      use pmeshpar, only:  ND, ndumcell
      use vpm_vars, only:  interf_iproj
      use vpm_size, only:  NN_bl, NN
      use pmgrid, only:    XMIN_pm, YMIN_pm, ZMIN_pm, &
                           XMAX_pm, YMAX_pm, ZMAX_pm, &
                           NXs_coarse_bl, NXf_coarse_bl, NYs_coarse_bl, &
                           NYf_coarse_bl, NZs_coarse_bl, NZf_coarse_bl, &
                           NXpm_coarse, NYpm_coarse, NZpm_coarse, DVPM, &
                           DXpm, DYpm, DZpm
      use pmlib, only: definepm
      use io, only: dummy_string, vpm_print, nocolor, blue, yellow
      use MPI

      Implicit None
      integer    :: nsiz(3), nsiz_bl(3)
      integer    :: i, j, k, np, my_rank, ierr, nb, istep, lev
      real(dp)   :: Xbound_tmp(6)
      integer    :: redifine_pm
      integer    :: NN_bl_tmp(6), NXbl, NYbl, NZbl 

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

      if(my_rank.eq.0) then
         st = MPI_WTIME()
         write (dummy_string, *) 'Defining Sizes'
         call vpm_print(dummy_string, blue, 2)
      end if
      tab_level = tab_level + 1

      if ((NTIME_pm .eq. 0) .or. idefine .eq. 0) then
         redifine_pm = 1
      else
         redifine_pm = 0
      end if


      if (my_rank .eq. 0) then
         if (redifine_pm.eq.1) then
            call get_domain_bounds_from_particles
            write (dummy_string, '(A)') 'The computational domain bounds are recalculated from the particle positions'
            call vpm_print(dummy_string,nocolor, 2)
         end if

         write (dummy_string, '(A)') 'The computational domain bounds are:'
         call vpm_print(dummy_string,nocolor, 2)
         write (dummy_string, '(A, F8.5, A, F8.5, A, F8.5)') &
                  achar(9)//'XMIN=', XMIN_pm, &
                  achar(9)//'YMIN=', YMIN_pm, &
                  achar(9)//'ZMIN=', ZMIN_pm
         call vpm_print(dummy_string,nocolor, 2)
         write (dummy_string, '(A, F8.5, A, F8.5, A, F8.5)') &
                  achar(9)//'XMAX=', XMAX_pm, &
                  achar(9)//'YMAX=', YMAX_pm, &
                  achar(9)//'ZMAX=', ZMAX_pm
         call vpm_print(dummy_string,nocolor, 2)
      end if

      ! BRADCAST THE MIN AND MAX OF THE COMPUTATIONAL DOMAIN
      call MPI_BCAST(XMIN_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(YMIN_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ZMIN_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call MPI_BCAST(XMAX_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(YMAX_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ZMAX_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      
      ! CREATE THE BOUNDING BOX
      Xbound(1) = XMIN_pm
      Xbound(2) = YMIN_pm
      Xbound(3) = ZMIN_pm
      Xbound(4) = XMAX_pm 
      Xbound(5) = YMAX_pm 
      Xbound(6) = ZMAX_pm
      
      Dpm(1) = DXpm
      Dpm(2) = DYpm
      Dpm(3) = DZpm
      
      nsiz(1) = NBI*ncoarse
      nsiz(2) = NBJ*ncoarse
      nsiz(3) = NBK*ncoarse

      if(my_rank.eq.0) then
         write (dummy_string, '(A)') 'The extended domain is defined'
         call vpm_print(dummy_string,nocolor, 2)
         write (dummy_string, '(A, I5, A, I1)') achar(9)//'NTIME_PM=', NTIME_pm, &
         achar(9)//'Redefine=', idefine
         call vpm_print(dummy_string,nocolor, 2)
         write (dummy_string, '(A)') achar(9)//'The number of cells (nodes-1) must be divisible by the processor subdivision:'
         call vpm_print(dummy_string,nocolor, 2)
         write (dummy_string, '(A, I5, A, I5, A, I5)') &
                  achar(9)//achar(9)//'X-dir: ', nsiz(1), &
                  achar(9)//achar(9)//'Y-dir: ', nsiz(2), &
                  achar(9)//achar(9)//'Z-dir: ', nsiz(3)
         call vpm_print(dummy_string,nocolor, 2)
      endif
      
      ! First Change Dpm so that the numbers of cells divides
      ! by nsize i.e with NBI,NBJ,NBK, ncoarse,levmax depending on the criterion
      ! thats why ndumcell=0
      ndumcell = 0
      if (redifine_pm .eq. 1) then
         call definepm(3, Xbound, Dpm, ND, ndumcell, nsiz, NN, NN_bl)
      endif

      ! THE new extended domain is defined
      XMIN_pm = Xbound(1)
      YMIN_pm = Xbound(2)
      ZMIN_pm = Xbound(3)
      XMAX_pm = Xbound(4)
      YMAX_PM = Xbound(5)
      ZMAX_pm = Xbound(6)

      NXpm_coarse = NN(1)
      NYpm_coarse = NN(2)
      NZpm_coarse = NN(3)

      NXs_coarse_bl = NN_bl(1)
      NYs_coarse_bl = NN_bl(2)
      NZs_coarse_bl = NN_bl(3)
      NXf_coarse_bl = NN_bl(4)
      NYf_coarse_bl = NN_bl(5)
      NZf_coarse_bl = NN_bl(6)

      DXpm = Dpm(1)
      DYpm = Dpm(2)
      DZpm = Dpm(3)

      DVpm = DXpm*DYpm
      if (ND .eq. 3) then
         DVpm = DVpm*DZpm
      end if

      if (my_rank .eq. 0) then
         write (dummy_string, '(A, F8.5, A, F8.5, A, F8.5)') &
                  achar(9)//'XMIN='//achar(9), XMIN_pm, &
                  achar(9)//'YMIN='//achar(9), YMIN_pm, &
                  achar(9)//'ZMIN='//achar(9), ZMIN_pm
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, F8.5, A, F8.5, A, F8.5)') &
                  achar(9)//'XMAX='//achar(9), XMAX_pm, &
                  achar(9)//'YMAX='//achar(9), YMAX_pm, &
                  achar(9)//'ZMAX='//achar(9), ZMAX_pm
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, F8.5, A, F8.5, A, F8.5)') &
                  achar(9)//'DX='//achar(9), DXpm, &
                  achar(9)//'DY='//achar(9), DYpm, &
                  achar(9)//'DZ='//achar(9), DZpm
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A, I5, A, I5, A, I5)') &
                  achar(9)//'Nodes X= ', NXpm_coarse, &
                  achar(9)//achar(9)//'Nodes Y= ', NYpm_coarse, &
                  achar(9)//achar(9)//'Nodes Z= ', NZpm_coarse
         call vpm_print(dummy_string,nocolor,2)

         write (dummy_string, '(A, F8.5)') achar(9) // 'DV=', DVpm
         call vpm_print(dummy_string,nocolor,2)
         write (dummy_string, '(A)') 'The extended coarse domain is defined everytime'
         call vpm_print(dummy_string,nocolor,2)
      end if

      ! define block grids
      ! so they are divided by ncoarse and ilevmax
      ! so as to have the coarse information at the boundaries exactly.
      BLOCKS = np
      if (.not. allocated(XBound_bl)) then
         allocate (Xbound_bl(6, BLOCKS), NNbl(3, BLOCKS), NNbl_bl(6, BLOCKS))
      end if
      
      ! DEVIDE THE DOMAIN INTO EQUAL BLOCKS 
      NXB = int(nint(((Xbound(4) - Xbound(1))/Dpm(1)))) ! NUMBER OF CELLS IN A BLOCK IN X DIRECTION 
      NYB = int(nint(((Xbound(5) - Xbound(2))/Dpm(2)))) ! NUMBER OF CELLS IN A BLOCK IN Y DIRECTION
      NZB = int(nint(((Xbound(6) - Xbound(3))/Dpm(3)))) ! NUMBER OF CELLS IN A BLOCK IN Z DIRECTION
      NXbl = NXB/NBI
      NYbl = NYB/NBJ
      NZbl = NZB/NBK
      
      ndumcell_bl = ncoarse
      nsiz_bl = ncoarse!ndumcell_coarse!*2*2**ilevmax

      if (my_rank.eq.0) then
         write(dummy_string, '(A)') achar(9) // 'The number of cells in each direction for the block grid'
         call vpm_print(dummy_string,nocolor,2)
         write(dummy_string, '(A, I5, A, I5, A, I5)') achar(9) // &
                  achar(9)//'NXB='//achar(9)//achar(9), NXB,      &
                  achar(9)//'NYB='//achar(9)//achar(9), NYB,      &
                  achar(9)//'NZB='//achar(9)//achar(9), NZB       
         call vpm_print(dummy_string,nocolor,2)                          
         write(dummy_string, '(A, I5, A, I5, A, I5)') achar(9) //  &
                  achar(9)//'NXbl='//achar(9)//achar(9), NXbl,    &
                  achar(9)//'NYbl='//achar(9)//achar(9), NYbl,    &
                  achar(9)//'NZbl='//achar(9)//achar(9), NZbl
         call vpm_print(dummy_string,nocolor,2)
         write(dummy_string, '(A)') achar(9)//achar(9)//'After calling definepm on each block '
         call vpm_print(dummy_string,nocolor,2)
      end if

      nb_i = -1
      nb_j = -1
      nb_k = -1
      do k = 1, NBK
         do j = 1, NBJ
            do i = 1, NBI
               nb = (k - 1)*NBJ*NBI + (j - 1)*NBI + i
               Xbound_bl(1, nb) = Xbound(1) + (i - 1)*(NXbl)*Dpm(1)
               Xbound_bl(4, nb) = Xbound(1) + (i)*(NXbl)*Dpm(1)
               
               Xbound_bl(2, nb) = Xbound(2) + (j - 1)*(NYbl)*Dpm(2)
               Xbound_bl(5, nb) = Xbound(2) + (j)*(NYbl)*Dpm(2)
               
               Xbound_bl(3, nb) = Xbound(3) + (k - 1)*(NZbl)*Dpm(3)
               Xbound_bl(6, nb) = Xbound(3) + (k)*(NZbl)*Dpm(3)
               Xbound_tmp(1:6) = Xbound_bl(1:6, nb)
               call definepm(1, Xbound_tmp, Dpm, ND, ndumcell_bl, nsiz_bl, NN_tmp, NN_bl_tmp)
               
               Xbound_bl(1:6, nb) = Xbound_tmp(1:6)
               NNbl(1:3, nb) = NN_tmp(1:3)         ! KOMVOI
               NNbl_bl(1:6, nb) = NN_bl_tmp(1:6)   ! KOMVOI START STOP in the local system that do not include the dummy cells
               
               if (my_rank .eq. 0) then
                  write (dummy_string, '(A,I3)') achar(9)//'Block ', nb
                  call vpm_print(dummy_string,blue,2)
                  write (dummy_string, '(A, F8.5, A, F8.5, A, F8.5)') achar(9)//&
                           achar(9)//'XMIN= '//achar(9), Xbound_bl(1, nb), &
                           achar(9)//'YMIN= '//achar(9), Xbound_bl(2, nb), &
                           achar(9)//'ZMIN= '//achar(9), Xbound_bl(3, nb)
                  call vpm_print(dummy_string,nocolor,2)
                  write (dummy_string, '(A, F8.5, A, F8.5, A, F8.5)') achar(9)// &
                           achar(9)//'XMAX= '//achar(9), Xbound_bl(4, nb), &
                           achar(9)//'YMAX= '//achar(9), Xbound_bl(5, nb), &
                           achar(9)//'ZMAX= '//achar(9), Xbound_bl(6, nb)
                  call vpm_print(dummy_string,nocolor,2)
                  write(dummy_string, '(A, I5)') achar(9) // achar(9)//&
                                       'X cells in block= '//achar(9)//achar(9), NNbl(1, nb)
                  call vpm_print(dummy_string,nocolor,2)
                  write(dummy_string, '(A, I5)') achar(9) // achar(9)//&
                                       'Y cells in block= '//achar(9)//achar(9), NNbl(2, nb)
                  call vpm_print(dummy_string,nocolor,2)
                  write(dummy_string, '(A, I5)') achar(9) // achar(9)//&
                                       'Z cells in block= '//achar(9)//achar(9), NNbl(3, nb)
                  call vpm_print(dummy_string,nocolor,2)
               end if

               if (nb .eq. my_rank + 1) then
                  nb_i = i
                  nb_j = j
                  nb_k = k
               end if
            end do
         end do
      end do
      if ((nb_i .eq. -1) .or. (nb_j .eq. -1) .or. (nb_k .eq. -1)) then
         print *, "Processor with rank ", my_rank, " does not have a block"
         stop
      end if

      !define coarse grid must cover block grids
      Xbound(1) = XMIN_pm!minval(Xbound_bl(1,:))
      Xbound(2) = YMIN_pm!minval(Xbound_bl(2,:))
      Xbound(3) = ZMIN_pm!minval(Xbound_bl(3,:))
      Xbound(4) = XMAX_pm!maxval(Xbound_bl(4,:))
      Xbound(5) = YMAX_pm!maxval(Xbound_bl(5,:))
      Xbound(6) = ZMAX_pm!maxval(Xbound_bl(6,:))
      Xbound_coarse = Xbound
      Dpm_coarse = ncoarse*Dpm
      ndumcell_coarse = 4!2**ilevmax
      nsiz_bl = 2**ilevmax
      call definepm(1, Xbound_coarse, Dpm_coarse, ND, ndumcell_coarse, nsiz_bl, NN_coarse, NN_bl_coarse)

      !add to dummy cells to the grid globally used for remeshing purposes mainly
      !ndumcell=4
      !Xbound(1)=XMIN_pm;Xbound(2)=YMIN_pm;Xbound(3)=ZMIN_pm
      !Xbound(4)=XMAX_pm;Xbound(5)=YMAX_pm;Xbound(6)=ZMAX_pm
      !Dpm   (1)=DXpm   ;Dpm(2)   =DYpm   ;Dpm(3)   =DZpm
      !call  definepm(1,Xbound,Dpm,ND,ndumcell,nsiz,NN,NN_bl)
      !XMIN_pm=Xbound(1);YMIN_pm=Xbound(2);ZMIN_pm=Xbound(3)
      !XMAX_pm=Xbound(4);YMAX_PM=Xbound(5);ZMAX_pm=Xbound(6)
      !NXpm=NN(1);NYpm=NN(2);NZpm=NN(3)
      !NXs_bl=NN_bl(1);NYs_bl=NN_bl(2);NZs_bl=NN_bl(3)
      !NXf_bl=NN_bl(4);NYf_bl=NN_bl(5);NZf_bl=NN_bl(6)
      !print *,'final mesh',NN
      !-----
      
      ! CALCULATE THE NUMBER OF LEVELS
      if (my_rank .eq. 0) then
         do lev = 0, ilevmax
            istep = 2**lev
            !!!!!!!!if not divided exactly dummy cell
            if (ND .eq. 2) then
               if (int((NNbl_bl(4, 1) - NNbl_bl(1, 1))/istep) .eq. 0 .or. &
               int((NNbl_bl(5, 1) - NNbl_bl(2, 1))/istep) .eq. 0) then
                  ilevmax = lev - 1
                  print *, 'Changing number of levels', ilevmax
                  exit
               end if
            else
               if (  int((NNbl_bl(4, 1) - NNbl_bl(1, 1))/istep) .eq. 0 .or. &
                     int((NNbl_bl(5, 1) - NNbl_bl(2, 1))/istep) .eq. 0 .or. &
                     int((NNbl_bl(6, 1) - NNbl_bl(3, 1))/istep) .eq. 0) then
                  ilevmax = lev - 1
                  print *, 'Changing number of levels', ilevmax
                  exit
               end if
            end if
         end do
      end if
      call MPI_BCAST(ilevmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      tab_level = tab_level - 1
      if (my_rank .eq. 0) then
         et = MPI_WTIME()
         write (dummy_string, *) 'Defining Sizes finished in:', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
         call vpm_print(dummy_string, yellow, 2)
      end if      
   end subroutine define_sizes

   subroutine get_domain_bounds_from_particles
      use parvar, only: NVR, XP
      use vpm_vars, only: interf_iproj
      use pmgrid, only: XMIN_pm, YMIN_pm, ZMIN_pm, DXpm, DYpm, DZpm, XMAX_pm, &
                        YMAX_pm, ZMAX_pm
      real(dp) :: XMIN_pm_old, YMIN_pm_old, ZMIN_pm_old, XMAX_pm_old, YMAX_pm_old, ZMAX_pm_old
      logical  :: bounds_changed

      XMIN_pm_old = XMIN_pm
      YMIN_pm_old = YMIN_pm
      ZMIN_pm_old = ZMIN_pm
      XMAX_pm_old = XMAX_pm
      YMAX_pm_old = YMAX_pm
      ZMAX_pm_old = ZMAX_pm

      XMIN_pm = minval(XP(1, 1:NVR)) - interf_iproj*DXpm
      YMIN_pm = minval(XP(2, 1:NVR)) - interf_iproj*DYpm
      ZMIN_pm = minval(XP(3, 1:NVR)) - interf_iproj*DZpm

      XMAX_pm = maxval(XP(1, 1:NVR)) + interf_iproj*DXpm
      YMAX_pm = maxval(XP(2, 1:NVR)) + interf_iproj*DYpm
      ZMAX_pm = maxval(XP(3, 1:NVR)) + interf_iproj*DZpm
      
      ! Check if the domain bounds have changed
      if (XMIN_pm .ne. XMIN_pm_old .or. YMIN_pm .ne. YMIN_pm_old .or. ZMIN_pm .ne. ZMIN_pm_old .or. &
          XMAX_pm .ne. XMAX_pm_old .or. YMAX_pm .ne. YMAX_pm_old .or. ZMAX_pm .ne. ZMAX_pm_old) then
         bounds_changed =  .TRUE.
      else
         bounds_changed = .FALSE.
      end if
   end subroutine get_domain_bounds_from_particles

   subroutine convect_first_order
      use vpm_vars, only: DT_c
      use parvar, only: NVR, XP, UP, GP, QP
      integer                      :: i

      do i = 1, NVR
         XP(1:3, i) = XP(1:3, i) + UP(1:3, i)*DT_c
         QP(1:neqpm, i) = QP(1:neqpm, i) + GP(1:neqpm, i)*DT_c
      end do
   end subroutine convect_first_order

   subroutine project_particles_parallel
      use pmgrid, only:    RHS_pm
      use pmeshpar, only:  SOL_pm, IDVPM, ND
      use parvar, only:    NVR
      use projlib, only:   projlibinit, project_particles_3D, project_vol3d
      
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
      if (my_rank .eq. 0) NVR_p = NVR_p + mod(NVR, np)

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

      call projlibinit(Xbound, Dpm, NN, NN_bl, IDVPM, ND)

      ! PROJECT PARTICLES TO PM
      st = MPI_WTIME()
      if (my_rank .eq. 0) then
         write (dummy_string, "(A)") "Projecting Particles on PM (RHS_pm)"
         call vpm_print(dummy_string, blue, 1) 
      endif
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call project_particles_3D(RHS_pm, QP_scatt, XP_scatt, NVR_projtype_scatt, NVR_p, neqpm, ieq, neqpm, QINF, NVR_p)
      
      if (my_rank .eq. 0) then
         write (dummy_string, "(A)") "Gathering RHS from all processors"
         call vpm_print(dummy_string, blue, 1) 
      endif

      call proj_gath(NN) ! RHS IS NOW FILLED

      tab_level = tab_level - 1
      if (my_rank .eq. 0) then
         write (dummy_string, "(A)") achar(9)//"Normalizing RHS by volume"
         call vpm_print(dummy_string, blue, 1)

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
      use parvar, only:    NVR
      use pmgrid, only:    velvrx_pm, velvry_pm, velvrz_pm, RHS_pm, &
                           deformx_pm, deformy_pm, deformz_pm
      use pmeshpar, only:  SOL_pm

      integer, intent(in)               :: itypeb
      integer, allocatable              :: ieq(:)
      real(dp), allocatable             :: QINF(:)
      integer                           :: ierr, my_rank, np, i

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_Size(MPI_COMM_WORLD, np, ierr)

      if (my_rank .eq. 0) st = MPI_WTIME()

      ! BROADCASTING
      call rhsbcast(RHS_pm, NN, neqpm)
      call rhsbcast(SOL_pm, NN, neqpm)
      if (itypeb .eq. 1) call velbcast_3d
      if (itypeb .eq. 2) call defbcast_3d
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
      ! WHEN ITYPEB = 2 WE GET THE GP FROM THE deformation(stored in SOL_pm)
      call back_to_particles_3D(SOL_pm, XP_scatt, QP_scatt, UP_scatt, GP_scatt,  &
                                velvrx_pm, velvry_pm, velvrz_pm,                 &
                                deformx_pm, deformy_pm, deformz_pm,              &
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

   subroutine print_velocity_stats
      use pmgrid, only: NXpm_coarse, NYpm_coarse, NZpm_coarse,&
                        velvrx_pm, velvry_pm, velvrz_pm 
      integer :: i, j, k 
      if (VERBOCITY >= 2) then 
         write (dummy_string, "(A)") 'Velocity info in the PM:'
         call vpm_print(dummy_string, red, 1)
         write (dummy_string, "(A,A, F8.3, A, F8.3, A, F8.3)") &
            achar(9)//"X:", &
            achar(9)//"Max:", maxval(velvrx_pm), &
            achar(9)//"Min:", minval(velvrx_pm), &
            achar(9)//"Mean:", sum(velvrx_pm)/NXpm_coarse/NYpm_coarse/NZpm_coarse
         call vpm_print(dummy_string, nocolor, 2)
         write (dummy_string, "(A,A, F8.3, A, F8.3, A, F8.3)") &
            achar(9)//"Y:", &
            achar(9)//"Max:", maxval(velvry_pm), &
            achar(9)//"Min:", minval(velvry_pm), &
            achar(9)//"Mean:", sum(velvry_pm)/NXpm_coarse/NYpm_coarse/NZpm_coarse
         call vpm_print(dummy_string, nocolor, 2)
         write (dummy_string, "(A,A, F8.3, A, F8.3, A, F8.3)") &
            achar(9)//"Z:", &
            achar(9)//"Max:", maxval(velvrz_pm), &
            achar(9)//"Min:", minval(velvrz_pm), &
            achar(9)//"Mean:", sum(velvrz_pm)/NXpm_coarse/NYpm_coarse/NZpm_coarse
         call vpm_print(dummy_string, nocolor, 2)
         ! CHeck for nan values
         if (any(isnan(velvrx_pm)) .or. any(isnan(velvry_pm)) .or. any(isnan(velvrz_pm))) then
            ! Print Velx
            do i = 1, NXpm_coarse
               do j = 1, NYpm_coarse
                  do k = 1, NZpm_coarse
                     if (isnan(velvrx_pm(i, j, k)) .or. isnan(velvry_pm(i, j, k)) .or. isnan(velvrz_pm(i, j, k))) then
                        write (*, "(A, I3, A, I3, A, I3)") &
                              achar(9)//"I:", i, achar(9)//"J:", j, achar(9)//"K:", k, achar(9)
                        ! write (*, "(A, 3F15.8)") &
                              ! achar(9)//achar(9)//"Velx:" , velvrx_pm(i, j, k), achar(9)//"Vely", velvry_pm(i, j, k), achar(9)//"Velz", velvrz_pm(i, j, k)
                        stop
                     end if
                  end do
               end do
            end do
            write (*, *) achar(9), 'VPM: NAN VALUES IN VELOCITY'
            stop
         end if
      endif
   end subroutine print_velocity_stats

   subroutine write_pm_solution(NTIME, NN_in, NNbl_in, RHS, SOL, velx, vely, velz, deformx_pm, deformy_pm, deformz_pm)
      use pmgrid, only: XMIN_pm, YMIN_pm, ZMIN_pm, DXpm, DYpm, DZpm
      
      integer, intent(in)              :: NTIME
      integer, intent(in)              :: NN_in(3), NNbl_in(6)
      real(dp), intent(in)             :: RHS(neqpm, NN_in(1), NN_in(2), NN_in(3)), &
                                          SOL(neqpm, NN_in(1), NN_in(2), NN_in(3))
      real(dp), intent(in)             :: velx(NN_in(1), NN_in(2), NN_in(3)),       &
                                          vely(NN_in(1), NN_in(2), NN_in(3)),       &
                                          velz(NN_in(1), NN_in(2), NN_in(3))
      real(dp), intent(in), optional   :: deformx_pm(NN_in(1), NN_in(2), NN_in(3)), &
                                          deformy_pm(NN_in(1), NN_in(2), NN_in(3)), &
                                          deformz_pm(NN_in(1), NN_in(2), NN_in(3))
      integer                          :: NXs, NYs, NZs, NXf, NYf, NZf

      character*50                     :: filout
      integer                          :: i, j, k, size_eq
      logical                          :: exist_flag
      real(dp)                         :: cellx, celly, cellz, velocx, velocy, velocz, &
                                          wmegax, wmegay, wmegaz, psi_1, psi_2, psi_3, &
                                          deformx, deformy, deformz

      
      write (filout, '(a,i5.5,a)') trim(vpm_write_folder)//"/", NTIME, trim(pm_output_file_suffix)

      write (dummy_string, "(A)") achar(9)//'Writing PM solution to file: '//trim(filout)
      call vpm_print(dummy_string, nocolor, 2)

      NXs = NNbl_in(1)
      NYs = NNbl_in(2)
      NZs = NNbl_in(3)
      NXf = NNbl_in(4)
      NYf = NNbl_in(5)
      NZf = NNbl_in(6)

      ! INQUIRE if file exists
      inquire (file=trim(filout), exist=exist_flag)
      ! if file exists open it with overwrite mode
      if (exist_flag) then
         open (1, file=trim(filout), status='replace')
      else
         open (1, file=trim(filout))
      end if

      if (present(deformx_pm)) then
         write (1, "(a)") 'VARIABLES = "X" "Y" "Z" "U" "V" "W" '//&
                          '"VORTX" "VORTY" "VORTZ" "PSI1" "PSI2" '//&
                          '"PSI3" "DEFORMX" "DEFORMY" "DEFORMZ"'
      else
         write (1, "(a)") 'VARIABLES = "X" "Y" "Z" "U" "V" "W" "VORTX" "VORTY" "VORTZ" "PSI1" "PSI2" "PSI3" '
      end if


      write (1, "(a,I5.5,a,I5.5,a,I5.5)") 'ZONE I=', NXf - NXs + 1, &
                                              ' J=', NYf - NYs + 1, &
                                              ' K=', NZf - NZs + 1, ' F=POINT'

      do k = NZs, NZf 
         do j = NYs, NYf
            do i = NXs, NXf
               ! Structured grid coordinates
               cellx =  XMIN_pm + (I - 1)*DXpm
               celly =  YMIN_pm + (J - 1)*DYpm
               cellz =  ZMIN_pm + (K - 1)*DZpm

               ! Calculated Velocity
               velocx = velx(i, j, k)
               velocy = vely(i, j, k)
               velocz = velz(i, j, k)

               ! Deformation of vorticity
               if (present(deformx_pm)) then
                  deformx = deformx_pm(i, j, k)
                  deformy = deformy_pm(i, j, k)
                  deformz = deformz_pm(i, j, k)
               else
                  deformx = 0.d0
                  deformy = 0.d0
                  deformz = 0.d0
               end if

               ! RHS
               wmegax = -RHS(1, I, J, K)
               wmegay = -RHS(2, I, J, K)
               wmegaz = -RHS(3, I, J, K)
               
               ! Solution
               psi_1 = SOL(1, I - NXs + 1, J - NYs + 1, K - NZs + 1) 
               psi_2 = SOL(2, I - NXs + 1, J - NYs + 1, K - NZs + 1)
               psi_3 = SOL(3, I - NXs + 1, J - NYs + 1, K - NZs + 1)

               write (1, '(15(E14.7,1x))') cellx, celly, cellz,      &
                                           velocx, velocy, velocz,   &
                                           wmegax, wmegay, wmegaz,   &
                                           psi_1, psi_2, psi_3,      &
                                           deformx, deformy, deformz
            end do
         end do
      end do
      close (1)
   end subroutine write_pm_solution

   subroutine write_particles(NTIME, XPR, UPR, QPR, NVR)
      Implicit None
      integer, intent(in) :: NTIME, NVR
      real(dp), intent(in):: XPR(3, NVR), QPR(3, NVR), UPR(3, NVR)
      integer ::i
      logical :: exist_flag
      character*80 :: filout1

      write (filout1, '(a,i5.5,a)') trim(vpm_write_folder)//"/", NTIME, trim(particle_output_file_suffix)
      write (dummy_string, "(A)") achar(9)//'Writing particles to file: '//trim(filout1)
      call vpm_print(dummy_string, nocolor, 2)
      
      ! INQUIRE if file exists
      inquire (file=trim(filout1), exist=exist_flag)
      ! if file exists open it with overwrite mode
      if (exist_flag) then
         open (10, file=trim(filout1), status='replace')
      else
         open (10, file=trim(filout1))
      end if
       write (10, "(A)") 'VARIABLES = "i" "X" "Y" "Z" "U" "V" "W" "VORTX" "VORTY" "VORTZ"'
      write (10, "(a,I5.5,a)") 'ZONE I=', NVR, ' F=POINT'
      do i = 1, NVR
         write (10, '(i5.5,1x,9(E15.7,1x))') i, &
            XPR(1, i), XPR(2, i), XPR(3, i), &
            UPR(1, i), UPR(2, i), UPR(3, i), &
            QPR(1, i), QPR(2, i), QPR(3, i)
      end do

      ! call system('~/bin/preplot '&
      !    //filout1//' >/dev/null')
      ! call system('rm '//filout1)
      close(10)
   end subroutine write_particles

end Module vpm_lib

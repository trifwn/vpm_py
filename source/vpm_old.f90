
Module openmpth
   integer                       ::OMPTHREADS
End Module openmpth

Module vpm_lib
   use base_types, only: dp
   use vpm_vars
   use vpm_size
   use openmpth

contains
   include "vpm_mpi.f90"
   include "vpm_remesh.f90"
   include "vpm_time.f90"
   include "vpm_gcalc.f90"

   Subroutine vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
                  RHS_pm_in, Velx, Vely, Velz, NTIME_in, NI_in, NVRM_in)
      !  -> XP : particle positions (3 * NVR)
      !  -> QP : particle quantities (neqpm + 1 ) * NVR)
      !  -> UP : particle velocities
      !  -> GP : particle defromation (wmega * \nabla)  u
      !  -> NVR : number of particles
      !  -> neqpm : number of equations
      !  -> WhatToDo : 0 - initialize, 1 - solve, 2 - convect, 3 - project, 4 - back, 5 - diffuse
      !  -> \nabla^2 u_i = RHS_i (NEQ)
      !  -> Velx, Vely, Velz : velocity field at grid points
      !  -> NTIME_IN : current time
      !  -> NI_in: viscocity -> DIFFUSION OF VORTICITY
      !  -> NVRM : NVR MAX

      ! -> FIRST INITIALIZE
      ! use vpm_vars
      ! use vpm_size

      use pmeshpar, only: PI, PI2, PI4, ND, SOL_PM, IDVPM
      use parvar, only: QP, XP, UP, GP, NVR, print_particle_info
      use pmgrid, only: VELVRX_PM, VELVRY_PM, VELVRZ_PM, RHS_PM, NXpm, NYpm, NZpm, EPSVOL, Nblocks, print_RHS_pm
      use pmlib, only: pmesh
      use projlib, only: projlibinit, project_particles_3D, project_vol3d
      use yapslib, only: yaps3d
      ! use openmpth, only: OMPTHREADS

      use MPI
      use mkl_service

      Implicit None

      ! ARGUEMENTS
      real(dp), intent(inout), target   :: XP_in(:, :), QP_in(:, :), UP_in(:, :), GP_in(:, :)
      real(dp), intent(inout), pointer  :: RHS_pm_in(:, :, :, :)
      real(dp), intent(inout), pointer  :: velx(:, :, :), vely(:, :, :), velz(:, :, :)
      integer, intent(inout)                    :: NVR_in
      integer, intent(in)                       :: neqpm_in, WhatToDo, NVRM_in, NTIME_in
      real(dp), intent(in)              :: NI_in

      ! LOCAL VARIABLES
      integer                                   :: nb
      integer                                   :: II, itypeb
      integer, allocatable                      :: ieq(:)
      real(dp), allocatable             :: QINF(:)
      real(dp), allocatable             :: SOL_pm_bl(:, :, :, :), RHS_pm_bl(:, :, :, :)

      integer                                   :: ierr, my_rank, np
      integer                                   :: i ,j, k 

      ! DEPRECATED
      ! real(dp)                          :: totmass, totvor, MACH, pr
      ! real(dp)                          :: errorCFD, errorf, error, XPM, YPM,
      ! character*50                              :: outfil
      !real(dp)                           :: Xbound_tmp(6), XO(3)
      !integer                                    :: NN_tmp(3),NN_bl_tmp(6),
      ! real(dp)                          :: xi, yi, ksi1, ksi2, th1, th2, w1, w2, xl, yl, r
      ! integer                                   :: omp_get_num_threads, omp_get_max_threads

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

      ! WRITE ARGUEMENTS
      if (my_rank .eq. 0) then
         write (*, *) achar(27)//'[1;31m'
         if (WhatToDo .eq. 0) then
            write (*, *) achar(27)//'VPM: Initialize '
         else if (WhatToDo .eq. 1) then
            write (*, *) 'VPM: Project Particles- Solve PM- Get Velocities on PM'
         else if (WhatToDo .eq. 2) then
            write (*, *) 'VPM: Project Particles- Solve PM- Get Velocities and Deformation on PM'
         else if (WhatToDo .eq. 3) then
            write (*, *) 'VPM: Project Particles- Solve PM '
         else if (WhatToDo .eq. 4) then
            write (*, *) 'VPM: Project Particles- Interpolate from PM to particles '
         else if (WhatToDo .eq. 5) then
            write (*, *) 'VPM: Project Particles- Diffuse Particle vorticity ', achar(27)//'[0m'
         end if
         write (*, *) achar(9), 'Input Arguments:', achar(27)//'[1;34m'
         write (*, *) achar(9), achar(9), 'NTIME_in = ', NTIME_in
         write (*, *) achar(9), achar(9), 'WhatToDo = ', WhatToDo
         write (*, *) achar(9), achar(9), 'NVR_in = ', NVR_in
         write (*, *) achar(9), achar(9), 'neqpm_in = ', neqpm_in
         write (*, *) achar(9), achar(9), 'NI_in = ', NI_in
         write (*, *) achar(9), achar(9), 'NVRM_in = ', NVRM_in, achar(27)//'[0m'
      end if

      ! II = 1; ! 0: Serial Pmesh 1: Yaps
      II = 1
      ND = 3
      NVR_size = NVRM_in
      NTIME_pm = NTIME_in
      neqpm = neqpm_in
      NI = NI_in

      if (II .eq. 1) then
         ! call mkl_set_num_threads(1)
         ! call mkl_domain_set_num_threads(1)
      end if

      if (my_rank .eq. 0) then
         NVR = NVR_in
         nullify (XP)
         XP => XP_in
         nullify (QP)
         QP => QP_in; 
         nullify (UP)
         UP => UP_in; 
         nullify (GP)
         GP => GP_in
      else
         nullify (QP)
         nullify (XP)
         nullify (UP)
         nullify (GP)
         !nullify(velvrx_pm); nullify(velvry_pm); nullify(velvrz_pm)
      end if

      nullify (RHS_pm_in)

      nb = my_rank + 1

      !  -> WhatToDo : 0 - initialize
      ! DEFINE SIZES + ALLOCATIONS = 0
      if (WhatToDo .eq. 0) then

         if (my_rank .eq. 0) st = MPI_WTIME()
         call define_sizes
         if (my_rank .eq. 0) then
            et = MPI_WTIME()
            write (*, *) achar(9), 'VPM: Defining Sizes:', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
         end if
         if (allocated(velvrx_pm)) then
            deallocate (velvrx_pm, velvry_pm, velvrz_pm)
            allocate (velvrx_pm(NXpm, NYpm, NZpm), velvry_pm(NXpm, NYpm, NZpm), velvrz_pm(Nxpm, NYpm, NZpm))
            velvrx_pm = 0.d0; 
            velvry_pm = 0.d0; 
            velvrz_pm = 0.d0
         else
            allocate (velvrx_pm(NXpm, NYpm, NZpm), velvry_pm(NXpm, NYpm, NZpm), velvrz_pm(Nxpm, NYpm, NZpm))
            velvrx_pm = 0.d0; 
            velvry_pm = 0.d0; 
            velvrz_pm = 0.d0
         end if

         nullify (velx); nullify (vely); nullify (velz)
         velx => velvrx_pm; vely => velvry_pm; velz => velvrz_pm
         iwrite = 0
         return
      end if

      ! BCAST NVR
      call MPI_BCAST(NVR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (NVR .eq. 0) return

      ! NN_tmp is the number of cells in each direction
      ! NN_bl_tmp is the start and finish of the cells in each direction
      ! Xbound_tmp is the boundary of the domain
      NN_tmp(1:3) = NNbl(1:3, nb)
      NN_bl_tmp(1:6) = NNbl_bl(1:6, nb)
      allocate (SOL_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3)), RHS_pm_bl(neqpm, NN_tmp(1), NN_tmp(2), NN_tmp(3)))
      SOL_pm_bl = 0.d0; RHS_pm_bl = 0.d0

      !allocate RHS_pm for parallel projection
      if (allocated(RHS_pm)) then
         deallocate (RHS_pm)
         allocate (RHS_pm(neqpm + 1, NXpm, NYpm, NZpm))
         RHS_pm = 0.d0
      else
         allocate (RHS_pm(neqpm + 1, NXpm, NYpm, NZpm))
         RHS_pm = 0.d0
      end if
      !allocate SOL_pm for parallel projection
      if (allocated(SOL_pm)) then
         deallocate (SOL_pm)
         allocate (SOL_pm(neqpm, NXpm, NYpm, NZpm))
         SOL_pm = 0.d0
      else
         allocate (SOL_pm(neqpm, NXpm, NYpm, NZpm))
         SOL_pm = 0.d0
      end if
      ! SOL_PM AND RHS_PM ARE 0 BOTH ON THE BL AND PM
      call project_particles ! PARTICLES TO -> PM MEANING FROM Qs We get RHS_pm

      !  -> WhatToDo: 4 - interpolate from PM to particles
      if (WhatToDo .eq. 4) then
         if (my_rank .eq. 0) then

            call calc_velocity_serial_3d(-1) ! ONLY DELATATION THETA STO PM
            ! WE LOST SOL_PM
            ! SOL_PM IS NOW FILLED WITH THETA (DEFORMATION)

            !    itypeb=1!normal back to particles
            !    call back_to_particles_3D(SOL_pm,XP,QP,UP,GP,&
            !                              velvrx_pm,velvry_pm,velvrz_pm,&
            !                              Xbound,Dpm,NN,NVR,neqpm,interf_iproj,itypeb,NVR_size)

            if (IPMWRITE .GT. 0) then
               do i = 1, IPMWRITE
                  if (NTIME_pm .ge. IPMWSTART(i) .and. NTIME_pm .le. (IPMWSTART(i) + IPMWSTEPS(i))) call write_pm_solution
               end do
            end if
            ! if (iynslice .ne. 0) call write_pm_solutionXavatar(iynslice)
         end if

         !
         itypeb = 1 ! MAYBE IT SHOULD BE 2?
         ! WHEN ITYPEB = 1 WE GET THE UP AND GP From the velocity
         call back_to_particles_par

         iwrite = 0
         !  RHS_pm_in=>RHS_pm
         !  velx=>velvrx_pm; vely=>velvry_pm;velz=>velvrz_pm

         ! DEALLOCATE THE MEMORY
         deallocate (SOL_pm)
         deallocate (RHS_pm)
         deallocate (SOL_pm_bl, RHS_pm_bl)
         return
      end if

      !  -> WhatToDo : 5 - diffuse
      if (WhatToDo .eq. 5) then
         if (allocated(velvrx_pm)) then
            deallocate (velvrx_pm, velvry_pm, velvrz_pm)
            allocate (velvrx_pm(NXpm, NYpm, NZpm), velvry_pm(NXpm, NYpm, NZpm), velvrz_pm(Nxpm, NYpm, NZpm))
            velvrx_pm = 0.d0; velvry_pm = 0.d0; velvrz_pm = 0.d0
         else
            allocate (velvrx_pm(NXpm, NYpm, NZpm), velvry_pm(NXpm, NYpm, NZpm), velvrz_pm(Nxpm, NYpm, NZpm))
            velvrx_pm = 0.d0; velvry_pm = 0.d0; velvrz_pm = 0.d0
         end if

         nullify (velx); nullify (vely); nullify (velz)
         velx => velvrx_pm; vely => velvry_pm; velz => velvrz_pm
         if (my_rank .eq. 0) then
            !diffusion stores -NI*grad^2 w * Vol in GP(1,:)

            call diffuse_vort_3d ! DIFFUSION OF VORTICITY
            ! SOL_pm = -VIS \nabla \cdot RHS_pm

            ! itypeb=2!back to particles the diffused vorticity
            ! call back_to_particles_3D(SOL_pm,XP,QP,UP,GP,&
            !                           velvrx_pm,velvry_pm,velvrz_pm,&
            !                           Xbound,Dpm,NN,NVR,neqpm,interf_iproj,itypeb,NVR_size)
         end if
         itypeb = 2
         ! WHEN ITYPEB = 2 WE GET THE GP FROM THE SOL_PM (DEFORMATION) and from QP
         call back_to_particles_par
         deallocate (SOL_pm, RHS_pm)
         deallocate (SOL_pm_bl, RHS_pm_bl)
         return

      end if

      !------------------------------ FOR WHATTODO = 1, 2, 3 --------------------------------

      if (my_rank .eq. 0) st = MPI_WTIME()
      !call diffuse_vort_3d

      ! SCATTER PROBLEM TO ALL PROCESSORS
      call rhsbcast(RHS_pm, NN, neqpm + 1) ! RHS PM -> TO ALL PROCESSORS -> RHS_PM_BL
      call print_RHS_pm()

      IF (II .ne. 0) then
         nb = my_rank + 1
         NN_tmp(1:3) = NNbl(1:3, nb)
         call rhsscat_3d(BLOCKS, NN_tmp, NNbl, NNbl_bl, NN_bl, nb_i, nb_j, nb_k, RHS_pm_bl) !
      END IF
      if (my_rank .eq. 0) then
         et = MPI_WTIME()
         write (*, *) achar(9), 'VPM:Broadcasting/Scattering RHS:', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
      end if


      ! SOLVE THE PROBLEM
      if (my_rank .eq. 0) st = MPI_WTIME()
      call pmesh_solve
      if (my_rank .eq. 0) then
         write (*, *) achar(9), '---> max SOL_PM', maxval(abs(SOL_pm))
         write (*, *) achar(9), '---> max RHSL_PM', maxval(abs(RHS_pm))
         write (*, *) ""
         et = MPI_WTIME()
         write (*, *) achar(9), 'VPM:Solving Particle Mesh', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
      end if

      ! CLEAR VELOCITIES
      if (allocated(velvrx_pm)) then
         deallocate (velvrx_pm, velvry_pm, velvrz_pm)
         allocate (velvrx_pm(NXpm, NYpm, NZpm), velvry_pm(NXpm, NYpm, NZpm), velvrz_pm(Nxpm, NYpm, NZpm))
         velvrx_pm = 0.d0; velvry_pm = 0.d0; velvrz_pm = 0.d0
      else
         allocate (velvrx_pm(NXpm, NYpm, NZpm), velvry_pm(NXpm, NYpm, NZpm), velvrz_pm(Nxpm, NYpm, NZpm))
         velvrx_pm = 0.d0; velvry_pm = 0.d0; velvrz_pm = 0.d0
      end if

      nullify (velx); nullify (vely); nullify (velz)
      velx => velvrx_pm; vely => velvry_pm; velz => velvrz_pm

      !  -> WhatToDo :  1 - GET VELOCITIES ON PM FROM PM SOLUTION
      if (WhatToDo .eq. 1) then
         if (my_rank .eq. 0) then
            velvrx_pm = 0.d0; velvry_pm = 0.d0; velvrz_pm = 0.d0

            !call convect_first_order(Xbound,Dpm,NN,NN_bl)

            ! FROM THE SOLUTION OF PM WE GET THE VELOCITIES ON THE GRID
            call calc_velocity_serial_3d(0) ! VELOCITY STO PM
            ! SOL_PM is stille the solution of vorticity

            !  if(mod(NTIME,NWRITE).eq.0) call write_pm_solution(NTIME)
            !if (ND.eq.3) then
            ! call hill_error(NN,NN_bl,Xbound,Dpm,SOL_pm,velvrx_pm,velvry_pm,velvrz_pm)
            ! call write_pm_solution(NTIME)
            !stop
            !endif

            !  RHS_pm_in=>RHS_pm
            deallocate (SOL_pm)
            deallocate (RHS_pm)
            ! Print Velx
            write (*, *) achar(9), 'MAX VELOCITY INSIDE PM:'
            write (*, *) achar(9), achar(9), maxval(velx)
            write (*, *) achar(9), achar(9), maxval(vely)
            write (*, *) achar(9), achar(9), maxval(velz)
            
         end if
         
         deallocate (SOL_pm_bl, RHS_pm_bl)
         return
      end if

      !  -> WhatToDo : 2 - GET VELOCITIES AND DEFORMATION ON PM FROM PM SOLUTION
      if (WhatToDo .eq. 2) then
         if (my_rank .eq. 0) then
            velvrx_pm = 0; velvry_pm = 0; velvrz_pm = 0
            !call convect_first_order(Xbound,Dpm,NN,NN_bl)
            
            st = MPI_WTIME()
            call calc_velocity_serial_3d(1) ! VELOCITY AND THETA STO PM
            write (*, *) achar(9), 'MAX VELOCITY INSIDE PM:'
            write (*, *) achar(9), achar(9), maxval(velvrx_pm)
            write (*, *) achar(9), achar(9), maxval(velvry_pm)
            write (*, *) achar(9), achar(9), maxval(velvrz_pm)
            write (*, *) achar(9), "Velocity shape: "
            write (*, *) achar(9), achar(9), shape(velvrx_pm)
            ! CHeck for nan values
            if (any(isnan(velvrx_pm)) .or. any(isnan(velvry_pm)) .or. any(isnan(velvrz_pm))) then
               ! Print Velx
               do i = 1, NXpm
                  do j = 1, NYpm
                     do k = 1, NZpm
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
            ! SOL_PM IS NOW FILLED WITH THETA (DEFORMATION)


            et = MPI_WTIME()
            write (*, *) achar(9), 'VPM: Velocity Calculation using FD', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
            !     call calc_antidiffusion
            ! itypeb=1
            ! call back_to_particles_3D(SOL_pm,XP,QP,UP,GP,&
            !                           velvrx_pm,velvry_pm,velvrz_pm,&
            !                           Xbound,Dpm,NN,NVR,neqpm,interf_iproj,itypeb,NVR_size)
            !      if(mod(NTIME_pm,20).eq.0.or.NTIME_pm.eq.1) call write_pm_solution
            !if(mod(NTIME_pm,100).eq.0.or.NTIME_pm.eq.1) call write_pm_solution
            if (IPMWRITE .GT. 0) then
               do i = 1, IPMWRITE
                  if (NTIME_pm .ge. IPMWSTART(i) .and. NTIME_pm .le. (IPMWSTART(i) + IPMWSTEPS(i))) then
                     st = MPI_WTIME()
                     call write_pm_solution
                     et = MPI_WTIME()
                     write (*, *) achar(9), 'VPM: Writing Solution:', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
                  end if
               end do
            end if
            ! if (iynslice .ne. 0) call write_pm_solutionXavatar(iynslice)
            ! call writeline

         end if
         itypeb = 1
         call back_to_particles_par ! INTERPOLATION FROM PM TO PARTICLES


         iwrite = 0
         !if(IPMWRITE.GT.0) then
         !  do i=1,IPMWRITE
         ! if(NTIME.ge.IPMWSTART(i).and.NTIME.le.(IPMWSTART(i)+IPMWSTEPS(i))) call write_pm_solution(NTIME)
         !  enddo
         !endif
         !if (ND.eq.3) then
         ! call hill_error(NN,NN_bl,Xbound,Dpm,SOL_pm,velvrx_pm,velvry_pm,velvrz_pm)
         ! call write_pm_solution(NTIME)
         !stop
         !endif

         !  RHS_pm_in=>RHS_pm
         deallocate (SOL_pm)
         deallocate (RHS_pm)
         ! deallocate (velvrx_pm, velvry_pm, velvrz_pm)
         deallocate (SOL_pm_bl, RHS_pm_bl)
         return
      end if

   contains
      ! CONTAINS

      Subroutine pmesh_solve !
         integer :: rank
         !Yaps or Serial Pmesh

         IF (II .eq. 0) then
            if (my_rank .eq. 0) write (*, *) achar(9), achar(9), 'Solving PM with Serial Pmesh'

            IF (my_rank .eq. 0) then
               write (*, *) 'Solving_pm'
               SOL_pm(1:neqpm, :, :, :) = 0.0
               itree = iyntree
               iynbc = 1!for infinite domain bc's
               Nblocks = 1
               call pmesh(SOL_pm, RHS_pm, QP, XP, Xbound, DPm, NN, NN_bl, ND, Nblocks, ibctyp, 1, neqpm, &
                          iynbc, NVR, itree, ilevmax)
               ! call calc_velocity_serial_3d(1)
            END IF
            !--------------------------------------------
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call velbcast_3d
         ELSE
            if (my_rank .eq. 0) write (*, *) achar(9), 'Solving PM with YAPS'
            iret = 0
            call yaps3d(SOL_pm_bl, RHS_pm_bl, Xbound_bl, Xbound_coarse, Dpm, Dpm_coarse, NNbl, NNbl_bl, &
                        NN_coarse, NN_bl_coarse, ND, BLOCKS, ibctyp, 1, neqpm, ncoarse, NBI, NBJ, NBK, nb_i, nb_j, nb_k, &
                        iret, iyntree, ilevmax, neqpm)

            if (my_rank .eq. 0) then
               write (*, *) achar(9), achar(27)//'[1;34m', 'Final PM block solution values'
            END IF
            do rank = 0, np-1
               if (my_rank .eq. rank) then
                  write (*, *) achar(9), achar(9), 'np=', my_rank, 'SOL', maxval(abs(SOL_pm_bl(neqpm, :, :, :)))
               end if
               call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            end do
            if (my_rank .eq. 0) write (*, *) achar(27)//'[0m'
            nb = my_rank + 1
            NN_tmp(1:3) = NNbl(1:3, nb)
            NN_bl_tmp(1:6) = NNbl_bl(1:6, nb)
            Xbound_tmp(1:6) = Xbound_bl(1:6, nb)
            call solget_3d(BLOCKS, NBI, NBJ, NBK, NN_tmp, NNbl, NNbl_bl, NN_bl, SOL_pm_bl) ! GATHER SOLUTION
            !if (my_rank.eq.0) call calc_velocity_serial_3d(1)
            ! call velbcast_3d
         END IF
         !--------------------------------------------

      End Subroutine pmesh_solve

      Subroutine project_particles
         ! BCAST NVR
         call MPI_BCAST(NVR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         ! SPLIT PARTICLES on PROCESSORS
         NVR_p = NVR/np
         ! ROUND
         if (my_rank .eq. 0) NVR_p = NVR_p + mod(NVR, np)
         ! SCATTER PARTICLES ON EACH PROCESSOR
         allocate (XP_scatt(3, NVR_p), QP_scatt(neqpm + 1, NVR_p), NVR_projscatt(NVR_p))
         NVR_projscatt = interf_iproj
         call particles_scat

         ! INITIALIZATION OF PROJECTION LIBRARY
         call projlibinit(Xbound, Dpm, NN, NN_bl, EPSVOL, IDVPM, ND)

         ! PROJECT PARTICLES TO PM
         st = MPI_WTIME()
         allocate (ieq(neqpm + 1), QINF(neqpm + 1))
         QINF = 0.d0
         do i = 1, neqpm + 1
            ieq(i) = i
         end do
         ! FILLS RHS_PM FROM PARTICLES
         if (my_rank .eq. 0) then 
            write (*, *) achar(9), achar(27)//'[1;34m', 'VPM: Projecting Particles to PM' , achar(27)//'[0m'
         endif
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         
         call project_particles_3D(RHS_pm, QP_scatt, XP_scatt, NVR_projscatt, NVR_p, neqpm + 1, ieq, neqpm + 1, QINF, NVR_p)
         call proj_gath(NN)
         ! RHS IS NOW FILLED

         if (my_rank .eq. 0) then
            !allocate(ieq(neqpm+1),QINF(neqpm+1))
            !QINF=0.d0
            !do i=1,neqpm+1
            !   ieq(i)=i
            !enddo
            !allocate (NVR_projscatt(NVR_size))
            !NVR_projscatt=interf_iproj
            !call  projlibinit(Xbound,Dpm,NN,NN_bl,EPSVOL,IDVPM,ND)
            !call project_particles_3D(RHS_pm,QP,XP,NVR_projscatt,NVR,neqpm+1,ieq,neqpm+1,QINF,NVR_size)

            call project_vol3d(RHS_pm, neqpm + 1, ieq, neqpm + 1, IDVPM)
            ! RHS_PM IS NOW NORMALIZED BY VOLUME (DENSITY)

            et = MPI_WTIME()
            write (*, *) achar(9), 'VPM: Projecting Particles to PM', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
            ! deallocate(ieq,QINF)
            ! deallocate(NVR_projscatt)
         end if
         deallocate (ieq, QINF)
         deallocate (XP_scatt, QP_scatt, NVR_projscatt)

      End Subroutine project_particles

      Subroutine back_to_particles_par
         if (my_rank .eq. 0) st = MPI_WTIME()

         ! BROADCASTING
         call rhsbcast(RHS_pm, NN, neqpm + 1)
         call rhsbcast(SOL_pm, NN, neqpm)
         if (itypeb .eq. 1) call velbcast_3d
         call MPI_BCAST(NVR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

         ! ALLOCATE
         NVR_p = NVR/np
         if (my_rank .eq. 0) NVR_p = NVR_p + mod(NVR, np)
         allocate (XP_scatt(3, NVR_p), QP_scatt(neqpm + 1, NVR_p), &
                   UP_scatt(3, NVR_p), GP_scatt(3, NVR_p))
         UP_scatt = 0.d0; GP_scatt = 0.d0

         ! SCATTERING XP AND QP
         call particles_scat
         if (my_rank .eq. 0) then
            et = MPI_WTIME()
            write (*, *) achar(9), 'VPM_LIB: Scattering Particles:', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
         end if
         if (my_rank .eq. 0) st = MPI_WTIME()

         ! WHEN ITYPEB = 1 WE GET THE UP AND GP From the velocity
         ! WHEN ITYPEB = 2 WE GET THE GP FROM THE SOL_PM (DEFORMATION) and from QP
         call back_to_particles_3D(SOL_pm, XP_scatt, QP_scatt, UP_scatt, GP_scatt, &
                                   velvrx_pm, velvry_pm, velvrz_pm, &
                                   Xbound, Dpm, NN, NVR_p, neqpm, interf_iproj, itypeb, NVR_p)
         ! GATHERS XP, QP, UP, GP
         call particles_gath

         deallocate (XP_scatt, QP_scatt, UP_scatt, GP_scatt)
         if (my_rank .eq. 0) then
            et = MPI_WTIME()
            write (*, *) achar(9), 'VPM_LIB: Interpolating Particles in: ', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
         end if

      End Subroutine back_to_particles_par

   End Subroutine vpm

   Subroutine define_sizes
      ! INITIALIZES GRID
      ! XBOUND : BOUNDARY OF THE DOMAIN (KATW ARISTERA ME PANW DEXIA)
      ! NN : NUMBER OF CELLS IN EACH DIRECTION
      ! DPM : CELL SIZE
      ! NN_BL(6) : (NXS,NYS,NZS,NXF,NYF,NZF) START FINISH
      ! call definepm(3, Xbound, Dpm, ND, ndumcell, nsiz, NN, NN_bl)

      use pmeshpar, only: nd, ndumcell
      use vpm_vars, only: interf_iproj
      use vpm_size, only: NN_bl, NN
      use parvar, only: XP, NVR
      use pmgrid, only: XMIN_pm, YMIN_pm, ZMIN_pm, DXpm, DYpm, DZpm, XMAX_pm, &
                        YMAX_pm, ZMAX_pm, NXs_bl, NXf_bl, NYS_bl, NYF_bl, &
                        NZs_bl, NZF_bl, NXPM, NYPM, NZPM, DVPM, DXPM2, DYPM2, &
                        DZPM2
      use pmlib, only: definepm
      use MPI

      Implicit None
      integer    :: nsiz(3), nsiz_bl(3)
      integer    :: i, j, k, np, my_rank, ierr, nb, istep, lev
      !real(dp)              :: Xbound_tmp(6)
      !integer                       :: NN_tmp(3),NN_bl_tmp(6)

      call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

      BLOCKS = np
      !-------First Change Dpm so that the numbers of cells divides
      !-------by nsize i.e with NBI,NBJ,ncoarse,levmax depending on the criterion
      !------that's why ndumcell=0
      if (my_rank .eq. 0) then
         if ((NTIME_pm .eq. 0 .and. idefine .eq. 1) .or. idefine .eq. 0) then
            XMIN_pm = minval(XP(1, 1:NVR)) - interf_iproj*DXpm
            YMIN_pm = minval(XP(2, 1:NVR)) - interf_iproj*DYpm
            ZMIN_pm = minval(XP(3, 1:NVR)) - interf_iproj*DZpm

            XMAX_pm = maxval(XP(1, 1:NVR)) + interf_iproj*DXpm
            YMAX_pm = maxval(XP(2, 1:NVR)) + interf_iproj*DYpm
            ZMAX_pm = maxval(XP(3, 1:NVR)) + interf_iproj*DZpm
            !else
            !  XMAX_pm=maxval(XP(1,1:NVR)) + interf_iproj*DXpm
         end if
      end if

      call MPI_BCAST(XMIN_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(YMIN_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ZMIN_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call MPI_BCAST(XMAX_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(YMAX_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ZMAX_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      Xbound(1) = XMIN_pm; Xbound(2) = YMIN_pm; Xbound(3) = ZMIN_pm
      Xbound(4) = XMAX_pm; Xbound(5) = YMAX_pm; Xbound(6) = ZMAX_pm
      Dpm(1) = DXpm; Dpm(2) = DYpm; Dpm(3) = DZpm
      ndumcell = 0
      nsiz(1) = NBI*ncoarse
      nsiz(2) = NBJ*ncoarse
      nsiz(3) = NBK*ncoarse
      if ((NTIME_pm .eq. 0 .and. idefine .eq. 1) .or. idefine .eq. 0) then
         ! DPM (DX , DY , DZ)  (CELL SIZE)
         ! ND
         ! NDUMCELL -> EPEKTASI SE CELLS TOU XWRIOU SE D_EXTENDED
         ! ITYPE :
         ! NSIZE : NN mod nsiz == 0
         !         ndum_new(1)  = nsize(1) - mod(NN(1)-1,nsize(1))
         if (my_rank .eq. 0) then
            write (*, *) achar(9), 'Calling definepm with idefine= ', idefine
         end if
         call definepm(3, Xbound, Dpm, ND, ndumcell, nsiz, NN, NN_bl)
         !else
         !  call  definepm(4,Xbound,Dpm,ND,ndumcell,nsiz,NN,NN_bl)
      end if
      XMIN_pm = Xbound(1); YMIN_pm = Xbound(2); ZMIN_pm = Xbound(3)
      XMAX_pm = Xbound(4); YMAX_PM = Xbound(5); ZMAX_pm = Xbound(6)
      NXpm = NN(1); NYpm = NN(2); NZpm = NN(3)
      NXs_bl(1) = NN_bl(1); NYs_bl(1) = NN_bl(2); NZs_bl(1) = NN_bl(3)
      NXf_bl(1) = NN_bl(4); NYf_bl(1) = NN_bl(5); NZf_bl(1) = NN_bl(6)
      DXpm = Dpm(1); DYpm = Dpm(2); DZpm = Dpm(3)
      if (ND .eq. 2) then
         DVpm = DXpm*DYpm
         DXpm2 = 2*DXpm
         DYpm2 = 2*DYpm
      else
         DVpm = DXpm*DYpm*DZpm
         DXpm2 = 2*DXpm
         DYpm2 = 2*DYpm
         DZpm2 = 2*DZpm
      end if
      if (my_rank .eq. 0) then
         write (*, *) achar(9), 'The PM Grid has:'
         write (*, *) achar(9), achar(9), "DPM=", Dpm(1), Dpm(2), Dpm(3)
         write (*, *) achar(9), achar(9), "XMIN_PM:", Xbound(1), "XMAX_PM:", Xbound(4)
         write (*, *) achar(9), achar(9), "YMIN_PM:", Xbound(2), "YMAX_PM:", Xbound(5)
         write (*, *) achar(9), achar(9), "ZMIN_PM:", Xbound(3), "ZMAX_PM:", Xbound(6)
         write (*, *) achar(9), achar(9), 'NN:', NN(1), NN(2), NN(3)
      end if
      !define block grid so the are divided by ncoarse and ilevmax
      !so as to have the coarse information at the boundaries exactly.
      ndumcell_bl = ncoarse
      nsiz_bl = ncoarse!ndumcell_coarse!*2*2**ilevmax
      if (.not. allocated(XBound_bl)) then
         allocate (Xbound_bl(6, BLOCKS), NNbl(3, BLOCKS), NNbl_bl(6, BLOCKS))
      end if
      NXB = int(nint(((Xbound(4) - Xbound(1))/Dpm(1))))
      NYB = int(nint(((Xbound(5) - Xbound(2))/Dpm(2))))
      NZB = int(nint(((Xbound(6) - Xbound(3))/Dpm(3))))
      NXbl = NXB/NBI
      NYbl = NYB/NBJ
      NZbl = NZB/NBK
      if (my_rank .eq. 0) then
         write (*, *) achar(9), "The size of the grid is ", NXB, NYB, NZB
         write (*, *) achar(9), "The subdivision is", NBI, "x", NBJ, "x", NBK
         write (*, *) achar(9), 'Check sizes:'
         write (*, *) achar(9), achar(9), 'NXbl->NXB/NBI', NXbl, " ->", float(NXB)/NBI
         write (*, *) achar(9), achar(9), 'NYbl->NYB/NBJ', NYbl, " ->", float(NYB)/NBJ
         write (*, *) achar(9), achar(9), 'NZbl->NZB/NBK', NZbl, " ->", float(NZB)/NBK
         write (*, *) achar(9), achar(9), 'ncourse', ncoarse
         write (*, *) achar(9), 'Check sizes coarse:'
         write (*, *) achar(9), achar(9), "NXbl/ncoarse", float(NXbl)/ncoarse
         write (*, *) achar(9), achar(9), "NYbl/ncoarse", float(NYbl)/ncoarse
         write (*, *) achar(9), achar(9), "NZbl/ncoarse", float(NZbl)/ncoarse
         write (*, *) achar(9), "Assigning a block grid to each processor"
      end if

      do k = 1, NBK
         do j = 1, NBJ
            do i = 1, NBI
               nb = (k - 1)*NBJ*NBI + (j - 1)*NBI + i
               if (my_rank .eq. 0) then
                  write (*, *) achar(9), achar(9), 'Block', nb
               end if
               Xbound_bl(1, nb) = Xbound(1) + (i - 1)*(NXbl)*Dpm(1)
               Xbound_bl(4, nb) = Xbound(1) + (i)*(NXbl)*Dpm(1)

               Xbound_bl(2, nb) = Xbound(2) + (j - 1)*(NYbl)*Dpm(2)
               Xbound_bl(5, nb) = Xbound(2) + (j)*(NYbl)*Dpm(2)

               Xbound_bl(3, nb) = Xbound(3) + (k - 1)*(NZbl)*Dpm(3)
               Xbound_bl(6, nb) = Xbound(3) + (k)*(NZbl)*Dpm(3)
               if (my_rank .eq. 0) then
                  ! write (*, *)  &
                  !             "XMIN=", Xbound_bl(1, nb),&
                  !             "YMIN=", Xbound_bl(2, nb),&
                  !             "ZMIN=", Xbound_bl(3, nb),&
                  !             "XMAX=", Xbound_bl(4, nb),&
                  !             "YMAX=", Xbound_bl(5, nb),&
                  !             "ZMAX=", Xbound_bl(6, nb)
               end if
               Xbound_tmp(1:6) = Xbound_bl(1:6, nb)
               call definepm(1, Xbound_tmp, Dpm, ND, ndumcell_bl, nsiz_bl, NN_tmp, NN_bl_tmp)
               Xbound_bl(1:6, nb) = Xbound_tmp(1:6)
               NNbl(1:3, nb) = NN_tmp(1:3) ! KOMVOI
               NNbl_bl(1:6, nb) = NN_bl_tmp(1:6) ! KOMVOI START STOP
               if (nb .eq. my_rank + 1) then
                  nb_i = i
                  nb_j = j
                  nb_k = k
               end if
            end do
         end do
      end do

      !--B
      nb = my_rank + 1
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
      !NXs_bl(1)=NN_bl(1);NYs_bl(1)=NN_bl(2);NZs_bl(1)=NN_bl(3)
      !NXf_bl(1)=NN_bl(4);NYf_bl(1)=NN_bl(5);NZf_bl(1)=NN_bl(6)
      !print *,'final mesh',NN
      !-----

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

               if (int((NNbl_bl(4, 1) - NNbl_bl(1, 1))/istep) .eq. 0 .or. &
                   int((NNbl_bl(5, 1) - NNbl_bl(2, 1))/istep) .eq. 0 .or. &
                   int((NNbl_bl(6, 1) - NNbl_bl(3, 1))/istep) .eq. 0) then
                  ilevmax = lev - 1
                  print *, 'Changing number of levels', ilevmax
                  exit
               end if

            end if
            !  print *, nn_lev(1:2,lev)
         end do
      end if
      call MPI_BCAST(ilevmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   End Subroutine define_sizes

   Subroutine write_pm_solution
      use pmgrid, only: NXf_bl , NXs_bl, NYf_bl, NYs_bl, NZf_bl, NZs_bl, &
                        DXpm, DYpm, DZpm, XMIN_pm, YMIN_pm, ZMIN_pm, &
                        velvrx_pm, velvry_pm, velvrz_pm, RHS_pm

      character*50        :: filout
      integer           :: i, j, k
      real(dp)  :: XPM, YPM, ZPM, velocx, velocy, velocz

      ! if(iwrite.ne.0) return
      write (filout, '(a,i5.5,a)') "sol/", NTIME_pm, 'solution.dat'
      write (*, *) achar(9), 'Writing solution to file: ', trim(filout)
      open (1, file=filout)
      WRITE (1, "(a)") 'VARIABLES = "X" "Y" "Z" "U" "V" "W" "VORTX" "VORTY" "VORTZ"'
      WRITE (1, "(a,I5.5,a,I5.5,a,I5.5)") 'ZONE I=', NXf_bl(1) - NXs_bl(1) + 1, ' J=', NYf_bl(1) - NYs_bl(1) + 1, &
         ' K=', NZf_bl(1) - NZs_bl(1) + 1, ' F=POINT'
      ! Write the size of the loops

      do k = NZs_bl(1), NZf_bl(1)
         do j = NYs_bl(1), NYf_bl(1)
            do i = NXs_bl(1), NXf_bl(1)
               ! WRITE(1,*)'ZONE I=',NXpm,' J=',NYpm,' F=POINT'
               ! do j=1,NYpm
               !   do i=1,NXpm
               XPM = XMIN_pm + (I - 1)*DXpm
               YPM = YMIN_pm + (J - 1)*DYpm
               ZPM = ZMIN_pm + (K - 1)*DZpm
               velocx = velvrx_pm(i, j, k)
               velocy = velvry_pm(i, j, k)
               velocz = velvrz_pm(i, j, k)

               WRITE (1, '(9(E20.10,1x))') XPM, YPM, ZPM, &
                                          velocx, velocy, velocz, &
                                          -RHS_pm(1, I, J, K), &
                                          -RHS_pm(2, I, J, K), &
                                          -RHS_pm(3, I, J, K)
                                          !,RHS_pm(4,I,J,K), &
                                          ! SOL_pm(1,I,J,K),SOL_pm(2,I,J,K), SOL_pm(3,I,J,K)

            end do
         end do
      end do
      close (1)
      iwrite = 1
   End Subroutine write_pm_solution

   Subroutine write_particles(NTIME, XPR, UPR, QPR, NVR)
      Implicit None
      integer, intent(in) :: NTIME, NVR
      real(dp), intent(in):: XPR(3, NVR), QPR(3, NVR), UPR(3, NVR)
      integer ::i
      character*80 :: filout1

      write (filout1, '(a,i5.5,a)') "sol/", NTIME, 'vr.dat'
      write (*, *) achar(9), 'Writing particles to file: ', trim(filout1)
      open (10, file=filout1)
      WRITE (10, *) 'VARIABLES = "i" "X" "Y" "Z" "U" "V" "W" "VORTX" "VORTY" "VORTZ"'
      do i = 1, NVR
         write (10, '(i5.5,9(E20.10,1x))') i, &
            XPR(1, i), XPR(2, i), XPR(3, i), &
            UPR(1, i), UPR(2, i), UPR(3, i), &
            QPR(1, i), QPR(2, i), QPR(3, i)
      end do
      ! call system('~/bin/preplot '&
      !    //filout1//' >/dev/null')
      ! call system('rm '//filout1)
      ! close(10)
   End Subroutine write_particles

End Module vpm_lib

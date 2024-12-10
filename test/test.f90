module test_mod
   use vpm_types, only: dp
   real(dp), allocatable, target    :: XPR(:, :), QPR(:, :), UPR(:, :), GPR(:, :), &
                                       XPO(:, :), QPO(:, :), &
                                       XP_in(:, :), QP_in(:, :)
   real(dp), pointer       :: RHS_ptr(:,:,:,:)
   real(dp), pointer       :: SOL_ptr(:,:,:,:)
   real(dp), pointer       :: vel_ptr(:,:,:,:)
   integer, allocatable    :: qflag(:)
   integer                 :: NVR_ext, NVR_ext_init
end module test_mod

Program test_pm
   use vpm_types, only: dp
   use pmgrid, only:     XMIN_pm, NXs_fine_bl, DXpm, DYpm, DZpm, set_RHS_pm, ncoarse, IDVPM, &
                         SOL_pm, velocity_pm, deform_pm
   use vpm_vars, only:   mrem, interf_iproj,   &
                         IPMWRITE, idefine, IPMWSTART, IPMWSTEPS, OMPTHREADS, nremesh, iyntree,&
                         ilevmax, ibctyp, NBI, NBJ, NBK
   use vpm_size, only:   st, et, fine_grid
   use test_mod, only:   XPR, QPR, UPR, GPR, NVR_ext, &
                         QPO, XPO, Qflag, &
                         QP_in, XP_in, &
                         RHS_ptr, vel_ptr, SOL_ptr
   use test_app, only:   hill_assign
   use parvar, only:     NVR
   use vpm_lib, only:    vpm 
   use vpm_remesh, only: remesh_particles_3d
   use file_io, only:    write_pm_solution_hdf5, write_particles_hdf5, case_folder 
   use console_io, only: vpm_print, red, green, blue, yellow, nocolor, dummy_string, tab_level, VERBOCITY
   use serial_vector_field_operators, only: divergence
   use MPI

   implicit none
   real(dp)             :: Vref, NI_in, DT_in, FACDEF, T, &
                           XMIN, XMAX, UINF(3)
   integer              :: NVR_size
   integer              :: my_rank, np, ierr, i,ieq, neq, j, max_iter, ncell_rem
   logical              :: pmfile_exists
   real(dp)             :: sphere_radius = 1.0_dp
   real(dp)             :: sphere_z_0 = 0.0_dp
   real(dp)             :: u_free_stream = 1.0_dp
   real(dp)             :: CFL_x, CFL_y, CFL_z, CFL
   real(dp)             :: CFL_treshold_up = 0.9_dp
   real(dp)             :: CFL_treshold_down = 0.5_dp
   real(dp)             :: CFL_target = 0.7_dp
   character(len=250)   :: cfl_file
   real(dp), allocatable :: divergence_HILL(:,:,:)

   call MPI_INIT(ierr)
   call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! READ SETTINGS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   case_folder = 'results_deform_dt=0.01/'
   inquire (file='pm.inp', exist=pmfile_exists)
   if(my_rank.eq.0) print "(A,L1)", 'pm.inp exists:', pmfile_exists
   if (pmfile_exists) then
      open (1, file='pm.inp')
      read (1, *) DXpm, DYpm, DZpm     ! CELL SIZES
      read (1, *) interf_iproj         ! INTERFACE PROJECTION FUNCTION 1 , 2, 3 , 4
      read (1, *) ibctyp               ! 1 - PERIODIC, 2 - INFLOW, 3 - OUTFLOW, 4 - WALL
      read (1, *) IDVPM                ! Variable/Constant Volume(0,1)
      read (1, *) ncoarse              ! NUMBER OF FINE CELLS PER COARSE CELL per dir
      read (1, *) NBI, NBJ, NBK        !  NBI x NBJ x NBK = NUM OF PROCESSORS (NP)
      read (1, *) nremesh, ncell_rem   ! 0: NO REMESHING, 1: REMESHING, ncell_rem: PARTICLE PER CELL
      read (1, *) iyntree, ilevmax     ! 1: TREE 0: NO TREE, 3: NUMB OF SUBDIVISION (2^3)
      read (1, *) OMPTHREADS           ! 1 - OPENMP THREADS
      read (1, *) idefine              ! 0: FREE GRID, 1: FIXED GRID
      read (1, *) IPMWRITE             ! SAVING PARAMETER
      if (IPMWRITE .gt. 10) stop       !maximume value of writes equal to 10
      if (IPMWRITE .GT. 0) then
         do i = 1, IPMWRITE            !max value 10
            read (1, *) IPMWSTART(i), IPMWSTEPS(i) ! START AND FREQ OF WRITING
         end do
      end if
      close (1)
      if (my_rank .eq. 0) then
         print *, 'Inputs read:'
         print *, achar(9), 'DXpm=', DXpm
         print *, achar(9), 'DYpm=', DYpm
         print *, achar(9), 'DZpm=', DZpm
      end if
   else
      DXpm = 0.1_dp
      DYpm = 0.1_dp
      DZpm = 0.1_dp
 
      NBI = 1
      NBJ = 2
      NBK = 3

      interf_iproj = 4
      ibctyp = 2
      IDVPM = 1
      ncoarse = 8
      ncell_rem = 1
      iyntree = 1
      ilevmax = 1
      OMPTHREADS = 2
      nremesh = 1
      idefine = 0

      IPMWRITE = 1; IPMWSTART(1) = 0; IPMWSTEPS(1) = 6000
   endif
   VERBOCITY = 0
   !--- END READ SETTINGS
   
   !--- Problem settings
   NI_in = -0.1_dp        ! Viscosity
   DT_in =  0.01_dp       ! =dx/U
   neq = 3                ! NUMBER OF EQUATIONS
   UINF = 0               ! INFLOW VELOCITY

   NVR = 100
   NVR_ext = NVR
   NVR_size = NVR_ext
   allocate (XPR(3, NVR_ext))
   allocate (QPR(4, NVR_ext))
   allocate (UPR(3, NVR_ext))
   allocate (GPR(3, NVR_ext))

   XPR = 0
   QPR = 0
   UPR = 0
   GPR = 0

   XPR(1:3, 1) = [-2*sphere_radius, -2*sphere_radius, -2*sphere_radius]
   XPR(1:3, 2) = [ 2*sphere_radius,  2*sphere_radius,  2*sphere_radius]
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !--- INITIALIZATION VPM
   call vpm(XPR, QPR, UPR, GPR, NVR_ext, neq, 0, RHS_ptr, vel_ptr, 0, NI_in, NVR_size)
   !--- END INITIALIZATION VPM
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!y
   tab_level = 0
   if (my_rank .eq. 0) then
      write (dummy_string, "(A)") 'Hill Vortex Initialization'
      call vpm_print(dummy_string, red, 1)
   end if
   allocate (RHS_ptr(neq, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
   call hill_assign(RHS_ptr,                                                               &
                    fine_grid%NN, fine_grid%NN_bl, fine_grid%Xbound, fine_grid%Dpm, neq,     &
                    sphere_radius, u_free_stream, sphere_z_0                                 )
   
                    
   if (my_rank.eq.0) then                    
      divergence_HILL = divergence(RHS_ptr, fine_grid%Dpm(1), fine_grid%Dpm(2), fine_grid%Dpm(3))
      write (*, *) "MAX DIV HILL", maxval(abs(divergence_HILL))
      write (*, *) "MIN DIV HILL", minval(abs(divergence_HILL))
      write (*, *) "MEAN DIV HILL", sum(abs(divergence_HILL))
   endif
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)

   call set_RHS_pm(RHS_ptr,size(RHS_ptr,1), size(RHS_ptr,2), size(RHS_ptr,3), size(RHS_ptr,4))
   !------------ Remeshing ----------------
   ! We remesh the particles in order to properly distribute them in the domain
   call remesh_particles_3d(-1,ncell_rem, XPR, QPR, GPR, UPR, NVR_EXT)
   ! Reinitalize the domain
   call vpm(XPR, QPR, UPR, GPR, NVR_ext, neq, 0, RHS_ptr, vel_ptr, 1, NI_in, NVR_size)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !------------ End Remeshing ----------------

   ! ------- Allocate memory for all particles
   if (my_rank .eq. 0) then
      allocate (QPO(1:neq + 1, NVR_ext))     ! OUTPUT STRENGTH
      allocate (XPO(1:3, NVR_ext))           ! OUTPUT POSITION
      allocate (QP_in(1:neq + 1, NVR_ext))   ! INPUT STRENGTH
      allocate (XP_in(1:3, NVR_ext))         ! INPUT POSITION
      allocate (Qflag(NVR_ext))
      Qflag = 0 ! 

      QPO = QPR
      XPO = XPR
      
      XMIN = XMIN_pm + (NXs_fine_bl + 4 - 1)*DXpm
      XMAX = maxval(XPO(1, :))
      XPO(1, :) = XPO(1, :) - (XMAX - XMIN)
      
      QP_in = QPO
      XP_in = XPO
      NVR = NVR_ext
      NVR_size = NVR_ext
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! ------- End Allocate memory for all particles

   !--- MAIN LOOP
   T = 0
   max_iter = 1000
   do i = 1, max_iter
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (my_rank .eq. 0) then
         write (*, *)
         write (*, *) "-----------------------------------------------------------------------------"
         write (*, *) achar(27)//'[1;92mITERATION= ', i, ' of', max_iter, achar(27)//'[0m'
         write (*, *) achar(27)//'[1;92mT=', T, achar(27)//'[0m'
         write (*, *) achar(27)//'[1;92mDT=', DT_in, achar(27)//'[0m'
         write (*, *) "-----------------------------------------------------------------------------"
      end if
      !get velocities and deformations
      T = T +  DT_in

      !--- ALLOCATIONS FOR ALL PARTICLES and sources
      if (my_rank .eq. 0) then
         if (allocated(UPR)) deallocate(UPR)
         if (allocated(GPR)) deallocate(GPR)
         allocate (UPR(3, NVR_EXT))
         allocate (GPR(3, NVR_EXT))
         UPR = 0; GPR = 0
      end if

      !--- VPM GETS VELOCITIES AND DEFORMATIONS FROM THE PM SOLUTION
      call vpm(XPR, QPR, UPR, GPR, NVR_EXT, neq, 2, RHS_ptr, vel_ptr, i, NI_in, NVR_size)
      tab_level = 0
      if (my_rank .eq. 0) then
         write (*, "(A)") ''
         write (*, "(A)") 'Velocity Field'
         write (*, "(A)") '---------------------------------'
         write (*, "(A)") 'X-Component'
         write (*, "(A, F8.3)") achar(9)//'min', minval(vel_ptr(1, :, :, :))
         write (*, "(A, F8.3)") achar(9)//'max', maxval(vel_ptr(1, :, :, :))
         write (*, "(A, F8.3)") achar(9)//'mean', sum(vel_ptr(1, :, :, :))/size(vel_ptr(1, :, :, :))
         write (*, "(A)") '---------------------------------'
         write (*, "(A)") 'Y-Component'
         write (*, "(A, F8.3)") achar(9)//'min', minval(vel_ptr(2, :, :, :))
         write (*, "(A, F8.3)") achar(9)//'max', maxval(vel_ptr(2, :, :, :))
         write (*, "(A, F8.3)") achar(9)//'mean', sum(vel_ptr(2, :, :, :))/size(vel_ptr(2, :, :, :))
         write (*, "(A)") '---------------------------------'
         write (*, "(A)") 'Z-Component'
         write (*, "(A, F8.3)") achar(9)//'min', minval(vel_ptr(3, :, :, :))
         write (*, "(A, F8.3)") achar(9)//'max', maxval(vel_ptr(3, :, :, :))
         write (*, "(A, F8.3)") achar(9)//'mean', sum(vel_ptr(3, :, :, :))/size(vel_ptr(3, :, :, :))
         write (*, "(A)") '---------------------------------'
         CFL_x = maxval(abs(vel_ptr(1, :, :, :)))*DT_in/DXpm
         CFL_y = maxval(abs(vel_ptr(2, :, :, :)))*DT_in/DYpm
         CFL_z = maxval(abs(vel_ptr(3, :, :, :)))*DT_in/DZpm
         CFL = CFL_x + CFL_y + CFL_z
         write (*, "(A, F8.3,A ,F8.3, A)")  achar(27)//'[1;33mCFL Criterion = ', CFL, &
                                            achar(9)//'with DT=', DT_IN, achar(27)//'[0m'
         write (*, "(A,F8.3, A, F8.3, A , F8.3)") &
                           achar(9)//'X-axis min', CFL_x,                                       &
                           achar(9)//'max', maxval(abs(vel_ptr(1, :, :, :)))*DT_in/DXpm, &
                           achar(9)//'mean', sum(abs(vel_ptr(1, :, :, :)))*DT_in/DXpm/size(vel_ptr(1, :, :, :))
         write (*, "(A,F8.3, A, F8.3, A , F8.3)") &
                           achar(9)//'Y-axis min', CFL_y,                                       &
                           achar(9)//'max', maxval(abs(vel_ptr(1, :, :, :)))*DT_in/DXpm, &
                           achar(9)//'mean', sum(abs(vel_ptr(1, :, :, :)))*DT_in/DXpm/size(vel_ptr(1, :, :, :))
         write (*, "(A,F8.3, A, F8.3, A , F8.3)") &
                           achar(9)//'Z-axis min', CFL_z,                                       &
                           achar(9)//'max', maxval(abs(vel_ptr(1, :, :, :)))*DT_in/DXpm, &
                           achar(9)//'mean', sum(abs(vel_ptr(1, :, :, :)))*DT_in/DXpm/size(vel_ptr(1, :, :, :))
         
         ! ! IF CFL > 0.9 adjust the time step so that the CFL is 0.9
         ! if (CFL .gt. CFL_treshold_up) then
         !    write (*, "(A)") achar(27)//'[1;31mCFL CRITERION EXCEEDED', achar(27)//'[0m'
         !    write (*, "(A,F8.3)") achar(9)//'Old Time Step = ', DT_in
         !    DT_in = (CFL_target/CFL) * DT_in
         !    write (*, "(A,F8.3)") achar(9)//'New Time Step = ', DT_in
         !    write (*, "(A,F8.3)") achar(9)//'New CFL = ', CFL_target
         ! end if
         ! ! IF CFL < 0.5 adjust the time step so that the CFL is 0.5
         ! if (CFL .lt. CFL_treshold_down) then
         !    DT_in = (CFL_target/CFL) * DT_in
         !    write (*, "(A)") achar(27)//'[1;91mCFL CRITERION EXCEEDED', achar(27)//'[0m'
         !    write (*, "(A,F8.3)") achar(9)//'New Time Step = ', DT_in
         !    write (*, "(A,F8.3)") achar(9)//'Old Time Step = ', DT_in
         !    write (*, "(A,F8.3)") achar(9)//'New CFL = ', CFL_target
         ! end if
         ! write (*, "(A)") ''

         !WRITE TO INFOFILE
         if (i.eq.1) then
            write (cfl_file, '(A)') trim(case_folder)//'cfl.csv'
            open (unit=10, file=cfl_file, status='replace', action='write')
            write (10, "(A)") 'Iteration, DT, DX, DY, DZ, CFL, CFL_x, CFL_y, CFL_z'
         else 
            open(unit=10, file=cfl_file, status='old', position='append', action='write')
         end if
         write (10, '(I0, 8(",",ES14.7))') i, DT_in, DXpm, DYpm, DZpm, CFL, CFL_x, CFL_y, CFL_z
         close (10)
      end if
      !--- END VPM GETS VELOCITIES AND DEFORMATIONS FROM THE PM SOLUTION
      
      ! CONVECTING PARTICLES
      if (my_rank .eq. 0) then
         ! Write particles to file
         write (dummy_string, "(A)") "Writing Particles and PM"
         call vpm_print(dummy_string, red, 1)
         st = MPI_WTIME()
         call write_particles_hdf5(i, XPR, UPR, QPR, GPR, neq, NVR, NVR_ext)
         call write_pm_solution_hdf5(i, fine_grid%NN, fine_grid%NN_bl, neq, RHS_ptr, SOL_pm, velocity_pm, deform_pm)
         et = MPI_WTIME()
         write (dummy_string, "(A,I3,A,F8.3,A)") &
               achar(9)//'finished in:', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
         call vpm_print(dummy_string, yellow, 1)

         write (dummy_string, "(A)") "Convection of Particles"
         call vpm_print(dummy_string, red, 1)
         st = MPI_WTIME()
         !$omp parallel private(j,FACDEF)
         !$omp shared(XPR,UPR,QPR,GPR)
         !$omp do
         do j = 1, NVR
            ! Move particles
            XPR(1:3, j) = XPR(1:3, j) + (UPR(1:3, j) + UINF(1:3))*DT_in

            ! Vortex Stretching 
            FACDEF = - 1._dp
            !OMET   = sqrt ( QPR(1,j)**2 + QPR(2,j)**2 + QPR(3,j)**2 )
            !OG     = sqrt ( GPR(1,j)**2 + GPR(2,j)**2 + GPR(3,j)**2 )
            ! if (OG.ne.0.) then                                                      !RMETM
            !    if (OMET.gt.0.001)  FACDEF = OMET*MIN(RMETM,DT_in*OG/OMET)/OG/DT_in  !RMETM
            ! endif                                                                   !RMETM
            QPR(1:3, j) = QPR(1:3, j) - FACDEF * GPR(1:3, j)*DT_in
         end do
         !$omp enddo
         !$omp end parallel
         et = MPI_WTIME()
         write (dummy_string, "(A,I3,A,F8.3,A)") &
               achar(9)//'finished in:', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
         call vpm_print(dummy_string, yellow, 1)
         
         write (dummy_string, "(A)") "Moving Particles in/out"
         call vpm_print(dummy_string, red, 1)
         st = MPI_WTIME()
         ! call find_par_in(T, UINF(1))
         ! call find_par_out
         et = MPI_WTIME()
         write (dummy_string, "(A,I3,A,F8.3,A)") &
               achar(9)//'finished in:', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
         call vpm_print(dummy_string, yellow, 1)
      end if
      ! END CONVECTING PARTICLES

      !--- VPM INITIALIZATION AND REMESHING
      ! Used to redefine the sizes
      call vpm(XPR, QPR, UPR, GPR, NVR_EXT, neq, 0, RHS_ptr, vel_ptr, i, NI_in, NVR_size)
      tab_level = 0

      if (mod(i, 1).eq.0) then         
         ! Remesh the particles
         ! call remesh_particles_3d(1, ncell_rem, XPR, QPR, GPR, UPR, NVR_EXT)
         ! ! BCAST NVR_EXT
         call MPI_BCAST(NVR_EXT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         NVR_size = NVR_EXT
      endif
      !--- END VPM INITIALIZATION AND REMESHING
      
      !--- VPM DIFFUSION
      ! call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,5,RHS_pm_in,vel,i,NI_in,NVR_ext)
      !if (my_rank.eq.0) then
      !    do j= 1,NVR_ext
      !        QPR(1:3,j) = QPR(1:3,j)  -GPR(1:3,j) * DT_in
      !    enddo
      !endif
   end do
   !--- END MAIN LOOP
   call MPI_FINALIZE(ierr)
end Program test_pm

subroutine find_par_out
   use vpm_types, only: dp
   use pmgrid, only: DXpm, DYpm, DZpm, XMIN_pm, YMIN_pm, ZMIN_pm, &
                     NXf_fine_bl, NXs_fine_bl, NYf_fine_bl, NYs_fine_bl, NZf_fine_bl, NZs_fine_bl, &
                     NXpm_fine, NYpm_fine, NZpm_fine
   use vpm_vars, only: interf_iproj, neqpm, mrem
   use test_mod, only: NVR_ext, XPR, QPR

   implicit none
   integer  :: neq
   integer  :: i, NVR_in, NVR_out, NVR_out_max
   integer, allocatable :: NVR_projout(:)
   real(dp)             :: XMAX, XMIN, YMAX, YMIN, ZMIN, ZMAX, EPSX, EPSY, EPSZ
   real(dp)             :: xc, yc, zc, dx, dy, dz
   real(dp), allocatable :: XP_out(:, :), QP_out(:, :), XP_tmp(:, :), QP_tmp(:, :)
   ! character *50                ::filout

   !an eps for arithmetic reasons
   EPSX = 0.01*DXpm
   EPSY = 0.01*DYpm
   EPSZ = 0.01*DZpm
   XMAX = XMIN_pm + (NXf_fine_bl - (interf_iproj/2 + 1) - 1)*DXpm - EPSX
   XMIN = XMIN_pm + (NXs_fine_bl + (interf_iproj/2 + 1) - 1)*DXpm + EPSX
   YMAX = YMIN_pm + (NYf_fine_bl - (interf_iproj/2 + 1) - 1)*DYpm - EPSY
   YMIN = YMIN_pm + (NYs_fine_bl + (interf_iproj/2 + 1) - 1)*DYpm + EPSY
   ZMAX = ZMIN_pm + (NZf_fine_bl - (interf_iproj/2 + 1) - 1)*DZpm - EPSZ
   ZMIN = ZMIN_pm + (NZs_fine_bl + (interf_iproj/2 + 1) - 1)*DZpm + EPSZ
   xc = 0.5d0*(XMAX + XMIN)
   yc = 0.5d0*(YMAX + YMIN)
   zc = 0.5d0*(ZMAX + ZMIN)
   dx = XMAX - XMIN
   dy = YMAX - YMIN
   dz = ZMAX - ZMIN
   NVR_out = 0

   ! do i=1,NVR_ext
   !    if (XPR(1,i).lt.XMIN.or.XPR(1,i).gt.XMAX) then
   !         XPR(1,i) = XPR(1,i)  - int( (XPR(1,i) - xc)/(0.5d0*dx))*dx
   !       !NVR_out = NVR_out + 1
   !       !XP(1,i) = XPR(1,i)  -2*(XPR(1,i)-xc)
   !    endif
   !    if (XPR(2,i).lt.YMIN.or.XPR(2,i).gt.YMAX) then
   !        XPR(2,i) = XPR(2,i)  - int( (XPR(2,i) - yc)/(0.5d0*dy))*dy
   !       !XPR(2,i) = XPR(2,i)  -2*(XPR(2,i)-yc)
   !       !NVR_out = NVR_out + 1
   !    endif
   !    if (XPR(3,i).lt.ZMIN.or.XPR(3,i).gt.ZMAX) then
   !        XPR(3,i) = XPR(3,i)  - int( (XPR(3,i) - zc)/(0.5d0*dz))*dz
   !       !XPR(3,i) = XPR(3,i)  -2*(XPR(3,i)-zc)
   !      !NVR_out = NVR_out + 1
   !    endif
   ! enddo

   neq = neqpm+1
   NVR_out_max = (2*NXpm_fine*NYpm_fine + 2*NYpm_fine*NZpm_fine + 2*NZpm_fine*NXpm_fine)*3*mrem**2
   allocate (XP_out(1:3, NVR_out_max), QP_out(1:neq, NVR_out_max), NVR_projout(NVR_out_max))
   allocate (XP_tmp(1:3, NVR_ext), QP_tmp(1:neq, NVR_ext))
   NVR_projout = 2!interf_iproj
   NVR_out = 0
   NVR_in = 0
   do i = 1, NVR_ext
      if (XPR(2, i) .lt. YMIN .or. XPR(2, i) .gt. YMAX &
            .or. XPR(3, i) .lt. ZMIN .or. XPR(3, i) .gt. ZMAX) then
         NVR_out = NVR_out + 1
         XP_out(1:3, NVR_out) = XPR(1:3, i)
         QP_out(1:neq, NVR_out) = QPR(1:neq, i)
      else
         NVR_in = NVR_in + 1
         XP_tmp(1:3, NVR_in) = XPR(1:3, i)
         QP_tmp(1:neq, NVR_in) = QPR(1:neq, i)
      end if
   end do
   !if (NVR_out.eq.0) RHS_pm_out=0.d0
   write (*, *) achar(9), 'Particles out', NVR_out
   deallocate (XPR, QPR)
   allocate (XPR(3, NVR_in), QPR(neq, NVR_in))
   XPR(1:3, 1:NVR_in) = XP_tmp(1:3, 1:NVR_in); QPR(1:neq, 1:NVR_in) = QP_tmp(1:neq, 1:NVR_in)
   NVR_ext = NVR_in
   !deallocate(XP_tmp,QP_tmp)
   ! allocate(ieq(neqpm+1),QINF(neqpm+1))
   ! QINF=0.d0
   ! do i=1,neqpm+1
   !   ieq(i)=i
   ! enddo
   ! call project_particles_3D(RHS_pm_out,QP_out,XP_out,NVR_projout,NVR_out,neqpm+1,ieq,neqpm+1,QINF,NVR_out_max)
   !RHS_pm_out(neqpm+1,:,:,:)=DVpm
   !write(*,*) maxval(abs(QP_out(1,:))),maxval(abs(QPR(1,:)))
   !call project_vol3d(RHS_pm_out,neqpm+1,ieq,neqpm+1,1)
   !write(*,*) 'out',NVR_out,NVR_in,maxval(abs(RHS_pm_out(1,:,:,:)))
   deallocate (XP_out, QP_out)
end subroutine find_par_out

subroutine find_par_in(T_in, U)
   use vpm_types, only: dp
   use vpm_vars, only: neqpm, mrem
   use test_mod, only: Qflag, QPO, QP_in, XPO, XP_in, NVR_ext, XPR, QPR, NVR_ext_init
   use pmgrid, only: XMIN_pm, NXs_fine_bl, DXpm, NXpm_fine, NYpm_fine, NZpm_fine

   !
   ! This subroutine moves the particles in and out of the domain
   ! It also resets the inflow particles
   !

   implicit none
   real(dp), intent(in)    :: T_in, U
   real(dp)                :: XO, XMIN
   integer                 :: NVR_in, neq
   integer                 :: i, NVR_in_max
   real(dp), allocatable   :: XP_tmp(:, :), QP_tmp(:, :)

   ! Find the inflow particles
   XMIN = XMIN_pm + (NXs_fine_bl + 4 - 1)*DXpm
   XO = XMIN
   
   ! Q flag is used to determine if a particle is inflow or not
   if (minval(Qflag) .eq. 1) then
      Qflag = 0
      QPO = QP_in; XPO = XP_in
      write (155, *) 'Reseting inflow Particles'
   end if

   ! Real number of equations
   neq = neqpm+1

   ! Max number of particles that can enter the domain
   NVR_in_max = (2*NXpm_fine*NYpm_fine + 2*NYpm_fine*NZpm_fine + 2*NZpm_fine*NXpm_fine)*3*mrem**2
   allocate (XP_tmp(3, NVR_ext + NVR_in_max), QP_tmp(neq, NVR_ext + NVR_in_max))
   
   XP_tmp(1:3, 1:NVR_ext) = XPR(1:3, 1:NVR_ext)
   QP_tmp(1:neq, 1:NVR_ext) = QPR(1:neq, 1:NVR_ext)
   NVR_in = 0
   do i = 1, NVR_ext_init
      XPO(1, i) = XPO(1, i) + U*T_in
      if (XPO(1, i) .gt. XO .and. qflag(i) .eq. 0) then
         NVR_in = NVR_in + 1
         Qflag(i) = 1
         XP_tmp(1:3, NVR_ext + NVR_in) = XPO(1:3, i)
         QP_tmp(1:neq, NVR_ext + NVR_in) = QPO(1:neq, i)
      end if
   end do

   deallocate (XPR, QPR)
   NVR_ext = NVR_ext + NVR_in
   allocate (XPR(1:3, NVR_ext), QPR(1:neq, NVR_ext))
   XPR = XP_tmp
   QPR = QP_tmp
   deallocate (XP_tmp, QP_tmp)
   write (*, *) achar(9), 'Particles in', NVR_in

   ! ADD SOURCES
   ! allocate (XP_tmp(3, NVR_sources + NVR_in_max))
   ! allocate (QP_tmp(neq, NVR_sources + NVR_in_max))
   ! XP_tmp(1:3, 1:NVR_sources) = XSOUR(1:3, 1:NVR_sources)
   ! QP_tmp(1:neq, 1:NVR_sources) = QSOUR(1:neq, 1:NVR_sources)
   NVR_in = 0
   ! do i = 1, NVR_sources_init
   !    XPO(1, NVR_ext_init + i) = XPO(1, NVR_ext_init + i) + U*T_in
   !    if (XPO(1, NVR_ext_init + i) .gt. XO .and. qflag(NVR_ext_init + i) .eq. 0) then
   !       NVR_in = NVR_in + 1
   !       Qflag(NVR_ext_init + i) = 1
   !       XP_tmp(1:3, NVR_sources + NVR_in) = XPO(1:3, NVR_ext_init + i)
   !       QP_tmp(1:neq, NVR_sources + NVR_in) = QPO(1:neq, NVR_ext_init + i)
   !    end if
   ! end do
   ! deallocate (XSOUR, QSOUR)
   ! NVR_sources = NVR_sources + NVR_in
   ! allocate (XSOUR(1:3, NVR_sources), QSOUR(1:neq, NVR_sources))
   ! XSOUR = XP_tmp
   ! QSOUR = QP_tmp
   ! deallocate (XP_tmp, QP_tmp)

   write (*, *) achar(9), 'Sources in', NVR_in
end subroutine find_par_in

subroutine move_par_out(DT)
   use vpm_types, only: dp
   use vpm_vars, only: interf_iproj
   use test_mod, only: NVR_ext, XPR, UPR
   use pmgrid, only: DXpm, DYpm, DZpm, XMIN_pm, NXf_fine_bl, NXs_fine_bl, &
                     YMIN_pm, NYf_fine_bl, NYs_fine_bl, NZf_fine_bl, NZs_fine_bl, ZMIN_pm

   ! subroutine to move particles out of the domain
   ! It moves the particles out of the domain if they are outside the domain
   
   implicit none
   real(dp), intent(in)::DT
   integer  :: i
   real(dp)             :: XMAX, XMIN, YMAX, YMIN, ZMIN, ZMAX, EPSX, EPSY, EPSZ
   real(dp)             :: xc, yc, zc, dx, dy, dz
   ! character *50                ::filout

   !an eps for arithmetic reasons
   EPSX = 0.01*DXpm
   EPSY = 0.01*DYpm
   EPSZ = 0.01*DZpm
   XMAX = XMIN_pm + (NXf_fine_bl - (interf_iproj/2 + 1) - 1)*DXpm - EPSX
   XMIN = XMIN_pm + (NXs_fine_bl + (interf_iproj/2 + 1) - 1)*DXpm + EPSX
   YMAX = YMIN_pm + (NYf_fine_bl - (interf_iproj/2 + 1) - 1)*DYpm - EPSY
   YMIN = YMIN_pm + (NYs_fine_bl + (interf_iproj/2 + 1) - 1)*DYpm + EPSY
   ZMAX = ZMIN_pm + (NZf_fine_bl - (interf_iproj/2 + 1) - 1)*DZpm - EPSZ
   ZMIN = ZMIN_pm + (NZs_fine_bl + (interf_iproj/2 + 1) - 1)*DZpm + EPSZ
   xc = 0.5d0*(XMAX + XMIN)
   yc = 0.5d0*(YMAX + YMIN)
   zc = 0.5d0*(ZMAX + ZMIN)
   dx = XMAX - XMIN
   dy = YMAX - YMIN
   dz = ZMAX - ZMIN

   do i = 1, NVR_ext
      if (XPR(2, i) .lt. YMIN .or. XPR(2, i) .gt. YMAX) then
         XPR(2, i) = XPR(2, i) - UPR(2, i)*DT

      end if
      if (XPR(3, i) .lt. ZMIN .or. XPR(3, i) .gt. ZMAX) then
         XPR(3, i) = XPR(3, i) - UPR(3, i)*DT
      end if
   end do

   return
end subroutine move_par_out

! subroutine create_sources(UINF)
!    use vpm_size, only: Dpm
!    use vpm_vars, only: interf_iproj, neqpm
!    use test_mod, only: NVR_sources, XSOUR, QSOUR
!    use pmgrid, only: NXf_bl, NXs_bl, NYf_bl, NYs_bl, NZf_bl, NZs_bl, &
!                      DXpm, DYpm, DZpm, XMIN_pm, YMIN_pm, ZMIN_pm

!    implicit none
!    real(dp), intent(in) :: UINF(3)
!    real(dp) :: X(8), Y(8), Z(8)
!    integer :: nxfin, nxstart, nyfin, nystart, nzfin, nzstart, NXpm1, NYpm1, NZpm1
!    integer :: i, j, k, npar
!    nxfin = NXf_bl - interf_iproj/2
!    nxstart = NXs_bl + interf_iproj/2
!    nyfin = NYf_bl - interf_iproj/2
!    nystart = NYs_bl + interf_iproj/2
!    nzfin = NZf_bl - interf_iproj/2
!    nzstart = NZs_bl + interf_iproj/2
!    NXpm1 = nxfin - nxstart
!    NYpm1 = nyfin - nystart
!    NZpm1 = nzfin - nzstart
!    NVR_sources = NXpm1*NYpm1*2 + NXpm1*NZpm1*2
!    Dpm(1) = DXpm; Dpm(2) = DYpm; Dpm(3) = DZpm
!    allocate (XSOUR(3, NVR_sources), QSOUR(neqpm+1, NVR_sources))
!    XSOUR = 0
!    QSOUR = 0
!    npar = 0
!    j = nystart
!    do k = nzstart, nzfin - 1
!       do i = nxstart, nxfin - 1
!          npar = npar + 1
!          X(1) = XMIN_pm + Dpm(1)*(i - 1)
!          X(2) = XMIN_pm + Dpm(1)*(i)

!          Y(1) = YMIN_pm + Dpm(2)*(j - 1)

!          Z(1) = ZMIN_pm + Dpm(3)*(k - 1)
!          Z(2) = ZMIN_pm + Dpm(3)*(k)

!          XSOUR(1, npar) = 0.5d0*(X(1) + X(2))
!          XSOUR(2, npar) = Y(1)
!          XSOUR(3, npar) = 0.5d0*(Z(1) + Z(2))

!          QSOUR(3, npar) = -UINF(1)*Dpm(1)*Dpm(3)
!       end do
!    end do

!    j = nyfin
!    do k = nzstart, nzfin - 1
!       do i = nxstart, nxfin - 1
!          npar = npar + 1
!          X(1) = XMIN_pm + Dpm(1)*(i - 1)
!          X(2) = XMIN_pm + Dpm(1)*(i)

!          Y(1) = YMIN_pm + Dpm(2)*(j - 1)

!          Z(1) = ZMIN_pm + Dpm(3)*(k - 1)
!          Z(2) = ZMIN_pm + Dpm(3)*(k)

!          XSOUR(1, npar) = 0.5d0*(X(1) + X(2))
!          XSOUR(2, npar) = Y(1)
!          XSOUR(3, npar) = 0.5d0*(Z(1) + Z(2))

!          QSOUR(3, npar) = UINF(1)*Dpm(1)*Dpm(3)
!       end do
!    end do

!    k = nzstart
!    do j = nystart, nyfin - 1
!       do i = nxstart, nxfin - 1
!          npar = npar + 1
!          X(1) = XMIN_pm + Dpm(1)*(i - 1)
!          X(2) = XMIN_pm + Dpm(1)*(i)

!          Y(1) = YMIN_pm + Dpm(2)*(j - 1)
!          Y(2) = YMIN_pm + Dpm(2)*(j)

!          Z(1) = ZMIN_pm + Dpm(3)*(k - 1)

!          XSOUR(1, npar) = 0.5d0*(X(1) + X(2))
!          XSOUR(2, npar) = 0.5d0*(Y(1) + Y(2))
!          XSOUR(3, npar) = Z(1)

!          QSOUR(2, npar) = UINF(1)*Dpm(1)*Dpm(2)
!       end do
!    end do

!    k = nzfin
!    do j = nystart, nyfin - 1
!       do i = nxstart, nxfin - 1
!          npar = npar + 1
!          X(1) = XMIN_pm + Dpm(1)*(i - 1)
!          X(2) = XMIN_pm + Dpm(1)*(i)

!          Y(1) = YMIN_pm + Dpm(2)*(j - 1)
!          Y(2) = YMIN_pm + Dpm(2)*(j)

!          Z(1) = ZMIN_pm + Dpm(3)*(k - 1)

!          XSOUR(1, npar) = 0.5d0*(X(1) + X(2))
!          XSOUR(2, npar) = 0.5d0*(Y(1) + Y(2))
!          XSOUR(3, npar) = Z(1)

!          QSOUR(2, npar) = -UINF(1)*Dpm(1)*Dpm(2)
!       end do
!    end do

!    QSOUR = -QSOUR
! end subroutine create_sources


Module api
   use vpm_lib
   use vpm_vars
   use vpm_size
   use openmpth
   use, intrinsic :: iso_c_binding, only: c_float, c_int, c_bool, c_null_ptr, c_double

   !private ::  starttime,endtime,st,et,ct
   !private ::  nb_i,nb_j,nb_k,NBB,NXbl,NYbl,NZbl,BLOCKS,NXB,NYB,NZB,ndumcell_coarse ,ndumcell_bl
   !private :: II,iynbc,iret,NBI,NBJ,NBK,NVR_out_thres,NREMESH,ntorder,&
   !                                   iyntree,ilevmax,itree,nsize_out,ibctyp

contains
   subroutine initialize(dx_pm, dy_pm, dz_pm, proj_type, bc_type, vol_type, eps_vol, &
                         num_coarse, num_nbi, num_nbj, num_nbk, remesh_type, num_remesh_cells, tree_type, &
                         max_level, omp_threads, grid_define, slice_type, write_type, write_start, write_steps) &
      bind(C, name='init')

      use pmgrid, only: DXpm, DYpm, DZpm, EPSVOL
      use vpm_vars, only: interf_iproj, ncoarse, ncell_rem, iynslice, IPMWRITE, idefine, IPMWSTART, IPMWSTEPS
      use vpm_size, only: NREMESH, iyntree, ilevmax, ibctyp, NBI, NBJ, NBK
      use pmeshpar, only: IDVPM
      use openmpth, only: OMPTHREADS
      use MPI

      ! Declare the parameters to be passed in
      implicit none
      real(c_double), intent(in) :: dx_pm, dy_pm, dz_pm
      integer(c_int), intent(in) :: proj_type, bc_type, vol_type, num_coarse
      integer(c_int), intent(in) :: num_nbi, num_nbj, num_nbk
      integer(c_int), intent(in) :: remesh_type, num_remesh_cells, tree_type, max_level
      integer(c_int), intent(in) :: omp_threads, grid_define, slice_type, write_type
      real(c_double), intent(in) :: eps_vol
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
      write (*, *) 'PM GRID: DXpm: ', DXpm, ' DYpm: ', DYpm, ' DZpm: ', DZpm
      EPSVOL = eps_vol

      ! VPM_VARS
      interf_iproj = proj_type
      ncoarse = num_coarse
      ncell_rem = num_remesh_cells
      iynslice = slice_type
      IPMWRITE = write_type
      idefine = grid_define
      ! Check the IPMWRITE parameter and process accordingly
      if (IPMWRITE .gt. 10) stop  ! maximum value of writes equal to 10
      if (IPMWRITE .GT. 0) then
         do i = 1, IPMWRITE ! maximum value 10
            IPMWSTART(i) = write_start(i)
            IPMWSTEPS(i) = write_steps(i)
         end do
      end if

      ! VPM_SIZE
      NREMESH = remesh_type
      iyntree = tree_type
      ilevmax = max_level
      ibctyp = bc_type
      NBI = num_nbi
      NBJ = num_nbj
      NBK = num_nbk

      ! PMESHPAR
      IDVPM = vol_type

      ! OPENMPTH
      OMPTHREADS = omp_threads

   end subroutine initialize

   Subroutine finalize() bind(C, name='finalize')
      use vpm_vars
      use vpm_size
      use pmeshpar
      use parvar
      use pmgrid
      use MPI
      use pmlib
      use projlib
      use yapslib
      use openmpth

      implicit none
      
      ! integer :: ierr
      ! print *, 'Finalizing MPI'
      ! call MPI_FINALIZE(ierr)
      ! if (ierr .ne. 0) then
      !    print *, 'Error finalizing MPI'
      !    stop
      ! end if
      ! Free the memory
      ! call free_vpm_vars()
      ! call free_vpm_size()
      ! call free_pmeshpar()
      ! call free_parvar()
      ! call free_pmgrid()
      ! call free_openmpth()
      ! call free_pmlib()
      ! call free_projlib()
      ! call free_yapslib()
   End Subroutine finalize

   Subroutine call_vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
                       RHS_pm_in, Velx, Vely, Velz, NTIME_in, NI_in, NVRM_in) bind(C, name='vpm')
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

      use vpm_vars
      use vpm_size
      use pmeshpar
      use parvar
      use pmgrid
      use MPI
      use pmlib
      use projlib
      use yapslib
      use openmpth

      ! Fortran to C bindings
      implicit none
      ! Interface for the arguments
      integer(c_int), intent(inout) :: NVR_in
      integer(c_int), intent(in)    :: neqpm_in, WhatToDo, NVRM_in, NTIME_in
      real(c_double), intent(in)    :: NI_in

      ! real(c_double), intent(inout), target :: XP_in(:,:), QP_in(:,:), UP_in(:,:), GP_in(:,:)
      real(c_double), intent(inout), target :: XP_in(3, NVR_in), QP_in(neqpm_in + 1, NVR_in)
      real(c_double), intent(inout), target :: UP_in(3, NVR_in), GP_in(3, NVR_in)
      real(c_double), intent(inout), pointer :: RHS_pm_in(:, :, :, :), velx(:, :, :), vely(:, :, :), velz(:, :, :)

      integer :: ierr, my_rank !, i, j, k, l

      ! Get the rank of the process
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      if (ierr .ne. 0) then
         print *, 'Error getting the rank of the process'
         stop
      end if

      ! print *, char(9), 'Calling VPM API. Rank: ', my_rank
      ! CALL THE VPM Subroutine
      call vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
               RHS_pm_in, velx, vely, velz, NTIME_in, NI_in, NVRM_in)
      ! print *, char(9), 'VPM API Returned. Rank: ', my_rank

   End Subroutine call_vpm

   Subroutine call_remesh_particles_3d(iflag) bind(C, name='remesh_particles_3d')
      ! use vpm_vars
      ! use vpm_size
      ! use pmeshpar
      ! use parvar
      ! use pmgrid
      ! use pmlib
      ! use projlib
      ! use yapslib
      ! use openmpth
      use vpm_lib, only: remesh_particles_3d
      use MPI

      implicit none
      integer(c_int), intent(in) :: iflag
      integer :: ierr, my_rank

      ! Print all the arguments one by one
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      if (ierr .ne. 0) then
         print *, 'Error getting the rank of the process'
         stop
      end if
      
      ! print *, 'Calling REMESH PARTICLES 3D. Entry point: API. Rank: ', my_rank
      call remesh_particles_3d(iflag)

      ! Convert the arguments back to the original type

      ! write (*, *) 'REMESH_PARTICLES_3D API FINISHED. Rank: ', my_rank
   End Subroutine call_remesh_particles_3d

   Subroutine print_pmeshpar() bind(C, name='print_pmeshpar')
      use pmeshpar, only: print_pmeshpar_info
      implicit none

      call print_pmeshpar_info()

   End Subroutine print_pmeshpar

   Subroutine print_projlib() bind(C, name='print_projlib')
      use projlib, only: print_projlib_info
      implicit none

      call print_projlib_info()

   End Subroutine print_projlib

   Subroutine print_pmgrid() bind(C, name='print_pmgrid')
      use pmgrid, only: print_pmgrid_info
      implicit none

      call print_pmgrid_info()
   End Subroutine print_pmgrid

   Subroutine print_vpm_vars() bind(C, name='print_vpm_vars')
      use vpm_vars, only: print_vpm_vars_info
      implicit none

      call print_vpm_vars_info()
   End Subroutine print_vpm_vars

   Subroutine print_vpm_size() bind(C, name='print_vpm_size')
      use vpm_size, only: print_vpm_size_info
      implicit none

      call print_vpm_size_info()
   End Subroutine print_vpm_size

   Subroutine print_parvar() bind(C, name='print_parvar')
      use parvar, only: print_parvar_info
      implicit none

      call print_parvar_info()
   End Subroutine print_parvar

End Module api

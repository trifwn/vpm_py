Module api
   use ND_Arrays
   use vpm_lib
   use vpm_vars
   use vpm_size
   use openmpth
   use, intrinsic :: iso_c_binding, only: c_float, c_int, c_bool, c_null_ptr, c_double, c_ptr

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
      use parvar, only: set_neq
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

      call set_neq(3)

   end subroutine initialize

   subroutine finalize() bind(C, name='finalize')
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
   End subroutine finalize

   subroutine call_vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
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
      real(c_double), intent(inout), pointer :: RHS_pm_in(:, :, :, :), Velx(:, :, :), Vely(:, :, :), Velz(:, :, :)
      integer :: ierr, my_rank 

      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      if (ierr .ne. 0) then
         print *, 'Error getting the rank of the process'
         stop
      end if
      call vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
               RHS_pm_in, Velx, Vely, Velz, NTIME_in, NI_in, NVRM_in)
   End subroutine call_vpm

   subroutine call_remesh_particles_3d(iflag, XP_arr, QP_arr, GP_arr, UP_arr, NVR_out) bind(C, name='remesh_particles_3d')
      use vpm_lib, only: remesh_particles_3d
      use base_types, only: dp
      use ND_Arrays

      implicit none
      integer(c_int), intent(in) :: iflag
      type(ND_Array), intent(out), target :: XP_arr, QP_arr, GP_arr, UP_arr
      integer(c_int), intent(out) :: NVR_out

      ! Local variables
      real(dp), allocatable, target :: XP_out(:, :), QP_out(:, :), GP_out(:, :), UP_out(:, :)

      call remesh_particles_3d(iflag, XP_out, QP_out, GP_out, UP_out, NVR_out)
      XP_arr = from_intrinsic(XP_out, shape(XP_out))
      QP_arr = from_intrinsic(QP_out, shape(QP_out))
      GP_arr = from_intrinsic(GP_out, shape(GP_out))
      UP_arr = from_intrinsic(UP_out, shape(UP_out))

   End subroutine call_remesh_particles_3d

   subroutine pmgrid_set_RHS_pm(RHS_pm_in,size1, size2, size3, size4) bind(C, name='set_RHS_pm')
      use pmgrid, only: set_RHS_pm
      implicit none
      integer(c_int), intent(in) :: size1, size2, size3, size4
      real(c_double), dimension(size1, size2, size3, size4), intent(in) :: RHS_pm_in
      
      call set_RHS_pm(RHS_pm_in)
      ! print *, 'RHS_PM SET'
   end subroutine pmgrid_set_RHS_pm

   subroutine get_neqpm(neqpm_out) bind(C, name='get_neqpm')
      use vpm_vars, only: neqpm
      implicit none
      integer(c_int) :: neqpm_out

      neqpm_out = neqpm
   end subroutine get_neqpm

   subroutine get_velocity_pm(velx_out, vely_out, velz_out) bind(C, name='get_velocity_pm')
      use pmgrid, only: velvrx_pm, velvry_pm, velvrz_pm
      implicit none
      real(c_double), dimension(:,:,:), pointer :: velx_out, vely_out, velz_out

      velx_out => velvrx_pm
      vely_out => velvry_pm
      velz_out => velvrz_pm
   End subroutine get_velocity_pm

   subroutine get_size_XP(size_out) bind(C, name='get_size_XP')
      use parvar, only: XP
      implicit none
      integer(c_int), dimension(2) :: size_out

      size_out = [size(XP, 1), size(XP, 2)]
   End subroutine get_size_XP

   subroutine pmgrid_get_NN(NN_out) bind(C, name='get_NN') 
      use vpm_size, only: get_NN
      implicit none
      integer(c_int), dimension(3) :: NN_out

      call get_NN(NN_out)
   End subroutine pmgrid_get_NN

   subroutine pmgrid_get_NN_bl(NN_bl_out) bind(C, name='get_NN_bl') 
      use vpm_size, only: get_NN_bl
      implicit none
      integer(c_int), dimension(6) :: NN_bl_out

      call get_NN_bl(NN_bl_out)
      
   End subroutine pmgrid_get_NN_bl

   subroutine pmgrid_get_NX_pm(NX_pm_out) bind(C, name='get_NX_pm') 
      use pmgrid, only: get_NXpm
      implicit none
      integer(c_int) :: NX_pm_out

      call get_NXpm(NX_pm_out)
   End subroutine pmgrid_get_NX_pm

   subroutine pmgrid_get_NY_pm(NY_pm_out) bind(C, name='get_NY_pm') 
      use pmgrid, only: get_NYpm
      implicit none
      integer(c_int) :: NY_pm_out

      call get_NYpm(NY_pm_out)
   End subroutine pmgrid_get_NY_pm

   subroutine pmgrid_get_NZ_pm(NZ_pm_out) bind(C, name='get_NZ_pm') 
      use pmgrid, only: get_NZpm
      implicit none
      integer(c_int) :: NZ_pm_out

      call get_NZpm(NZ_pm_out)
   End subroutine pmgrid_get_NZ_pm

   subroutine pmgrid_get_Xbound(Xbound_out) bind(C, name='get_Xbound') 
      use vpm_size, only: get_Xbound
      implicit none
      real(c_double), dimension(6) :: Xbound_out
      call get_Xbound(Xbound_out)
   End subroutine pmgrid_get_Xbound


   !! LIBRARY PRINTS
   subroutine print_pmeshpar() bind(C, name='print_pmeshpar')
      use pmeshpar, only: print_pmeshpar_info
      implicit none

      call print_pmeshpar_info()

   End subroutine print_pmeshpar

   subroutine print_projlib() bind(C, name='print_projlib')
      use projlib, only: print_projlib_info
      implicit none

      call print_projlib_info()

   End subroutine print_projlib

   subroutine print_pmgrid() bind(C, name='print_pmgrid')
      use pmgrid, only: print_pmgrid_info
      implicit none

      call print_pmgrid_info()
   End subroutine print_pmgrid

   subroutine print_vpm_vars() bind(C, name='print_vpm_vars')
      use vpm_vars, only: print_vpm_vars_info
      implicit none

      call print_vpm_vars_info()
   End subroutine print_vpm_vars

   subroutine print_vpm_size() bind(C, name='print_vpm_size')
      use vpm_size, only: print_vpm_size_info
      implicit none

      call print_vpm_size_info()
   End subroutine print_vpm_size

   subroutine print_parvar() bind(C, name='print_parvar')
      use parvar, only: print_parvar_info
      implicit none

      call print_parvar_info()
   End subroutine print_parvar

End Module api

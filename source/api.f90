module nd_array_mod
  use ISO_C_BINDING
  implicit none

  type, bind(C) :: nd_array_t
      integer :: ndims   ! Number of dimensions
      integer :: total_size  ! Total size of the data array
      type(C_PTR) :: shape_ptr  !  shape array 
      type(C_PTR) :: data_ptr   ! Pointer to the actual data
  end type nd_array_t

contains

   function create_nd_array(ndims, arr_shape) result(arr) 
      implicit none
      integer(c_int) :: ndims
      integer(c_int), allocatable, dimension(:), intent(in), target :: arr_shape
      integer, dimension(:), contiguous, pointer :: temp_data
      type(nd_array_t) :: arr
      integer :: total_size, i

      print *, 'Creating ND Array'
      print *, 'Number of dimensions:', ndims
      arr%ndims = ndims
      
      print *, 'Array shape:', arr_shape(1:ndims)
      arr%shape_ptr = C_LOC(arr_shape)  ! Get the C pointer to the shape array
      
      total_size = 1
      do i = 1, ndims
         total_size = total_size * arr_shape(i)
      end do
      print *, 'Total size:', total_size
      arr%total_size = total_size

      ! Allocate the data array
      call c_f_pointer(arr%data_ptr, temp_data, [total_size])
      temp_data = 0  ! Initialize the data array to zero

  end function create_nd_array

  subroutine destroy_nd_array(arr)
      implicit none
      type(nd_array_t), intent(inout), target :: arr
      type(C_PTR), pointer :: data_ptr
      integer, pointer :: shape_ptr(:)

      data_ptr => arr%data_ptr
      shape_ptr => arr%shape_ptr

      if (associated(data_ptr)) then
         deallocate(data_ptr)
      end if
      if (associated(shape_ptr)) then
         deallocate(shape_ptr)
      end if
  end subroutine destroy_nd_array

  function get_nd_array_data_ptr(arr) result(data_) 
      !  use ISO_C_BINDING
      implicit none
      type(nd_array_t), intent(in) :: arr
      type(C_PTR) :: data_
      
      data_ = arr%data_ptr
  end function get_nd_array_data_ptr

  subroutine ndarray_print(arr)
      use ISO_C_BINDING
      implicit none
      type(nd_array_t), intent(in) :: arr
      real(C_DOUBLE), dimension(:), pointer :: data_
      integer(C_INT), dimension(:), pointer :: shape_
      integer :: i, j  ! Adjust for higher dimensions as needed

      call c_f_pointer(arr%shape_ptr, shape_, [arr%ndims])
      call c_f_pointer(arr%data_ptr, data_, [arr%total_size])
      
      print *, 'Array in Fortran (shape:', shape_, '):'
      print *, "Total size:", arr%total_size
      
      write(*,*) "Address pointed to by data:", C_LOC(data_), "(decimal)"
      write(*,*) "Address pointed to by shape:", C_LOC(shape_), "(hexadecimal)"
      ! Implement a while loop that loops the dimensions 
      ! and prints the data in the correct order
      do i = 1, arr%ndims
         write(*,*) 'Dimension:', i, 'Size:', shape_(i)
         do j = i, i + shape_(i)
            write(*,*) data_(j)
         end do
      end do



      write (*, *) 'End of array'
      write (*, *) ''
  end subroutine ndarray_print  
end module nd_array_mod

Module api
   use nd_array_mod
   use vpm_lib
   use vpm_vars
   use vpm_size
   use openmpth
   use, intrinsic :: iso_c_binding, only: c_float, c_int, c_bool, c_null_ptr, c_double, c_ptr

   private
   real(c_double), dimension(:), allocatable, target :: test_array
   type(nd_array_t), target :: test_ndarr

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
      double precision, allocatable, target :: XP_out(:, :), QP_out(:, :)
      integer :: NVR_out
      

      ! Print all the arguments one by one
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      if (ierr .ne. 0) then
         print *, 'Error getting the rank of the process'
         stop
      end if
      
      ! print *, 'Calling REMESH PARTICLES 3D. Entry point: API. Rank: ', my_rank
      call remesh_particles_3d(iflag, XP_out, QP_out, NVR_out)

      ! Convert the arguments back to the original type

      ! write (*, *) 'REMESH_PARTICLES_3D API FINISHED. Rank: ', my_rank
   End Subroutine call_remesh_particles_3d

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

   subroutine get_nvr(NVR_out) bind(C, name='get_num_particles')
      use parvar, only: NVR
      implicit none
      integer(c_int) :: NVR_out

      NVR_out = NVR
   end subroutine get_nvr

   subroutine get_particle_positions(XP_out) bind(C, name='get_particle_positions')
      use parvar, only: XP, NVR
      use vpm_vars, only: neqpm
      implicit none
      real(c_double), dimension(3, NVR) :: XP_out
      XP_out = XP
   end subroutine get_particle_positions

   subroutine set_particle_positions(XP_in) bind(C, name='set_particle_positions')
      use parvar, only: XP, NVR
      implicit none
      real(c_double), dimension(3, NVR) :: XP_in
      XP = XP_in
   end subroutine set_particle_positions

   subroutine get_particle_strengths(QP_out) bind(C, name='get_particle_strengths')
      use parvar, only: QP, NVR
      use vpm_vars, only: neqpm
      implicit none
      real(c_double), dimension(neqpm + 1, NVR) :: QP_out
      QP_out = QP
   end subroutine get_particle_strengths

   subroutine set_particle_strengths(QP_in) bind(C, name='set_particle_strengths')
      use parvar, only: QP, NVR
      use vpm_vars, only: neqpm
      implicit none
      real(c_double), dimension(neqpm + 1, NVR) :: QP_in
      QP = QP_in
   end subroutine set_particle_strengths

   subroutine get_particle_velocities(UP_out) bind(C, name='get_particle_velocities')
      use parvar, only: UP, NVR
      use vpm_vars, only: neqpm
      implicit none
      real(c_double), dimension(3, NVR) :: UP_out
      UP_out = UP
   end subroutine get_particle_velocities

   subroutine set_particle_velocities(UP_in) bind(C, name='set_particle_velocities')
      use parvar, only: UP, NVR
      use vpm_vars, only: neqpm
      implicit none
      real(c_double), dimension(3, NVR) :: UP_in
      UP = UP_in
   end subroutine set_particle_velocities

   subroutine get_particle_deformation(GP_out) bind(C, name='get_particle_deformation')
      use parvar, only: GP, NVR
      use vpm_vars, only: neqpm
      implicit none
      real(c_double), dimension(3, NVR) :: GP_out
      GP_out = GP
   end subroutine get_particle_deformation

   subroutine set_particle_deformation(GP_in) bind(C, name='set_particle_deformation')
      use parvar, only: GP, NVR
      use vpm_vars, only: neqpm
      implicit none
      real(c_double), dimension(3, NVR) :: GP_in
      GP = GP_in
   end subroutine set_particle_deformation

   Subroutine get_velocity_pm(velx_out, vely_out, velz_out) bind(C, name='get_velocity_pm')
      use pmgrid, only: velvrx_pm, velvry_pm, velvrz_pm
      implicit none
      real(c_double), dimension(:,:,:), pointer :: velx_out, vely_out, velz_out

      velx_out => velvrx_pm
      vely_out => velvry_pm
      velz_out => velvrz_pm
   End Subroutine get_velocity_pm

   Subroutine get_size_XP(size_out) bind(C, name='get_size_XP')
      use parvar, only: XP
      implicit none
      integer(c_int), dimension(2) :: size_out

      size_out = [size(XP, 1), size(XP, 2)]
   End Subroutine get_size_XP

   Subroutine pmgrid_get_NN(NN_out) bind(C, name='get_NN') 
      use vpm_size, only: get_NN
      implicit none
      integer(c_int), dimension(3) :: NN_out

      call get_NN(NN_out)
   End Subroutine pmgrid_get_NN

   Subroutine pmgrid_get_NN_bl(NN_bl_out) bind(C, name='get_NN_bl') 
      use vpm_size, only: get_NN_bl
      implicit none
      integer(c_int), dimension(6) :: NN_bl_out

      call get_NN_bl(NN_bl_out)
      
   End Subroutine pmgrid_get_NN_bl

   Subroutine pmgrid_get_NX_pm(NX_pm_out) bind(C, name='get_NX_pm') 
      use pmgrid, only: get_NXpm
      implicit none
      integer(c_int) :: NX_pm_out

      call get_NXpm(NX_pm_out)
   End Subroutine pmgrid_get_NX_pm

   Subroutine pmgrid_get_NY_pm(NY_pm_out) bind(C, name='get_NY_pm') 
      use pmgrid, only: get_NYpm
      implicit none
      integer(c_int) :: NY_pm_out

      call get_NYpm(NY_pm_out)
   End Subroutine pmgrid_get_NY_pm

   Subroutine pmgrid_get_NZ_pm(NZ_pm_out) bind(C, name='get_NZ_pm') 
      use pmgrid, only: get_NZpm
      implicit none
      integer(c_int) :: NZ_pm_out

      call get_NZpm(NZ_pm_out)
   End Subroutine pmgrid_get_NZ_pm

   Subroutine pmgrid_get_Xbound(Xbound_out) bind(C, name='get_Xbound') 
      use vpm_size, only: get_Xbound
      implicit none
      real(c_double), dimension(6) :: Xbound_out

      call get_Xbound(Xbound_out)
   End Subroutine pmgrid_get_Xbound

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

!!!!!!!!!!!!!!!!!!!!!!!!!! TESTS
   subroutine print_array() bind(C, name='print_array')
      use ISO_C_BINDING
      implicit none
      integer :: i

      write(*,*) 'Array in Fortran:'
      do i = 1, size(test_array)
         write(*,*) test_array(i)
      end do
   end subroutine

   subroutine init_array(n) bind(C, name='init_array')
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value :: n
      integer :: i

      if (allocated(test_array)) deallocate(test_array)  ! Deallocate if already allocated
      allocate(test_array(n))

      do i = 1, n
         test_array(i) = i * 1.0 
      end do
   end subroutine

   subroutine change_array() bind(C, name='change_array')
      use ISO_C_BINDING
      implicit none
      integer :: i

      do i = 1, size(test_array)
         test_array(i) = test_array(i) * 2.0
      end do
   end subroutine

   function get_module_array() result(array) bind(C, name='get_module_array')
      use ISO_C_BINDING
      implicit none
      type(C_PTR) :: array

      array = C_LOC(test_array)  ! Get the C pointer to the module array
   end function

   subroutine free_array() bind(C, name='free_array')
      use ISO_C_BINDING
      implicit none

      if (allocated(test_array)) deallocate(test_array)
   end subroutine

   Subroutine init_ndarray(arr_shape, ndims) bind(C, name='init_ndarray')
      use ISO_C_BINDING
      use nd_array_mod
      implicit none
      integer(c_int) :: ndims
      integer(c_int), dimension(ndims) , intent(in) :: arr_shape
      test_ndarr = create_nd_array(ndims, arr_shape)
   end Subroutine init_ndarray

   Subroutine free_ndarray() bind(C, name='free_ndarray')
      use ISO_C_BINDING
      use nd_array_mod
      implicit none

      call destroy_nd_array(test_ndarr)
   end Subroutine free_ndarray

   function get_ndarray_data_ptr() result(array) bind(C, name='get_ndarray_data_ptr')
      use ISO_C_BINDING
      use nd_array_mod
      implicit none
      type(C_PTR) :: array
      real(C_DOUBLE), dimension(:), pointer :: data_
      real(C_DOUBLE), dimension(:), allocatable, target :: temp_data
      
      print *, '1 Getting ND Array Data Pointer'
      data_ = get_nd_array_data_ptr(test_ndarr)
      print *, '2 Got Something'

      temp_data = data_
      array = C_LOC(temp_data)  ! Get the C pointer to the module array

   end function get_ndarray_data_ptr

   Subroutine print_ndarray() bind(C, name='print_ndarray')
      use ISO_C_BINDING
      use nd_array_mod
      implicit none

      call ndarray_print(test_ndarr)
   end Subroutine print_ndarray

End Module api

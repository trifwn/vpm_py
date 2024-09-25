Module api
   use ND_Arrays
   use vpm_lib
   use vpm_vars
   use vpm_size
   use openmpth
   use, intrinsic :: iso_c_binding, only: c_float, c_int, c_bool, c_null_ptr, c_double, c_ptr

   integer, parameter :: MAX_STRING_LENGTH = 256
   character(kind=c_char, len=MAX_STRING_LENGTH) :: test_string = 'This is a test string'
contains
   subroutine initialize(dx_pm, dy_pm, dz_pm, proj_type, bc_type, vol_type,                                     &
                         num_coarse, num_nbi, num_nbj, num_nbk, remesh_type, tree_type,                         &
                         max_level, omp_threads, grid_define, write_type, write_start, write_steps              &
                         ) &
      bind(C, name='init')

      use pmgrid, only: DXpm, DYpm, DZpm, IDVPM, ncoarse
      use vpm_vars, only: interf_iproj, IPMWRITE, idefine, IPMWSTART, IPMWSTEPS, OMPTHREADS, &
                          nremesh, iyntree, ilevmax, ibctyp, NBI, NBJ, NBK
      use parvar, only: set_neq
      use MPI

      ! Declare the parameters to be passed in
      implicit none
      real(c_double), intent(in) :: dx_pm, dy_pm, dz_pm
      integer(c_int), intent(in) :: proj_type, bc_type, vol_type, num_coarse
      integer(c_int), intent(in) :: num_nbi, num_nbj, num_nbk
      integer(c_int), intent(in) :: remesh_type, tree_type, max_level
      integer(c_int), intent(in) :: omp_threads, grid_define, write_type
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
      IDVPM = vol_type
      ncoarse = num_coarse

      ! VPM_VARS
      interf_iproj = proj_type
      idefine = grid_define
      IPMWRITE = write_type
      ! Check the IPMWRITE parameter and process accordingly
      if (IPMWRITE .GT. 0) then
         do i = 1, IPMWRITE ! maximum value 10
            IPMWSTART(i) = write_start(i)
            IPMWSTEPS(i) = write_steps(i)
            if (IPMWRITE .gt. 10) stop  ! maximum value of writes equal to 10
         end do
      end if

      ! VPM_SIZE
      nremesh = remesh_type
      iyntree = tree_type
      ilevmax = max_level
      ibctyp = bc_type
      NBI = num_nbi
      NBJ = num_nbj
      NBK = num_nbk


      ! OPENMPTH
      OMPTHREADS = omp_threads

      call set_neq(3)

   end subroutine initialize

   subroutine set_verbose_level(verbocity_in) bind(C, name='set_verbosity')
      use console_io, only: verbocity
      implicit none
      integer(c_int), intent(in) :: verbocity_in
      verbocity = verbocity_in
   end subroutine set_verbose_level

   subroutine finalize() bind(C, name='finalize')
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
      ! call free_parvar()
      ! call free_pmgrid()
      ! call free_openmpth()
      ! call free_pmlib()
      ! call free_projlib()
      ! call free_yapslib()
   End subroutine finalize

   subroutine call_vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
                       RHS_pm_out, Velx_out, Vely_out, Velz_out, NTIME_in, NI_in,&
                       NVR_size_in, deformx_out, deformy_out, deformz_out) bind(C, name='vpm')
      !  -> XP : particle positions (3 * NVR)
      !  -> QP : particle quantities (neqpm * NVR)
      !  -> UP : particle velocities (3 * NVR)
      !  -> GP : particle defromation (3*NVR) (wmega * \nabla)  u
      !  -> NVR : number of particles
      !  -> neqpm : number of equations
      !  -> WhatToDo : 0 - initialize, 1 - solve, 2 - convect, 3 - project, 4 - back, 5 - diffuse
      !  -> \nabla^2 u_i = RHS_i (NEQ)
      !  -> Velx, Vely, Velz : velocity field at grid points
      !  -> NTIME_IN : current time
      !  -> NI_in: viscocity -> DIFFUSION OF VORTICITY
      !  -> NVRM : NVR MAX
      use vpm_lib, only: vpm
      use ND_Arrays
      use, intrinsic :: ieee_arithmetic

      ! Fortran to C bindings
      implicit none
      ! Interface for the arguments
      integer(c_int), intent(inout)          :: NVR_in
      integer(c_int), intent(in)             :: neqpm_in, WhatToDo, NVR_size_in, NTIME_in
      real(c_double), intent(in)             :: NI_in
      real(c_double), intent(inout), target  :: XP_in(3, NVR_in), QP_in(neqpm_in + 1, NVR_in)
      real(c_double), intent(inout), target  :: UP_in(3, NVR_in), GP_in(3, NVR_in)
      type(ND_Array), intent(out)            :: RHS_pm_out, Velx_out, Vely_out, Velz_out
      type(ND_Array), intent(out), optional  :: deformx_out, deformy_out, deformz_out

      real(c_double), pointer          :: RHS_pm_ptr(:, :, :, :), &
                                          Velx_ptr(:, :, :), Vely_ptr(:, :, :), Velz_ptr(:, :, :)
      real(c_double), pointer          :: deformx_ptr(:, :, :), deformy_ptr(:, :, :), deformz_ptr(:, :, :)
      integer :: ierr, my_rank

      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      if (ierr .ne. 0) then
         print *, 'Error getting the rank of the process'
         stop
      end if
      
      call vpm(XP_in, QP_in, UP_in, GP_in, NVR_in, neqpm_in, WhatToDo, &
               RHS_pm_ptr, Velx_ptr, Vely_ptr, Velz_ptr, NTIME_in, NI_in, NVR_size_in, &
               deformx_ptr, deformy_ptr, deformz_ptr)

      
      ! Copy the data back to the arrays
      ! Assign the pointers to the arrays
      if (associated(RHS_pm_ptr)) then 
         RHS_pm_out = from_intrinsic(RHS_pm_ptr, shape(RHS_pm_ptr))
      end if

      if (associated(Velx_ptr)) then
         Velx_out = from_intrinsic(Velx_ptr, shape(Velx_ptr))
      end if

      if (associated(Vely_ptr)) then
         Vely_out = from_intrinsic(Vely_ptr, shape(Vely_ptr))
      end if

      if (associated(Velz_ptr)) then
         Velz_out = from_intrinsic(Velz_ptr, shape(Velz_ptr))
      end if

   
      ! if ((associated(deformx_ptr)).and.(present(deformx_out))) then
      !    deformx_out = from_intrinsic(deformx_ptr, shape(deformx_ptr))
      ! end if

      ! if ((associated(deformy_ptr)).and.(present(deformy_out))) then
      !    deformy_out = from_intrinsic(deformy_ptr, shape(deformy_ptr))
      ! end if

      ! if ((associated(deformz_ptr)).and.(present(deformz_out))) then
      !    deformz_out = from_intrinsic(deformz_ptr, shape(deformz_ptr))
      ! end if
   End subroutine call_vpm

   subroutine call_remesh_particles_3d(iflag, npar_per_cell, XP_arr, QP_arr, GP_arr, UP_arr, NVR_out, cuttof_value) &
         bind(C, name='remesh_particles_3d')
      use vpm_lib, only: remesh_particles_3d
      use base_types, only: dp
      use ND_Arrays

      implicit none
      integer(c_int), intent(in) :: iflag, npar_per_cell
      type(ND_Array), intent(out), target :: XP_arr, QP_arr, GP_arr, UP_arr
      integer(c_int), intent(out) :: NVR_out
      real(c_double), intent(in), optional :: cuttof_value

      ! Local variables
      real(dp), allocatable, target, save :: XP_out(:, :), QP_out(:, :), GP_out(:, :), UP_out(:, :)

      if (present(cuttof_value)) then
         call remesh_particles_3d(iflag, npar_per_cell, XP_out, QP_out, GP_out, UP_out, NVR_out, cuttof_value)
      else
         call remesh_particles_3d(iflag, npar_per_cell, XP_out, QP_out, GP_out, UP_out, NVR_out)
      end if

      XP_arr = from_intrinsic(XP_out, shape(XP_out))
      QP_arr = from_intrinsic(QP_out, shape(QP_out))
      GP_arr = from_intrinsic(GP_out, shape(GP_out))
      UP_arr = from_intrinsic(UP_out, shape(UP_out))

   End subroutine call_remesh_particles_3d

   subroutine write_particle_mesh_solution(folder, filename) bind(C, name='write_particle_mesh_solution')
      use console_io, only: vpm_write_folder, pm_output_file_suffix
      use file_io, only: write_pm_solution
      use pmgrid, only: velvrx_pm, velvry_pm, velvrz_pm, deformx_pm, deformy_pm, deformz_pm, RHS_pm, SOL_pm
      use vpm_vars, only: NTIME_pm, neqpm
      use vpm_size, only: fine_grid 
      implicit none
      character(kind = c_char), intent(in), optional :: folder(*), filename(*)
      integer :: NN(3), NN_bl(6)

      if (present(folder)) then
         call set_string_f_c(vpm_write_folder, folder)
      endif

      if (present(filename)) then
         call set_string_f_c(pm_output_file_suffix, filename)
      endif

      NN = fine_grid%NN
      NN_bl = fine_grid%NN_bl
      if ((allocated(deformx_pm)).and.(allocated(deformy_pm)).and.(allocated(deformz_pm))) then
         call write_pm_solution(NTIME_pm, NN, NN_bl, neqpm, RHS_pm, SOL_pm, velvrx_pm, velvry_pm, velvrz_pm, &
                                deformx_pm, deformy_pm, deformz_pm)
      else
         call write_pm_solution(NTIME_pm, NN, NN_bl, neqpm, RHS_pm, SOL_pm, velvrx_pm, velvry_pm, velvrz_pm)
      end if
   End subroutine write_particle_mesh_solution

   subroutine write_particles_stored(folder, filename) bind(C, name='write_particles')
      use file_io, only: write_particles
      use parvar, only: XP, QP, UP, GP, NVR, NVR_size
      use vpm_vars, only: neqpm, NTIME_pm
      use console_io, only: vpm_write_folder, particle_output_file_suffix
      implicit none
      character(kind=c_char), intent(in), optional :: folder(*), filename(*)

      if (present(folder)) then
         call set_string_f_c(vpm_write_folder, folder)
      endif
      if (present(filename)) then
         call set_string_f_c(particle_output_file_suffix, filename)
      endif

      call write_particles(NTIME_pm, XP, UP, QP, GP, neqpm, NVR, NVR_size)
   End subroutine write_particles_stored

   subroutine write_particles_stored_hdf5(folder, filename) bind(C, name='write_particles_hdf5')
      use file_io, only: write_particles_hdf5
      use parvar, only: XP, QP, UP, GP, NVR, NVR_size
      use vpm_vars, only: neqpm, NTIME_pm
      use console_io, only: vpm_write_folder, particle_output_file_suffix
      implicit none
      character(kind=c_char), intent(in), optional :: folder(*), filename(*)
      if (present(folder)) then
         call set_string_f_c(vpm_write_folder, folder)
      endif
      if (present(filename)) then
         call set_string_f_c(particle_output_file_suffix, filename)
      endif
      call write_particles_hdf5(NTIME_pm, XP, UP, QP, GP, neqpm, NVR, NVR_size)
   End subroutine write_particles_stored_hdf5

   subroutine write_particle_mesh_solution_hdf5(folder, filename) bind(C, name='write_particle_mesh_solution_hdf5')
      use file_io, only: write_pm_solution_hdf5
      use pmgrid, only: velvrx_pm, velvry_pm, velvrz_pm, deformx_pm, deformy_pm, deformz_pm, RHS_pm, SOL_pm
      use vpm_vars, only: NTIME_pm, neqpm
      use vpm_size, only: fine_grid
      use console_io, only: vpm_write_folder, pm_output_file_suffix
      implicit none
      character(kind = c_char), intent(in), optional :: folder(*), filename(*)
      integer :: NN(3), NN_bl(6)

      if (present(folder)) then
         call set_string_f_c(vpm_write_folder, folder)
      endif

      if (present(filename)) then
         call set_string_f_c(pm_output_file_suffix, filename)
      endif

      NN = fine_grid%NN
      NN_bl = fine_grid%NN_bl
      if ((allocated(deformx_pm)).and.(allocated(deformy_pm)).and.(allocated(deformz_pm))) then
         call write_pm_solution_hdf5(NTIME_pm, NN, NN_bl,neqpm, RHS_pm, SOL_pm, velvrx_pm, velvry_pm, velvrz_pm, &
                                deformx_pm, deformy_pm, deformz_pm)
      else
         call write_pm_solution_hdf5(NTIME_pm, NN, NN_bl,neqpm, RHS_pm, SOL_pm, velvrx_pm, velvry_pm, velvrz_pm)
      end if
   End subroutine write_particle_mesh_solution_hdf5


   subroutine get_particle_positions(XP_out) bind(C, name='get_particle_positions')
      use ND_Arrays
      use parvar, only: XP
      implicit none
      type(ND_Array), intent(out) :: XP_out
      XP_out = from_intrinsic(XP, shape(XP))
   end subroutine get_particle_positions

   subroutine get_particle_strengths(QP_out) bind(C, name='get_particle_strengths')
      use iso_c_binding
      use ND_Arrays
      use parvar, only: QP
      implicit none
      type(ND_Array), intent(out) :: QP_out
      QP_out = from_intrinsic(QP, shape(QP))
   end subroutine get_particle_strengths

   subroutine get_particle_deformation(GP_out) bind(C, name='get_particle_deformation')
      use iso_c_binding
      use ND_Arrays
      use parvar, only: GP
      implicit none
      type(ND_Array), intent(out) :: GP_out
      GP_out = from_intrinsic(GP, shape(GP))
   end subroutine get_particle_deformation

   subroutine get_particle_velocities(UP_out) bind(C, name='get_particle_velocities')
      use iso_c_binding
      use ND_Arrays
      use parvar, only: UP
      implicit none
      type(ND_Array), intent(out) :: UP_out
      UP_out = from_intrinsic(UP, shape(UP))
   end subroutine get_particle_velocities

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

   subroutine get_deformation_pm(deformx_out, deformy_out, deformz_out) bind(C, name='get_deformation_pm')
      use pmgrid, only: deformx_pm, deformy_pm, deformz_pm
      implicit none
      real(c_double), dimension(:,:,:), pointer :: deformx_out, deformy_out, deformz_out

      deformx_out => deformx_pm
      deformy_out => deformy_pm
      deformz_out => deformz_pm
   End subroutine get_deformation_pm

!! VPM SIZE

!! LIBRARY PRINTS 
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

   !! STRING MANIPULATION
   subroutine get_string_f_c(fstring, cptr)
      use iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: cptr
      character(len=MAX_STRING_LENGTH), intent(in) :: fstring

      character(kind=c_char, len=MAX_STRING_LENGTH), target :: fstring_t
      integer :: i, trimmed_length, full_length

      fstring_t = trim(fstring) // c_null_char
      trimmed_length = len_trim(fstring_t)
      full_length = len(fstring_t)
      do i = trimmed_length + 1, 256
         fstring_t(i:i) = c_null_char
      end do
      cptr = c_loc(fstring_t)
      ! print *, "Fortran string: ", trim(fstring)
      ! print *, "Fortran string full length: ", full_length
      ! print *, "Fortran string trimmed length: ", trimmed_length
      ! print *, "Fortran string bytes: ", iachar(fstring(1:trimmed_length))
   end subroutine get_string_f_c

   subroutine set_string_f_c(fstring, new_string)
      use iso_c_binding
      implicit none
      character(kind=c_char), intent(in) :: new_string(*)
      character(len=MAX_STRING_LENGTH), intent(out) :: fstring
      integer :: i
      
      fstring = ""
      do i = 1, 256
         if (new_string(i) == c_null_char) exit
         fstring(i:i) = new_string(i)
      end do
      ! print *, "Set test_string to: ", trim(test_string)
      ! print *, "Length of test_string: ", len_trim(test_string)
   end subroutine set_string_f_c

   !!! FILE I/O !!!
   subroutine set_particle_file_suffix(suffix)  bind(C, name='set_particle_file_suffix')
      use iso_c_binding
      use console_io, only: particle_output_file_suffix
      character(kind = c_char), dimension(*), intent(in) :: suffix
      call set_string_f_c(particle_output_file_suffix, suffix)
   end subroutine set_particle_file_suffix

   subroutine set_pm_file_suffix(suffix)  bind(C, name='set_pm_file_suffix')
      use iso_c_binding
      use console_io, only: pm_output_file_suffix
      character(kind = c_char), dimension(*), intent(in) :: suffix
      call set_string_f_c(pm_output_file_suffix, suffix)
   end subroutine set_pm_file_suffix

   subroutine set_vpm_write_folder(folder)  bind(C, name='set_vpm_write_folder')
      use iso_c_binding
      use console_io, only: vpm_write_folder
      character(kind = c_char), dimension(*), intent(in) :: folder
      call set_string_f_c(vpm_write_folder, folder)
   end subroutine set_vpm_write_folder

   subroutine print_write_settings() bind(C, name = 'print_write_settings')
      use console_io, only: vpm_print, particle_output_file_suffix, pm_output_file_suffix, vpm_write_folder
      call vpm_print("Particle output file suffix: " // trim(particle_output_file_suffix), NOCOLOR, 1)
      call vpm_print("PM output file suffix: " // trim(pm_output_file_suffix), NOCOLOR, 1)
      call vpm_print("VPM write folder: " // trim(vpm_write_folder), NOCOLOR, 1)
   end subroutine print_write_settings
End Module api

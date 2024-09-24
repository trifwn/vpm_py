module pmgrid
   use base_types, only: dp
   real(dp), allocatable, target, save :: velvrx_pm(:, :, :), velvry_pm(:, :, :), velvrz_pm(:, :, :)
   real(dp), allocatable, target, save :: deformx_pm(:, :, :), deformy_pm(:, :, :), deformz_pm(:, :, :)
   real(dp), allocatable, target, save :: RHS_pm(:, :, :, :)
   real(dp), allocatable, target, save :: SOL_pm(:, :, :, :) 

   real(dp), save                      :: XMIN_pm, XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
   real(dp), save                      :: DXpm, DYpm, DZpm, DVpm
   
   integer, save                       :: NXpm_coarse, NYpm_coarse, NZpm_coarse
   integer, save                       :: NXs_coarse_bl, NYs_coarse_bl, NXf_coarse_bl,&
                                          NYf_coarse_bl, NZs_coarse_bl, NZf_coarse_bl
   integer, save                       :: NBlocks

   !! PMESH PARAMETERS
   integer                             :: ND
   integer, save                       :: nbound, ndumcell, ncoarse, IDVPM


   ! GETTERS
   public :: get_NXpm_coarse, get_NYpm_coarse, get_NZpm_coarse
   ! SETTERS
   public :: set_RHS_pm 
   ! Printers
   public :: print_pmgrid_info

contains

   subroutine associate_velocities(velx_ptr, vely_ptr, velz_ptr)
      real(dp), pointer, intent(out) :: velx_ptr(:,:,:), vely_ptr(:,:,:), velz_ptr(:,:,:)
      if (associated(velx_ptr))     nullify (velx_ptr)
      if (associated(vely_ptr))     nullify (vely_ptr)
      if (associated(velz_ptr))     nullify (velz_ptr) 
      velx_ptr => velvrx_pm
      vely_ptr => velvry_pm
      velz_ptr => velvrz_pm
   end subroutine associate_velocities

   subroutine associate_deformations(deformx_ptr, deformy_ptr, deformz_ptr)
      real(dp), pointer, intent(out) :: deformx_ptr(:,:,:), deformy_ptr(:,:,:), deformz_ptr(:,:,:)
      if (associated(deformx_ptr))  nullify (deformx_ptr)
      if (associated(deformy_ptr))  nullify (deformy_ptr)
      if (associated(deformz_ptr))  nullify (deformz_ptr) 
      deformx_ptr => deformx_pm
      deformy_ptr => deformy_pm
      deformz_ptr => deformz_pm
   end subroutine associate_deformations

   subroutine set_pm_velocities_zero
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


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! GETTERS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine get_NXpm_coarse(NX_pm_out) bind(C, name='get_NX_pm_coarse') 
      use, intrinsic :: iso_c_binding, only: c_int
      implicit none
      integer(c_int), intent(out) :: NX_pm_out
      NX_pm_out = NXpm_coarse
   End subroutine get_NXpm_coarse

   subroutine get_NYpm_coarse(NY_pm_out) bind(C, name='get_NY_pm_coarse') 
      use, intrinsic :: iso_c_binding, only: c_int
      implicit none
      integer(c_int), intent(out) :: NY_pm_out
      NY_pm_out = NYpm_coarse
   End subroutine get_NYpm_coarse

   subroutine get_NZpm_coarse(NZ_pm_out) bind(C, name='get_NZ_pm_coarse') 
      use, intrinsic :: iso_c_binding, only: c_int
      implicit none
      integer(c_int), intent(out) :: NZ_pm_out
      NZ_pm_out = NZpm_coarse
   End subroutine get_NZpm_coarse


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! SETTERS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine set_RHS_pm(RHS_pm_in,size1, size2, size3, size4) bind(C, name='set_RHS_pm')
      use, intrinsic :: iso_c_binding, only: c_int, c_double
      implicit none
      integer(c_int), intent(in) :: size1, size2, size3, size4
      real(c_double), dimension(size1, size2, size3, size4), intent(in) :: RHS_pm_in
      ! Deallocate previous RHS_pm if it is allocated
      if (allocated(RHS_pm)) deallocate (RHS_pm)
      ! Allocate RHS_pm to match the shape of RHS_pm_in
      allocate (RHS_pm(size(RHS_pm_in,1), size(RHS_pm_in,2), size(RHS_pm_in,3), size(RHS_pm_in,4)))
      ! Copy data from RHS_pm_in to RHS_pm
      RHS_pm = RHS_pm_in
   end subroutine set_RHS_pm
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! PRINTERS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine print_pmgrid_info()
      use console_io, only: dp_3d_alloc_info, dp_4d_alloc_info, i_1d_array_info
      print *, "PM GRID INFO"
      print *, "============"
      print *, ""
      call dp_3d_alloc_info("velvrx_pm", velvrx_pm)
      call dp_3d_alloc_info("velvry_pm", velvry_pm)
      call dp_3d_alloc_info("velvrz_pm", velvrz_pm)
      call dp_3d_alloc_info("deformx_pm", deformx_pm)
      call dp_3d_alloc_info("deformy_pm", deformy_pm)
      call dp_3d_alloc_info("deformz_pm", deformz_pm)
      call dp_4d_alloc_info("RHS_pm", RHS_pm)
      print *, achar(9), "XMIN_pm: ", XMIN_pm
      print *, achar(9), "XMAX_pm: ", XMAX_pm
      print *, achar(9), "YMIN_pm: ", YMIN_pm
      print *, achar(9), "YMAX_pm: ", YMAX_pm
      print *, achar(9), "ZMIN_pm: ", ZMIN_pm
      print *, achar(9), "ZMAX_pm: ", ZMAX_pm
      print *, achar(9), "DXpm: ", DXpm
      print *, achar(9), "DYpm: ", DYpm
      print *, achar(9), "DZpm: ", DZpm
      print *, achar(9), "DVpm: ", DVpm
      print *, achar(9), "NXpm: ", NXpm_coarse
      print *, achar(9), "NYpm: ", NYpm_coarse
      print *, achar(9), "NZpm: ", NZpm_coarse
      print *, achar(9), "NXs_bl", NXs_coarse_bl
      print *, achar(9), "NYs_bl", NYs_coarse_bl
      print *, achar(9), "NXf_bl", NXf_coarse_bl
      print *, achar(9), "NYf_bl", NYf_coarse_bl
      print *, achar(9), "NZs_bl", NZs_coarse_bl
      print *, achar(9), "NZf_bl", NZf_coarse_bl
      print *, achar(9), "NBlocks: ", NBlocks
      print '(A,I6)',    achar(9)//"ND        =", ND
      print '(A,I6)',    achar(9)//"nbound    =", nbound
      print '(A,I6)',    achar(9)//"ndumcell  =", ndumcell
      print '(A,I6)',    achar(9)//"IDVPM     =", IDVPM
      call dp_4d_alloc_info("SOL_pm", SOL_pm)
   end subroutine print_pmgrid_info

   subroutine print_RHS_pm()
      use console_io, only: dp_4d_alloc_info 
      call dp_4d_alloc_info("RHS_pm", RHS_pm) 
   end subroutine print_RHS_pm

   subroutine print_velocity_stats
      use console_io, only: vpm_print, nocolor, red, VERBOCITY, dummy_string
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

end module pmgrid
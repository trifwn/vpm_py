module pmgrid
   use base_types, only: dp
   real(dp), allocatable, target, save :: velvrx_pm(:, :, :), velvry_pm(:, :, :), velvrz_pm(:, :, :)
   real(dp), allocatable, target, save :: RHS_pm(:, :, :, :)

   real(dp), save                      :: XMIN_pm, XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
   real(dp), save                      :: DXpm, DYpm, DZpm, DVpm, EPSVOL
   integer, save                       :: NXpm_coarse, NYpm_coarse, NZpm_coarse
   integer, save                       :: NXs_coarse_bl, NYs_coarse_bl, NXf_coarse_bl,&
                                          NYf_coarse_bl, NZs_coarse_bl, NZf_coarse_bl
   integer, save                       :: NBlocks

   ! GETTERS
   public :: get_NXpm_coarse, get_NYpm_coarse, get_NZpm_coarse
   ! SETTERS
   public :: set_RHS_pm 
   ! Printers
   public :: print_pmgrid_info

contains

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
      use io, only: dp_3d_alloc_info, dp_4d_alloc_info, i_1d_array_info
      print *, "PM GRID INFO"
      print *, "============"
      print *, ""
      call dp_3d_alloc_info("velvrx_pm", velvrx_pm)
      call dp_3d_alloc_info("velvry_pm", velvry_pm)
      call dp_3d_alloc_info("velvrz_pm", velvrz_pm)
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
      print *, achar(9), "EPSVOL: ", EPSVOL
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
   end subroutine print_pmgrid_info

   subroutine print_RHS_pm()
      use io, only: dp_4d_alloc_info 
      call dp_4d_alloc_info("RHS_pm", RHS_pm) 
   end subroutine print_RHS_pm

end module pmgrid
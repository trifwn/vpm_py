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
   public :: get_NXpm, get_NYpm, get_NZpm

   ! SETTERS
   public :: set_RHS_pm 

   ! Printers
   public :: print_pmgrid_info

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! GETTERS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine get_NXpm(NXpm_)
      integer, intent(out) :: NXpm_
      NXpm_ = NXpm_coarse
   end subroutine get_NXpm

   subroutine get_NYpm(NYpm_)
      integer, intent(out) :: NYpm_
      NYpm_ = NYpm_coarse
   end subroutine get_NYpm

   subroutine get_NZpm(NZpm_)
      integer, intent(out) :: NZpm_
      NZpm_ = NZpm_coarse
   end subroutine get_NZpm

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! SETTERS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine set_RHS_pm(RHS_pm_in)
      real(dp), intent(in) :: RHS_pm_in(:,:,:,:)
      
      ! Deallocate previous RHS_pm if it is allocated
      if (allocated(RHS_pm)) deallocate(RHS_pm)

      ! Allocate RHS_pm to match the shape of RHS_pm_in
      allocate(RHS_pm(size(RHS_pm_in,1), size(RHS_pm_in,2), size(RHS_pm_in,3), size(RHS_pm_in,4)))

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
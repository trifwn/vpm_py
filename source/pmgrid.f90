module pmgrid
   use base_types, only: dp
   real(dp), allocatable, target, save :: velvrx_pm(:, :, :), velvry_pm(:, :, :), velvrz_pm(:, :, :)
   real(dp), allocatable, target, save :: RHS_pm(:, :, :, :)

   real(dp), save                      :: XMIN_pm, XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
   real(dp), save                      :: DD, DXpm, DYpm, DZpm, DXpm2, DYpm2, DZpm2, DVpm, EPSVOL
   integer, save                       :: NXpm, NYpm, NZpm, NXpm_par, NYpm_par, NZpm_par
   integer, save                       :: NXs_bl(10), NYs_bl(10), NXf_bl(10), NYf_bl(10), NZs_bl(10), NZf_bl(10), NBlocks

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

   Subroutine get_NXpm(NXpm_)
      integer, intent(out) :: NXpm_
      NXpm_ = NXpm
   end Subroutine get_NXpm

   Subroutine get_NYpm(NYpm_)
      integer, intent(out) :: NYpm_
      NYpm_ = NYpm
   end Subroutine get_NYpm

   Subroutine get_NZpm(NZpm_)
      integer, intent(out) :: NZpm_
      NZpm_ = NZpm
   end Subroutine get_NZpm

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

   Subroutine print_pmgrid_info()
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
      print *, achar(9), "DD: ", DD
      print *, achar(9), "DXpm: ", DXpm
      print *, achar(9), "DYpm: ", DYpm
      print *, achar(9), "DZpm: ", DZpm
      print *, achar(9), "DXpm2: ", DXpm2
      print *, achar(9), "DYpm2: ", DYpm2
      print *, achar(9), "DZpm2: ", DZpm2
      print *, achar(9), "DVpm: ", DVpm
      print *, achar(9), "EPSVOL: ", EPSVOL
      print *, achar(9), "NXpm: ", NXpm
      print *, achar(9), "NYpm: ", NYpm
      print *, achar(9), "NZpm: ", NZpm
      print *, achar(9), "NXpm_par: ", NXpm_par
      print *, achar(9), "NYpm_par: ", NYpm_par
      print *, achar(9), "NZpm_par: ", NZpm_par
      call i_1d_array_info("NXs_bl", NXs_bl, 10)
      call i_1d_array_info("NYs_bl", NYs_bl, 10)
      call i_1d_array_info("NXf_bl", NXf_bl, 10)
      call i_1d_array_info("NYf_bl", NYf_bl, 10)
      call i_1d_array_info("NZs_bl", NZs_bl, 10)
      call i_1d_array_info("NZf_bl", NZf_bl, 10)
      print *, achar(9), "NBlocks: ", NBlocks
   end Subroutine print_pmgrid_info

   Subroutine print_RHS_pm()
      use io, only: dp_4d_alloc_info 
      call dp_4d_alloc_info("RHS_pm", RHS_pm) 
   end Subroutine print_RHS_pm

end module pmgrid
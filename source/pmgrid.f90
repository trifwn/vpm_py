module pmgrid
   double precision, allocatable, target, save        :: velvrx_pm(:, :, :), velvry_pm(:, :, :), velvrz_pm(:, :, :)
   double precision, allocatable, target, save       :: RHS_pm(:, :, :, :)

   double precision, save               :: XMIN_pm, XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
   double precision, save               :: DD, DXpm, DYpm, DZpm, DXpm2, DYpm2, DZpm2, DVpm, EPSVOL
   integer, save                        :: NXpm, NYpm, NZpm, NXpm_par, NYpm_par, NZpm_par
   integer, save                        :: NXs_bl(10), NYs_bl(10), NXf_bl(10), NYf_bl(10), NZs_bl(10), NZf_bl(10), NBlocks

   ! GETTERS
   public :: get_velvrx_pm, get_velvry_pm, get_velvrz_pm, get_RHS_pm
   public :: get_XMIN_pm, get_XMAX_pm, get_YMIN_pm, get_YMAX_pm, get_ZMIN_pm, get_ZMAX_pm
   public :: get_DD, get_DXpm, get_DYpm, get_DZpm, get_DXpm2, get_DYpm2, get_DZpm2, get_DVpm, get_EPSVOL
   public :: get_NXpm, get_NYpm, get_NZpm, get_NXpm_par, get_NYpm_par, get_NZpm_par
   public :: get_NXs_bl, get_NYs_bl, get_NXf_bl, get_NYf_bl, get_NZs_bl, get_NZf_bl, get_NBlocks

   ! SETTERS
   public :: set_velvrx_pm, set_velvry_pm, set_velvrz_pm, set_RHS_pm
   public :: set_XMIN_pm, set_XMAX_pm, set_YMIN_pm, set_YMAX_pm, set_ZMIN_pm, set_ZMAX_pm
   public :: set_DD, set_DXpm, set_DYpm, set_DZpm, set_DXpm2, set_DYpm2, set_DZpm2, set_DVpm, set_EPSVOL
   public :: set_NXpm, set_NYpm, set_NZpm, set_NXpm_par, set_NYpm_par, set_NZpm_par
   public :: set_NXs_bl, set_NYs_bl, set_NXf_bl, set_NYf_bl, set_NZs_bl, set_NZf_bl, set_NBlocks

   ! Printers
   public :: print_pmgrid_info

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! GETTERS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   Subroutine get_velvrx_pm(velvrx_pm_)
      double precision, allocatable, intent(out),dimension(:,:,:) :: velvrx_pm_
      velvrx_pm_ = velvrx_pm
   end Subroutine get_velvrx_pm

   Subroutine get_velvry_pm(velvry_pm_)
      double precision, allocatable, intent(out),dimension(:,:,:) :: velvry_pm_
      velvry_pm_ = velvry_pm 
   end Subroutine get_velvry_pm

   Subroutine get_velvrz_pm(velvrz_pm_)
      double precision, allocatable, intent(out),dimension(:,:,:) :: velvrz_pm_
      velvrz_pm_ = velvrz_pm
   end Subroutine get_velvrz_pm

   Subroutine get_RHS_pm(RHS_pm_)
      double precision, allocatable, intent(out), dimension(:,:,:,:) :: RHS_pm_
      RHS_pm_ = RHS_pm
   end Subroutine get_RHS_pm

   Subroutine get_XMIN_pm(XMIN_pm_)
      double precision, intent(out) :: XMIN_pm_
      XMIN_pm_ = XMIN_pm
   end Subroutine get_XMIN_pm

   Subroutine get_XMAX_pm(XMAX_pm_)
      double precision, intent(out) :: XMAX_pm_
      XMAX_pm_ = XMAX_pm
   end Subroutine get_XMAX_pm

   Subroutine get_YMIN_pm(YMIN_pm_)
      double precision, intent(out) :: YMIN_pm_
      YMIN_pm_ = YMIN_pm
   end Subroutine get_YMIN_pm

   Subroutine get_YMAX_pm(YMAX_pm_)
      double precision, intent(out) :: YMAX_pm_
      YMAX_pm_ = YMAX_pm
   end Subroutine get_YMAX_pm

   Subroutine get_ZMIN_pm(ZMIN_pm_)
      double precision, intent(out) :: ZMIN_pm_
      ZMIN_pm_ = ZMIN_pm
   end Subroutine get_ZMIN_pm

   Subroutine get_ZMAX_pm(ZMAX_pm_)
      double precision, intent(out) :: ZMAX_pm_
      ZMAX_pm_ = ZMAX_pm
   end Subroutine get_ZMAX_pm

   Subroutine get_DD(DD_)
      double precision, intent(out) :: DD_
      DD_ = DD
   end Subroutine get_DD

   Subroutine get_DXpm(DXpm_)
      double precision, intent(out) :: DXpm_
      DXpm_ = DXpm
   end Subroutine get_DXpm

   Subroutine get_DYpm(DYpm_)
      double precision, intent(out) :: DYpm_
      DYpm_ = DYpm
   end Subroutine get_DYpm

   Subroutine get_DZpm(DZpm_)
      double precision, intent(out) :: DZpm_
      DZpm_ = DZpm
   end Subroutine get_DZpm

   Subroutine get_DXpm2(DXpm2_)
      double precision, intent(out) :: DXpm2_
      DXpm2_ = DXpm2
   end Subroutine get_DXpm2

   Subroutine get_DYpm2(DYpm2_)
      double precision, intent(out) :: DYpm2_
      DYpm2_ = DYpm2
   end Subroutine get_DYpm2

   Subroutine get_DZpm2(DZpm2_)
      double precision, intent(out) :: DZpm2_
      DZpm2_ = DZpm2
   end Subroutine get_DZpm2

   Subroutine get_DVpm(DVpm_)
      double precision, intent(out) :: DVpm_
      DVpm_ = DVpm
   end Subroutine get_DVpm

   Subroutine get_EPSVOL(EPSVOL_)
      double precision, intent(out) :: EPSVOL_
      EPSVOL_ = EPSVOL
   end Subroutine get_EPSVOL

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

   Subroutine get_NXpm_par(NXpm_par_)
      integer, intent(out) :: NXpm_par_
      NXpm_par_ = NXpm_par
   end Subroutine get_NXpm_par

   Subroutine get_NYpm_par(NYpm_par_)
      integer, intent(out) :: NYpm_par_
      NYpm_par_ = NYpm_par
   end Subroutine get_NYpm_par

   Subroutine get_NZpm_par(NZpm_par_)
      integer, intent(out) :: NZpm_par_
      NZpm_par_ = NZpm_par
   end Subroutine get_NZpm_par

   Subroutine get_NXs_bl(NXs_bl_)
      integer, dimension(10), intent(out) :: NXs_bl_
      NXs_bl_ = NXs_bl
   end Subroutine get_NXs_bl

   Subroutine get_NYs_bl(NYs_bl_)
      integer, dimension(10), intent(out) :: NYs_bl_
      NYs_bl_ = NYs_bl
   end Subroutine get_NYs_bl

   Subroutine get_NXf_bl(NXf_bl_)
      integer, dimension(10), intent(out) :: NXf_bl_
      NXf_bl_ = NXf_bl
   end Subroutine get_NXf_bl

   Subroutine get_NYf_bl(NYf_bl_)
      integer, dimension(10), intent(out) :: NYf_bl_
      NYf_bl_ = NYf_bl
   end Subroutine get_NYf_bl

   Subroutine get_NZs_bl(NZs_bl_)
      integer, dimension(10), intent(out) :: NZs_bl_
      NZs_bl_ = NZs_bl
   end Subroutine get_NZs_bl

   Subroutine get_NZf_bl(NZf_bl_)
      integer, dimension(10), intent(out) :: NZf_bl_
      NZf_bl_ = NZf_bl
   end Subroutine get_NZf_bl

   Subroutine get_NBlocks(NBlocks_)
      integer, intent(out) :: NBlocks_
      NBlocks_ = NBlocks
   end Subroutine get_NBlocks

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! SETTERS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   Subroutine set_velvrx_pm(velvrx_pm_)
      double precision, allocatable, intent(in) :: velvrx_pm_
      velvrx_pm = velvrx_pm_
   end Subroutine set_velvrx_pm

   Subroutine set_velvry_pm(velvry_pm_)
      double precision, allocatable, intent(in) :: velvry_pm_
      velvry_pm = velvry_pm_
   end Subroutine set_velvry_pm

   Subroutine set_velvrz_pm(velvrz_pm_)
      double precision, allocatable, intent(in) :: velvrz_pm_
      velvrz_pm = velvrz_pm_
   end Subroutine set_velvrz_pm

   subroutine set_RHS_pm(RHS_pm_in)
      double precision, intent(in) :: RHS_pm_in(:,:,:,:)
      
      ! Deallocate previous RHS_pm if it is allocated
      if (allocated(RHS_pm)) deallocate(RHS_pm)

      ! Allocate RHS_pm to match the shape of RHS_pm_in
      allocate(RHS_pm(size(RHS_pm_in,1), size(RHS_pm_in,2), size(RHS_pm_in,3), size(RHS_pm_in,4)))

      ! Copy data from RHS_pm_in to RHS_pm
      RHS_pm = RHS_pm_in
   end subroutine set_RHS_pm

   Subroutine set_XMIN_pm(XMIN_pm_)
      double precision, intent(in) :: XMIN_pm_
      XMIN_pm = XMIN_pm_
   end Subroutine set_XMIN_pm

   Subroutine set_XMAX_pm(XMAX_pm_)
      double precision, intent(in) :: XMAX_pm_
      XMAX_pm = XMAX_pm_
   end Subroutine set_XMAX_pm

   Subroutine set_YMIN_pm(YMIN_pm_)
      double precision, intent(in) :: YMIN_pm_
      YMIN_pm = YMIN_pm_
   end Subroutine set_YMIN_pm

   Subroutine set_YMAX_pm(YMAX_pm_)
      double precision, intent(in) :: YMAX_pm_
      YMAX_pm = YMAX_pm_
   end Subroutine set_YMAX_pm

   Subroutine set_ZMIN_pm(ZMIN_pm_)
      double precision, intent(in) :: ZMIN_pm_
      ZMIN_pm = ZMIN_pm_
   end Subroutine set_ZMIN_pm

   Subroutine set_ZMAX_pm(ZMAX_pm_)
      double precision, intent(in) :: ZMAX_pm_
      ZMAX_pm = ZMAX_pm_
   end Subroutine set_ZMAX_pm

   Subroutine set_DD(DD_)
      double precision, intent(in) :: DD_
      DD = DD_
   end Subroutine set_DD

   Subroutine set_DXpm(DXpm_)
      double precision, intent(in) :: DXpm_
      DXpm = DXpm_
   end Subroutine set_DXpm

   Subroutine set_DYpm(DYpm_)
      double precision, intent(in) :: DYpm_
      DYpm = DYpm_
   end Subroutine set_DYpm

   Subroutine set_DZpm(DZpm_)
      double precision, intent(in) :: DZpm_
      DZpm = DZpm_
   end Subroutine set_DZpm

   Subroutine set_DXpm2(DXpm2_)
      double precision, intent(in) :: DXpm2_
      DXpm2 = DXpm2_
   end Subroutine set_DXpm2

   Subroutine set_DYpm2(DYpm2_)
      double precision, intent(in) :: DYpm2_
      DYpm2 = DYpm2_
   end Subroutine set_DYpm2

   Subroutine set_DZpm2(DZpm2_)
      double precision, intent(in) :: DZpm2_
      DZpm2 = DZpm2_
   end Subroutine set_DZpm2

   Subroutine set_DVpm(DVpm_)
      double precision, intent(in) :: DVpm_
      DVpm = DVpm_
   end Subroutine set_DVpm

   Subroutine set_EPSVOL(EPSVOL_)
      double precision, intent(in) :: EPSVOL_
      EPSVOL = EPSVOL_
   end Subroutine set_EPSVOL

   Subroutine set_NXpm(NXpm_)
      integer, intent(in) :: NXpm_
      NXpm = NXpm_
   end Subroutine set_NXpm

   Subroutine set_NYpm(NYpm_)
      integer, intent(in) :: NYpm_
      NYpm = NYpm_
   end Subroutine set_NYpm

   Subroutine set_NZpm(NZpm_)
      integer, intent(in) :: NZpm_
      NZpm = NZpm_
   end Subroutine set_NZpm

   Subroutine set_NXpm_par(NXpm_par_)
      integer, intent(in) :: NXpm_par_
      NXpm_par = NXpm_par_
   end Subroutine set_NXpm_par

   Subroutine set_NYpm_par(NYpm_par_)
      integer, intent(in) :: NYpm_par_
      NYpm_par = NYpm_par_
   end Subroutine set_NYpm_par

   Subroutine set_NZpm_par(NZpm_par_)
      integer, intent(in) :: NZpm_par_
      NZpm_par = NZpm_par_
   end Subroutine set_NZpm_par

   Subroutine set_NXs_bl(NXs_bl_)
      integer, dimension(10), intent(in) :: NXs_bl_
      NXs_bl = NXs_bl_
   end Subroutine set_NXs_bl

   Subroutine set_NYs_bl(NYs_bl_)
      integer, dimension(10), intent(in) :: NYs_bl_
      NYs_bl = NYs_bl_
   end Subroutine set_NYs_bl

   Subroutine set_NXf_bl(NXf_bl_)
      integer, dimension(10), intent(in) :: NXf_bl_
      NXf_bl = NXf_bl_
   end Subroutine set_NXf_bl

   Subroutine set_NYf_bl(NYf_bl_)
      integer, dimension(10), intent(in) :: NYf_bl_
      NYf_bl = NYf_bl_
   end Subroutine set_NYf_bl

   Subroutine set_NZs_bl(NZs_bl_)
      integer, dimension(10), intent(in) :: NZs_bl_
      NZs_bl = NZs_bl_
   end Subroutine set_NZs_bl

   Subroutine set_NZf_bl(NZf_bl_)
      integer, dimension(10), intent(in) :: NZf_bl_
      NZf_bl = NZf_bl_
   end Subroutine set_NZf_bl

   Subroutine set_NBlocks(NBlocks_)
      integer, intent(in) :: NBlocks_
      NBlocks = NBlocks_
   end Subroutine set_NBlocks

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


end module pmgrid
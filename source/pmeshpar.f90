module pmeshpar
   double precision, save              :: PI, PI2, PI4, DT

   double precision, allocatable, save   :: velx_pm(:, :, :), vely_pm(:, :, :), velz_pm(:, :, :), qx_pm(:, :, :), qy_pm(:, :, :)
   double precision, allocatable, save   :: velphix_pm(:, :, :), velphiy_pm(:, :, :), velphiz_pm(:, :, :)
   double precision, allocatable, save   :: SOL_pm(:, :, :, :), SOL_0_pm(:, :, :, :)
   double precision, allocatable, save    :: source_bound(:, :), x_s(:, :), y_s(:, :), z_s(:, :), d_s(:), cos_s(:), sin_s(:)
   double precision,allocatable,save    :: source_bound_lev(:,:,:),xs_lev(:,:,:),ys_lev(:,:,:),zs_lev(:,:,:),ds_lev(:,:,:)
   integer, allocatable, save    :: nbound_lev(:)
   ! double precision, allocatable, save    :: Psiz_pm_0(:, :, :), Psiz_pm_f(:, :, :)
   ! double precision, allocatable, save   :: Cont_pm(:, :)
   integer                              :: levmax, npar_cell, ND

   integer, save                        :: nbound, ndumcell, NVR_CFD_sa, IDVPM

   ! Getters
   public :: get_PI, get_PI2, get_PI4, get_DT
   public :: get_velx_pm, get_vely_pm, get_velz_pm, get_qx_pm, get_qy_pm
   public :: get_velphix_pm, get_velphiy_pm, get_velphiz_pm
   public :: get_SOL_pm, get_SOL_0_pm
   public :: get_source_bound, get_x_s, get_y_s, get_z_s, get_d_s, get_cos_s, get_sin_s
   public :: get_source_bound_lev, get_xs_lev, get_ys_lev, get_zs_lev, get_ds_lev
   public :: get_nbound_lev
   public :: get_levmax, get_npar_cell, get_ND
   public :: get_nbound, get_ndumcell, get_NVR_CFD_sa, get_IDVPM

   ! Setters
   public :: set_PI, set_PI2, set_PI4, set_DT
   public :: set_velx_pm, set_vely_pm, set_velz_pm, set_qx_pm, set_qy_pm
   public :: set_velphix_pm, set_velphiy_pm, set_velphiz_pm
   public :: set_SOL_pm, set_SOL_0_pm
   public :: set_source_bound, set_x_s, set_y_s, set_z_s, set_d_s, set_cos_s, set_sin_s
   public :: set_source_bound_lev, set_xs_lev, set_ys_lev, set_zs_lev, set_ds_lev
   public :: set_nbound_lev
   public :: set_levmax, set_npar_cell, set_ND
   public :: set_nbound, set_ndumcell, set_NVR_CFD_sa, set_IDVPM

   ! Printers
   public :: print_pmeshpar_info

contains 


    !!!!!!!!!!!!!!!!!!!
    ! Getters
    !!!!!!!!!!!!!!!!!!!
   function get_PI() result(val)
      double precision :: val
      val = PI
   end function get_PI

    function get_PI2() result(val)
        double precision :: val
        val = PI2
    end function get_PI2

    function get_PI4() result(val)
        double precision :: val
        val = PI4
    end function get_PI4

    function get_DT() result(val)
        double precision :: val
        val = DT
    end function get_DT

    function get_velx_pm() result(arr)
      double precision, dimension(:,:,:), allocatable :: arr
      arr = velx_pm
   end function get_velx_pm

    function get_vely_pm() result(arr)
        double precision, dimension(:,:,:), allocatable :: arr
        arr = vely_pm
    end function get_vely_pm

    function get_velz_pm() result(arr)
        double precision, dimension(:,:,:), allocatable :: arr
        arr = velz_pm
    end function get_velz_pm

    function get_qx_pm() result(arr)
        double precision, dimension(:,:,:), allocatable :: arr
        arr = qx_pm
    end function get_qx_pm

    function get_qy_pm() result(arr)
        double precision, dimension(:,:,:), allocatable :: arr
        arr = qy_pm
    end function get_qy_pm

    function get_velphix_pm() result(arr)
        double precision, dimension(:,:,:), allocatable :: arr
        arr = velphix_pm
    end function get_velphix_pm

    function get_velphiy_pm() result(arr)
        double precision, dimension(:,:,:), allocatable :: arr
        arr = velphiy_pm
    end function get_velphiy_pm

    function get_velphiz_pm() result(arr)
        double precision, dimension(:,:,:), allocatable :: arr
        arr = velphiz_pm
    end function get_velphiz_pm

    function get_SOL_pm() result(arr)
        double precision, dimension(:,:,:,:), allocatable :: arr
        arr = SOL_pm
    end function get_SOL_pm

    function get_SOL_0_pm() result(arr)
        double precision, dimension(:,:,:,:), allocatable :: arr
        arr = SOL_0_pm
    end function get_SOL_0_pm

    function get_source_bound() result(arr)
        double precision, dimension(:,:), allocatable :: arr
        arr = source_bound
    end function get_source_bound

    function get_x_s() result(arr)
        double precision, dimension(:,:), allocatable :: arr
        arr = x_s
    end function get_x_s

    function get_y_s() result(arr)
        double precision, dimension(:,:), allocatable :: arr
        arr = y_s
    end function get_y_s

    function get_z_s() result(arr)
        double precision, dimension(:,:), allocatable :: arr
        arr = z_s
    end function get_z_s

    function get_d_s() result(arr)
        double precision, dimension(:), allocatable :: arr
        arr = d_s
    end function get_d_s

    function get_cos_s() result(arr)
        double precision, dimension(:), allocatable :: arr
        arr = cos_s
    end function get_cos_s

    function get_sin_s() result(arr)
        double precision, dimension(:), allocatable :: arr
        arr = sin_s
    end function get_sin_s

    function get_source_bound_lev() result(arr)
        double precision, dimension(:,:,:), allocatable :: arr
        arr = source_bound_lev
    end function get_source_bound_lev

    function get_xs_lev() result(arr)
        double precision, dimension(:,:,:), allocatable :: arr
        arr = xs_lev
    end function get_xs_lev

    function get_ys_lev() result(arr)
        double precision, dimension(:,:,:), allocatable :: arr
        arr = ys_lev
    end function get_ys_lev

    function get_zs_lev() result(arr)
        double precision, dimension(:,:,:), allocatable :: arr
        arr = zs_lev
    end function get_zs_lev

    function get_ds_lev() result(arr)
        double precision, dimension(:,:,:), allocatable :: arr
        arr = ds_lev
    end function get_ds_lev

    function get_nbound_lev() result(arr)
        integer, dimension(:), allocatable :: arr
        arr = nbound_lev
    end function get_nbound_lev

    function get_levmax() result(val)
        integer :: val
        val = levmax
    end function get_levmax

    function get_npar_cell() result(val)
        integer :: val
        val = npar_cell
    end function get_npar_cell

    function get_ND() result(val)
        integer :: val
        val = ND
    end function get_ND

    function get_nbound() result(val)
        integer :: val
        val = nbound
    end function get_nbound

    function get_ndumcell() result(val)
        integer :: val
        val = ndumcell
    end function get_ndumcell

    function get_NVR_CFD_sa() result(val)
        integer :: val
        val = NVR_CFD_sa
    end function get_NVR_CFD_sa

    function get_IDVPM() result(val)
        integer :: val
        val = IDVPM
    end function get_IDVPM

    !!!!!!!!!!!!!!!!!!!
    ! Setters
    !!!!!!!!!!!!!!!!!!!

    subroutine set_PI(val)
        double precision, intent(in) :: val
        PI = val
    end subroutine set_PI

    subroutine set_PI2(val)
        double precision, intent(in) :: val
        PI2 = val
    end subroutine set_PI2

    subroutine set_PI4(val)
        double precision, intent(in) :: val
        PI4 = val
    end subroutine set_PI4

    subroutine set_DT(val)
        double precision, intent(in) :: val
        DT = val
    end subroutine set_DT

    subroutine set_velx_pm(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(velx_pm)) deallocate(velx_pm)
        allocate(velx_pm, source=arr)
    end subroutine set_velx_pm

    subroutine set_vely_pm(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(vely_pm)) deallocate(vely_pm)
        allocate(vely_pm, source=arr)
    end subroutine set_vely_pm

    subroutine set_velz_pm(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(velz_pm)) deallocate(velz_pm)
        allocate(velz_pm, source=arr)
    end subroutine set_velz_pm

    subroutine set_qx_pm(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(qx_pm)) deallocate(qx_pm)
        allocate(qx_pm, source=arr)
    end subroutine set_qx_pm

    subroutine set_qy_pm(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(qy_pm)) deallocate(qy_pm)
        allocate(qy_pm, source=arr)
    end subroutine set_qy_pm

    subroutine set_velphix_pm(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(velphix_pm)) deallocate(velphix_pm)
        allocate(velphix_pm, source=arr)
    end subroutine set_velphix_pm

    subroutine set_velphiy_pm(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(velphiy_pm)) deallocate(velphiy_pm)
        allocate(velphiy_pm, source=arr)
    end subroutine set_velphiy_pm

    subroutine set_velphiz_pm(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(velphiz_pm)) deallocate(velphiz_pm)
        allocate(velphiz_pm, source=arr)
    end subroutine set_velphiz_pm

    subroutine set_SOL_pm(arr)
        double precision, dimension(:,:,:,:), intent(in) :: arr
        if (allocated(SOL_pm)) deallocate(SOL_pm)
        allocate(SOL_pm, source=arr)
    end subroutine set_SOL_pm

    subroutine set_SOL_0_pm(arr)
        double precision, dimension(:,:,:,:), intent(in) :: arr
        if (allocated(SOL_0_pm)) deallocate(SOL_0_pm)
        allocate(SOL_0_pm, source=arr)
    end subroutine set_SOL_0_pm

    subroutine set_source_bound(arr)
        double precision, dimension(:,:), intent(in) :: arr
        if (allocated(source_bound)) deallocate(source_bound)
        allocate(source_bound, source=arr)
    end subroutine set_source_bound

    subroutine set_x_s(arr)
        double precision, dimension(:,:), intent(in) :: arr
        if (allocated(x_s)) deallocate(x_s)
        allocate(x_s, source=arr)
    end subroutine set_x_s

    subroutine set_y_s(arr)
        double precision, dimension(:,:), intent(in) :: arr
        if (allocated(y_s)) deallocate(y_s)
        allocate(y_s, source=arr)
    end subroutine set_y_s

    subroutine set_z_s(arr)
        double precision, dimension(:,:), intent(in) :: arr
        if (allocated(z_s)) deallocate(z_s)
        allocate(z_s, source=arr)
    end subroutine set_z_s

    subroutine set_d_s(arr)
        double precision, dimension(:), intent(in) :: arr
        if (allocated(d_s)) deallocate(d_s)
        allocate(d_s, source=arr)
    end subroutine set_d_s

    subroutine set_cos_s(arr)
        double precision, dimension(:), intent(in) :: arr
        if (allocated(cos_s)) deallocate(cos_s)
        allocate(cos_s, source=arr)
    end subroutine set_cos_s

    subroutine set_sin_s(arr)
        double precision, dimension(:), intent(in) :: arr
        if (allocated(sin_s)) deallocate(sin_s)
        allocate(sin_s, source=arr)
    end subroutine set_sin_s

    subroutine set_source_bound_lev(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(source_bound_lev)) deallocate(source_bound_lev)
        allocate(source_bound_lev, source=arr)
    end subroutine set_source_bound_lev

    subroutine set_xs_lev(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(xs_lev)) deallocate(xs_lev)
        allocate(xs_lev, source=arr)
    end subroutine set_xs_lev

    subroutine set_ys_lev(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(ys_lev)) deallocate(ys_lev)
        allocate(ys_lev, source=arr)
    end subroutine set_ys_lev

    subroutine set_zs_lev(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(zs_lev)) deallocate(zs_lev)
        allocate(zs_lev, source=arr)
    end subroutine set_zs_lev

    subroutine set_ds_lev(arr)
        double precision, dimension(:,:,:), intent(in) :: arr
        if (allocated(ds_lev)) deallocate(ds_lev)
        allocate(ds_lev, source=arr)
    end subroutine set_ds_lev

    subroutine set_nbound_lev(arr)
        integer, dimension(:), intent(in) :: arr
        if (allocated(nbound_lev)) deallocate(nbound_lev)
        allocate(nbound_lev, source=arr)
    end subroutine set_nbound_lev

    subroutine set_levmax(val)
        integer, intent(in) :: val
        levmax = val
    end subroutine set_levmax

    subroutine set_npar_cell(val)
        integer, intent(in) :: val
        npar_cell = val
    end subroutine set_npar_cell

    subroutine set_ND(val)
        integer, intent(in) :: val
        ND = val
    end subroutine set_ND

    subroutine set_nbound(val)
        integer, intent(in) :: val
        nbound = val
    end subroutine set_nbound

    subroutine set_ndumcell(val)
        integer, intent(in) :: val
        ndumcell = val
    end subroutine set_ndumcell

    subroutine set_NVR_CFD_sa(val)
        integer, intent(in) :: val
        NVR_CFD_sa = val
    end subroutine set_NVR_CFD_sa

    subroutine set_IDVPM(val)
        integer, intent(in) :: val
        IDVPM = val
    end subroutine set_IDVPM

    !!!!!!!!!!!!!!!!!!!
    ! Printers
    !!!!!!!!!!!!!!!!!!!

    subroutine print_pmeshpar_info()
        use io, only: dp_3d_alloc_info, dp_4d_alloc_info, dp_2d_alloc_info, i_1d_alloc_info, &
                      dp_1d_alloc_info
        print *, "Pmeshpar INFO"
        print *, "============"
        print *, ""
        print '(A,F12.6)', achar(9)//"PI  =", PI
        print '(A,F12.6)', achar(9)//"PI2 =", PI2
        print '(A,F12.6)', achar(9)//"PI4 =", PI4
        print '(A,F12.6)', achar(9)//"DT  =", DT
        print '(A,I6)',    achar(9)//"levmax    =", levmax
        print '(A,I6)',    achar(9)//"npar_cell =", npar_cell
        print '(A,I6)',    achar(9)//"ND        =", ND
        print '(A,I6)',    achar(9)//"nbound    =", nbound
        print '(A,I6)',    achar(9)//"ndumcell  =", ndumcell
        print '(A,I6)',    achar(9)//"NVR_CFD_sa=", NVR_CFD_sa
        print '(A,I6)',    achar(9)//"IDVPM     =", IDVPM
        
        print '(A)', achar(9)//""
        print '(A)', achar(9)//"Array Information:"
        
        ! Velocity arrays
        call dp_3d_alloc_info("velx_pm", velx_pm)
        call dp_3d_alloc_info("vely_pm", vely_pm)
        call dp_3d_alloc_info("velz_pm", velz_pm)
        call dp_3d_alloc_info("qx_pm", qx_pm)
        call dp_3d_alloc_info("qy_pm", qy_pm)
        
        ! Velocity phi arrays
        call dp_3d_alloc_info("velphix_pm", velphix_pm)
        call dp_3d_alloc_info("velphiy_pm", velphiy_pm)
        call dp_3d_alloc_info("velphiz_pm", velphiz_pm)
        
        ! SOL arrays
        call dp_4d_alloc_info("SOL_pm", SOL_pm)
        call dp_4d_alloc_info("SOL_0_pm", SOL_0_pm)
        
        ! Source and coordinate arrays
        call dp_2d_alloc_info("source_bound", source_bound)
        call dp_2d_alloc_info("x_s", x_s)
        call dp_2d_alloc_info("y_s", y_s)
        call dp_2d_alloc_info("z_s", z_s)
        call dp_1d_alloc_info("d_s", d_s)
        call dp_1d_alloc_info("cos_s", cos_s)
        call dp_1d_alloc_info("sin_s", sin_s)
        
        ! Level-specific arrays
        call dp_3d_alloc_info("source_bound_lev", source_bound_lev)
        call dp_3d_alloc_info("xs_lev", xs_lev)
        call dp_3d_alloc_info("ys_lev", ys_lev)
        call dp_3d_alloc_info("zs_lev", zs_lev)
        call dp_3d_alloc_info("ds_lev", ds_lev)
        call i_1d_alloc_info("nbound_lev", nbound_lev)

    end subroutine print_pmeshpar_info

end module pmeshpar
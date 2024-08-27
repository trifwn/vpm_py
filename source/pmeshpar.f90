module pmeshpar
    use base_types, only: dp
    use constants, only: PI, PI2, PI4
    real(dp), save        :: DT

    real(dp), allocatable, save      :: velx_pm(:, :, :), vely_pm(:, :, :), velz_pm(:, :, :)
    real(dp), allocatable, save      :: velphix_pm(:, :, :), velphiy_pm(:, :, :), velphiz_pm(:, :, :)
    real(dp), allocatable, save      :: qx_pm(:, :, :), qy_pm(:, :, :)
    real(dp), allocatable, save      :: SOL_pm(:, :, :, :), SOL_0_pm(:, :, :, :)
    real(dp), allocatable, save      :: deformation_pm(:,:,:,:)
    real(dp), allocatable, save      :: source_bound(:, :), source_bound_lev(:,:,:)
    real(dp), allocatable, save      :: x_s(:, :), y_s(:, :), z_s(:, :), d_s(:), cos_s(:), sin_s(:)
    real(dp),allocatable,save        :: xs_lev(:,:,:),ys_lev(:,:,:),zs_lev(:,:,:),ds_lev(:,:,:)
    integer, allocatable, save       :: nbound_lev(:)
    integer                          :: levmax, npar_cell, ND
    integer, save                    :: nbound, ndumcell, NVR_CFD_sa, IDVPM
    ! real(dp), allocatable, save    :: Psiz_pm_0(:, :, :), Psiz_pm_f(:, :, :)
    ! real(dp), allocatable, save    :: Cont_pm(:, :)


    ! Printers
    public :: print_pmeshpar_info

contains 

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
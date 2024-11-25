module pmgrid
    use base_types, only: dp

    type  :: cartesian_grid
        real(dp) :: Xbound(6) ! Domain boundaries
        real(dp) :: Dpm(3)    ! Grid spacing
        integer  :: NN(3)     ! Number of nodes in each direction
        integer  :: NN_bl(6)  ! Indices of the domain that do not include the dummy cells (start finish)
    end type cartesian_grid

    real(dp), allocatable, target, save :: velocity_pm(:, :, :, :)
    real(dp), allocatable, target, save :: deform_pm(:, :, :, :)
    real(dp), allocatable, target, save :: RHS_pm(:, :, :, :)
    real(dp), allocatable, target, save :: SOL_pm(:, :, :, :)
    real(dp), allocatable, target       :: SOL_pm_bl(:, :, :, :)
    real(dp), allocatable, target       :: RHS_pm_bl(:, :, :, :)

    ! FINE GRID CHARACTERISTICS
   real(dp), save                      :: XMIN_pm,     YMIN_pm,     ZMIN_pm,           &
                                          XMAX_pm,     YMAX_pm,     ZMAX_pm
   real(dp), save                      :: DXpm,        DYpm,        DZpm,        DVpm
   integer, save                       :: NXpm_fine,   NYpm_fine,   NZpm_fine
   integer, save                       :: NXs_fine_bl, NYs_fine_bl, NXf_fine_bl,       &
                                          NYf_fine_bl, NZs_fine_bl, NZf_fine_bl

   !! PMESH PARAMETERS
    integer, save                       :: nbound, ndumcell, ncoarse, IDVPM

    ! GETTERS
    public :: get_NXpm, get_NYpm, get_NZpm
    ! SETTERS
    public :: set_RHS_pm
    ! Printers
    public :: print_pmgrid_info

contains

    subroutine associate_velocities(vel_ptr)
        real(dp), pointer, intent(out) :: vel_ptr(:, :, :, :)
        if (associated(vel_ptr)) nullify (vel_ptr)
        vel_ptr => velocity_pm
    end subroutine associate_velocities

    subroutine associate_deformations(deform_ptr)
        real(dp), pointer, intent(out) :: deform_ptr(:, :, :, :)
        if (associated(deform_ptr)) nullify (deform_ptr)
        deform_ptr => deform_pm
    end subroutine associate_deformations

    subroutine set_pm_velocities_zero
        if (allocated(velocity_pm)) deallocate (velocity_pm)
        allocate (velocity_pm(3, NXpm_fine, NYpm_fine, NZpm_fine))
        velocity_pm = 0.d0
    end subroutine set_pm_velocities_zero

    subroutine set_pm_deformations_zero
        if (allocated(deform_pm)) deallocate (deform_pm)
        allocate (deform_pm(3, NXpm_fine, NYpm_fine, NZpm_fine))
        deform_pm = 0.d0
    end subroutine set_pm_deformations_zero

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! GETTERS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine get_NXpm(NX_pm_out) bind(C, name='get_NX_pm')
        use, intrinsic :: iso_c_binding, only: c_int
        implicit none
        integer(c_int), intent(out) :: NX_pm_out
        NX_pm_out = NXpm_fine
    end subroutine get_NXpm

    subroutine get_NYpm(NY_pm_out) bind(C, name='get_NY_pm')
        use, intrinsic :: iso_c_binding, only: c_int
        implicit none
        integer(c_int), intent(out) :: NY_pm_out
        NY_pm_out = NYpm_fine
    end subroutine get_NYpm

    subroutine get_NZpm(NZ_pm_out) bind(C, name='get_NZ_pm')
        use, intrinsic :: iso_c_binding, only: c_int
        implicit none
        integer(c_int), intent(out) :: NZ_pm_out
        NZ_pm_out = NZpm_fine
    end subroutine get_NZpm

    subroutine get_dpm(DXpm_out, DYpm_out, DZpm_out) bind(C, name='get_dpm')
        use, intrinsic :: iso_c_binding, only: c_double
        implicit none
        real(c_double), intent(out) :: DXpm_out, DYpm_out, DZpm_out
        DXpm_out = DXpm
        DYpm_out = DYpm
        DZpm_out = DZpm
    end subroutine get_dpm

    subroutine get_Xbounds(XMIN_pm_out, XMAX_pm_out) bind(C, name='get_Xbounds')
        use, intrinsic :: iso_c_binding, only: c_double
        implicit none
        real(c_double), intent(out) :: XMIN_pm_out, XMAX_pm_out
        XMIN_pm_out = XMIN_pm
        XMAX_pm_out = XMAX_pm
    end subroutine get_Xbounds

    subroutine get_Ybounds(YMIN_pm_out, YMAX_pm_out) bind(C, name='get_Ybounds')
        use, intrinsic :: iso_c_binding, only: c_double
        implicit none
        real(c_double), intent(out) :: YMIN_pm_out, YMAX_pm_out
        YMIN_pm_out = YMIN_pm
        YMAX_pm_out = YMAX_pm
    end subroutine get_Ybounds

    subroutine get_Zbounds(ZMIN_pm_out, ZMAX_pm_out) bind(C, name='get_Zbounds')
        use, intrinsic :: iso_c_binding, only: c_double
        implicit none
        real(c_double), intent(out) :: ZMIN_pm_out, ZMAX_pm_out
        ZMIN_pm_out = ZMIN_pm
        ZMAX_pm_out = ZMAX_pm
    end subroutine get_Zbounds

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! SETTERS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine set_RHS_pm(RHS_pm_in, size1, size2, size3, size4) bind(C, name='set_RHS_pm')
        use, intrinsic :: iso_c_binding, only: c_int, c_double
        implicit none
        integer(c_int), intent(in) :: size1, size2, size3, size4
        real(c_double), dimension(size1, size2, size3, size4), intent(in) :: RHS_pm_in
        ! Deallocate previous RHS_pm if it is allocated
        if (allocated(RHS_pm)) deallocate (RHS_pm)
        ! Allocate RHS_pm to match the shape of RHS_pm_in
        allocate (RHS_pm(size(RHS_pm_in, 1), size(RHS_pm_in, 2), size(RHS_pm_in, 3), size(RHS_pm_in, 4)))
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
        call dp_4d_alloc_info("velocity_pm", velocity_pm)
        call dp_4d_alloc_info("deform_pm", deform_pm)
        call dp_4d_alloc_info("RHS_pm", RHS_pm)
        print *, achar(9), "XMIN_pm: ", XMIN_pm
        print *, achar(9), "XMAX_pm: ", XMAX_pm
        print *, achar(9), "YMIN_pm: ", YMIN_pm
        print *, achar(9), "YMAX_pm: ", YMAX_pm
        print *, achar(9), "ZMIN_pm: ", ZMIN_pm
        print *, achar(9), "ZMAX_pm: ", ZMAX_pm
        print *, achar(9), "DXpm (fine): ", DXpm
        print *, achar(9), "DYpm (fine): ", DYpm
        print *, achar(9), "DZpm (fine): ", DZpm
        print *, achar(9), "DVpm (fine): ", DVpm
        print *, achar(9), "NXpm (fine): ", NXpm_fine
        print *, achar(9), "NYpm (fine): ", NYpm_fine
        print *, achar(9), "NZpm (fine): ", NZpm_fine
        print *, achar(9), "NXs_bl (fine)", NXs_fine_bl
        print *, achar(9), "NYs_bl (fine)", NYs_fine_bl
        print *, achar(9), "NXf_bl (fine)", NXf_fine_bl
        print *, achar(9), "NYf_bl (fine)", NYf_fine_bl
        print *, achar(9), "NZs_bl (fine)", NZs_fine_bl
        print *, achar(9), "NZf_bl (fine)", NZf_fine_bl
        print '(A,I6)', achar(9)//"nbound    =", nbound
        print '(A,I6)', achar(9)//"ndumcell  =", ndumcell
        print '(A,I6)', achar(9)//"IDVPM     =", IDVPM
        call dp_4d_alloc_info("SOL_pm", SOL_pm)
    end subroutine print_pmgrid_info

    subroutine print_RHS_pm(RHS)
        use console_io, only: dp_4d_alloc_info
        implicit none
        real(dp), intent(in), allocatable :: RHS(:, :, :, :)
        call dp_4d_alloc_info("RHS_pm", RHS)
    end subroutine print_RHS_pm

    subroutine print_velocity_stats
        use console_io, only: vpm_print, nocolor, red, VERBOCITY, dummy_string
        integer :: i, j, k
        if (VERBOCITY >= 2) then
            write (dummy_string, "(A)") 'Velocity U=∇xψ+∇ϕ info in the PM:'
            call vpm_print(dummy_string, red, 1)
            write (dummy_string, "(A,A, F8.3, A, F8.3, A, F8.3)") &
                achar(9)//"X:", &
                achar(9)//"Max:", maxval(velocity_pm(1, :, :, :)), &
                achar(9)//"Min:", minval(velocity_pm(1, :, :, :)), &
                achar(9)//"Mean:", sum(velocity_pm(1, :, :, :))/NXpm_fine/NYpm_fine/NZpm_fine
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, "(A,A, F8.3, A, F8.3, A, F8.3)") &
                achar(9)//"Y:", &
                achar(9)//"Max:", maxval(velocity_pm(2, :, :, :)), &
                achar(9)//"Min:", minval(velocity_pm(2, :, :, :)), &
                achar(9)//"Mean:", sum(velocity_pm(2, :, :, :))/NXpm_fine/NYpm_fine/NZpm_fine
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, "(A,A, F8.3, A, F8.3, A, F8.3)") &
                achar(9)//"Z:", &
                achar(9)//"Max:", maxval(velocity_pm(3, :, :, :)), &
                achar(9)//"Min:", minval(velocity_pm(3, :, :, :)), &
                achar(9)//"Mean:", sum(velocity_pm(3, :, :, :))/NXpm_fine/NYpm_fine/NZpm_fine
            call vpm_print(dummy_string, nocolor, 2)
            ! CHeck for nan values
            if (any(isnan(velocity_pm))) then
                ! Print Velx
                do i = 1, NXpm_fine
                    do j = 1, NYpm_fine
                        do k = 1, NZpm_fine
                            if (isnan(velocity_pm(1, i, j, k)) .or. &
                                isnan(velocity_pm(2, i, j, k)) .or. &
                                isnan(velocity_pm(3, i, j, k))) then
                                write (*, "(A, I3, A, I3, A, I3)") &
                                    achar(9)//"I:", i, achar(9)//"J:", j, achar(9)//"K:", k, achar(9)
                                stop
                            end if
                        end do
                    end do
                end do
                write (*, *) achar(9), 'VPM: NAN VALUES IN VELOCITY'
                stop
            end if
        end if
    end subroutine print_velocity_stats

    subroutine print_vortex_stretching_stats
        use console_io, only: vpm_print, nocolor, red, VERBOCITY, dummy_string
        integer :: i, j, k
        if (VERBOCITY >= 2) then
            write (dummy_string, "(A)") 'Vortex Stretching (w∇)·u info in the PM:'
            call vpm_print(dummy_string, red, 1)
            write (dummy_string, "(A,A, F8.3, A, F8.3, A, F8.3)") &
                achar(9)//"X:", &
                achar(9)//"Max:", maxval(deform_pm(1, :, :, :)), &
                achar(9)//"Min:", minval(deform_pm(1, :, :, :)), &
                achar(9)//"Mean:", sum(deform_pm(1, :, :, :))/NXpm_fine/NYpm_fine/NZpm_fine
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, "(A,A, F8.3, A, F8.3, A, F8.3)") &
                achar(9)//"Y:", &
                achar(9)//"Max:", maxval(deform_pm(2, :, :, :)), &
                achar(9)//"Min:", minval(deform_pm(2, :, :, :)), &
                achar(9)//"Mean:", sum(deform_pm(2, :, :, :))/NXpm_fine/NYpm_fine/NZpm_fine
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, "(A,A, F8.3, A, F8.3, A, F8.3)") &
                achar(9)//"Z:", &
                achar(9)//"Max:", maxval(deform_pm(3, :, :, :)), &
                achar(9)//"Min:", minval(deform_pm(3, :, :, :)), &
                achar(9)//"Mean:", sum(deform_pm(3, :, :, :))/NXpm_fine/NYpm_fine/NZpm_fine
            call vpm_print(dummy_string, nocolor, 2)
            ! CHeck for nan values
            if (any(isnan(deform_pm))) then
                ! Print Velx
                do i = 1, NXpm_fine
                    do j = 1, NYpm_fine
                        do k = 1, NZpm_fine
                            if (isnan(deform_pm(1, i, j, k)) .or. &
                                isnan(deform_pm(2, i, j, k)) .or. &
                                isnan(deform_pm(3, i, j, k))) then
                                write (*, "(A, I3, A, I3, A, I3)") &
                                    achar(9)//"I:", i, achar(9)//"J:", j, achar(9)//"K:", k, achar(9)
                                stop
                            end if
                        end do
                    end do
                end do
                write (*, *) achar(9), 'VPM: NAN VALUES IN VELOCITY'
                stop
            end if
        end if
    end subroutine print_vortex_stretching_stats
end module pmgrid

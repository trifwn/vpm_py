module console_io
    use vpm_types, only: dp, cartesian_grid
    integer, save :: VERBOCITY = 2
    ! Verbose level
    ! 0: No output
    ! 1: Basic output
    ! 2: Detailed output
    ! 3: Debug output

    integer, save :: tab_level = 0
    character(len=400), save :: dummy_string = ""
    public :: vpm_print, set_verbose_level
    public :: tab_level, dummy_string, VERBOCITY

    ! Print allocatable integer arrays
    public:: i_1d_alloc_info, i_2d_alloc_info
    ! Print fixed-size integer arrays
    public:: i_1d_array_info, i_2d_array_info

    ! Print allocatable double precision arrays
    public:: dp_1d_alloc_info, dp_2d_alloc_info, dp_3d_alloc_info, dp_4d_alloc_info
    ! Print pointer double precision arrays
    public:: dp_1d_ptr_info, dp_2d_ptr_info
    ! Print fixed-size double precision arrays
    public:: dp_1d_array_info, dp_2d_array_info

    integer, parameter :: RED = 31
    integer, parameter :: GREEN = 32
    integer, parameter :: BLUE = 34
    integer, parameter :: YELLOW = 33
    integer, parameter :: MAGENTA = 35
    integer, parameter :: CYAN = 36
    integer, parameter :: NOCOLOR = -1
contains
    subroutine vpm_print(msg, color, importance)
        character(len=*), intent(in) :: msg
        integer, intent(in) :: color
        integer, intent(in) :: importance
        character(len=256) :: tabbed_msg
        integer :: i

        ! Prepend tabs based on tab_level
        if (importance > VERBOCITY) then
            return
        end if
        tabbed_msg = ''
        do i = 1, tab_level
            tabbed_msg = trim(tabbed_msg)//char(9)
        end do
        tabbed_msg = trim(tabbed_msg)//trim(msg)
        ! Print the message with the specified color
        if (color == 31) then
            write (*, "(A)") achar(27)//'[1;31m'//trim(tabbed_msg)//achar(27)//'[0m'
        else if (color == 32) then
            write (*, "(A)") achar(27)//'[1;32m'//trim(tabbed_msg)//achar(27)//'[0m'
        else if (color == 33) then
            write (*, "(A)") achar(27)//'[1;33m'//trim(tabbed_msg)//achar(27)//'[0m'
        else if (color == 34) then
            write (*, "(A)") achar(27)//'[1;34m'//trim(tabbed_msg)//achar(27)//'[0m'
        else if (color == 35) then
            write (*, "(A)") achar(27)//'[1;35m'//trim(tabbed_msg)//achar(27)//'[0m'
        else if (color == 36) then
            write (*, "(A)") achar(27)//'[1;36m'//trim(tabbed_msg)//achar(27)//'[0m'
        else if (color == -1) then
            write (*, "(A)") trim(tabbed_msg)
        else
            write (*, "(A)") trim(tabbed_msg)
        end if
    end subroutine vpm_print

    subroutine set_verbose_level(level)
        integer, intent(in) :: level
        VERBOCITY = level
    end subroutine set_verbose_level
    
    !! FIELD INFORMATION !!
    subroutine print_stats_rank4(grid, field, field_name)
        implicit none
        type(cartesian_grid), intent(in) :: grid
        real(dp), target, intent(in)     :: field(3, grid%NN(1), grid%NN(2), grid%NN(3))
        character(len=*), intent(in)     :: field_name
        real(dp)                         :: min_val(3), max_val(3), mean_val(3), std_val(3)
        integer                          :: i, N, nx, ny, nz

        N = grid%NN(1)*grid%NN(2)*grid%NN(3)
        nx = grid%NN(1)
        ny = grid%NN(2)
        nz = grid%NN(3)
        ! Calculate the min, max and mean values for each component
        do i = 1, 3
            min_val(i) = minval(   field(i, 3:nx-2, 3:ny-2, 3:nz-2))
            max_val(i) = maxval(   field(i, 3:nx-2, 3:ny-2, 3:nz-2))
            mean_val(i) = sum(     field(i, 3:nx-2, 3:ny-2, 3:nz-2)) / N
            std_val(i) = sqrt(sum((field(i, 3:nx-2, 3:ny-2, 3:nz-2) - mean_val(i))**2) / N)
        end do

        print *, achar(27)//'[1;33m'//field_name//achar(27)//'[0m'
        print *, '------------------------------------------------'
        print *, 'Component  |   Min Value   |   Max Value   |   Mean Value  |   Std Dev'
        print *, '------------------------------------------------'
        print '(A10,4F15.8)', 'X', min_val(1), max_val(1), mean_val(1), std_val(1)
        print '(A10,4F15.8)', 'Y', min_val(2), max_val(2), mean_val(2), std_val(2)
        print '(A10,4F15.8)', 'Z', min_val(3), max_val(3), mean_val(3), std_val(3)
        print *, '------------------------------------------------'
    end subroutine print_stats_rank4 

    subroutine print_stats_rank3(grid, field, field_name)
        implicit none
        type(cartesian_grid), intent(in) :: grid
        real(dp), target, intent(in)     :: field(grid%NN(1), grid%NN(2), grid%NN(3))
        character(len=*), intent(in)     :: field_name
        real(dp)                         :: min_val, max_val, mean_val, std_val 
        integer                          :: N
        ! Calculate the min, max and mean values for each component
        N = grid%NN(1)*grid%NN(2)*grid%NN(3)

        min_val = minval(field)
        max_val = maxval(field)
        mean_val = sum(field)/ N
        std_val = sqrt(sum((field - mean_val)**2)/N)

        print *, achar(27)//'[1;33m'//field_name//achar(27)//'[0m' 
        print *, '------------------------------------------------'
        print *, 'Min Value   |   Max Value   |   Mean Value   |   Std Dev'
        print *, '------------------------------------------------'
        print '(4E15.6)', min_val, max_val, mean_val, std_val
        print *, '------------------------------------------------'
    end subroutine print_stats_rank3


    !! TIMESTEP INFORMATION !!
    subroutine print_timestep_information(timestep_info)
        use vpm_types, only: timestepInformation, solveInformation
        implicit none
        type(timestepInformation), intent(in) :: timestep_info
        real(dp)           :: mean_div_w, max_div_w, min_div_w
        real(dp)           :: mean_div_u, max_div_u, min_div_u
        real(dp)           :: total_kinetic_energy, total_enstrophy
        real(dp)           :: total_momentum_x, total_momentum_y, total_momentum_z
        integer            :: VERBOCITY_PRINT
        
        mean_div_w = timestep_info%mean_div_w
        max_div_w  = timestep_info%max_div_w
        min_div_w  = timestep_info%min_div_w
        
        mean_div_u = timestep_info%mean_div_u
        max_div_u  = timestep_info%max_div_u
        min_div_u  = timestep_info%min_div_u
        
        total_momentum_x = timestep_info%total_momentum_x_pm 
        total_momentum_y = timestep_info%total_momentum_y_pm 
        total_momentum_z = timestep_info%total_momentum_z_pm 
        total_kinetic_energy = timestep_info%total_kinetic_energy_pm
        total_enstrophy = timestep_info%total_enstrophy_pm

        VERBOCITY_PRINT = 0

        write(dummy_string, *) ""
        call vpm_print(dummy_string, nocolor, VERBOCITY_PRINT)

        write(dummy_string, "(A)") 'Divergence of the Ψ field'
        call vpm_print(dummy_string, yellow, VERBOCITY_PRINT)
        
        write(dummy_string, "(A,E11.4,A,E11.4,A,E11.4)") achar(9)//"div(ω)"//achar(9)// &
            achar(9)//" min : ", timestep_info%min_div_w, &
            achar(9)//" max : ", timestep_info%max_div_w, &
            achar(9)//" mean: ", timestep_info%mean_div_w
        call vpm_print(dummy_string, blue, VERBOCITY_PRINT)

        write(dummy_string, "(A)") 'Divergence of the velocity field'
        call vpm_print(dummy_string, yellow, VERBOCITY_PRINT)
        write(dummy_string, "(A,E11.4,A,E11.4,A,E11.4)") achar(9)//"div(u)"//achar(9)// &
            achar(9)//" min : ", min_div_u, &
            achar(9)//" max : ", max_div_u, &
            achar(9)//" mean: ", mean_div_u
        call vpm_print(dummy_string, blue, VERBOCITY_PRINT)

        write(dummy_string, "(A)") 'Total Momentum in the domain'
        call vpm_print(dummy_string, yellow, VERBOCITY_PRINT)
        write(dummy_string, "(A,E11.4)") achar(9)//'Total Momentum x : ', total_momentum_x
        call vpm_print(dummy_string, blue, VERBOCITY_PRINT)
        write(dummy_string, "(A,E11.4)") achar(9)//'Total Momentum y : ', total_momentum_y
        call vpm_print(dummy_string, blue, VERBOCITY_PRINT)
        write(dummy_string, "(A,E11.4)") achar(9)//'Total Momentum z : ', total_momentum_z
        call vpm_print(dummy_string, blue, VERBOCITY_PRINT)

        write(dummy_string, "(A)") 'Total Enstrophy in the domain'
        call vpm_print(dummy_string, yellow, VERBOCITY_PRINT)

        write(dummy_string, "(A,E11.4)") achar(9)//'sum(Enstrophy) : ', total_enstrophy
        call vpm_print(dummy_string, blue, VERBOCITY_PRINT)

        write (dummy_string, *) ""
        call vpm_print(dummy_string, nocolor, VERBOCITY_PRINT)
    end subroutine print_timestep_information

    subroutine print_solve_information(solve_info)
        use vpm_types, only: solveInformation
        implicit none
        type(solveInformation), intent(in) :: solve_info
        integer :: VERBOCITY_PRINT, i

        VERBOCITY_PRINT = 0
        write (dummy_string, "(A)") 'Residuals of the solution'
        call vpm_print(dummy_string, yellow, VERBOCITY)
        do i = 1, size(solve_info%f_min)
            ! For each equation write the laplacian - RHS_pm
            write (dummy_string, "(A, I3, A)") '   Equation =', i, ":   Δf = RHS"
            call vpm_print(dummy_string, blue, VERBOCITY_PRINT)
            write (dummy_string, "(A, E11.4, A, E11.4, A, E11.4)") achar(9)//'Forcing (RHS)'// &
                achar(9)//'min : ', solve_info%f_min(i), & 
                achar(9)//'max : ', solve_info%f_max(i), & 
                achar(9)//'mean: ', solve_info%f_mean(i) 
            call vpm_print(dummy_string, nocolor, VERBOCITY_PRINT)
            write (dummy_string, "(A, E11.4, A, E11.4, A, E11.4)") achar(9)//"Solution"// &
                achar(9)//'min : ', solve_info%sol_min(i), & 
                achar(9)//'max : ', solve_info%sol_max(i), & 
                achar(9)//'mean: ', solve_info%sol_mean(i) 
            call vpm_print(dummy_string, nocolor, VERBOCITY_PRINT)
            write (dummy_string, "(A, E11.4, A, E11.4, A, E11.4)") achar(9)//'Res:=Δf-RHS'// &
                achar(9)//'min : ', solve_info%residual_min(i), & 
                achar(9)//'max : ', solve_info%residual_max(i), & 
                achar(9)//'mean: ', solve_info%residual_mean(i) 
            call vpm_print(dummy_string, nocolor, VERBOCITY_PRINT)
        end do
    end subroutine print_solve_information

    !!! INTEGER ARRAYS !!!
    subroutine i_1d_alloc_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        integer, dimension(:), allocatable, intent(in) :: arr
        integer, dimension(4) :: sample_values
        character(len=100) :: name

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        if (allocated(arr)) then
            print *, achar(9)//trim(name), " (1D): Size = (", size(arr), ")"
            sample_values = arr(1:min(4, size(arr)))
            print '(A,4I12)', achar(9)//"Sample values: ", sample_values
            ! Print the number of non-zero elements
            print *, achar(9)//"Number of non-zero elements: ", count(arr /= 0)
            ! Print the mean/max/min values
            print *, achar(9)//"Mean value: ", sum(arr)/size(arr)
            print *, achar(9)//"Max value: ", maxval(arr)
            print *, achar(9)//"Min value: ", minval(arr)
        else
            print *, achar(9)//trim(name), " NOT ALLOCATED"
        end if
        print '(A)', ""
    end subroutine i_1d_alloc_info

    subroutine i_1d_array_info(name_in, arr, size_arr)
        integer, intent(in) :: size_arr
        character(len=*), intent(in) :: name_in
        character(len=100) :: name
        integer, dimension(size_arr), intent(in) :: arr

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        print *, achar(9)//trim(name), " (1D): Size = (", size(arr), ")"
        print '(A,4I12)', achar(9)//"Sample values: ", arr(1:min(4, size(arr)))
        ! Print the number of non-zero elements
        print *, achar(9)//"Number of non-zero elements: ", count(arr /= 0)
        ! Print the mean/max/min values
        print *, achar(9)//"Mean value: ", sum(arr)/size(arr)
        print *, achar(9)//"Max value: ", maxval(arr)
        print *, achar(9)//"Min value: ", minval(arr)
        print '(A)', ""
    end subroutine i_1d_array_info

    subroutine i_2d_alloc_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        character(len=100) :: name
        integer, dimension(:, :), allocatable, intent(in) :: arr
        integer, dimension(4) :: sample_values

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        if (allocated(arr)) then
            print *, achar(9)//trim(name), " (2D): Size = (", size(arr, 1), ",", size(arr, 2), ")"
            sample_values = arr(1, 1:min(4, size(arr, 2)))
            print "(A,4I12)", achar(9)//"Sample values: ", sample_values
            ! Print the number of non-zero elements
            print *, achar(9)//"Number of non-zero elements: ", count(arr /= 0)
            ! Print the mean/max/min values
            print *, achar(9)//"Mean value: ", sum(arr)/size(arr)
            print *, achar(9)//"Max value: ", maxval(arr)
            print *, achar(9)//"Min value: ", minval(arr)
        else
            print '(A,A)', achar(9)//trim(name), ": Not allocated"
        end if
        print '(A)', ""
    end subroutine i_2d_alloc_info

    subroutine i_2d_array_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        integer, dimension(:, :), intent(in) :: arr
        character(len=100) :: name

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        print *, achar(9)//trim(name), " (2D): Size = (", size(arr, 1), ",", size(arr, 2), ")"
        print '(A,4I12)', achar(9)//"Sample values: ", arr(1, 1:min(4, size(arr, 2)))
        ! Print the number of non-zero elements
        print *, achar(9)//"Number of non-zero elements: ", count(arr /= 0)
        ! Print the mean/max/min values
        print *, achar(9)//"Mean value: ", sum(arr)/size(arr)
        print *, achar(9)//"Max value: ", maxval(arr)
        print *, achar(9)//"Min value: ", minval(arr)
        print '(A)', ""
    end subroutine i_2d_array_info

    ! DOUBLE PRECISION ARRAYS
    subroutine dp_1d_alloc_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        character(len=100) :: name
        real(dp), dimension(:), allocatable, intent(in) :: arr
        real(dp), dimension(4) :: sample_values

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        if (allocated(arr)) then
            print '(A,A,A,I0,A)', achar(9)//trim(name), " (1D): Size = (", size(arr), ")"
            sample_values = arr(1:min(4, size(arr)))
            print '(A,4F12.6)', achar(9)//"Sample values: ", sample_values
            ! Print the number of non-zero elements
            print *, achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
            ! Print the mean/max/min values
            print *, achar(9)//"Mean value: ", sum(arr)/size(arr)
            print *, achar(9)//"Max value: ", maxval(arr)
            print *, achar(9)//"Min value: ", minval(arr)
        else
            print '(A,A)', achar(9)//trim(name), ": Not allocated"
        end if
        print '(A)', ""
    end subroutine dp_1d_alloc_info

    subroutine dp_1d_array_info(name_in, arr, size_arr)
        integer, intent(in) :: size_arr
        character(len=*), intent(in) :: name_in
        real(dp), dimension(size_arr), intent(in) :: arr
        character(len=100) :: name

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        print *, achar(9)//trim(name), " (1D): Size = (", size(arr), ")"
        print '(A,4E12.4)', achar(9)//"Sample values: ", arr(1:min(4, size(arr)))
        print *, achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
        print *, achar(9)//"Mean value: ", sum(arr)/size(arr)
        print *, achar(9)//"Max value: ", maxval(arr)
        print '(A)', ""
    end subroutine dp_1d_array_info

    subroutine dp_1d_ptr_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        character(len=100) :: name
        real(dp), pointer, intent(in) :: arr(:)
        real(dp), dimension(4) :: sample_values

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        print *, achar(9)//trim(name), " (1D): Size = (", size(arr), ")"
        sample_values = arr(1:min(4, size(arr)))
        print '(A,4E12.4)', achar(9)//"Sample values: ", sample_values
        ! Print the number of non-zero elements
        print *, achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
        ! Print the mean/max/min values
        print *, achar(9)//"Mean value: ", sum(arr)/size(arr)
        print *, achar(9)//"Max value: ", maxval(arr)
        print *, achar(9)//"Min value: ", minval(arr)
        print '(A)', ""
    end subroutine dp_1d_ptr_info

    subroutine dp_2d_alloc_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        character(len=100) :: name
        real(dp), dimension(:, :), allocatable, intent(in) :: arr
        real(dp), dimension(4) :: sample_values

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        if (allocated(arr)) then
            print *, achar(9)//trim(name), " (2D): Size = (", size(arr, 1), ",", size(arr, 2), ")"
            sample_values = arr(1, 1:min(4, size(arr, 2)))
            print '(A,4F12.6)', achar(9)//"Sample values: ", sample_values
            ! Print the number of non-zero elements
            print *, achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
            ! Print the mean/max/min values
            print *, achar(9)//"Mean value: ", sum(arr)/size(arr)
            print *, achar(9)//"Max value: ", maxval(arr)
            print *, achar(9)//"Min value: ", minval(arr)
        else
            print '(A,A)', achar(9)//trim(name), ": Not allocated"
        end if
        print '(A)', ""
    end subroutine dp_2d_alloc_info

    subroutine dp_2d_ptr_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        character(len=100) :: name
        real(dp), pointer, intent(in) :: arr(:, :)
        real(dp), dimension(4) :: sample_values

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        print *, achar(9)//trim(name), " (2D): Size = (", size(arr, 1), ", ", size(arr, 2), ")"
        sample_values = arr(1, 1:min(4, size(arr, 2)))
        print '(A,4E12.4)', achar(9)//"Sample values: ", sample_values
        ! Print the number of non-zero elements
        print *, achar(9)//"Number of zero elements: ", count(arr == 0.0d0)
        print *, achar(9)//"Number of NaN elements: ", count(isnan(arr))
        print *, achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
        ! Print the mean/max/min values
        print *, achar(9)//"Mean value: ", sum(arr)/size(arr)
        print *, achar(9)//"Max value: ", maxval(arr)
        print *, achar(9)//"Min value: ", minval(arr)
        print '(A)', ""
    end subroutine dp_2d_ptr_info

    subroutine dp_2d_array_info(name_in, arr, size1, size2)
        character(len=*), intent(in) :: name_in
        character(len=100) :: name

        integer, intent(in) :: size1, size2
        real(dp), dimension(size1, size2), intent(in) :: arr

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        print *, achar(9)//trim(name), " (2D): Size = (", size(arr, 1), ",", size(arr, 2), ")"
        print '(A,4E12.4)', achar(9)//"Sample values: ", arr(1, 1:min(4, size(arr, 2)))
        ! Print the number of non-zero elements
        print *, achar(9)//"Number of zero elements: ", count(arr == 0.0d0)
        print *, achar(9)//"Number of NaN elements: ", count(isnan(arr))
        print *, achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
        ! Print the mean/max/min values
        print *, achar(9)//"Mean value: ", sum(arr)/size(arr)
        print *, achar(9)//"Max value: ", maxval(arr)
        print *, achar(9)//"Min value: ", minval(arr)
        print '(A)', ""
    end subroutine dp_2d_array_info

    subroutine dp_3d_alloc_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        character(len=100) :: name
        real(dp), dimension(:, :, :), allocatable, intent(in) :: arr
        real(dp) :: sample_value
        integer :: i, j, k

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        if (allocated(arr)) then
            print *, achar(9), trim(name), " (3D): Size = (", size(arr, 1), ",", size(arr, 2), ",", size(arr, 3), ")"
            print *, achar(9), "Sample values:"
            do i = 1, min(2, size(arr, 1))
            do j = 1, min(2, size(arr, 2))
                write (*, '(A,I0,A,I0,A)', advance='no') achar(9)//"arr(", i, ",", j, ",1:4) = "
                do k = 1, min(4, size(arr, 3))
                    sample_value = arr(i, j, k)
                    write (*, '(A,F12.6,A)', advance='no') achar(9), sample_value, ", "
                end do
                print *
            end do
            end do
            ! Print the number of non-zero elements
            print *, achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
            ! Print the mean/max/min values
            print *, achar(9)//"Mean value: ", sum(arr)/size(arr)
            print *, achar(9)//"Max value: ", maxval(arr)
            print *, achar(9)//"Min value: ", minval(arr)
        else
            print '(A,A)', achar(9)//trim(name), ": Not allocated"
        end if
        print '(A)', ""
    end subroutine dp_3d_alloc_info

    subroutine dp_4d_alloc_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        real(dp), dimension(:, :, :, :), allocatable, intent(in) :: arr

        if (allocated(arr)) then
            write (dummy_string, "(A)") trim(name_in)//" (4D):"
            call vpm_print(dummy_string, nocolor, 2)

            write (dummy_string, '(A,I3,A,I3,A,I3,A,I3,A)') achar(9)//"Size = (", size(arr, 1), &
                ",", size(arr, 2), ",", size(arr, 3), ",", size(arr, 4), ")"
            call vpm_print(dummy_string, nocolor, 2)

            ! Print the number of non-zero elements
            write (dummy_string, *) achar(9)//"Number of elements: ", size(arr)
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9)//"Number of zero elements: ", count(arr == 0.0d0)
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9)//"Number of NaN elements: ", count(isnan(arr))
            call vpm_print(dummy_string, nocolor, 2)
            ! Print the mean/max/min values
            write (dummy_string, *) achar(9)//"Mean value: ", sum(arr)/size(arr)
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9)//"Max value: ", maxval(arr)
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9)//"Min value: ", minval(arr)
            call vpm_print(dummy_string, nocolor, 2)
        else
            write (*, '(A,A)') achar(9)//achar(9)//trim(name_in), ": Not allocated"
        end if
        print '(A)', ""
    end subroutine dp_4d_alloc_info
end module console_io

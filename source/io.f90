module io
    use base_types, only: dp
    integer, save :: VERBOCITY = 2
    ! Verbose level
    ! 0: No output
    ! 1: Basic output
    ! 2: Detailed output
    ! 3: Debug output
    character (len=100), save :: particle_output_file_suffix = "particles.dat"
    character (len=100), save :: pm_output_file_suffix = "pm_output.dat"
    character (len=100), save :: vpm_write_folder = "results/"
    character (len=100), save :: vpm_speed_output_file = "speed_profile.dat"
    integer, save :: tab_level = 0
    character (len=400), save :: dummy_string = " "
    
    public :: vpm_print,set_verbose_level
    public :: particle_output_file_suffix, pm_output_file_suffix, vpm_write_folder, vpm_speed_output_file
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
        character(len=256) :: formatted_msg
        character(len=256) :: tabbed_msg
        integer :: i

        ! Prepend tabs based on tab_level
        if (importance > VERBOCITY) then
            return
        endif
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

    !!! INTEGER ARRAYS !!!
    subroutine i_1d_alloc_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        integer, dimension(:),allocatable, intent(in) :: arr
        integer, dimension(4) :: sample_values
        character(len=100) :: name

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        if (allocated(arr)) then
            print  *, achar(9)//trim(name), " (1D): Size = (", size(arr), ")"
            sample_values = arr(1:min(4,size(arr)))
            print '(A,4I12)', achar(9)//"Sample values: ", sample_values
            ! Print the number of non-zero elements
            print  *, achar(9)//"Number of non-zero elements: ", count(arr /= 0)
            ! Print the mean/max/min values
            print  *, achar(9)//"Mean value: ", sum(arr)/size(arr)
            print  *, achar(9)//"Max value: ", maxval(arr)
            print  *, achar(9)//"Min value: ", minval(arr)
        else
            print  *, achar(9)//trim(name), " NOT ALLOCATED"
        end if
        print '(A)', ""
    end subroutine i_1d_alloc_info
    
    subroutine i_1d_array_info(name_in, arr, size_arr)
        integer, intent(in) :: size_arr
        character(len=*), intent(in) :: name_in
        character(len=100) :: name
        integer, dimension(size_arr), intent(in) :: arr
        
        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        print  *, achar(9)//trim(name), " (1D): Size = (", size(arr), ")"
        print '(A,4I12)', achar(9)//"Sample values: ", arr(1:min(4,size(arr)))
        ! Print the number of non-zero elements
        print  *, achar(9)//"Number of non-zero elements: ", count(arr /= 0)
        ! Print the mean/max/min values
        print  *, achar(9)//"Mean value: ", sum(arr)/size(arr)
        print  *, achar(9)//"Max value: ", maxval(arr)
        print  *, achar(9)//"Min value: ", minval(arr)
        print '(A)', ""
    end subroutine i_1d_array_info

    subroutine i_2d_alloc_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        character(len=100) :: name
        integer, dimension(:,:), allocatable, intent(in) :: arr
        integer, dimension(4) :: sample_values
        
        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        if (allocated(arr)) then
            print *, achar(9)//trim(name), " (2D): Size = (", size(arr,1), ",", size(arr,2), ")"
            sample_values = arr(1,1:min(4,size(arr,2)))
            print "(A,4I12)", achar(9)//"Sample values: ", sample_values
            ! Print the number of non-zero elements
            print  *, achar(9)//"Number of non-zero elements: ", count(arr /= 0)
            ! Print the mean/max/min values
            print  *, achar(9)//"Mean value: ", sum(arr)/size(arr)
            print  *, achar(9)//"Max value: ", maxval(arr)
            print  *, achar(9)//"Min value: ", minval(arr)
        else
            print '(A,A)', achar(9)//trim(name), ": Not allocated"
        end if
        print '(A)', ""
    end subroutine i_2d_alloc_info 

    subroutine i_2d_array_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        integer, dimension(:,:), intent(in) :: arr
        character(len=100) :: name

        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        print  *, achar(9)//trim(name), " (2D): Size = (", size(arr,1), ",", size(arr,2), ")"
        print '(A,4I12)', achar(9)//"Sample values: ", arr(1,1:min(4,size(arr,2)))
        ! Print the number of non-zero elements
        print  *, achar(9)//"Number of non-zero elements: ", count(arr /= 0)
        ! Print the mean/max/min values
        print  *, achar(9)//"Mean value: ", sum(arr)/size(arr)
        print  *, achar(9)//"Max value: ", maxval(arr)
        print  *, achar(9)//"Min value: ", minval(arr)
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
            sample_values = arr(1:min(4,size(arr)))
            print '(A,4F12.6)', achar(9)//"Sample values: ", sample_values
            ! Print the number of non-zero elements
            print  *, achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
            ! Print the mean/max/min values
            print  *, achar(9)//"Mean value: ", sum(arr)/size(arr)
            print  *, achar(9)//"Max value: ", maxval(arr)
            print  *, achar(9)//"Min value: ", minval(arr)
        else
            print '(A,A)', achar(9)//trim(name), ": Not allocated"
        end if
        print '(A)', ""
    end subroutine dp_1d_alloc_info
    
    subroutine dp_1d_array_info(name_in,arr, size_arr)
        integer, intent(in) :: size_arr
        character(len=*), intent(in) :: name_in
        real(dp), dimension(size_arr), intent(in) :: arr
        character(len=100) :: name
        
        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        print  *, achar(9)//trim(name), " (1D): Size = (", size(arr), ")"
        print '(A,4E12.4)', achar(9)//"Sample values: ", arr(1:min(4,size(arr)))
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
        print  *, achar(9)//trim(name), " (1D): Size = (", size(arr), ")"
        sample_values =  arr(1:min(4,size(arr)))
        print '(A,4E12.4)', achar(9)//"Sample values: ", sample_values
        ! Print the number of non-zero elements
        print  *, achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
        ! Print the mean/max/min values
        print  *, achar(9)//"Mean value: ", sum(arr)/size(arr)
        print  *, achar(9)//"Max value: ", maxval(arr)
        print  *, achar(9)//"Min value: ", minval(arr)
        print '(A)', ""
    end subroutine dp_1d_ptr_info
    
    subroutine dp_2d_alloc_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        character(len=100) :: name
        real(dp), dimension(:,:), allocatable, intent(in) :: arr
        real(dp), dimension(4) :: sample_values
        
        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        if (allocated(arr)) then
            print *, achar(9)//trim(name), " (2D): Size = (", size(arr,1), ",", size(arr,2), ")"
            sample_values = arr(1,1:min(4,size(arr,2)))
            print '(A,4F12.6)', achar(9)//"Sample values: ", sample_values
            ! Print the number of non-zero elements
            print  *, achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
            ! Print the mean/max/min values
            print  *, achar(9)//"Mean value: ", sum(arr)/size(arr)
            print  *, achar(9)//"Max value: ", maxval(arr)
            print  *, achar(9)//"Min value: ", minval(arr)
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
        print  *, achar(9)//trim(name), " (2D): Size = (", size(arr, 1), ", ", size(arr, 2), ")"
        sample_values =  arr(1,1:min(4,size(arr,2)))
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
        print  *, achar(9)//trim(name), " (2D): Size = (", size(arr,1), ",", size(arr,2), ")"
        print '(A,4E12.4)', achar(9)//"Sample values: ", arr(1,1:min(4,size(arr,2)))
        ! Print the number of non-zero elements
        print *, achar(9)//"Number of zero elements: ", count(arr == 0.0d0)
        print *, achar(9)//"Number of NaN elements: ", count(isnan(arr))
        print *, achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
        ! Print the mean/max/min values
        print  *, achar(9)//"Mean value: ", sum(arr)/size(arr)
        print  *, achar(9)//"Max value: ", maxval(arr)
        print  *, achar(9)//"Min value: ", minval(arr)
        print '(A)', ""
    end subroutine dp_2d_array_info

    subroutine dp_3d_alloc_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        character(len=100) :: name
        real(dp), dimension(:,:,:), allocatable, intent(in) :: arr
        real(dp) :: sample_value
        integer :: i, j, k
        
        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        if (allocated(arr)) then
            print *, achar(9), trim(name), " (3D): Size = (", size(arr,1), ",", size(arr,2), ",", size(arr,3), ")"
            print *, achar(9), "Sample values:"
            do i = 1, min(2, size(arr,1))
            do j = 1, min(2, size(arr,2))
                write(*, '(A,I0,A,I0,A)', advance='no') achar(9)//"arr(", i, ",", j, ",1:4) = "
                do k = 1, min(4, size(arr,3))
                    sample_value = arr(i,j,k)
                    write(*, '(A,F12.6,A)', advance='no') achar(9), sample_value, ", "
                end do
                print *
            end do
            end do
            ! Print the number of non-zero elements
            print  *, achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
            ! Print the mean/max/min values
            print  *, achar(9)//"Mean value: ", sum(arr)/size(arr)
            print  *, achar(9)//"Max value: ", maxval(arr)
            print  *, achar(9)//"Min value: ", minval(arr)
        else
            print '(A,A)', achar(9)//trim(name), ": Not allocated"
        end if
        print '(A)', ""
    end subroutine dp_3d_alloc_info

    subroutine dp_4d_alloc_info(name_in, arr)
        character(len=*), intent(in) :: name_in
        real(dp), dimension(:,:,:,:), allocatable, intent(in) :: arr
        real(dp), dimension(4) :: sample_value
        
        if (allocated(arr)) then
            write(dummy_string, "(A)") trim(name_in)//" (4D):"
            call vpm_print(dummy_string, nocolor, 2)

            write(dummy_string, '(A,I3,A,I3,A,I3,A,I3,A)') achar(9)//"Size = (", size(arr,1),&
                     ",", size(arr,2), ",", size(arr,3), ",", size(arr,4), ")"
            call vpm_print(dummy_string, nocolor, 2)

            ! Print the number of non-zero elements
            write (dummy_string,*) achar(9)//"Number of elements: ", size(arr)
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string,*) achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string,*) achar(9)//"Number of zero elements: ", count(arr == 0.0d0)
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string,*) achar(9)//"Number of NaN elements: ", count(isnan(arr))
            call vpm_print(dummy_string, nocolor, 2)
            ! Print the mean/max/min values
            write (dummy_string,*) achar(9)//"Mean value: ", sum(arr)/size(arr)
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string,*) achar(9)//"Max value: ", maxval(arr)
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string,*) achar(9)//"Min value: ", minval(arr)
            call vpm_print(dummy_string, nocolor, 2)
        else
            write (*,'(A,A)') achar(9)//achar(9)//trim(name_in), ": Not allocated"
        end if
        print '(A)', ""
    end subroutine dp_4d_alloc_info
end module IO
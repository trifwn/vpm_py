
module IO
    use base_types, only: dp
    integer, save :: VERBOCITY = 1
    ! Verbose level
    ! 0: No output
    ! 1: Basic output
    ! 2: Detailed output
    ! 3: Debug output
    character (len=100), save :: particle_output_file_suffix = "particles"
    character (len=100), save :: pm_output_file_suffix = "pm_output"
    character (len=100), save :: vpm_write_folder = "results"
    character (len=100), save :: vpm_speed_output_file = "profile_speed.dat"
    integer, save :: tab_level = 0
    character (len=400), save :: dummy_string = " "

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
    integer, parameter :: YELLOW = 33
    integer, parameter :: BLUE = 34
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
            print *, "Importance: ", importance, " Verbocity: ", VERBOCITY
            return
        endif
        tabbed_msg = ''
        do i = 1, tab_level
            tabbed_msg = tabbed_msg//achar(9)
        end do
        tabbed_msg = tabbed_msg//msg

        tabbed_msg = trim(tabbed_msg)
        print *, "PRINTING"
        print *, "Color: ", color
        print *, "Message: ", tabbed_msg
        ! Print the message with the specified color
        if (color == 31) then
            write (*, "(A)") achar(27)//'[1;31m'//tabbed_msg//achar(27)//'[0m'
        else if (color == 32) then
            write (*, "(A)") achar(27)//'[1;32m'//tabbed_msg//achar(27)//'[0m'
        else if (color == 33) then
            write (*, "(A)") achar(27)//'[1;33m'//tabbed_msg//achar(27)//'[0m'
        else if (color == 34) then
            write (*, "(A)") achar(27)//'[1;34m'//tabbed_msg//achar(27)//'[0m'
        else if (color == 35) then
            write (*, "(A)") achar(27)//'[1;35m'//tabbed_msg//achar(27)//'[0m'
        else if (color == 36) then
            write (*, "(A)") achar(27)//'[1;36m'//tabbed_msg//achar(27)//'[0m'
        else if (color == -1) then
            write (*, "(A)") tabbed_msg
        else
            write (*, "(A)") tabbed_msg
        end if
        print *, "PRINTED"
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
        character(len=100) :: name
        real(dp), dimension(:,:,:,:), allocatable, intent(in) :: arr
        real(dp), dimension(4) :: sample_value
        
        name = achar(27)//'[1;33m'//name_in//achar(27)//'[0m'
        if (allocated(arr)) then
            write(*, "(A)") achar(9)//trim(name)//" (4D):"
            write(*, '(A,I3,A,I3,A,I3,A,I3,A)') achar(9)//achar(9)//"Size = (", size(arr,1),&
                     ",", size(arr,2), ",", size(arr,3), ",", size(arr,4), ")"
            sample_value(1:min(4,size(arr,4))) = arr(1,1,1,1:min(4,size(arr,4)))
            print '(A,4F12.6)', achar(9)//achar(9)//"Sample values: ", sample_value(1: min(4, size(arr, 4)))
            ! Print the number of non-zero elements
            print *, achar(9)//achar(9)//"Number of elements: ", size(arr)
            print *, achar(9)//achar(9)//"Number of non-zero elements: ", count(arr /= 0.0d0)
            print *, achar(9)//achar(9)//"Number of zero elements: ", count(arr == 0.0d0)
            print *, achar(9)//achar(9)//"Number of NaN elements: ", count(isnan(arr))
            ! Print the mean/max/min values
            print *, achar(9)//achar(9)//"Mean value: ", sum(arr)/size(arr)
            print *, achar(9)//achar(9)//"Max value: ", maxval(arr)
            print *, achar(9)//achar(9)//"Min value: ", minval(arr)
        else
            print '(A,A)', achar(9)//achar(9)//trim(name), ": Not allocated"
        end if
        print '(A)', ""
    end subroutine dp_4d_alloc_info
end module IO
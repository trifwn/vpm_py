module ND_Arrays
    use ISO_C_BINDING
    use vpm_types, only: dp
    implicit none

    type, bind(C) :: ND_Array
        integer(c_int)  :: ndims                  ! Number of dimensions
        integer(c_int)  :: total_size             ! Total size of the data array
        type(C_PTR)     :: shape_ptr              ! Pointer to the shape array
        type(C_PTR)     :: data_ptr               ! Pointer to the actual data
        integer(c_int)  :: own_data               ! Does the array own the data?
        ! Optional members
        ! integer :: ref_count                    ! Reference count
        ! logical :: is_view                      ! Is the array a view?
        ! type(ND_Array), pointer :: parent_arr   ! Pointer to the parent array
        ! integer :: data_location                ! 0: Fortran, 1: C
        ! real(dp) , allocatable  :: data_arr(:)
    end type ND_Array

    ! Create a new type that includes the data array it should inherit from the ND_Array type
    ! type, extends(ND_Array) :: F_Array
    !     type(C_PTR) :: shape_ptr
    !     type(C_PTR) :: data_ptr
    !     integer :: ndims
    !     integer :: total_size
    !     real(dp), dimension(:) :: data
    ! end type F_Array

    type, bind(C) :: Pointers
        type(C_PTR) :: shape_ptr
        type(C_PTR) :: data_ptr
    end type Pointers

    interface size
        module procedure get_size
    end interface size

    interface shape
        module procedure get_shape
    end interface shape

    interface operator(+)
        module procedure add__op
    end interface operator(+)

    interface operator(-)
        module procedure sub__op
    end interface operator(-)

    interface operator(*)
        module procedure mul__op
    end interface operator(*)

    interface operator(/)
        module procedure div__op
    end interface operator(/)

    interface print
        module procedure print_arr
    end interface print

    interface assignment(=)
        module procedure assign_array
    end interface assignment(=)

contains

    !!!!!!!!!!!!!!!!!! INITIALIZATION AND DEALLOCATION !!!!!!!!!!!!!!!!!!!!!!!!!!
    function allocate_and_create(ndims, shape_ptr) result(arr)
        implicit none
        integer, intent(in) :: ndims
        type(C_PTR), intent(in) :: shape_ptr

        ! Local variables
        type(ND_Array) :: arr
        integer(c_int), dimension(:), pointer :: shape_arr
        real(c_double), dimension(:), pointer, contiguous :: data_arr
        integer :: total_size, i

        ! Get the shape array
        call C_F_POINTER(shape_ptr, shape_arr, [ndims])
        total_size = 1
        do i = 1, ndims
            total_size = total_size*shape_arr(i)
        end do

        allocate (data_arr(total_size))
        ! Assign the C pointers
        arr%ndims = ndims
        arr%total_size = total_size
        arr%shape_ptr = C_LOC(data_arr)
        arr%data_ptr = C_LOC(shape_arr)
        arr%own_data = 1 
    end function allocate_and_create

    function create_dtype(ndims, total_size, arr_shape, arr_data, own_data) result(arr) bind(C, name='create_dtype')
        implicit none
        integer(c_int), intent(in) :: ndims
        integer(c_int), intent(in) :: total_size
        integer(c_int), dimension(ndims), intent(in), target :: arr_shape
        real(c_double), dimension(total_size), intent(in), target :: arr_data
        integer(c_int), intent(in) :: own_data
        type(ND_Array) :: arr

        ! Assign the C pointers
        arr%ndims = ndims
        arr%total_size = total_size
        arr%shape_ptr = C_LOC(arr_shape)
        arr%data_ptr = C_LOC(arr_data)
        arr%own_data = own_data
    end function create_dtype

    subroutine free_array(arr) bind(C, name="free_array")
        implicit none
        type(ND_Array), intent(inout) :: arr
        real(c_double), pointer, contiguous :: data_arr(:)

        ! Free data if allocated
        if (arr%own_data.eq.1) then
            ! Check for null pointer
            call C_F_POINTER(arr%data_ptr, data_arr, [arr%total_size])
            if (associated(data_arr)) then
                deallocate(data_arr)
                print * , 'Data deallocated'
            endif
        end if
            
        ! Reset pointers to NULL
        arr%data_ptr = C_NULL_PTR
        arr%shape_ptr = C_NULL_PTR
    end subroutine free_array

    !!!!!!!!!!!!!!!!!!END INITIALIZATION AND DEALLOCATION !!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!! ACCESSING AND MANIPULATING ARRAYS !!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_data_ptr(arr) result(data_) bind(C, name='get_data_ptr')
        implicit none
        type(ND_Array), intent(in) :: arr
        type(C_PTR) :: data_

        data_ = arr%data_ptr
    end function get_data_ptr

    function get_shape_ptr(arr) result(shape_) bind(C, name='get_shape_ptr')
        implicit none
        type(ND_Array), intent(in) :: arr
        type(C_PTR) :: shape_

        shape_ = arr%shape_ptr
    end function get_shape_ptr

    subroutine set_element(arr, subscripts, elem)
        integer, dimension(:), intent(in) :: subscripts
        type(ND_Array), intent(inout) :: arr
        real(dp), intent(in) :: elem
        integer :: offset, i
        real(c_double), dimension(:), pointer :: data_
        integer(c_int), dimension(:), pointer :: f_shape

        call C_F_POINTER(arr%shape_ptr, f_shape, [arr%ndims])
        call C_F_POINTER(arr%data_ptr, data_, [arr%total_size])

        ! Calculate the offset
        offset = 0
        do i = 1, arr%ndims
            if (subscripts(i) < 1 .or. subscripts(i) > f_shape(i)) then
                error stop 'Subscript out of bounds'
            end if
            offset = offset*f_shape(i) + (subscripts(i) - 1)
        end do

        data_(offset + 1) = elem
    end subroutine set_element

    function slice_array(arr, start, end_) result(slice)
        implicit none
        type(ND_Array), intent(in) :: arr
        type(ND_Array) :: slice
        integer, intent(in) :: start, end_
        real(c_double), dimension(:), pointer :: data_, slice_data
        integer(c_int), dimension(:), pointer :: shape_, slice_shape
        integer :: slice_size

        call C_F_POINTER(arr%shape_ptr, shape_, [arr%ndims])
        call C_F_POINTER(arr%data_ptr, data_, [arr%total_size])

        slice_size = end_ - start + 1
        slice%ndims = 1
        slice%total_size = slice_size
        allocate (slice_shape(1))
        slice_shape(1) = slice_size
        slice%shape_ptr = C_LOC(slice_shape)

        allocate (slice_data(slice_size))
        slice%data_ptr = C_LOC(slice_data)

        slice_data = data_(start:end_)

    end function slice_array

    function get_element(arr, subscripts) result(elem)
        implicit none
        type(ND_Array), intent(in) :: arr
        integer, dimension(:), intent(in) :: subscripts
        real(c_double) :: elem
        integer :: offset, i
        real(c_double), dimension(:), pointer :: data_
        integer(c_int), dimension(:), pointer :: shape_

        call C_F_POINTER(arr%shape_ptr, shape_, [arr%ndims])
        call C_F_POINTER(arr%data_ptr, data_, [arr%total_size])

        ! Calculate the offset for linear index
        offset = 0
        do i = 1, arr%ndims
            offset = offset*shape_(i) + (subscripts(i) - 1)
        end do

        elem = data_(offset + 1)
    end function get_element

    function get_size(arr) result(s)
        integer :: s
        type(ND_Array), intent(in) :: arr

        s = arr%total_size

    end function get_size

    function get_shape(arr) result(s)
        integer, dimension(:), pointer :: s
        type(ND_Array), intent(in) :: arr
        integer(c_int), dimension(:), pointer :: f_shape

        call C_F_POINTER(arr%shape_ptr, f_shape, [arr%ndims])
        s => f_shape  ! Associate with the temporary Fortran array

    end function get_shape

    subroutine assign_array(arr1, arr2)
        ! subroutine to assign the array1 to array2. It does not create a new array.
        ! It justs changes the pointers of array2 to the pointers of array1.
        type(ND_Array), intent(out) :: arr1
        type(ND_Array), intent(in) :: arr2

        ! Deallocate the data array if allocated
        call free_array(arr1)

        ! Assign the C pointers
        arr1%ndims = arr2%ndims
        arr1%total_size = arr2%total_size
        arr1%shape_ptr = arr2%shape_ptr
        arr1%data_ptr = arr2%data_ptr

    end subroutine assign_array
    !!!!!!!!!!!!!!!!!!END ACCESSING AND MANIPULATING ARRAYS !!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!! OPERATIONS ON ARRAYS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine add__(arr1, arr2, arr_res) bind(C, name='add')
        implicit none
        type(ND_Array), intent(in) :: arr1, arr2
        type(ND_Array), intent(out) :: arr_res
        real(c_double), dimension(:), pointer :: data1, data2, data_res
        integer(c_int), dimension(:), pointer :: shape1, shape2
        integer :: i

        call C_F_POINTER(arr1%shape_ptr, shape1, [arr1%ndims])
        call C_F_POINTER(arr2%shape_ptr, shape2, [arr2%ndims])

        if (arr1%ndims /= arr2%ndims) error stop 'Arrays must have the same number of dimensions.'
        do i = 1, arr1%ndims
            if (shape1(i) /= shape2(i)) error stop 'Arrays must have the same shape.'
        end do

        call C_F_POINTER(arr1%data_ptr, data1, [arr1%total_size])
        call C_F_POINTER(arr2%data_ptr, data2, [arr2%total_size])
        call C_F_POINTER(arr_res%data_ptr, data_res, [arr_res%total_size])

        !$OMP PARALLEL SHARED(data1, data2, data_res) PRIVATE(i)
        !$OMP DO
        do i = 1, arr1%total_size
            data_res(i) = data1(i) + data2(i)
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine add__

    subroutine sub__(arr1, arr2, arr_res) bind(C, name='sub')
        implicit none
        type(ND_Array), intent(in) :: arr1, arr2
        type(ND_Array), intent(out) :: arr_res
        real(c_double), dimension(:), pointer :: data1, data2, data_res
        integer(c_int), dimension(:), pointer :: shape1, shape2
        integer :: i

        call C_F_POINTER(arr1%shape_ptr, shape1, [arr1%ndims])
        call C_F_POINTER(arr2%shape_ptr, shape2, [arr2%ndims])

        if (arr1%ndims /= arr2%ndims) error stop 'Arrays must have the same number of dimensions.'
        do i = 1, arr1%ndims
            if (shape1(i) /= shape2(i)) error stop 'Arrays must have the same shape.'
        end do

        call C_F_POINTER(arr1%data_ptr, data1, [arr1%total_size])
        call C_F_POINTER(arr2%data_ptr, data2, [arr2%total_size])
        call C_F_POINTER(arr_res%data_ptr, data_res, [arr_res%total_size])

        !$OMP PARALLEL SHARED(data1, data2, data_res) PRIVATE(i)
        !$OMP DO
        do i = 1, arr1%total_size
            data_res(i) = data1(i) - data2(i)
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine sub__

    subroutine mul__(arr1, arr2, arr_res) bind(C, name='mul')
        implicit none
        type(ND_Array), intent(in) :: arr1, arr2
        type(ND_Array), intent(inout) :: arr_res
        real(c_double), dimension(:), pointer :: data1, data2, data_res
        integer(c_int), dimension(:), pointer :: shape1, shape2
        integer :: i

        call C_F_POINTER(arr1%shape_ptr, shape1, [arr1%ndims])
        call C_F_POINTER(arr2%shape_ptr, shape2, [arr2%ndims])

        if (arr1%ndims /= arr2%ndims) error stop 'Arrays must have the same number of dimensions.'
        do i = 1, arr1%ndims
            if (shape1(i) /= shape2(i)) error stop 'Arrays must have the same shape.'
        end do

        call C_F_POINTER(arr_res%data_ptr, data_res, [arr_res%total_size])
        call C_F_POINTER(arr1%data_ptr, data1, [arr1%total_size])
        call C_F_POINTER(arr2%data_ptr, data2, [arr2%total_size])

        !$OMP PARALLEL SHARED(data1, data2, data_res) PRIVATE(i)
        !$OMP DO
        do i = 1, arr1%total_size
            data_res(i) = data1(i)*data2(i)
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine mul__

    subroutine div__(arr1, arr2, arr_res) bind(C, name='div')
        implicit none
        type(ND_Array), intent(in) :: arr1, arr2
        type(ND_Array), intent(inout) :: arr_res
        real(c_double), dimension(:), pointer :: data1, data2, data_res
        integer(c_int), dimension(:), pointer :: shape1, shape2
        integer :: i

        call C_F_POINTER(arr1%shape_ptr, shape1, [arr1%ndims])
        call C_F_POINTER(arr2%shape_ptr, shape2, [arr2%ndims])

        if (arr1%ndims /= arr2%ndims) error stop 'Arrays must have the same number of dimensions.'
        do i = 1, arr1%ndims
            if (shape1(i) /= shape2(i)) error stop 'Arrays must have the same shape.'
        end do

        call C_F_POINTER(arr_res%data_ptr, data_res, [arr_res%total_size])
        call C_F_POINTER(arr1%data_ptr, data1, [arr1%total_size])
        call C_F_POINTER(arr2%data_ptr, data2, [arr2%total_size])

        !$OMP PARALLEL SHARED(data1, data2, data_res) PRIVATE(i)
        !$OMP DO
        do i = 1, arr1%total_size
            data_res(i) = data1(i)/data2(i)
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine div__

    function add__op(arr1, arr2) result(arr_res)
        type(ND_Array), intent(in) :: arr1, arr2
        type(ND_Array) :: arr_res

        arr_res = allocate_and_create(arr1%ndims, arr1%shape_ptr)
        call add__(arr1, arr2, arr_res)
    end function add__op

    function sub__op(arr1, arr2) result(arr_res)
        type(ND_Array), intent(in) :: arr1, arr2
        type(ND_Array) :: arr_res

        arr_res = allocate_and_create(arr1%ndims, arr1%shape_ptr)
        call sub__(arr1, arr2, arr_res)

    end function sub__op

    function mul__op(arr1, arr2) result(arr_res)
        type(ND_Array), intent(in) :: arr1, arr2
        type(ND_Array) :: arr_res

        arr_res = allocate_and_create(arr1%ndims, arr1%shape_ptr)
        call mul__(arr1, arr2, arr_res)
    end function mul__op

    function div__op(arr1, arr2) result(arr_res)
        type(ND_Array), intent(in) :: arr1, arr2
        type(ND_Array) :: arr_res

        arr_res = allocate_and_create(arr1%ndims, arr1%shape_ptr)
        call div__(arr1, arr2, arr_res)
    end function div__op

    !!!!!!!!!!!!!!!!!!!!!END OPERATIONS ON ARRAYS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!! CONVERSION FUNCTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function from_intrinsic(data_array, shape_arr) result(arr)
        implicit none
        real(c_double), dimension(..), intent(in), target :: data_array ! SHOULD BE ASSUMED-RANK INSTEAD OF ASSUMED-SIZE
        integer(c_int), dimension(:), intent(in), target :: shape_arr
        ! real(c_double), pointer :: first_element
        type(ND_Array) :: arr

        ! Local variables
        integer :: i, total_size
        integer, dimension(:), pointer :: shape_ptr
        real(c_double), dimension(:), pointer :: data_ptr
        integer(c_int) :: own_data = 0

        total_size = 1
        do i = 1, size(shape_arr)
            total_size = total_size*shape_arr(i)
        end do

        allocate (shape_ptr(size(shape_arr)))
        do i = 1, size(shape_arr)
            shape_ptr(i) = shape_arr(i)
        end do

        ! I want to create a pointer that points to the data_array
        ! Since the data_array is assumed to be contiguous, I can just point to the first element
        ! of the data_array
        ! first_element => data_array(1)
        call C_F_POINTER(C_LOC(data_array), data_ptr, [total_size])
        
        arr%ndims = size(shape_arr)
        arr%total_size = total_size
        arr%shape_ptr = C_LOC(shape_ptr)
        arr%data_ptr = C_LOC(data_ptr)
        arr%own_data = own_data 
    end function from_intrinsic

    subroutine convert_to_1D_array(arr, f_1D)
        implicit none
        type(ND_Array), intent(in) :: arr
        real(c_double), dimension(:), pointer :: f_1D
        if (arr%ndims /= 1) error stop 'Array must have 4 dimensions.'
        call C_F_POINTER(arr%data_ptr, f_1D, [arr%total_size])
    end subroutine convert_to_1D_array

    subroutine convert_to_2D_array(arr, f_2D)
        implicit none
        type(ND_Array), intent(in) :: arr
        real(c_double), dimension(:, :), pointer :: f_2D
        integer(c_int), dimension(:), pointer :: f_shape
        if (arr%ndims /= 2) error stop 'Array must have 4 dimensions.'

        call C_F_POINTER(arr%shape_ptr, f_shape, [arr%ndims])
        call C_F_POINTER(arr%data_ptr, f_2D, f_shape)
    end subroutine convert_to_2D_array

    subroutine convert_to_3D_array(arr, f_3D)
        implicit none
        type(ND_Array), intent(in) :: arr
        real(c_double), dimension(:, :, :), pointer :: f_3D
        integer(c_int), dimension(:), pointer :: f_shape
        if (arr%ndims /= 3) error stop 'Array must have 4 dimensions.'

        call C_F_POINTER(arr%shape_ptr, f_shape, [arr%ndims])
        call C_F_POINTER(arr%data_ptr, f_3D, f_shape)
    end subroutine convert_to_3D_array

    subroutine convert_to_4D_array(arr, f_4D)
        implicit none
        type(ND_Array), intent(in) :: arr
        real(dp), pointer :: f_4D(:, :, :, :)
        integer(c_int), pointer :: f_shape(:)
        call C_F_POINTER(arr%shape_ptr, f_shape, [arr%ndims])
        call C_F_POINTER(arr%data_ptr, f_4D, f_shape)
    end subroutine convert_to_4D_array

    function to_intrinsic_flattened(arr) result(f_array)
        implicit none
        type(ND_Array), intent(in) :: arr
        real(c_double), dimension(:), pointer :: f_array
        call C_F_POINTER(arr%data_ptr, f_array, [arr%total_size])
    end function to_intrinsic_flattened

    !!!!!!!!!!!!!!!!!!!!!END CONVERSION FUNCTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!! PRINTING AND DEBUGGING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine print_arr(arr) bind(C, name='print_array')
        implicit none
        type(ND_Array), intent(in) :: arr
        real(c_double), dimension(:), pointer :: data_
        integer(c_int), dimension(:), pointer :: shape_
        integer :: i, j, k
        print *, 'Fortran:'
        print *, "Number of dimensions:", arr%ndims
        print *, "Total size:", arr%total_size
        if (arr%ndims .gt. 100) error stop 'Printing arrays with more than 100 dimensions is not supported.'
        if (arr%total_size .gt. 1000) error stop 'Printing arrays with more than 1000 elements is not supported.'
        if (arr%ndims < 1) error stop 'Array must have at least one dimension.'
        if (arr%total_size < 1) error stop 'Array must have at least one element.'
        ! Map the shape and data pointers
        call c_f_pointer(arr%shape_ptr, shape_, [arr%ndims])
        call c_f_pointer(arr%data_ptr, data_, [arr%total_size])

        print *, 'Data (shape:', shape_, '):'

        ! Handling printing based on dimensions
        if (arr%ndims == 1) then
            do i = 1, shape_(1)
                write (*, *) data_(i)
            end do
        else if (arr%ndims == 2) then
            do i = 1, shape_(1)
                write (*, '(1000F8.2)') (data_((i - 1)*shape_(2) + j), j=1, shape_(2))
            end do
        else if (arr%ndims == 3) then
            do i = 1, shape_(1)
                do j = 1, shape_(2)
                    do k = 1, shape_(3)
                        write (*, *) "(", i, j, k, "):", data_((i - 1)*shape_(2)*shape_(3) + (j - 1)*shape_(3) + k)
                    end do
                end do
            end do
        else
            print *, 'Printing arrays with more than 3 dimensions is not yet supported.'
        end if

        write (*, *) 'End of array'
        write (*, *) ''
    end subroutine print_arr
    !!!!!!!!!!!!!!!!!!!!!END PRINTING AND DEBUGGING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ND_Arrays

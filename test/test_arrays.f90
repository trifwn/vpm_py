program test_arrays
    use ND_Arrays
    use vpm_types, only: dp
    implicit none

    type(ND_Array) :: arr1, arr2, arr_add, arr_sub, arr_mul, arr_div, arr_slice!, test_ndarr
    real(dp) , dimension(:, :) , allocatable :: data1, data2
    integer, dimension(2) , target:: shape_arr
    integer :: i, j

    shape_arr = [3, 4]
    allocate(data1(shape_arr(1), shape_arr(2)))
    allocate(data2(shape_arr(1), shape_arr(2)))
    data1 = 0.0d0
    data2 = 0.0d0


    ! Create arrays
    arr1 = from_intrinsic(data1, shape(data1))
    arr2 = from_intrinsic(data2, shape(data2))
    
    ! Initialize arrays with some values
    call set_element(arr1, [1, 1], 1.0d0)
    call set_element(arr1, [1, 2], 2.0d0)
    call set_element(arr1, [1, 3], 3.0d0)
    call set_element(arr1, [1, 4], 4.0d0)
    call set_element(arr1, [2, 1], 5.0d0)
    call set_element(arr1, [2, 2], 6.0d0)
    call set_element(arr1, [2, 3], 7.0d0)
    call set_element(arr1, [2, 4], 8.0d0)
    call set_element(arr1, [3, 1], 9.0d0)
    call set_element(arr1, [3, 2], 10.0d0)
    call set_element(arr1, [3, 3], 11.0d0)
    call set_element(arr1, [3, 4], 12.0d0)

    print *, 'Array 1:'
    call print_arr(arr1)
    print *, "Internal Data:"
    do j = 1, shape_arr(2)
        do i = 1, shape_arr(1)
            print *, data1(i, j)
        end do
    end do

    call set_element(arr2, [1, 1], 12.0d0)
    call set_element(arr2, [1, 2], 11.0d0)
    call set_element(arr2, [1, 3], 10.0d0)
    call set_element(arr2, [1, 4], 9.0d0)
    call set_element(arr2, [2, 1], 8.0d0)
    call set_element(arr2, [2, 2], 7.0d0)
    call set_element(arr2, [2, 3], 6.0d0)
    call set_element(arr2, [2, 4], 5.0d0)
    call set_element(arr2, [3, 1], 4.0d0)
    call set_element(arr2, [3, 2], 3.0d0)
    call set_element(arr2, [3, 3], 2.0d0)
    call set_element(arr2, [3, 4], 1.0d0) 

    write (*,*) achar(27)//'[1;31m'//"PRINTING ARRAY 1"//achar(27)//'[0m'
    call print_arr(arr1)
    write (*,*) achar(27)//'[1;31m'//"PRINTING ARRAY 2"//achar(27)//'[0m'
    call print_arr(arr2)


    ! Perform operations
    write (*,*) achar(27)//'[1;31m'//"ARR1 + ARR2"//achar(27)//'[0m'
    arr_add = arr1 + arr2
    print *, 'Array Addition (arr1 + arr2):'
    call print_arr(arr_add)
    
    write (*,*) achar(27)//'[1;31m'//"ARR1 - ARR2"//achar(27)//'[0m'
    arr_sub = arr1 - arr2
    print *, 'Array Subtraction (arr1 - arr2):'
    call print_arr(arr_sub)
    
    write (*,*) achar(27)//'[1;31m'//"ARR1 * ARR2"//achar(27)//'[0m'
    arr_mul = arr1 * arr2
    print *, 'Array Multiplication (arr1 * arr2):'
    call print_arr(arr_mul)

    write (*,*) achar(27)//'[1;31m'//"ARR1 / ARR2"//achar(27)//'[0m'
    arr_div = arr1 / arr2
    print *, 'Array Division (arr1 / arr2):'
    call print_arr(arr_div)

    print *, 'Did the arrays change?'
    call print_arr(arr_add)
    
    ! Access and print individual elements
    print *, ""
    print *, ""
    print *, 'Element at (2, 2) in arr1: ', get_element(arr1, [2, 2])
    print *, 'Element at (2, 2) in arr2: ', get_element(arr2, [2, 2])
    
    ! Perform and print slicing operation
    arr_slice = slice_array(arr1, 1, 3)
    print *, 'Slice of arr1 from 1 to 2:'
    call print_arr(arr_slice)
    
    ! Cleanup
    print *, ''
    print *, ''
    call free_array(arr1)
    call free_array(arr2)
    call free_array(arr_add)
    call free_array(arr_sub)
    call free_array(arr_mul)
    call free_array(arr_div)
    call free_array(arr_slice)

end program test_arrays
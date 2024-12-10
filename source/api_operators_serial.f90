module vector_field_operators_bindings
    use, intrinsic :: iso_c_binding
    use serial_vector_field_operators
    use ND_Arrays
    use vpm_types, only: dp
    implicit none

contains

   !! OPERATORS SERIAL
    subroutine calc_derivative(f, dimensions, ds, order, direction, result) bind(C, name='derivative')
        use serial_vector_field_operators, only: calculate_derivative
        implicit none
        type(ND_Array), intent(in), target :: f
        integer(c_int), intent(in) :: dimensions
        real(c_double), intent(in), target :: ds(dimensions)
        integer(c_int), intent(in) :: order
        integer(c_int), intent(in) :: direction
        type(ND_Array), intent(out), target :: result
        ! Local variables
        integer :: ndim_f, ndim, order_f, direction_f
        real(dp) :: dx, dy, dz
        real(dp), pointer :: f_1D_ptr(:), f_2D_ptr(:, :), f_3D_ptr(:, :, :), f_4D_ptr(:, :, :, :)
        real(dp), allocatable, save :: df_1D(:), df_2D(:, :), df_3D(:, :, :), df_4D(:, :, :, :)

        ! Find the ND_Arrays dimensions
        ndim_f = f%ndims
        ndim = size(ds)
        order_f = order
        direction_f = direction
        select case (ndim_f)
        case (1)
            call convert_to_1D_array(f, f_1D_ptr)
            dx = ds(1)
            if (ndim /= 1) then
                print *, 'Error: 1D arrays require 1 element in the ds array'
                stop
            end if
            df_1D = calculate_derivative(f_1D_ptr, dx, order_f, direction_f)
            result = from_intrinsic(df_1D, shape(df_1D))
        case (2)
            call convert_to_2D_array(f, f_2D_ptr)
            if (ndim == 1) then
                dx = ds(1)
                df_2D = calculate_derivative(f_2D_ptr, dx, order_f, direction_f)
            else if (ndim == 2) then
                dx = ds(1)
                dy = ds(2)
                df_2D = calculate_derivative(f_2D_ptr, dx, dy, order_f, direction_f)
            else
                print *, 'Error: 2D arrays require 1 or 2 elements in the ds array'
                stop
            end if
            result = from_intrinsic(df_2D, shape(df_2D))
        case (3)
            call convert_to_3D_array(f, f_3D_ptr)
            if (ndim == 2) then
                dx = ds(1)
                dy = ds(2)
                df_3D = calculate_derivative(f_3D_ptr, dx, dy, order_f, direction_f)
            else if (ndim == 3) then
                dx = ds(1)
                dy = ds(2)
                dz = ds(3)
                df_3D = calculate_derivative(f_3D_ptr, dx, dy, dz, order_f, direction_f)
            else
                print *, 'Error: 3D arrays require 2 or 3 elements in the ds array'
                stop
            end if
            result = from_intrinsic(df_3D, shape(df_3D))
        case (4)
            call convert_to_4D_array(f, f_4D_ptr)
            if (ndim /= 3) then
                print *, 'Error: 4D arrays require 3 elements in the ds array'
                stop
            end if
            dx = ds(1)
            dy = ds(2)
            dz = ds(3)
            df_4D = calculate_derivative(f_4D_ptr, dx, dy, dz, order_f, direction_f)
            result = from_intrinsic(df_4D, shape(df_4D))
        case default
            print *, 'Error: ND_Arrays with more than 4 dimensions are not supported'
            stop
        end select

    end subroutine calc_derivative

    subroutine calc_divergence(f, dimensions, ds, result) bind(C, name='divergence')
        type(ND_Array), intent(in), target :: f
        integer(c_int), intent(in) :: dimensions
        real(c_double), intent(in), target :: ds(dimensions)
        type(ND_Array), intent(out), target :: result
        ! Local variables
        real(dp), pointer :: f_4D_ptr(:, :, :, :)
        real(dp), allocatable, save :: div_result(:, :, :)
        real(dp) :: dx, dy, dz

        if (dimensions /= 3 .or. f%ndims /= 4) then
            print *, 'Error: Divergence requires a 3D vector field (4D array) and 3 elements in ds'
            stop
        end if

        call convert_to_4D_array(f, f_4D_ptr)
        dx = ds(1)
        dy = ds(2)
        dz = ds(3)

        div_result = divergence(f_4D_ptr, dx, dy, dz)
        result = from_intrinsic(div_result, shape(div_result))
    end subroutine calc_divergence

    subroutine calc_curl(f, dimensions, ds, result) bind(C, name='curl')
        type(ND_Array), intent(in), target :: f
        integer(c_int), intent(in) :: dimensions
        real(c_double), intent(in), target :: ds(dimensions)
        type(ND_Array), intent(out), target :: result
        ! Local variables
        real(dp), pointer :: f_4D_ptr(:, :, :, :)
        real(dp), allocatable, save :: curl_result(:, :, :, :)
        real(dp) :: dx, dy, dz

        if (dimensions /= 3 .or. f%ndims /= 4) then
            print *, 'Error: Curl requires a 3D vector field (4D array) and 3 elements in ds'
            stop
        end if

        call convert_to_4D_array(f, f_4D_ptr)
        dx = ds(1)
        dy = ds(2)
        dz = ds(3)

        curl_result = curl(f_4D_ptr, dx, dy, dz)
        result = from_intrinsic(curl_result, shape(curl_result))
    end subroutine calc_curl

    subroutine calc_laplacian(f, dimensions, ds, result) bind(C, name='laplacian')
        type(ND_Array), intent(in), target :: f
        integer(c_int), intent(in) :: dimensions
        real(c_double), intent(in), target :: ds(dimensions)
        type(ND_Array), intent(out), target :: result
        ! Local variables
        real(dp), pointer :: f_4D_ptr(:, :, :, :)
        real(dp), pointer :: f_3D_ptr(:, :, :)
        real(dp), allocatable, save :: lap_result4D(:, :, :, :)
        real(dp), allocatable, save :: lap_result3D(:, :, :)
        real(dp) :: dx, dy, dz

        if (dimensions /= 3) then
            print *, 'Error: Laplacian requires a 3D field (4D array) and 3 elements in ds'
            stop
        end if
        dx = ds(1)
        dy = ds(2)
        dz = ds(3)

        if (f%ndims == 4) then
            call convert_to_4D_array(f, f_4D_ptr)
            lap_result4D = laplacian(f_4D_ptr, dx, dy, dz)
            result = from_intrinsic(lap_result4D, shape(lap_result4D))
        else if (f%ndims == 3) then
            call convert_to_3D_array(f, f_3D_ptr)
            lap_result3D = laplacian(f_3D_ptr, dx, dy, dz)
            result = from_intrinsic(lap_result3D, shape(lap_result3D))
        else
            print *, 'Error: Laplacian requires a 3D field (4D array) or a 3D scalar field (3D array)'
            stop
        end if

    end subroutine calc_laplacian

    subroutine calc_gradient(f, dimensions, ds, result) bind(C, name='gradient')
        type(ND_Array), intent(in), target :: f
        integer(c_int), intent(in) :: dimensions
        real(c_double), intent(in), target :: ds(dimensions)
        type(ND_Array), intent(out), target :: result
        ! Local variables
        real(dp), pointer :: f_3D_ptr(:, :, :)
        real(dp), allocatable, save :: grad_result(:, :, :, :)
        real(dp) :: dx, dy, dz

        if (dimensions /= 3 .or. f%ndims /= 3) then
            print *, 'Error: Gradient requires a 3D scalar field and 3 elements in ds'
            stop
        end if

        call convert_to_3D_array(f, f_3D_ptr)
        dx = ds(1)
        dy = ds(2)
        dz = ds(3)

        grad_result = gradient(f_3D_ptr, dx, dy, dz)
        result = from_intrinsic(grad_result, shape(grad_result))
    end subroutine calc_gradient

    subroutine calc_vector_laplacian(f, dimensions, ds, result) bind(C, name='vector_laplacian')
        type(ND_Array), intent(in), target :: f
        integer(c_int), intent(in) :: dimensions
        real(c_double), intent(in), target :: ds(dimensions)
        type(ND_Array), intent(out), target :: result
        ! Local variables
        real(dp), pointer :: f_4D_ptr(:, :, :, :)
        real(dp), allocatable, save :: vec_lap_result(:, :, :, :)
        real(dp) :: dx, dy, dz

        if (dimensions /= 3 .or. f%ndims /= 4) then
            print *, 'Error: Vector Laplacian requires a 3D vector field (4D array) and 3 elements in ds'
            stop
        end if

        call convert_to_4D_array(f, f_4D_ptr)
        dx = ds(1)
        dy = ds(2)
        dz = ds(3)

        vec_lap_result = calc_vector_laplacian_expansion(f_4D_ptr, dx, dy, dz)
        result = from_intrinsic(vec_lap_result, shape(vec_lap_result))
    end subroutine calc_vector_laplacian

    subroutine calc_hessian(f, dimensions, ds, result) bind(C, name='hessian')
        type(ND_Array), intent(in), target :: f
        integer(c_int), intent(in) :: dimensions
        real(c_double), intent(in), target :: ds(dimensions)
        type(ND_Array), intent(out), target :: result
        ! Local variables
        real(dp), pointer :: f_4D_ptr(:, :, :, :)
        real(dp), pointer :: f_3D_ptr(:, :, :)
        real(dp), allocatable, save :: hes_result_4D(:, :, :, :, :, :)
        real(dp), allocatable, save :: hes_result_3D(:, :, :, :, :)
        real(dp) :: dx, dy, dz

        if (dimensions /= 3) then
            print *, 'Error: Hessian requires a 3D scalar field (4D array with first dimension 1) and 3 elements in ds'
            stop
        end if
        dx = ds(1)
        dy = ds(2)
        dz = ds(3)

        if (f%ndims == 4) then
            call convert_to_4D_array(f, f_4D_ptr)
            hes_result_4D = hessian(f_4D_ptr, dx, dy, dz)
            result = from_intrinsic(hes_result_4D, shape(hes_result_4D))
        else if (f%ndims == 3) then
            call convert_to_3D_array(f, f_3D_ptr)
            hes_result_3D = hessian(f_3D_ptr, dx, dy, dz)
            result = from_intrinsic(hes_result_3D, shape(hes_result_3D))
        else
            print *, 'Error: Hessian requires a 3D scalar field (4D array with first dimension 1) or a 3D scalar field'// &
                '(3D array and returns the gradient)'
            stop
        end if
    end subroutine calc_hessian

    subroutine calc_jacobian(f, dimensions, ds, result) bind(C, name='jacobian')
        type(ND_Array), intent(in), target :: f
        integer(c_int), intent(in) :: dimensions
        real(c_double), intent(in), target :: ds(dimensions)
        type(ND_Array), intent(out), target :: result
        ! Local variables
        real(dp), pointer :: f_4D_ptr(:, :, :, :)
        real(dp), pointer :: f_3D_ptr(:, :, :)
        real(dp), allocatable, save :: jac_result_4D(:, :, :, :, :)
        real(dp), allocatable, save :: jac_result_3D(:, :, :, :)
        real(dp) :: dx, dy, dz

        if (dimensions /= 3) then
            print *, 'Error: Jacobian requires and 3 elements in ds'
            stop
        end if
        dx = ds(1)
        dy = ds(2)
        dz = ds(3)

        if (f%ndims == 4) then
            call convert_to_4D_array(f, f_4D_ptr)
            jac_result_4D = jacobian(f_4D_ptr, dx, dy, dz)
            result = from_intrinsic(jac_result_4D, shape(jac_result_4D))
        else if (f%ndims == 3) then
            call convert_to_3D_array(f, f_3D_ptr)
            jac_result_3D = jacobian(f_3D_ptr, dx, dy, dz)
            result = from_intrinsic(jac_result_3D, shape(jac_result_3D))
        else
            print *, 'Error: Jacobian requires a 3D field (4D array) or a 3D scalar field'// &
                '(3D array and returns the gradient)'
            stop
        end if

    end subroutine calc_jacobian

end module vector_field_operators_bindings

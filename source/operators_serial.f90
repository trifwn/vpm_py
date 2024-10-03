module serial_vector_field_operators
    use base_types, only: dp
    use omp_lib
    implicit none

    private
    public :: divergence, curl, laplacian, gradient, vector_laplacian
    public :: apply_differential_operator
    public :: calculate_derivative, hessian, jacobian

    interface calculate_derivative
        module function calculate_derivative_1d(f, dx, dim, order) result(df)
            implicit none
            real(dp), intent(in) :: f(:), dx
            integer, intent(in) :: dim, order
            real(dp), allocatable :: df(:)
        end function

        module function calculate_vector_derivative_1d(f, dx, dim, order) result(df)
            implicit none
            real(dp), intent(in) :: f(:,:), dx
            integer, intent(in) :: dim, order
            real(dp), allocatable :: df(:,:)
        end function

        module function calculate_derivative_2d(f, dx, dy, dim, order) result(df)
            implicit none
            real(dp), intent(in) :: f(:,:), dx, dy
            integer, intent(in) :: dim, order
            real(dp), allocatable :: df(:,:)
        end function

        module function calculate_vector_derivative_2d(f, dx, dy, dim, order) result(df)
            implicit none
            real(dp), intent(in) :: f(:,:,:), dx, dy
            integer, intent(in) :: dim, order
            real(dp), allocatable :: df(:,:,:)
        end function

        module function calculate_derivative_3d(f, dx, dy, dz, dim, order) result(df)
            implicit none
            real(dp), intent(in) :: f(:,:,:), dx, dy, dz
            integer, intent(in) :: dim, order
            real(dp), allocatable :: df(:,:,:)
        end function

        module function calculate_vector_derivative_3d(f, dx, dy, dz, dim, order) result(df)
            implicit none
            real(dp), intent(in) :: f(:,:,:,:), dx, dy, dz
            integer, intent(in) :: dim, order
            real(dp), allocatable :: df(:,:,:,:)
        end function
    end interface calculate_derivative

    interface hessian
        module function calculate_scalar_hessian(f, dx, dy, dz) result(hes)
            real(dp), intent(in) :: f(:,:,:), dx, dy, dz
            real(dp), allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
            real(dp), allocatable :: hes(:,:,:,:,:)
        end function calculate_scalar_hessian

        module function calculate_vector_hessian(f, dx, dy, dz) result(hes)
            real(dp), intent(in) :: f(:,:,:,:), dx, dy, dz
            real(dp), allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
            real(dp), allocatable :: hes(:,:,:,:,:,:)
        end function calculate_vector_hessian
    end interface hessian

    interface jacobian
        module function calculate_scalar_jacobian(f, dx, dy, dz) result(jac)
            real(dp), intent(in) :: f(:,:,:), dx, dy, dz
            real(dp), allocatable :: jac(:,:,:,:)
            integer :: nx, ny, nz
        end function calculate_scalar_jacobian
        
        module function calculate_vector_jacobian(f, dx, dy, dz) result(jac)
            real(dp), intent(in) :: f(:,:,:,:), dx, dy, dz
            real(dp), allocatable :: jac(:,:,:,:,:)
            integer :: nx, ny, nz, i, neq
        end function calculate_vector_jacobian
    end interface jacobian
contains
    module function calculate_derivative_1d(f, dx, dim, order) result(df)
        real(dp), intent(in) :: f(:), dx
        integer, intent(in) :: dim, order
        real(dp), allocatable :: df(:)
        integer :: n, i

        n = size(f)
        allocate(df(n))

        select case(order)
        case(1)  ! First-order derivative (second-order accurate)
            df(1) = (-3*f(1) + 4*f(2) - f(3)) / (2*dx)
            df(2:n-1) = (f(3:n) - f(1:n-2)) / (2*dx)
            df(n) = (3*f(n) - 4*f(n-1) + f(n-2)) / (2*dx)
        case(2)  ! Second-order derivative
            df(1) = (f(3) - 2*f(2) + f(1)) / (dx**2)
            df(2:n-1) = (f(3:n) - 2*f(2:n-1) + f(1:n-2)) / (dx**2)
            df(n) = (f(n) - 2*f(n-1) + f(n-2)) / (dx**2)
        end select
    end function calculate_derivative_1d

    module function calculate_derivative_2d(f, dx, dy, dim, order) result(df)
        real(dp), intent(in) :: f(:,:), dx, dy
        integer, intent(in) :: dim, order
        real(dp), allocatable :: df(:,:)
        integer :: nx, ny, i, j

        nx = size(f, 1)
        ny = size(f, 2)
        allocate(df(nx, ny))

        select case(dim)
        case(1)  ! x-derivative
            do j = 1, ny
                df(:,j) = calculate_derivative_1d(f(:,j), dx, 1, order)
            end do
        case(2)  ! y-derivative
            do i = 1, nx
                df(i,:) = calculate_derivative_1d(f(i,:), dy, 1, order)
            end do
        end select
    end function calculate_derivative_2d

    module function calculate_derivative_3d(f, dx, dy, dz, dim, order) result(df)
         real(dp), intent(in) :: f(:,:,:), dx, dy, dz
         integer, intent(in) :: dim, order
         real(dp), allocatable :: df(:,:,:)
         integer :: nx, ny, nz, i, j, k

         nx = size(f, 1)
         ny = size(f, 2)
         nz = size(f, 3)

         allocate(df(nx, ny, nz))

         !$omp parallel do collapse(3)
         select case(dim)
         case(1)  ! x-derivative
             do k = 1, nz
                 do j = 1, ny
                     df(:,j,k) = calculate_derivative_1d(f(:,j,k), dx, 1, order)
                 end do
             end do
         case(2)  ! y-derivative
               do k = 1, nz
                  do i = 1, nx
                        df(i,:,k) = calculate_derivative_1d(f(i,:,k), dy, 1, order)
                  end do
               end do
         case(3)  ! z-derivative
               do j = 1, ny
                  do i = 1, nx
                        df(i,j,:) = calculate_derivative_1d(f(i,j,:), dz, 1, order)
                  end do
               end do
         end select
         !$omp end parallel do
    end function calculate_derivative_3d

    module function calculate_vector_derivative_1d(f,dx, dim, order) result(df)
        real(dp), intent(in) :: f(:,:), dx
        integer, intent(in) :: dim, order
        real(dp), allocatable :: df(:,:)
        integer :: neq, n, i

        neq = size(f, 1)
        n = size(f, 2)
        allocate(df(neq, n))

        !$omp parallel do
        do i = 1, neq
            df(i,:) = calculate_derivative_1d(f(i,:), dx, dim, order)
        end do
        !$omp end parallel do
    end function calculate_vector_derivative_1d

    module function calculate_vector_derivative_2d(f, dx, dy, dim, order) result(df)
        real(dp), intent(in) :: f(:,:,:), dx, dy
        integer, intent(in) :: dim, order
        real(dp), allocatable :: df(:,:,:)
        integer :: neq, nx, ny, i

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        allocate(df(neq, nx, ny))

        !$omp parallel do
        do i = 1, neq
            select case(dim)
            case(1)  ! x-derivative
                df(i,:,:) = calculate_derivative_2d(f(i,:,:), dx, dy, 1, order)
            case(2)  ! y-derivative
                df(i,:,:) = calculate_derivative_2d(f(i,:,:), dx, dy, 2, order)
            end select
        end do
        !$omp end parallel do
    end function calculate_vector_derivative_2d

    module function calculate_vector_derivative_3d(f, dx, dy, dz, dim, order) result(df)
        real(dp), intent(in) :: f(:,:,:,:), dx, dy, dz
        integer, intent(in) :: dim, order
        real(dp), allocatable :: df(:,:,:,:)
        integer :: neq, nx, ny, nz, i, j, k, eq

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)
        allocate(df(neq, nx, ny, nz))
         
        !$omp parallel do
        do eq = 1, neq
            select case(dim)
            case(1)  ! x-derivative
                df(eq,:,:,:) = calculate_derivative_3d(f(eq,:,:,:), dx, dy, dz, 1, order)
            case(2)  ! y-derivative
                df(eq,:,:,:) = calculate_derivative_3d(f(eq,:,:,:), dx, dy, dz, 2, order)
            case(3)  ! z-derivative
                df(eq,:,:,:) = calculate_derivative_3d(f(eq,:,:,:), dx, dy, dz, 3, order)
            end select
        end do
        !$omp end parallel do

    end function calculate_vector_derivative_3d

    module function calculate_scalar_jacobian(f, dx,dy, dz) result(jac)
        real(dp), intent(in) :: f(:,:,:), dx, dy, dz
        real(dp), allocatable :: jac(:,:,:,:)
        integer :: nx, ny, nz

        nx = size(f, 1)
        ny = size(f, 2)
        nz = size(f, 3)

        allocate(jac(3, nx, ny, nz))
        jac(1,:,:,:) = calculate_derivative(f, dx, dy, dz, 1, 1)
        jac(2,:,:,:) = calculate_derivative(f, dx, dy, dz, 2, 1)
        jac(3,:,:,:) = calculate_derivative(f, dx, dy, dz, 3, 1)
    end function calculate_scalar_jacobian

    module function calculate_vector_jacobian(f, dx, dy, dz) result(jac)
        real(dp), intent(in) :: f(:,:,:,:), dx, dy, dz
        real(dp), allocatable :: jac(:,:,:,:,:)
        integer :: nx, ny, nz, i, neq

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)

        allocate(jac(neq, 3, nx, ny, nz))

        !$omp parallel do
        do i = 1, neq
            jac(i,:,:,:,:) = calculate_scalar_jacobian(f(i,:,:,:), dx, dy, dz)
        end do
        !$omp end parallel do
    end function calculate_vector_jacobian

    module function calculate_scalar_hessian(f, dx, dy, dz) result(hes)
        real(dp), intent(in) :: f(:,:,:), dx, dy, dz
        real(dp), allocatable :: hes(:,:,:,:,:)
        real(dp), allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)

        allocate(hes(3,3, size(f, 1), size(f, 2), size(f, 3)))

        fx = calculate_derivative(f, dx, dy, dz, 1, 1)
        fy = calculate_derivative(f, dx, dy, dz, 2, 1)
        fz = calculate_derivative(f, dx, dy, dz, 3, 1)

        hes(1,1,:,:,:) = calculate_derivative(fx, dx, dy, dz, 1, 1)
        hes(1,2,:,:,:) = calculate_derivative(fx, dx, dy, dz, 2, 1)
        hes(1,3,:,:,:) = calculate_derivative(fx, dx, dy, dz, 3, 1)
        hes(2,1,:,:,:) = calculate_derivative(fy, dx, dy, dz, 1, 1)
        hes(2,2,:,:,:) = calculate_derivative(fy, dx, dy, dz, 2, 1)
        hes(2,3,:,:,:) = calculate_derivative(fy, dx, dy, dz, 3, 1)
        hes(3,1,:,:,:) = calculate_derivative(fz, dx, dy, dz, 1, 1)
        hes(3,2,:,:,:) = calculate_derivative(fz, dx, dy, dz, 2, 1)
        hes(3,3,:,:,:) = calculate_derivative(fz, dx, dy, dz, 3, 1)
    end function calculate_scalar_hessian

    module function calculate_vector_hessian(f, dx, dy, dz) result(hes)
        real(dp), intent(in) :: f(:,:,:,:), dx, dy, dz
        real(dp), allocatable :: hes(:,:,:,:,:,:)
        real(dp), allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
        integer :: nx, ny, nz, i, neq

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)

        allocate(hes(neq, 3, 3, nx, ny, nz))

        !$omp parallel do
        do i = 1,neq
            hes(i,:,:,:,:,:) = calculate_scalar_hessian(f(i,:,:,:), dx, dy, dz)
        enddo
        !$omp end parallel do
    end function calculate_vector_hessian

    subroutine divergence(f, dx, dy, dz, result)
        real(dp), intent(in) :: f(:,:,:,:), dx, dy, dz
        real(dp), allocatable, intent(out) :: result(:,:,:,:)
        integer :: nx, ny, nz

        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)

        if (size(f, 1) /= 3) then
            print *, "Error: Divergence requires a 3D vector field"
            return
        end if

        allocate(result(1, nx, ny, nz))
        result(1,:,:,:) = calculate_derivative(f(1,:,:,:), dx, dy, dz, 1, 1) + &
                        calculate_derivative(f(2,:,:,:), dx, dy, dz, 2, 1) + &
                        calculate_derivative(f(3,:,:,:), dx, dy, dz, 3, 1)
    end subroutine divergence

    subroutine curl(f, dx, dy, dz, result)
        real(dp), intent(in) :: f(:,:,:,:), dx, dy, dz
        real(dp), allocatable, intent(out) :: result(:,:,:,:)
        real(dp), allocatable              :: Jtmp1(:,:,:), Jtmp2(:,:,:)
        integer :: nx, ny, nz

        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)

        if (size(f, 1) /= 3) then
            print *, "Error: Curl requires a 3D vector field"
            return
        end if

        allocate(result(3, nx, ny, nz))
        
        ! We need to calculate the curl of a vector field F = (F1, F2, F3)
        ! That means calculating the offdiagnoal components of the Jacobian matrix
        ! of the vector field F. The Jacobian matrix is given by:
        ! J = | dF1/dx  dF1/dy  dF1/dz |
        !     | dF2/dx  dF2/dy  dF2/dz |
        !     | dF3/dx  dF3/dy  dF3/dz |
        
        ! The curl of a vector field is given by:
        ! curl(F) = (dF3/dy - dF2/dz, dF1/dz - dF3/dx, dF2/dx - dF1/dy)
        ! This is equivalent to the following:
        ! curl(F) = (J32 - J23, J13 - J31, J21 - J12)
        Jtmp1 = calculate_derivative(f(3,:,:,:), dx, dy, dz, 2, 1) ! dF3/dy = J32
        Jtmp2 = calculate_derivative(f(2,:,:,:), dx, dy, dz, 3, 1) ! dF2/dz = J23
        result(1,:,:,:) = Jtmp1 - Jtmp2                            ! dF3/dy - dF2/dz = J32 - J23

        Jtmp1 = calculate_derivative(f(1,:,:,:), dx, dy, dz, 3, 1) ! dF1/dz = J13
        Jtmp2 = calculate_derivative(f(3,:,:,:), dx, dy, dz, 1, 1) ! dF3/dx = J31
        result(2,:,:,:) = Jtmp1 - Jtmp2                            ! dF1/dz - dF3/dx = J13 - J31
        
        Jtmp1 = calculate_derivative(f(2,:,:,:), dx, dy, dz, 1, 1) ! dF2/dx = J21
        Jtmp2 = calculate_derivative(f(1,:,:,:), dx, dy, dz, 2, 1) ! dF1/dy = J12
        result(3,:,:,:) = Jtmp1 - Jtmp2                            ! dF2/dx - dF1/dy = J21 - J12

    end subroutine curl

    subroutine laplacian(f, dx, dy, dz, result)
        real(dp), intent(in) :: f(:,:,:,:), dx, dy, dz
        real(dp), allocatable, intent(out) :: result(:,:,:,:)
        integer :: neq, nx, ny, nz, i

        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)

        allocate(result(neq, nx, ny, nz))
        !$omp parallel do
        do i = 1, neq
            result(i,:,:,:) = calculate_derivative(f(i,:,:,:), dx, dy, dz, 1, 2) + &
                            calculate_derivative(f(i,:,:,:), dx, dy, dz, 2, 2) + &
                            calculate_derivative(f(i,:,:,:), dx, dy, dz, 3, 2)
        end do
        !$omp end parallel do
    end subroutine laplacian

    subroutine gradient(f, dx, dy, dz, result)
        real(dp), intent(in) :: f(:,:,:,:), dx, dy, dz
        real(dp), allocatable, intent(out) :: result(:,:,:,:)
        integer :: nx, ny, nz

        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)

        allocate(result(3, nx, ny, nz))
        result(1,:,:,:) = calculate_derivative(f(1,:,:,:), dx, dy, dz, 1, 1)
        result(2,:,:,:) = calculate_derivative(f(1,:,:,:), dx, dy, dz, 2, 1)
        result(3,:,:,:) = calculate_derivative(f(1,:,:,:), dx, dy, dz, 3, 1)
    end subroutine gradient

    subroutine vector_laplacian(f, dx, dy, dz, result)
        real(dp), intent(in) :: f(:,:,:,:), dx, dy, dz
        real(dp), allocatable, intent(out) :: result(:,:,:,:)
        real(dp), allocatable :: temp(:,:,:,:), grad_div(:,:,:,:)

        call curl(f, dx, dy, dz, temp)
        call curl(temp, dx, dy, dz, result)

        call divergence(f, dx, dy, dz, temp)
        call gradient(temp, dx, dy, dz, grad_div)

        result = result - grad_div
    end subroutine vector_laplacian

    function apply_differential_operator(f, dx, dy, dz, operator_type) result(result)
        real(dp), intent(in) :: f(:,:,:,:), dx, dy, dz
        character(len=*), intent(in) :: operator_type
        real(dp), allocatable :: result(:,:,:,:)

        ! Size information
        integer :: neq, nx, ny, nz
        neq = size(f, 1)
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 4)

        ! Select the appropriate operator based on the operator_type
        select case(trim(operator_type))
        case ("divergence")
            call divergence(f, dx, dy, dz, result)
            
        case ("curl")
            call curl(f, dx, dy, dz, result)
            
        case ("laplacian")
            call laplacian(f, dx, dy, dz, result)
            
        case ("gradient")
            call gradient(f, dx, dy, dz, result)
            
        case ("vector_laplacian")
            call vector_laplacian(f, dx, dy, dz, result)

        case default
            ! Unknown operator type handling
            print *, "Error: Unknown operator type: ", trim(operator_type)
            stop
        end select
    end function apply_differential_operator
end module serial_vector_field_operators
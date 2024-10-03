module test_operators
   use iso_fortran_env, only: real64
   use MPI
   use data_communication
   use serial_vector_field_operators
   implicit none

   real(real64), parameter :: tol = 1.0e-6_real64
   integer :: my_rank, num_p, ierror
contains

   subroutine run_all_tests()
      integer :: total_tests, passed_tests
      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
      call MPI_Comm_size(MPI_COMM_WORLD, num_p, ierror)
      total_tests = 0
      passed_tests = 0
      
      if(my_rank == 0) then   
         print *, "===================================="
         print *, "Starting Operators Test Suite"
         print *, "===================================="
         print *
      end if
      call test_distribute_and_collect(total_tests, passed_tests)
      call test_scalar_derivative(total_tests, passed_tests)
      call test_vector_derivative(total_tests, passed_tests)
      call test_scalar_jacobian(total_tests, passed_tests)
      call test_vector_jacobian(total_tests, passed_tests)
      call test_scalar_hessian(total_tests, passed_tests)
      call test_vector_hessian(total_tests, passed_tests)
      call test_differential_operators(total_tests, passed_tests)
      if(my_rank == 0) then
         print *
         print *, "===================================="
         print *, "Test Suite Summary:"
         print *, "Total tests:  ", total_tests
         print *, "Passed tests: ", passed_tests
         print *, "Failed tests: ", total_tests - passed_tests
         print *, "===================================="
      end if
   end subroutine run_all_tests

   subroutine test_scalar_derivative(total_tests, passed_tests)
      integer, intent(inout) :: total_tests, passed_tests
      integer :: i, j, k, l 
      logical :: test_passed
      real(real64), allocatable :: local_f(:, :, :), local_df(:, :, :), collected_data(:,:,:)
      real(real64), allocatable :: scalar_f(:,:,:), scalar_df_analytic(:,:,:), scalar_df(:,:,:)
      real(real64) :: dx, dy, dz, x, y, z
      integer      :: neq, nx, ny, nz, local_nx, local_ny, local_nz
      character(len=1024) :: filename
      ! SCALAR F
      nx = 50
      ny = 50
      nz = 50
      dx = 0.1_real64
      dy = 1.2_real64
      dz = 0.3_real64
 
      ! if (my_rank == 0) then
         allocate(scalar_f(nx, ny, nz))
         allocate(scalar_df_analytic(nx, ny, nz))
         do l = 1, nz
            do k = 1, ny
               do j = 1, nx
                  x = real(j, real64) * dx
                  y = real(k, real64) * dy
                  z = real(l, real64) * dz
                  scalar_f(j,k,l) = x*x + y*y + z*z
                  scalar_df_analytic(j,k,l) = 2.0_real64 * x
               end do
            end do
         end do
      ! endif 

      if (my_rank == 0) then
         total_tests = total_tests + 1
         print *, "Test 1: Derivative Calculation"
         print *, "  - Calculating x-derivative..."
         scalar_df = calculate_derivative(scalar_f, dx, dy, dz, 1, 1)
         ! Check if the derivative is calculated correctly
         test_passed = .true.
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  if (abs(scalar_df(i,j,k) - scalar_df_analytic(i,j,k)) > tol) then
                     test_passed = .false.
                     exit
                  end if
               end do
               if (.not. test_passed) exit
            end do
            if (.not. test_passed) exit
         end do

         if (test_passed) then
            print *, "  - Serial Test passed: x-derivative calculated correctly."
            passed_tests = passed_tests + 1
         else
            print *, "  - Serial Test failed: Incorrect x-derivative calculation."
         end if
         deallocate(scalar_df)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierror)

      ! Parallel derivative calculation
      if (my_rank == 0) print *, "  - Calculating x-derivative in parallel..."
     
      call commit_array(scalar_f, [3, 2, 1], 1)
      call distribute_data(scalar_f, local_f)
      call collect_data(collected_data, local_f)
      scalar_df = calculate_derivative(local_f, dx, dy, dz, 1, 1)
      call collect_data(collected_data, scalar_df)
      call free_vars()

      if (my_rank == 0) then
         total_tests = total_tests + 1
         ! Check if the derivative is calculated correctly
         test_passed = .true.
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  if (abs(collected_data(i,j,k) - scalar_df_analytic(i,j,k)) > tol) then
                     test_passed = .false.
                     exit
                  end if
               end do
               if (.not. test_passed) exit
            end do
            if (.not. test_passed) exit
         end do

         if (test_passed) then
            print *, "  - Parallel Test passed: x-derivative calculated correctly."
            print *
            passed_tests = passed_tests + 1
         else
            print *, "  - Parallel Test failed: Incorrect x-derivative calculation."
            print *
         end if
      end if
      
      ! Deallocate arrays
      deallocate(scalar_df)
      deallocate(scalar_f)
      deallocate(scalar_df_analytic)
   end subroutine test_scalar_derivative

   subroutine test_vector_derivative(total_tests, passed_tests)
      integer, intent(inout) :: total_tests, passed_tests
      real(real64), allocatable :: f(:,:,:,:), df(:,:,:,:), real_df(:,:,:,:)
      real(real64) :: dx, dy, dz, x, y, z
      integer :: neq, nx, ny, nz, i, j, k, l
      logical :: test_passed

      if (my_rank == 0) then
         print *, "Test 2: Derivative Calculation (Vector Field)"
         
         total_tests = total_tests + 1

         neq = 3
         nx  = 6
         ny  = 6
         nz  = 6
         dx  = 0.1_real64
         dy  = 0.3_real64
         dz  = 0.5_real64

         allocate(f(3, nx, ny, nz))
         allocate(real_df(3, nx, ny, nz))

         do l = 1, nz
            do k = 1, ny
               do j = 1, nx
                  x = real(j, real64) * dx
                  y = real(k, real64) * dy
                  z = real(l, real64) * dz
                  ! Function assignments
                  f(1,j,k,l) = x*x + y*y + z*z
                  f(2,j,k,l) = x*y + y*z + z*x
                  f(3,j,k,l) = x + y + z
                  ! Gradient (partial derivatives) assignments
                  real_df(1,j,k,l) = 2.0 * x  ! ∂f(1)/∂x for f(1) = x² + y² + z²
                  real_df(2,j,k,l) = y + z    ! ∂f(2)/∂x for f(2) = xy + yz + zx
                  real_df(3,j,k,l) = 1.0      ! ∂f(3)/∂x for f(3) = x + y + z
               end do
            end do
         end do

         print *, "  - Calculating x-derivative of the vector field..."
         df = calculate_derivative(f, dx, dy, dz, 1, 1)

         test_passed = .true.
         do l = 1, nz
            do k = 1, ny
               do j = 1, nx
                  do i = 1, neq
                     if (abs(df(i,j,k,l) - real_df(i,j,k,l)) > tol) then
                        write(*,*) "Mismatch at (", i, j, k, l, "):", df(i,j,k,l), real_df(i,j,k,l)
                        test_passed = .false.
                        exit
                     end if
                  end do
                  ! if (.not. test_passed) exit
               end do
               ! if (.not. test_passed) exit
            end do
            ! if (.not. test_passed) exit
         end do

         if (test_passed) then
            print *, "  - Test passed: x-derivative calculated correctly."
            passed_tests = passed_tests + 1
         else
            print *, "  - Test failed: Incorrect x-derivative calculation."
         end if

         print *
         ! Deallocate arrays
         deallocate(f)
         deallocate(real_df)
         deallocate(df)
      end if
   end subroutine test_vector_derivative

   subroutine test_scalar_jacobian(total_tests, passed_tests)
      integer, intent(inout) :: total_tests, passed_tests
      real(real64), allocatable :: f(:,:,:), jac(:,:,:,:), analytic_jac(:,:,:,:)
      real(real64) :: dx, dy, dz, x, y, z
      integer :: nx, ny, nz, i, j, k, l
      logical :: test_passed

      if (my_rank == 0) then
         print *, "Test 3: Jacobian Calculation"

         total_tests = total_tests + 1

         nx = 5
         ny = 5
         nz = 5
         dx = 0.1_real64 
         dy = 0.1_real64
         dz = 0.1_real64

         allocate(f(nx, ny, nz))
         allocate(analytic_jac(3, nx, ny, nz))

         do l = 1, nz
            do k = 1, ny
               do j = 1, nx
                  x = real(j, real64) * dx
                  y = real(k, real64) * dy
                  z = real(l, real64) * dz
                  f(j,k,l) = x*x + y*y + z*z
                  analytic_jac(1,j,k,l) = 2.0_real64 * x
                  analytic_jac(2,j,k,l) = 2.0_real64 * y
                  analytic_jac(3,j,k,l) = 2.0_real64 * z
               end do
            end do
         end do

         print *, "  - Calculating Jacobian..."
         jac = jacobian(f, dx, dy, dz)

         test_passed = .true.
         do i = 1, 3
            do l = 1, nz
               do k = 1, ny
                  do j = 1, nx
                     if (abs(jac(i,j,k,l) - analytic_jac(i,j,k,l)) > tol) then
                        print *, "  - Mismatch at (", j, k, l, "):", jac(i,j,k,l), analytic_jac(i,j,k,l)
                        test_passed = .false.
                        exit
                     end if
                  end do
                  if (.not. test_passed) exit
               end do
               if (.not. test_passed) exit
            end do
         end do

         if (test_passed) then
            print *, "  - Test passed: Jacobian calculated correctly."
            passed_tests = passed_tests + 1
         else
            print *, "  - Test failed: Incorrect Jacobian calculation."
         end if

         print *
         ! Deallocate arrays
         deallocate(f)
         deallocate(analytic_jac)
         deallocate(jac)
      endif
   end subroutine test_scalar_jacobian

   subroutine test_vector_jacobian(total_tests, passed_tests)
      integer, intent(inout) :: total_tests, passed_tests
      real(real64), allocatable :: f(:,:,:,:), jac(:,:,:,:,:), analytic_jac(:,:,:,:,:)
      real(real64) :: dx, dy, dz, x, y, z
      integer :: nx, ny, nz, i, j, k, l, neq
      logical :: test_passed

      if (my_rank == 0) then
         print *, "Test 4: Jacobian Calculation (Vector Field)"

         total_tests = total_tests + 1

         neq = 3
         nx = 5
         ny = 5
         nz = 5
         dx = 1.0_real64
         dy = 1.0_real64
         dz = 1.0_real64

         allocate(f(3, nx, ny, nz))
         allocate(analytic_jac(3, 3, nx, ny, nz))

         do l = 1, nz
            do k = 1, ny
               do j = 1, nx
                  x = real(j, real64) * dx
                  y = real(k, real64) * dy
                  z = real(l, real64) * dz

                  ! EQ 1: f₁ = x² + y² + z²
                  f(1,j,k,l) = x*x + y*y + z*z
                  analytic_jac(1,1,j,k,l) = 2.0_real64 * x   ! ∂f₁/∂x
                  analytic_jac(1,2,j,k,l) = 2.0_real64 * y   ! ∂f₁/∂y
                  analytic_jac(1,3,j,k,l) = 2.0_real64 * z   ! ∂f₁/∂z

                  ! EQ 2: f₂ = xy + yz + zx
                  f(2,j,k,l) = x*y + y*z + z*x
                  analytic_jac(2,1,j,k,l) = y + z           ! ∂f₂/∂x
                  analytic_jac(2,2,j,k,l) = x + z           ! ∂f₂/∂y
                  analytic_jac(2,3,j,k,l) = x + y           ! ∂f₂/∂z

                  ! EQ 3: f₃ = x + y + z
                  f(3,j,k,l) = x + y + z
                  analytic_jac(3,1,j,k,l) = 1.0_real64      ! ∂f₃/∂x
                  analytic_jac(3,2,j,k,l) = 1.0_real64      ! ∂f₃/∂y
                  analytic_jac(3,3,j,k,l) = 1.0_real64      ! ∂f₃/∂z

               end do
            end do
         end do

         print *, "  - Calculating Jacobian of the vector field..."
         jac = jacobian(f, dx, dy, dz)

         test_passed = .true.
         do i = 1, 3
            do l = 1, nz
               do k = 1, ny
                  do j = 1, nx
                     if (abs(jac(i,i,j,k,l) - analytic_jac(i,i,j,k,l)) > tol) then
                        test_passed = .false.
                        exit
                     end if
                  end do
                  if (.not. test_passed) exit
               end do
               if (.not. test_passed) exit
            end do
         end do

         if (test_passed) then
            print *, "  - Test passed: Jacobian calculated correctly."
            passed_tests = passed_tests + 1
         else
            print *, "  - Test failed: Incorrect Jacobian calculation."
         end if
         print *
         ! Deallocate arrays
         deallocate(f)
         deallocate(analytic_jac)
         deallocate(jac)
      endif
   end subroutine test_vector_jacobian
   
   subroutine test_scalar_hessian(total_tests, passed_tests)
      integer, intent(inout) :: total_tests, passed_tests
      real(real64), allocatable :: f(:,:,:), hes(:,:,:,:,:), analytic_hes(:,:,:,:,:)
      real(real64) :: dx, dy, dz, x, y, z
      integer :: nx, ny, nz, i, j, k, d1, d2
      logical :: test_passed

      if (my_rank == 0) then
         print *, "Test 5: Hessian Calculation"

         total_tests = total_tests + 1

         nx = 5
         ny = 5
         nz = 5
         dx = 1.0_real64
         dy = 1.0_real64
         dz = 1.0_real64

         allocate(f(nx, ny, nz))
         allocate(analytic_hes(3, 3,nx, ny, nz))

         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  x = real(i, real64) * dx
                  y = real(j, real64) * dy
                  z = real(k, real64) * dz

                  f(i,j,k) = x*x + y*y + z*z + x*y + y*z + z*x
                  analytic_hes(1,1,i,j,k) = 2.0 * dx
                  analytic_hes(2,2,i,j,k) = 2.0 * dy
                  analytic_hes(3,3,i,j,k) = 2.0 * dz

                  ! Diagonal terms of the Hessian
                  analytic_hes(1,1,i,j,k) = 2.0
                  analytic_hes(2,2,i,j,k) = 2.0
                  analytic_hes(3,3,i,j,k) = 2.0

                  ! Off-diagonal terms of the Hessian
                  analytic_hes(1,2,i,j,k) = 1.0
                  analytic_hes(2,1,i,j,k) = 1.0
                  analytic_hes(1,3,i,j,k) = 1.0
                  analytic_hes(3,1,i,j,k) = 1.0
                  analytic_hes(2,3,i,j,k) = 1.0
                  analytic_hes(3,2,i,j,k) = 1.0
               end do
            end do
         end do

         print *, "  - Calculating Hessian..."
         hes = hessian(f, dx, dy, dz)

         test_passed = .true.

         do k = 1, nz
            do j = 1, ny
               do i = 1, nx

                  ! Check the Hessian matrix
                  do d1 = 1, 3
                     do d2 = 1, 3
                        if (abs(hes(d1,d2,i,j,k) - analytic_hes(d1,d2,i,j,k)) > tol) then
                           print *, "  - Mismatch at (", i, j, k, "): (",d1, d2, "):", hes(d1,d2,i,j,k), analytic_hes(d1,d2,i,j,k)
                           test_passed = .false.
                           exit
                        end if
                     end do
                  end do
                  if (.not. test_passed) exit
                  ! Exit if the test failed
               end do
               if (.not. test_passed) exit
            end do
            if (.not. test_passed) exit
         end do


         if (test_passed) then
            print *, "  - Test passed: Hessian calculated correctly."
            passed_tests = passed_tests + 1
         else
            print *, "  - Test failed: Incorrect Hessian calculation."
         end if

         print *
         ! Deallocate arrays
         deallocate(f)
         deallocate(analytic_hes)
         deallocate(hes)
      endif
   end subroutine test_scalar_hessian
   
   subroutine test_vector_hessian(total_tests, passed_tests)
      integer, intent(inout) :: total_tests, passed_tests
      real(real64), allocatable :: f(:,:,:,:), hes(:,:,:,:,:,:), analytic_hes(:,:,:,:,:,:)
      real(real64) :: dx, dy, dz, x, y, z, diff, sin_tolerance, other_tolerance, tolerance
      integer :: nx, ny, nz, i, j, k, l, d1, d2, neq
      logical :: test_passed

      if (my_rank == 0) then
         print *, "Test 6: Hessian Calculation (Vector Field)"

         total_tests = total_tests + 1

         neq = 3
         nx = 5
         ny = 5
         nz = 5
         dx = 1.0e-5_real64 
         dy = 2.0e-5_real64
         dz = 3.0e-5_real64

         ! Allocate for a 3-component vector function
         allocate(f(neq, nx, ny, nz))
         allocate(analytic_hes(neq, 3, 3, nx, ny, nz))

         ! Initialize the analytic Hessian for each component
         do l = 1, nz
            do k = 1, ny
               do j = 1, nx
                  x = real(j, real64) * dx
                  y = real(k, real64) * dy
                  z = real(l, real64) * dz

                  ! Component 1: f1 = x^2 + y^2 + z^2
                  ! Diagonal terms for f1
                  f(1,j,k,l) = x*x + y*y + z*z
                  analytic_hes(1, 1, 1, j, k, l) = 2.0
                  analytic_hes(1, 2, 2, j, k, l) = 2.0
                  analytic_hes(1, 3, 3, j, k, l) = 2.0
                  ! Mixed derivatives for f1
                  analytic_hes(1, 1, 2, j, k, l) = 0.0
                  analytic_hes(1, 1, 3, j, k, l) = 0.0
                  analytic_hes(1, 2, 3, j, k, l) = 0.0
                  analytic_hes(1, 2, 1, j, k, l) = 0.0
                  analytic_hes(1, 3, 1, j, k, l) = 0.0
                  analytic_hes(1, 3, 2, j, k, l) = 0.0

                  ! Component 2: f2 = x*y + y*z + z*x
                  f(2,j,k,l) = x*y + y*z + z*x
                  ! Diagonal terms for f2
                  analytic_hes(2, 1, 1, j, k, l) = 0.0
                  analytic_hes(2, 2, 2, j, k, l) = 0.0
                  analytic_hes(2, 3, 3, j, k, l) = 0.0
                  ! Mixed derivatives for f2
                  analytic_hes(2, 1, 2, j, k, l) = 1.0
                  analytic_hes(2, 1, 3, j, k, l) = 1.0
                  analytic_hes(2, 2, 3, j, k, l) = 1.0
                  analytic_hes(2, 2, 1, j, k, l) = 1.0
                  analytic_hes(2, 3, 1, j, k, l) = 1.0
                  analytic_hes(2, 3, 2, j, k, l) = 1.0

                  ! Component 3: f3 = sin(x) * cos(y) * z^2
                  f(3,j,k,l) = sin(x) * cos(y) * z**2
                  ! Diagonal terms for f4
                  analytic_hes(3, 1, 1, j, k, l) = cos(x) * cos(y) * z**2
                  analytic_hes(3, 2, 2, j, k, l) = -sin(x) * sin(y) * z**2
                  analytic_hes(3, 3, 3, j, k, l) = 2.0 * sin(x) * cos(y)
                  ! Mixed derivatives for f4
                  analytic_hes(3, 1, 2, j, k, l) = -cos(x) * sin(y) * z**2
                  analytic_hes(3, 1, 3, j, k, l) = 2.0 * cos(x) * cos(y) * z
                  analytic_hes(3, 2, 3, j, k, l) = -2.0 * sin(x) * sin(y) * z
                  analytic_hes(3, 2, 1, j, k, l) = -cos(x) * sin(y) * z**2
                  analytic_hes(3, 3, 1, j, k, l) = 2.0 * cos(x) * cos(y) * z
                  analytic_hes(3, 3, 2, j, k, l) = -2.0 * sin(x) * sin(y) * z
               end do
            end do
         end do
         print *, "  - Calculating Hessian for vector function..."
         hes = hessian(f, dx, dy, dz)
         test_passed = .true.
         
         sin_tolerance = 1.0e-6
         other_tolerance = 1.0e-12

         do l = 1, nz
            do k = 1, ny
               do j = 1, nx
                  ! Check the Hessian matrix for each component
                  do d1 = 1, 3
                     do d2 = 1, 3
                        do i = 1, neq
                           if (i == 3) then
                              tolerance = sin_tolerance
                           else
                              tolerance = other_tolerance
                           end if

                           if (abs(hes(i,d1,d2,j,k,l) - analytic_hes(i,d1,d2,j,k,l)) > tolerance) then
                              print "(A,I3, A, 3I3, A, 2I3, A, 2F10.6, 2F10.6)",          &
                                    "  - Mismatch at eq: ", i,                            &
                                    " idx: (", j, k, l, "): (dx_i*dx_j):(",d1, d2, "):",  &
                                    hes(i, d1, d2, j,k,l),                                &
                                    analytic_hes(i, d1, d2, j,k,l)
                              diff = abs(hes(i, d1, d2, j,k,l) - analytic_hes(i, d1, d2, j,k,l))
                              print "(A, E12.6)", "  - Difference: ", diff
                              print "(A, E12.6)", "  - Tolerance: ", tol
                              test_passed = .false.
                              exit
                           end if
                        end do
                     end do
                  end do

                  if (.not. test_passed) exit
                  ! Exit if the test failed
               end do
               if (.not. test_passed) exit
            end do
            if (.not. test_passed) exit
         end do

         if (test_passed) then
            print *, "  - Test passed: Vector Hessian calculated correctly."
            passed_tests = passed_tests + 1
         else
            print *, "  - Test failed: Incorrect Vector Hessian calculation."
         end if

         print *
         ! Deallocate arrays
         deallocate(f)
         deallocate(analytic_hes)
         deallocate(hes)
      endif
   end subroutine test_vector_hessian

   subroutine test_differential_operators(total_tests, passed_tests)
      integer, intent(inout) :: total_tests, passed_tests
      real(real64), allocatable :: f(:,:,:,:), result(:,:,:,:)
      real(real64) :: dx, dy, dz, x, y, z, tol
      integer :: neq, nx, ny, nz, j, k, l
      logical :: test_passed

      if (my_rank == 0) then
         print *, "Test 7: Differential Operators"

         total_tests = total_tests + 3  ! We're testing 3 operators (divergence, curl, Laplacian)

         neq = 3
         nx = 10
         ny = 10
         nz = 10
         dx = 0.1_real64
         dy = 0.1_real64
         dz = 0.1_real64
         tol = 1e-6_real64  ! Tolerance for checking the results

         allocate(f(neq, nx, ny, nz))

         ! Define a simple vector field f(x, y, z)
         ! Let's define f as:
         ! f(1) = x^2 + y^2
         ! f(2) = y^2 + z^2
         ! f(3) = z^2 + x^2
         do l = 1, nz
            do k = 1, ny
                  do j = 1, nx
                     x = real(j, real64) * dx
                     y = real(k, real64) * dy
                     z = real(l, real64) * dz

                     f(1,j,k,l) = x**2 + y**2
                     f(2,j,k,l) = y**2 + z**2
                     f(3,j,k,l) = z**2 + x**2
                  end do
            end do
         end do

         ! Test divergence
         print *, "  - Testing divergence..."
         result = apply_differential_operator(f, dx, dy, dz, "divergence")
         test_passed = .true.
         ! Analytical divergence: div(f) = 2x + 2y + 2z
         do l = 1, nz
            do k = 1, ny
                  do j = 1, nx
                     x = real(j, real64) * dx
                     y = real(k, real64) * dy
                     z = real(l, real64) * dz
                     if (abs(result(1,j,k,l) - (2*x + 2*y + 2*z)) > tol) then
                        print *, "    - Divergence test failed at (", j, k, l, "):", result(1,j,k,l), "expected:", 2*x + 2*y + 2*z
                        test_passed = .false.
                     end if
                  end do
            end do
         end do

         if (test_passed) then
            print *, "    - Divergence test passed."
            passed_tests = passed_tests + 1
         else
            print *, "    - Divergence test failed."
         end if

         ! Test curl
         print *, "  - Testing curl..."
         result = apply_differential_operator(f, dx, dy, dz, "curl")
         test_passed = .true.
         ! Analytical curl: curl(f) = [-2z, 2x, -2y]
         do l = 1, nz
            do k = 1, ny
                  do j = 1, nx
                     x = real(j, real64) * dx
                     y = real(k, real64) * dy
                     z = real(l, real64) * dz
                     if ( abs(result(1,j,k,l) + 2*z) > tol .or. &
                        abs(result(2,j,k,l) + 2*x) > tol .or. &
                        abs(result(3,j,k,l) + 2*y) > tol ) then
                        print *, "    - Curl test failed at (", j, k, l, ")"
                        print *, "      - Result: ", result(:,j,k,l)
                        print *, "      - Expected: ", [-2*z, -2*x, -2*y]
                        test_passed = .false.
                     end if
                  end do
            end do
         end do

         if (test_passed) then
            print *, "    - Curl test passed."
            passed_tests = passed_tests + 1
         else
            print *, "    - Curl test failed."
         end if

         ! Test Laplacian
         print *, "  - Testing Laplacian..."
         result = apply_differential_operator(f, dx, dy, dz, "laplacian")
         test_passed = .true.
         ! Analytical Laplacian: lap(f) = [2 + 2, 2 + 2, 2 + 2] = [4, 4, 4]
         do l = 1, nz
            do k = 1, ny
                  do j = 1, nx
                     if (abs(result(1,j,k,l) - 4.0) > tol .or. &
                        abs(result(2,j,k,l) - 4.0) > tol .or. &
                        abs(result(3,j,k,l) - 4.0) > tol) then
                        print *, "    - Laplacian test failed at (", j, k, l, "):", result(:,j,k,l), "expected: 4.0"
                        test_passed = .false.
                     end if
                  end do
            end do
         end do

         if (test_passed) then
            print *, "    - Laplacian test passed."
            passed_tests = passed_tests + 1
         else
            print *, "    - Laplacian test failed."
         end if
         print *

         ! Deallocate arrays
         deallocate(f)
         deallocate(result)
      endif
   end subroutine test_differential_operators

   subroutine test_distribute_and_collect(total_tests, passed_tests)
      integer, intent(inout) :: total_tests, passed_tests
      real(real64), allocatable :: global_data(:,:,:,:), local_data(:,:,:,:)
      real(real64), allocatable :: collected_data(:,:,:,:)
      integer :: dim, nx, ny, nz, neq, elements, i, j, k, local_nx, local_ny, local_nz,eq
      real(real64) :: st, et
      character(len=1024) :: filename

      ! Set dimensions
      if (my_rank == 0) then
         print *, "Test 0: Distribute and Collect Data"
      end if
      ! Repeat 5 time for different dimensions from 100 to 10.000.000 elements
      do dim = 3, 7
         !  Split the elements into 3 dimensions
         neq = 2
         elements = 10**dim
         nx = 10
         ny = 10
         nz = elements / (nx * ny)
         if (nz < 10) then
            nz = 10
         end if

         if (my_rank == 0) then 
            print "(A,I10,A,I3,A,I3,A,I3,A,I10, A)",                 &
               "  - Testing with ", neq *nx *ny *nz,                 &
               ", elements (", neq, "x",  nx, "x", ny, "x", nz, ")"
            st = MPI_Wtime()
         end if

         ! Allocate and initialize global data on root process
         if (my_rank == 0) then
            allocate(global_data(neq,nx, ny, nz))
            global_data = reshape([(real(i, real64), i = 1, neq *nx * ny * nz)], [neq, nx, ny, nz])
         end if

         ! Distribute data
         call commit_array(global_data, [1, 1 , 2 , 3], 1)
         call distribute_data(global_data, local_data)
         ! Gather data
         call collect_data(collected_data, local_data)
         call free_vars()
         
         ! Validate that collected data matches original global data
         if (my_rank == 0) then
            total_tests = total_tests + 1
            et = MPI_Wtime()
            print "(A,F9.5,A)", achar(9)//"- Time taken: ", et - st, " s"
            if (all(collected_data == global_data)) then
                  passed_tests = passed_tests + 1
                  print "(A)", achar(9)//"- Test passed: Data distributed and collected correctly."
            else 
                  print "(A)", achar(9)//"- Test failed: Data distributed and collected incorrectly."
            end if
         end if
         ! Cleanup
         deallocate(local_data)
         if (my_rank == 0) then
            deallocate(global_data)
            deallocate(collected_data)
         end if
      enddo
      if (my_rank == 0) print *
   end subroutine test_distribute_and_collect
end module test_operators

program run_tests
   use test_operators
   use MPI
   implicit none
   integer :: ierr

   call MPI_INIT(ierr)
   call run_all_tests()
   call MPI_FINALIZE(ierr)
end program run_tests
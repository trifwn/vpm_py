module vpm_gcalc
    use vpm_types, only: dp
    use serial_vector_field_operators
    implicit none
contains
    subroutine calc_velocity
        use MPI
        use pmgrid, only: velocity_pm, SOL_pm
        use vpm_vars, only: neqpm
        use console_io, only: vpm_print, blue, yellow, dummy_string
        use vpm_types, only: dp
        use vpm_size, only: fine_grid
        implicit none
        real(dp) ::  st, et
        real(dp) ::  DXpm, DYpm, DZpm

        DXpm = fine_grid%Dpm(1)
        DYpm = fine_grid%Dpm(2)
        DZpm = fine_grid%Dpm(3)

        st = MPI_WTIME()
        write (dummy_string, "(A)") "Calculating Velocities on PM using FD"
        call vpm_print(dummy_string, blue, 1)

        ! Velocity due to vorticity
        velocity_pm = curl(SOL_pm(1:3, :, :, :), DXpm, DYpm, DZpm)
        ! Velocity due to irrotational part
        if (neqpm .eq. 4) then
            velocity_pm = velocity_pm + gradient(SOL_pm(4, :, :, :), DXpm, DYpm, DZpm)
        end if

        et = MPI_WTIME()
        write (dummy_string, "(A,I5,A,F8.2,A)") &
            achar(9)//'finished in:', int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
        call vpm_print(dummy_string, yellow, 1)
    end subroutine calc_velocity

    subroutine calc_vortex_stretching_conservative(velocity, deformation)
        use MPI
        use pmgrid, only: deform_pm, RHS_pm
        use vpm_vars, only: neqpm
        use vpm_size, only: fine_grid
        use console_io, only: vpm_print, blue, yellow, dummy_string
        implicit none
        real(dp) ::  DXpm2, DYpm2, DZpm2
        real(dp) ::  st, et
        real(dp), intent(in) :: velocity(:, :, :, :)
        real(dp), allocatable, intent(out):: deformation(:, :, :, :)
        real(dp) ::  upi, upj, upk, vpi, vpj, vpk, wpi, wpj, wpk, &
                    umi, umj, umk, vmi, vmj, vmk, wmi, wmj, wmk, &
                    wdudx, wdvdy, wdwdz
        integer  :: i, j, k
        integer  :: NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, NXf_fine_bl, NYf_fine_bl, NZf_fine_bl

        st = MPI_WTIME()
        write (dummy_string, "(A)") 'Calculating Deformations on PM using FD'
        call vpm_print(dummy_string, blue, 1)

        if (allocated(deformation)) deallocate (deformation)
        allocate (deformation, mold=velocity)

        DXpm2 = 2*fine_grid%Dpm(1)
        DYpm2 = 2*fine_grid%Dpm(2)
        DZpm2 = 2*fine_grid%Dpm(3)

        NXs_fine_bl = fine_grid%NN_bl(1)
        NYs_fine_bl = fine_grid%NN_bl(2)
        NZs_fine_bl = fine_grid%NN_bl(3)
        NXf_fine_bl = fine_grid%NN_bl(4)
        NYf_fine_bl = fine_grid%NN_bl(5)
        NZf_fine_bl = fine_grid%NN_bl(6)

        deform_pm = 0.d0
        !$omp parallel private(upi,upj,upk,vpi,vpj,vpk,wpi,wpj,wpk,umi,umj,umk,vmi,vmj,vmk,&
        !$omp                  wmi,wmj,wmk,          &
        !$omp                  i,j,k,                &
        !$omp                  wdudx,wdvdy,wdwdz )   &
        !$omp shared(velocity,RHS_pm,SOL_pm,deform_pm, deformation, total_error, max_error)
        !$omp num_threads(OMPTHREADS)
        !$omp do
        do k = NZs_fine_bl + 2, NZf_fine_bl - 2
            do j = NYs_fine_bl + 2, NYf_fine_bl - 2
                do i = NXs_fine_bl + 2, NXf_fine_bl - 2
                    !REMEMBER VORTICITY CARRIED IS -OMEGA and the quantity transfered is -OMEGA thus
                    ! deformation = - (omega*\nabla) u =  \div(w (x) u) - u * (\grad w) =
                    ! We will calculate the deformation by using \div(w (x) u) as it is a conservative form
                    upi = velocity(1, i + 1, j, k)
                    upj = velocity(1, i, j + 1, k)
                    upk = velocity(1, i, j, k + 1)

                    vpi = velocity(2, i + 1, j, k)
                    vpj = velocity(2, i, j + 1, k)
                    vpk = velocity(2, i, j, k + 1)

                    wpi = velocity(3, i + 1, j, k)
                    wpj = velocity(3, i, j + 1, k)
                    wpk = velocity(3, i, j, k + 1)

                    umi = velocity(1, i - 1, j, k)
                    umj = velocity(1, i, j - 1, k)
                    umk = velocity(1, i, j, k - 1)

                    vmi = velocity(2, i - 1, j, k)
                    vmj = velocity(2, i, j - 1, k)
                    vmk = velocity(2, i, j, k - 1)

                    wmi = velocity(3, i - 1, j, k)
                    wmj = velocity(3, i, j - 1, k)
                    wmk = velocity(3, i, j, k - 1)
                    !DEFORMATION WITH A MINUS BECAUSE WE HAVE STORED MINUS VORTICITY
                    wdudx = -(RHS_pm(1, i + 1, j, k)*upi - (RHS_pm(1, i - 1, j, k))*umi)/DXpm2
                    wdvdy = -(RHS_pm(2, i, j + 1, k)*upj - (RHS_pm(2, i, j - 1, k))*umj)/DYpm2
                    wdwdz = -(RHS_pm(3, i, j, k + 1)*upk - (RHS_pm(3, i, j, k - 1))*umk)/DZpm2
                    deform_pm(1, i, j, k) = wdudx + wdvdy + wdwdz

                    wdudx = -(RHS_pm(1, i + 1, j, k)*vpi - (RHS_pm(1, i - 1, j, k))*vmi)/DXpm2
                    wdvdy = -(RHS_pm(2, i, j + 1, k)*vpj - (RHS_pm(2, i, j - 1, k))*vmj)/DYpm2
                    wdwdz = -(RHS_pm(3, i, j, k + 1)*vpk - (RHS_pm(3, i, j, k - 1))*vmk)/DZpm2
                    deform_pm(2, i, j, k) = wdudx + wdvdy + wdwdz

                    wdudx = -(RHS_pm(1, i + 1, j, k)*wpi - (RHS_pm(1, i - 1, j, k))*wmi)/DXpm2
                    wdvdy = -(RHS_pm(2, i, j + 1, k)*wpj - (RHS_pm(2, i, j - 1, k))*wmj)/DYpm2
                    wdwdz = -(RHS_pm(3, i, j, k + 1)*wpk - (RHS_pm(3, i, j, k - 1))*wmk)/DZpm2
                    deform_pm(3, i, j, k) = wdudx + wdvdy + wdwdz
                end do
            end do
        end do
        !$omp enddo
        !$omp endparallel

        et = MPI_WTIME()
        write (dummy_string, "(A,I5,A,F8.2,A)") &
            achar(9)//'finished in:', int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
        call vpm_print(dummy_string, yellow, 1)
    end subroutine calc_vortex_stretching_conservative

    subroutine calc_vortex_stretching_definition(velocity, deformation)
        use MPI
        use pmgrid, only: RHS_pm
        use vpm_size, only: fine_grid
        implicit none
        real(dp), intent(in) :: velocity(:, :, :, :)
        real(dp), allocatable, intent(out):: deformation(:, :, :, :)
        real(dp), allocatable :: jac_u(:, :, :, :, :)
        integer  :: NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, NXf_fine_bl, NYf_fine_bl, NZf_fine_bl
        integer  :: i, j, k

        ! We will calculate the deformation using the exact expression deformation = - (omega*\nabla) u

        NXs_fine_bl = fine_grid%NN_bl(1)
        NYs_fine_bl = fine_grid%NN_bl(2)
        NZs_fine_bl = fine_grid%NN_bl(3)
        NXf_fine_bl = fine_grid%NN_bl(4)
        NYf_fine_bl = fine_grid%NN_bl(5)
        NZf_fine_bl = fine_grid%NN_bl(6)

        jac_u = jacobian(velocity, fine_grid%Dpm(1), fine_grid%Dpm(2), fine_grid%Dpm(3))

        if (allocated(deformation)) deallocate (deformation)
        allocate (deformation(3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
        deformation = 0.d0
        !$omp parallel private(i,j,k )
        !$omp shared(velocity_pm,RHS_pm,SOL_pm, deformation)
        !$omp num_threads(OMPTHREADS)
        !$omp do
        do k = NZs_fine_bl + 2, NZf_fine_bl - 2
            do j = NYs_fine_bl + 2, NYf_fine_bl - 2
                do i = NXs_fine_bl + 2, NXf_fine_bl - 2
                    deformation(1, i, j, k) = dot_product(-RHS_pm(1:3, i, j, k), jac_u(1, :, i, j, k))
                    deformation(2, i, j, k) = dot_product(-RHS_pm(1:3, i, j, k), jac_u(2, :, i, j, k))
                    deformation(3, i, j, k) = dot_product(-RHS_pm(1:3, i, j, k), jac_u(3, :, i, j, k))
                end do
            end do
        end do
        !$omp enddo
        !$omp endparallel
        deallocate (jac_u)
    end subroutine calc_vortex_stretching_definition

    subroutine calc_vortex_stretching_conservative2(velocity, deformation)
        use MPI
        use pmgrid, only: RHS_pm
        use vpm_size, only: fine_grid
        use vpm_types, only: dp
        implicit none
        real(dp), intent(in) :: velocity(:, :, :, :)
        real(dp), allocatable, intent(out):: deformation(:, :, :, :)
        real(dp), allocatable :: jac_u(:, :, :, :, :)
        real(dp), allocatable :: div_w(:, :, :)
        integer  :: i, j, k
        integer  :: NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, NXf_fine_bl, NYf_fine_bl, NZf_fine_bl
        real(dp) :: DXpm, DYpm, DZpm

        NXs_fine_bl = fine_grid%NN_bl(1)
        NYs_fine_bl = fine_grid%NN_bl(2)
        NZs_fine_bl = fine_grid%NN_bl(3)
        NXf_fine_bl = fine_grid%NN_bl(4)
        NYf_fine_bl = fine_grid%NN_bl(5)
        NZf_fine_bl = fine_grid%NN_bl(6)

        DXpm = fine_grid%Dpm(1)
        DYpm = fine_grid%Dpm(2)
        DZpm = fine_grid%Dpm(3)

        if (allocated(deformation)) deallocate (deformation)
        allocate (deformation(3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3)))
        deformation = 0.d0

        ! We will calculate the deformation using the exact expression
        ! deformation = - (omega*\nabla) u = \div(w (x) u) - u * (\grad w) =
        ! We will calculate the deformation by using \div(w (x) u) as it is a conservative form
        ! \div(w (x) u) = \div(w) * u + w \div(u)
        jac_u = jacobian(velocity, DXpm, DYpm, DZpm)
        div_w = divergence(-RHS_PM(1:3, :, :, :), DXpm, DYpm, DZpm)

        !$omp parallel private(i,j,k )
        !$omp shared(velocity_pm,RHS_pm, deformation)
        !$omp num_threads(OMPTHREADS)
        !$omp do
        do k = NZs_fine_bl + 2, NZf_fine_bl - 2
            do j = NYs_fine_bl + 2, NYf_fine_bl - 2
                do i = NXs_fine_bl + 2, NXf_fine_bl - 2
                    deformation(1, i, j, k) = dot_product(-RHS_pm(1:3, i, j, k), jac_u(1, :, i, j, k))
                    deformation(2, i, j, k) = dot_product(-RHS_pm(1:3, i, j, k), jac_u(2, :, i, j, k))
                    deformation(3, i, j, k) = dot_product(-RHS_pm(1:3, i, j, k), jac_u(3, :, i, j, k))

                    deformation(1:3, i, j, k) = deformation(1:3, i, j, k) + velocity(1:3, i, j, k)*div_w(i, j, k)
                end do
            end do
        end do
        !$omp enddo
        !$omp endparallel
        deallocate (jac_u)
        deallocate (div_w)
    end subroutine calc_vortex_stretching_conservative2

    module subroutine diffuse_vort_3d(NI)
        use MPI
        use vpm_vars, only: neqpm
        use vpm_size, only: fine_grid
        use pmgrid, only: RHS_pm, deform_pm
        use serial_vector_field_operators, only: laplacian 
        implicit none
        real(dp), intent(in) :: NI
        real(dp) :: viscosity
        real(dp), allocatable :: laplace_vort(:, :, :, :)
        integer  :: i, j, k
        integer  :: NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, NXf_fine_bl, NYf_fine_bl, NZf_fine_bl
        real(dp) :: DXpm, DYpm, DZpm
        logical  :: variable_volume

        NXs_fine_bl = fine_grid%NN_bl(1)
        NYs_fine_bl = fine_grid%NN_bl(2)
        NZs_fine_bl = fine_grid%NN_bl(3)
        NXf_fine_bl = fine_grid%NN_bl(4)
        NYf_fine_bl = fine_grid%NN_bl(5)
        NZf_fine_bl = fine_grid%NN_bl(6)

        DXpm = fine_grid%Dpm(1)
        DYpm = fine_grid%Dpm(2)
        DZpm = fine_grid%Dpm(3)
        laplace_vort = laplacian( RHS_pm(1:3, :, :, :), DXpm, DYpm, DZpm) 

        variable_volume = .false.
        if (neqpm .eq. 3) then
            viscosity = NI
        elseif (neqpm .ge. 4) then
            variable_volume = .true.
        else
            print *, "Error: Invalid number of equations"
            stop
        endif
        deform_pm = 0.0d0
        !$omp parallel private(i,j,k)
        !$omp shared(laplace_vort, RHS_pm, viscosity, deform_pm)
        !$omp do
        do k = NZs_fine_bl + 1, NZf_fine_bl - 1
            do j = NYs_fine_bl + 1, NYf_fine_bl - 1
                do i = NXs_fine_bl + 1, NXf_fine_bl - 1
                    if (variable_volume) viscosity = RHS_pm(4, i, j, k) + NI
                    ! We use minus sign because RHS = -w
                    deform_pm(1, i, j, k) = - viscosity * laplace_vort(1, i, j, k) 
                    deform_pm(2, i, j, k) = - viscosity * laplace_vort(2, i, j, k) 
                    deform_pm(3, i, j, k) = - viscosity * laplace_vort(3, i, j, k) 
                end do
            end do
        end do
        !$omp enddo
        !$omp endparallel
        deallocate (laplace_vort)
    end subroutine diffuse_vort_3d

    module subroutine calc_antidiffusion
        use pmgrid, only: RHS_pm, deform_pm
        use vpm_size, only: fine_grid
        implicit none
        real(dp), allocatable   :: laplace_vort(:, :, :, :)
        real(dp)                :: dwxdx, dwydy, dwzdz, Ct
        real(dp)                :: DXpm_sq, DYpm_sq, DZpm_sq
        integer                 :: i, j, k
        integer  :: NXs_fine_bl, NYs_fine_bl, NZs_fine_bl, NXf_fine_bl, NYf_fine_bl, NZf_fine_bl
        real(dp) :: DXpm, DYpm, DZpm

        NXs_fine_bl = fine_grid%NN_bl(1)
        NYs_fine_bl = fine_grid%NN_bl(2)
        NZs_fine_bl = fine_grid%NN_bl(3)
        NXf_fine_bl = fine_grid%NN_bl(4)
        NYf_fine_bl = fine_grid%NN_bl(5)
        NZf_fine_bl = fine_grid%NN_bl(6)

        DXpm = fine_grid%Dpm(1)
        DYpm = fine_grid%Dpm(2)
        DZpm = fine_grid%Dpm(3)

        DXpm_sq = DXpm**2
        DYpm_sq = DYpm**2
        DZpm_sq = DZpm**2
        Ct = 6.8*DXpm**2/4

        laplace_vort = -laplacian(RHS_pm(1:3, :, :, :), DXpm, DYpm, DZpm)
        !$omp parallel private(i,j,k, dwxdx, dwydy, dwzdz) num_threads(OMPTHREADS)
        !$omp shared(laplace_vort, deform_pm)
        !$omp do
        do k = NZs_fine_bl + 1, NZf_fine_bl - 1
            do j = NYs_fine_bl + 1, NYf_fine_bl - 1
                do i = NXs_fine_bl + 1, NXf_fine_bl - 1
                    !Minus because of (-w) has been included in laplvort
                    dwxdx = (laplace_vort(1, i + 1, j, k) - 2*laplace_vort(1, i, j, k) &
                             + laplace_vort(1, i - 1, j, k))/DXpm_sq
                    dwydy = (laplace_vort(1, i, j + 1, k) - 2*laplace_vort(1, i, j, k) &
                             + laplace_vort(1, i, j - 1, k))/DYpm_sq
                    dwzdz = (laplace_vort(1, i, j, k + 1) - 2*laplace_vort(1, i, j, k) &
                             + laplace_vort(1, i, j, k - 1))/DZpm_sq
                    ! U = grad x psi
                    deform_pm(1, i, j, k) = deform_pm(1, i, j, k) + Ct*dwxdx + dwydy + dwzdz

                    dwxdx = (laplace_vort(2, i + 1, j, k) - 2*laplace_vort(2, i, j, k) &
                             + laplace_vort(2, i - 1, j, k))/DXpm_sq
                    dwydy = (laplace_vort(2, i, j + 1, k) - 2*laplace_vort(2, i, j, k) &
                             + laplace_vort(2, i, j - 1, k))/DYpm_sq
                    dwzdz = (laplace_vort(2, i, j, k + 1) - 2*laplace_vort(2, i, j, k) &
                             + laplace_vort(2, i, j, k - 1))/DZpm_sq
                    ! U = grad x psi
                    deform_pm(2, i, j, k) = deform_pm(2, i, j, k) + Ct*dwxdx + dwydy + dwzdz

                    dwxdx = (laplace_vort(3, i + 1, j, k) - 2*laplace_vort(3, i, j, k) &
                             + laplace_vort(3, i - 1, j, k))/DXpm_sq
                    dwydy = (laplace_vort(3, i, j + 1, k) - 2*laplace_vort(3, i, j, k) &
                             + laplace_vort(3, i, j - 1, k))/DYpm_sq
                    dwzdz = (laplace_vort(3, i, j, k + 1) - 2*laplace_vort(3, i, j, k) &
                             + RHS_pm(3, i, j, k - 1))/DZpm_sq
                    ! U = grad x psi
                    deform_pm(3, i, j, k) = deform_pm(3, i, j, k) + Ct*dwxdx + dwydy + dwzdz
                end do
            end do
        end do
        !$omp enddo
        !$omp endparallel
        deallocate (laplace_vort)
    end subroutine calc_antidiffusion

    subroutine compute_cross_product(vec1, vec2, cross_product)
        real(dp), intent(in)  :: vec1(:, :, :, :)  
        real(dp), intent(in)  :: vec2(:, :, :, :) 
        real(dp), intent(out) :: cross_product(:, :, :, :)
        integer :: i, j, k

        ! ! Validate input dimensions
        ! if (size(vec1,1) /= 3 .or. size(vec2,1) /= 3 .or. size(cross_product,1) /= 3) then
        !     print *, "Error: Input vectors must have 3 components"
        !     return
        ! end if
        ! if (size(vec1,2) /= size(vec2,2) .or. size(vec1,3) /= size(vec2,3) .or. size(vec1,4) /= size(vec2,4)) then
        !     print *, "Error: Input vectors must have same dimensions"
        !     return
        ! end if
        cross_product = 0.0d0
        !$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(cross_product, vec1, vec2)
        do i = 1, size(vec1,2)
            do j = 1, size(vec1,3)
                do k = 1, size(vec1,4)
                    cross_product(1, i, j, k) = vec1(2, i, j, k)*vec2(3, i, j, k) - &
                                                vec1(3, i, j, k)*vec2(2, i, j, k)
                    cross_product(2, i, j, k) = vec1(3, i, j, k)*vec2(1, i, j, k) - &
                                                vec1(1, i, j, k)*vec2(3, i, j, k)
                    cross_product(3, i, j, k) = vec1(1, i, j, k)*vec2(2, i, j, k) - &
                                                vec1(2, i, j, k)*vec2(1, i, j, k)
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine compute_cross_product

    !> \brief Computes the divergence of a vector field on a grid using the Gauss divergence theorem.
    !>        For each cell, computes the net flux through all faces by considering:
    !>        1. Right face minus left face flux in x-direction
    !>        2. Top face minus bottom face flux in y-direction
    !>        3. Front face minus back face flux in z-direction
    !>        The face fluxes are computed using averaging of adjacent cell centers.
    !>
    !> \param[in] grid The grid on which the divergence is computed
    !> \param[in] Q The vector field (3 components) for which divergence is computed
    !> \param[out] divQ The resulting divergence at cell centers
    subroutine compute_gauss_divergence(grid, Q_nodes, divQ_center)
        use vpm_types, only: cartesian_grid
        implicit none

        ! Input arguments
        type(cartesian_grid), intent(in)  :: grid
        real(dp), intent(in)    :: Q_nodes(3, grid%NN(1), grid%NN(2), grid%NN(3))
        real(dp), intent(out)   :: divQ_center(1, grid%NN(1), grid%NN(2), grid%NN(3))

        ! Local variables for grid spacing and fluxes
        real(dp) :: DXpm, DYpm, DZpm, DVpm
        real(dp) :: flux_x_right, flux_x_left
        real(dp) :: flux_y_top, flux_y_bottom
        real(dp) :: flux_z_front, flux_z_back
        integer :: i, j, k, nx, ny, nz

        ! Initialize divergence field
        divQ_center = 0.0d0

        ! Get grid spacings and cell volume
        DXpm = grid%Dpm(1)
        DYpm = grid%Dpm(2)
        DZpm = grid%Dpm(3)
        DVpm = DXpm*DYpm*DZpm

        ! Get grid dimensions
        nx = grid%NN(1)
        ny = grid%NN(2)
        nz = grid%NN(3)

        ! Compute divergence using the Divergence Theorem
        !$omp parallel do default(none) &
        !$omp private(i, j, k, flux_x_right, flux_x_left, flux_y_top, flux_y_bottom, flux_z_front, flux_z_back) &
        !$omp shared(divQ_center, Q_nodes, DXpm, DYpm, DZpm, DVpm, nx, ny, nz)
        do k = 1, nz - 1
            do j = 1, ny - 1
                do i = 1, nx - 1
                    ! X-direction fluxes (normal component through x-faces)
                    flux_x_left = (Q_nodes(1, i, j, k) + Q_nodes(1, i, j, k + 1) + &
                                   Q_nodes(1, i, j + 1, k) + Q_nodes(1, i, j + 1, k + 1))/4.0d0*DYpm*DZpm

                    flux_x_right = (Q_nodes(1, i + 1, j, k) + Q_nodes(1, i + 1, j, k + 1) + &
                                    Q_nodes(1, i + 1, j + 1, k) + Q_nodes(1, i + 1, j + 1, k + 1))/4.0d0*DYpm*DZpm

                    ! Y-direction fluxes (normal component through y-faces)
                    flux_y_bottom = (Q_nodes(2, i, j, k) + Q_nodes(2, i, j, k + 1) + &
                                     Q_nodes(2, i + 1, j, k) + Q_nodes(2, i + 1, j, k + 1))/4.0d0*DXpm*DZpm

                    flux_y_top = (Q_nodes(2, i, j + 1, k) + Q_nodes(2, i, j + 1, k + 1) + &
                                  Q_nodes(2, i + 1, j + 1, k) + Q_nodes(2, i + 1, j + 1, k + 1))/4.0d0*DXpm*DZpm

                    ! Z-direction fluxes (normal component through z-faces)
                    flux_z_back = (Q_nodes(3, i, j, k) + Q_nodes(3, i, j + 1, k) + &
                                   Q_nodes(3, i + 1, j, k) + Q_nodes(3, i + 1, j + 1, k))/4.0d0*DXpm*DYpm

                    flux_z_front = (Q_nodes(3, i, j, k + 1) + Q_nodes(3, i, j + 1, k + 1) + &
                                    Q_nodes(3, i + 1, j, k + 1) + Q_nodes(3, i + 1, j + 1, k + 1))/4.0d0*DXpm*DYpm

                    ! Compute divergence as net flux through all faces divided by cell volume
                    divQ_center(1, i, j, k) = ((flux_x_right - flux_x_left) + &
                                               (flux_y_top - flux_y_bottom) + &
                                               (flux_z_front - flux_z_back))/DVpm
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine compute_gauss_divergence

    !> \brief Computes the gradient of a scalar field at grid nodes using a flux-based averaging approach.
    !>     - Computes fluxes at cell faces by averaging scalar field values at adjacent nodes.
    !>     - Uses flux differences to calculate the gradient at grid nodes.
    !
    !> \param[in] grid Grid information.
    !> \param[in] Q_center Cell-centered scalar field values.
    !> \param[out] gradQ_nodes Gradient at nodes in x, y, z directions.
    subroutine compute_gradient_at_nodes(grid, Q_center, gradQ_nodes)
        use vpm_types, only: cartesian_grid
        implicit none

        ! Input arguments
        type(cartesian_grid), intent(in)  :: grid
        real(dp), intent(in)    :: Q_center(:, :, :)                                           ! Cell-centered scalar values
        real(dp), intent(out)   :: gradQ_nodes(3, grid%NN(1), grid%NN(2), grid%NN(3))  ! Gradient at nodes

        ! Local variables
        integer :: i, j, k
        integer :: nx, ny, nz
        real(dp) :: DXpm, DYpm, DZpm

        real(dp) :: flux_x_right, flux_x_left ! Fluxes through x-faces in the x-direction (i)
        real(dp) :: flux_y_top, flux_y_bottom ! Fluxes through y-faces in the y-direction (j)
        real(dp) :: flux_z_front, flux_z_back ! Fluxes through z-faces in the z-direction (k)

        ! Grid dimensions and spacings
        DXpm = grid%Dpm(1)
        DYpm = grid%Dpm(2)
        DZpm = grid%Dpm(3)

        nx = grid%NN(1)
        ny = grid%NN(2)
        nz = grid%NN(3)

        gradQ_nodes = 0.0d0

        ! Compute divergence using the Divergence Theorem
        !$omp parallel do default(none) &
        !$omp private(i, j, k, flux_x_right, flux_x_left, flux_y_top, flux_y_bottom, flux_z_front, flux_z_back) &
        !$omp shared(gradQ_nodes, Q_center, DXpm, DYpm, DZpm, nx, ny, nz)
        do k = 2, nz
            do j = 2, ny
                do i = 2, nx
                    flux_x_left = (Q_center(i - 1, j - 1, k - 1) + Q_center(i - 1, j, k - 1) + &
                                   Q_center(i - 1, j - 1, k) + Q_center(i - 1, j, k) &
                                   )/4.0d0

                    flux_x_right = (Q_center(i, j - 1, k - 1) + Q_center(i, j, k - 1) + &
                                    Q_center(i, j - 1, k) + Q_center(i, j, k) &
                                    )/4.0d0

                    flux_y_bottom = (Q_center(i - 1, j - 1, k - 1) + Q_center(i, j - 1, k - 1) + &
                                     Q_center(i - 1, j - 1, k) + Q_center(i, j - 1, k) &
                                     )/4.0d0

                    flux_y_top = (Q_center(i - 1, j, k - 1) + Q_center(i, j, k - 1) + &
                                  Q_center(i - 1, j, k) + Q_center(i, j, k) &
                                  )/4.0d0

                    flux_z_back = (Q_center(i - 1, j - 1, k - 1) + Q_center(i, j - 1, k - 1) + &
                                   Q_center(i - 1, j, k - 1) + Q_center(i, j, k - 1) &
                                   )/4.0d0

                    flux_z_front = (Q_center(i - 1, j - 1, k) + Q_center(i, j - 1, k) + &
                                    Q_center(i - 1, j, k) + Q_center(i, j, k) &
                                    )/4.0d0

                    ! Compute divergence as net flux through all faces divided by cell volume
                    gradQ_nodes(1, i, j, k) = (flux_x_right - flux_x_left)/DXpm
                    gradQ_nodes(2, i, j, k) = (flux_y_top - flux_y_bottom)/DYpm
                    gradQ_nodes(3, i, j, k) = (flux_z_front - flux_z_back)/DZpm
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine compute_gradient_at_nodes

    !> \brief For the operator nabla _{3D}^{2}=\partial _{x}^{2}+\partial _{y}^{2}+\partial _{z}^{2}
    !>        the fundemental solution is the Green's function G(x,y,z;x',y',z')=1/(4\pi r) where r
    !>        is the distance between the two points (x,y,z) and (x',y',z'). This subroutine computes
    !>        the Green's function for the 3D Laplacian operator.
    subroutine laplace_3D_green(x, y, z, x_prime, y_prime, z_prime, G)
        use constants, only: pi
        implicit none
        real(dp), intent(in)  :: x, y, z, x_prime, y_prime, z_prime
        real(dp), intent(out) :: G
        real(dp) :: r
        r = sqrt((x - x_prime)**2 + (y - y_prime)**2 + (z - z_prime)**2)
        if (r .lt. 1.0d-14) then
            G = 0.0d0
        else
            G = 1.0d0/(4.0d0*pi*r)
        end if
    end subroutine laplace_3D_green

    !> \brief Solves the poisson equation using the Green's represantation theorem.
    !> NOTE: THE BOUNDARY INTEGRALS ARE NOT IMPLEMENTED YET
    !> \param[in] grid The computational grid (cartesian)
    !> \param[in] f The source term
    !> \param[out] solution The solution to the poisson equation
    subroutine solve_3D_poisson_Biot(grid, f, solution)
        use vpm_types, only: cartesian_grid
        implicit none

        ! Input parameters
        type(cartesian_grid), intent(in) :: grid
        real(dp), dimension(:, :, :), intent(in) :: f
        real(dp), dimension(:, :, :), intent(out) :: solution

        ! Local variables
        integer :: i, j, k, i_prime, j_prime, k_prime
        real(dp) :: x, y, z, x_prime, y_prime, z_prime
        real(dp) :: G
        real(dp) :: DXpm, DYpm, DZpm

        !  r = |x-y|. For the poisson equation \nabla^{2} u = f  the solution is given by:
        !
        !  u(x) = \int_{\Gamma} g(x}, y)\frac{\partial u}{\partial n(y)}    dS(y)
        !       - \int_{\Gamma} \frac{\partial g(x,y)}{\partial n(y)} u(y)  dS(y)
        !       + \int_{\Omega} g(x, y)f(y)                                 dV(y)
        ! where:
        ! - \Gamma is the boundary of the domain
        ! - \Omega is the domain
        ! - n(y) is the normal vector at y
        ! - dS(y) is the surface element at y
        ! - dV(y) is the volume element at y
        ! - g(x,y) is the Green's function for the 3D Laplacian operator
        ! The two boundary integrals specify the Newmann and Dirichlet boundary conditions and the
        ! volume integral specifies the source term f(y).
        ! The Green's function for the 3D Laplacian operator is given by G(x;y)=1/(4\pi r) where
        ! r = |x-y| is the distance between the two points x and y (see laplace_3D_green subroutine).

        DXpm = grid%Dpm(1)
        DYpm = grid%Dpm(2)
        DZpm = grid%Dpm(3)

        !$OMP PARALLEL DO COLLAPSE(3) &
        !$OMP PRIVATE(x, y, z, i_prime, j_prime, k_prime, x_prime, y_prime, z_prime, G) &
        !$OMP SHARED(grid, f, solution, DXpm, DYpm, DZpm)
        do k = 1, grid%NN(3)
            do j = 1, grid%NN(2)
                do i = 1, grid%NN(1)
                    x = i*DXpm
                    y = j*DYpm
                    z = k*DZpm
                    solution(i, j, k) = 0.0d0

                    ! Volume integral computation
                    do k_prime = 1, grid%NN(3)
                        do j_prime = 1, grid%NN(2)
                            do i_prime = 1, grid%NN(1)
                                x_prime = i_prime*DXpm
                                y_prime = j_prime*DYpm
                                z_prime = k_prime*DZpm

                                call laplace_3D_green(x, y, z, x_prime, y_prime, z_prime, G)
                                solution(i, j, k) = solution(i, j, k) + G*f(i_prime, j_prime, k_prime)*DXpm*DYpm*DZpm
                            end do
                        end do
                    end do

                end do
            end do
        end do
        !$OMP END PARALLEL DO

        ! Note: Boundary integral computation is not implemented yet

    end subroutine solve_3D_poisson_Biot
end module vpm_gcalc

module vpm_mpi
    use vpm_types, only: dp
    implicit none
    ! Create interfaces for all module procedures
contains
    !> @brief Broadcasts the RHS_pm array to all processes.
    !>
    !> This subroutine is responsible for broadcasting the RHS_pm array, which
    !> contains the right-hand side values.
    !>
    !> @param RHS_pm The array containing the right-hand side values to be broadcasted.
    !> @param NN The grid size of the RHS_pm array.
    !> @param neq The number of equations.
    subroutine rhsbcast(RHS_pm, NN, neq)
        use MPI
        use mpi_matrices, only: mpimat4
        implicit none

        integer, intent(in)        :: NN(3), neq
        real(dp), intent(inout)    :: RHS_pm(neq, NN(1), NN(2), NN(3))
        integer                    :: my_rank, np, ierr, mat4

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        !---------------------------------------------
        call mpimat4(mat4, neq, NN(1), NN(2), NN(3))
        call MPI_BCAST(RHS_pm, 1, mat4, 0, MPI_COMM_WORLD, ierr)
        call MPI_TYPE_FREE(mat4, ierr)
        !-----------------------------------
    end subroutine rhsbcast

    subroutine rhsscat(fine_grid, RHS_pm, block_grid, RHS_pm_bl, nb_i, nb_j, nb_k)
        use vpm_vars, only: neqpm
        use vpm_types, only: cartesian_grid
        use MPI
        implicit none
        type(cartesian_grid), intent(in)  :: fine_grid, block_grid
        real(dp), intent(out)   :: RHS_pm_bl(neqpm, block_grid%NN(1), block_grid%NN(2), block_grid%NN(3))
        real(dp), intent(in)    :: RHS_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3))
        integer, intent(in)     :: nb_i, nb_j, nb_k

        integer :: my_rank, ierr, i, np
        integer :: ixs, jxs, kxs, ixf, jxf, kxf, nb, NXs, NXf, NYs, NYf, NZs, NZf

        RHS_pm_bl = 0
        NXs = block_grid%NN_bl(1)
        NYs = block_grid%NN_bl(2)
        NZs = block_grid%NN_bl(3)
        NXf = block_grid%NN_bl(4)
        NYf = block_grid%NN_bl(5)
        NZf = block_grid%NN_bl(6)

        ixs = (nb_i - 1)*(NXf - NXs) + fine_grid%NN_bl(1)
        jxs = (nb_j - 1)*(NYf - NYs) + fine_grid%NN_bl(2)
        kxs = (nb_k - 1)*(NZf - NZs) + fine_grid%NN_bl(3)

        ! ixf= NN_bl(4) - (NBI - nb_i)*(NXf-NXs)
        ! jxf= NN_bl(5) - (NBJ - nb_j)*(NYf-NYs)
        ! kxf= NN_bl(6) - (NBK - nb_k)*(NZf-NZs)
        ixf = fine_grid%NN_bl(1) + (nb_i)*(NXf - NXs)
        jxf = fine_grid%NN_bl(2) + (nb_j)*(NYf - NYs)
        kxf = fine_grid%NN_bl(3) + (nb_k)*(NZf - NZs)

        RHS_pm_bl(1:neqpm, NXs:NXf, NYs:NYf, NZs:NZf) = RHS_pm(1:neqpm, ixs:ixf, jxs:jxf, kxs:kxf)
        if (nb_i .gt. 1) RHS_pm_bl(:, NXs, :, :) = 0.d0
        if (nb_j .gt. 1) RHS_pm_bl(:, :, NYs, :) = 0.d0
        if (nb_k .gt. 1) RHS_pm_bl(:, :, :, NZs) = 0.d0

    end subroutine rhsscat

    subroutine solget(BLOCKS, NBI, NBJ, NBK, my_block, all_blocks, fine_grid, SOL_pm, SOL_pm_bl)
        use vpm_vars, only: neqpm
        use vpm_types, only: cartesian_grid
        use mpi_matrices, only: mpimat4
        use console_io, only: vpm_print, blue
        use MPI

        implicit none
        integer, intent(in)     :: BLOCKS, NBI, NBJ, NBK
        type(cartesian_grid), intent(in)  :: my_block, all_blocks(BLOCKS), fine_grid
        real(dp), intent(in)    :: SOL_pm_bl(neqpm, my_block%NN(1), my_block%NN(2), my_block%NN(3))
        real(dp), intent(out)   :: SOL_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3))

        ! Local variables
        real(dp), allocatable   :: SOL_pm_tmp(:, :, :, :)
        integer                 :: my_rank, ierr, source, dest, status(MPI_STATUS_SIZE), mat4
        integer                 :: ixs, jxs, kxs, ixf, jxf, kxf, nb, NXs, NXf, NYs, NYf, NZs, NZf, j, k, NN_block(3), i, nbs

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)

        if (my_rank .eq. 0) call vpm_print('Assembling Solution', blue, 1)
        nb = my_rank + 1

        if (my_rank .eq. 0) then
            j = 1
            i = 1
            k = 1

            NXs = my_block%NN_bl(1)
            NYs = my_block%NN_bl(2)
            NZs = my_block%NN_bl(3)
            NXf = my_block%NN_bl(4)
            NYf = my_block%NN_bl(5)
            NZf = my_block%NN_bl(6)

            ixs = (i - 1)*(NXf - NXs) + fine_grid%NN_bl(1)
            jxs = (j - 1)*(NYf - NYs) + fine_grid%NN_bl(2)
            kxs = (k - 1)*(NZf - NZs) + fine_grid%NN_bl(3)

            ixf = ixs + (NXf - NXs + 1) - 1
            jxf = jxs + (NYf - NYs + 1) - 1
            kxf = kxs + (NZf - NZs + 1) - 1

            SOL_pm(1:neqpm, ixs:ixf, jxs:jxf, kxs:kxf) = SOL_pm_bl(1:neqpm, NXs:NXf, NYs:NYf, NZs:NZf)
            !-->Assign
            do k = 1, NBK
                do j = 1, NBJ
                    do i = 1, NBI
                        nbs = (k - 1)*NBI*NBJ + (j - 1)*NBI + i
                        if (nbs .eq. 1) cycle

                        NN_block(1:3) = all_blocks(nbs)%NN(1:3)
                        NXs = all_blocks(nbs)%NN_bl(1)
                        NYs = all_blocks(nbs)%NN_bl(2)
                        NZs = all_blocks(nbs)%NN_bl(3)
                        NXf = all_blocks(nbs)%NN_bl(4)
                        NYf = all_blocks(nbs)%NN_bl(5)
                        NZf = all_blocks(nbs)%NN_bl(6)

                        ixs = (i - 1)*(NXf - NXs) + fine_grid%NN_bl(1)
                        jxs = (j - 1)*(NYf - NYs) + fine_grid%NN_bl(2)
                        kxs = (k - 1)*(NZf - NZs) + fine_grid%NN_bl(3)

                        ixs = (i - 1)*(NXf - NXs) + fine_grid%NN_bl(1)
                        jxs = (j - 1)*(NYf - NYs) + fine_grid%NN_bl(2)
                        kxs = (k - 1)*(NZf - NZs) + fine_grid%NN_bl(3)

                        ixf = ixs + (NXf - NXs + 1) - 1
                        jxf = jxs + (NYf - NYs + 1) - 1
                        kxf = kxs + (NZf - NZs + 1) - 1

                        ! Allocate and create Dtype
                        allocate (SOL_pm_tmp(neqpm, NN_block(1), NN_block(2), NN_block(3)))
                        call mpimat4(mat4, neqpm, NN_block(1), NN_block(2), NN_block(3))

                        source = nbs - 1
                        call MPI_RECV(SOL_pm_tmp, 1, mat4, source, 1, MPI_COMM_WORLD, status, ierr)
                        SOL_pm(1:neqpm, ixs:ixf, jxs:jxf, kxs:kxf) = SOL_pm_tmp(1:neqpm, NXs:NXf, NYs:NYf, NZs:NZf)
                        call MPI_TYPE_FREE(mat4, ierr)
                        deallocate (SOL_pm_tmp)
                    end do
                end do
            end do
        else
            dest = 0
            call mpimat4(mat4, neqpm, my_block%NN(1), my_block%NN(2), my_block%NN(3))
            call MPI_SEND(SOL_pm_bl, 1, mat4, dest, 1, MPI_COMM_WORLD, ierr)
            call MPI_TYPE_FREE(mat4, ierr)
        end if

    end subroutine solget

    subroutine velbcast
        use pmgrid, only: velocity_pm
        use vpm_size, only: fine_grid
        use mpi_matrices, only: mpimat4
        use MPI
        implicit none
        integer :: my_rank, np, ierr, mat4

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        !---------------------------------------------
        call mpimat4(mat4, 3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3))
        call MPI_BCAST(velocity_pm, 1, mat4, 0, MPI_COMM_WORLD, ierr)
        call MPI_TYPE_FREE(mat4, ierr)
        !--------------------------------------------
    end subroutine velbcast

    subroutine defbcast
        use pmgrid, only: deform_pm
        use vpm_size, only: fine_grid
        use mpi_matrices, only: mpimat4
        use MPI
        implicit none
        integer :: my_rank, np, ierr, mat4

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        !---------------------------------------------
        call mpimat4(mat4, 3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3))
        call MPI_BCAST(deform_pm, 1, mat4, 0, MPI_COMM_WORLD, ierr)
        call MPI_TYPE_FREE(mat4, ierr)
        !--------------------------------------------
    end subroutine defbcast

    subroutine particles_scat
        use vpm_vars, only: neqpm
        use parvar, only: XP, QP, NVR, NVR_size, XP_scatt, QP_scatt, NVR_p
        use mpi_matrices, only: mpimat2_pm
        use MPI
        use console_io, only: vpm_print, blue, tab_level, dummy_string, nocolor, yellow

        implicit none
        integer :: my_rank, np, ierr, i
        integer :: dest, NVR_pr, NVR_r, mat2
        integer :: status(MPI_STATUS_SIZE)
        real(dp):: st, et
        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        if (my_rank .eq. 0) then
            st = MPI_WTIME()
            write (dummy_string, "(A)") 'Scattering Particles to all Processes'
            call vpm_print(dummy_string, blue, 1)
        end if
        tab_level = tab_level + 1

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        !---------------------------------------------
        if (my_rank .eq. 0) then
            XP_scatt(1:3, 1:NVR_p) = XP(1:3, 1:NVR_p)
            QP_scatt(1:neqpm + 1, 1:NVR_p) = QP(1:neqpm + 1, 1:NVR_p)
            NVR_pr = NVR_p
            NVR_r = NVR/np

            write (dummy_string, "(A,I5,A,I5, A, I5, A, I5)") &
                "Processor ", my_rank, " got : NVR_p = ", NVR_p, " and NVR_r = ", NVR_r, "  NVR message size = ", NVR_size
            call vpm_print(dummy_string, nocolor, 2)
            if (NVR_r .gt. 0) then
                do i = 2, np
                    dest = i - 1
                    call mpimat2_pm(mat2, 3, NVR_size, 3, NVR_r, NVR_pr)!NVR_pr because counting starts from 0
                    call MPI_SEND(XP, 1, mat2, dest, 1, MPI_COMM_WORLD, ierr)
                    if (ierr .ne. 0) then
                        write (*, "(A,I5)") "Error in receiving XP_scatt on processor ", my_rank
                        write (*, "(A,I5)") "Got error", ierr
                    end if
                    call MPI_TYPE_FREE(mat2, ierr)
                    if (ierr .ne. 0) then
                        write (*, "(A,I5)") "Error in receiving XP_scatt on processor ", my_rank
                        write (*, "(A,I5)") "Got error", ierr
                    end if

                    call mpimat2_pm(mat2, neqpm + 1, NVR_size, neqpm + 1, NVR_r, NVR_pr)
                    call MPI_SEND(QP, 1, mat2, dest, 1, MPI_COMM_WORLD, ierr)
                    if (ierr .ne. 0) then
                        write (*, "(A,I5)") "Error in receiving XP_scatt on processor ", my_rank
                        write (*, "(A,I5)") "Got error", ierr
                    end if
                    call MPI_TYPE_FREE(mat2, ierr)
                    if (ierr .ne. 0) then
                        write (*, "(A,I5)") "Error in receiving XP_scatt on processor ", my_rank
                        write (*, "(A,I5)") "Got error", ierr
                    end if

                    NVR_pr = NVR_pr + NVR_r
                end do
            end if
        else
            if (NVR_p .gt. 0) then
                call mpimat2_pm(mat2, 3, NVR_p, 3, NVR_p, 0)
                call MPI_RECV(XP_scatt, 1, mat2, 0, 1, MPI_COMM_WORLD, status, ierr)
                if (ierr .ne. 0) then
                    write (*, "(A,I5)") "Error in receiving XP_scatt on processor ", my_rank
                    write (*, "(A,I5)") "Got error", ierr
                end if
                call MPI_TYPE_FREE(mat2, ierr)
                if (ierr .ne. 0) then
                    write (*, "(A,I5)") "Error in freeing XP_scatt on processor ", my_rank
                    write (*, "(A,I5)") "Got error", ierr
                end if

                call mpimat2_pm(mat2, neqpm + 1, NVR_p, neqpm + 1, NVR_p, 0)
                call MPI_RECV(QP_scatt, 1, mat2, 0, 1, MPI_COMM_WORLD, status, ierr)
                if (ierr .ne. 0) then
                    write (*, "(A,I5)") "Error in receiving QP_scatt on processor ", my_rank
                    write (*, "(A,I5)") "Got error", ierr
                end if
                call MPI_TYPE_FREE(mat2, ierr)
                if (ierr .ne. 0) then
                    write (*, "(A,I5)") "Error in freeing QP_scatt on processor ", my_rank
                    write (*, "(A,I5)") "Got error", ierr
                end if

                write (dummy_string, "(A,I5,A,I5)") "Processor ", my_rank, " got : NVR_p = ", NVR_p
                call vpm_print(dummy_string, nocolor, 2)
            else
                write (dummy_string, "(A,I5)") "Got 0 particles on processor ", my_rank
                call vpm_print(dummy_string, nocolor, 2)
            end if
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (my_rank .eq. 0) then
            write (dummy_string, "(A,I5)") "Total number of particles distributed = ", NVR_pr
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, "(A,I5)") "Total number of particles = ", NVR
            call vpm_print(dummy_string, nocolor, 2)
            et = MPI_WTIME()
            write (dummy_string, "(A,I5,A,F8.2,A)") &
                'finished in:', int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
            call vpm_print(dummy_string, yellow, 1)
        end if
        tab_level = tab_level - 1

    end subroutine particles_scat

    subroutine particles_gath
        use vpm_vars, only: neqpm
        use parvar, only: XP, QP, UP, GP, NVR, NVR_p, XP_scatt, QP_scatt, UP_scatt, GP_scatt
        use mpi_matrices, only: mpimat2_pm
        use MPI

        implicit none
        integer :: my_rank, np, ierr, i
        integer :: dest, NVR_pr, NVR_r, mat2
        integer :: status(MPI_STATUS_SIZE)
        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        !---------------------------------------------
        if (my_rank .eq. 0) then
            XP(1:3, 1:NVR_p) = XP_scatt(1:3, 1:NVR_p)
            QP(1:neqpm + 1, 1:NVR_p) = QP_scatt(1:neqpm + 1, 1:NVR_p)
            UP(1:3, 1:NVR_p) = UP_scatt(1:3, 1:NVR_p)
            GP(1:3, 1:NVR_p) = GP_scatt(1:3, 1:NVR_p)
            deallocate (QP_scatt, XP_scatt, UP_scatt, GP_scatt)
            NVR_pr = NVR_p
            NVR_r = NVR/np
            allocate (XP_scatt(3, NVR_r), QP_scatt(neqpm + 1, NVR_r), UP_scatt(3, NVR_r), GP_scatt(3, NVR_r))
            if (NVR_r .gt. 0) then
                do i = 2, np
                    dest = i - 1
                    call mpimat2_pm(mat2, 3, NVR_r, 3, NVR_r, 0)!NVR_pr because counting starts from 0
                    call MPI_RECV(XP_scatt, 1, mat2, dest, 1, MPI_COMM_WORLD, status, ierr)
                    call MPI_TYPE_FREE(mat2, ierr)

                    call mpimat2_pm(mat2, neqpm + 1, NVR_r, neqpm + 1, NVR_r, 0)
                    call MPI_RECV(QP_scatt, 1, mat2, dest, 1, MPI_COMM_WORLD, status, ierr)
                    call MPI_TYPE_FREE(mat2, ierr)

                    call mpimat2_pm(mat2, 3, NVR_r, 3, NVR_r, 0)!NVR_pr because counting starts from 0
                    call MPI_RECV(UP_scatt, 1, mat2, dest, 1, MPI_COMM_WORLD, status, ierr)
                    call MPI_TYPE_FREE(mat2, ierr)

                    call mpimat2_pm(mat2, 3, NVR_r, 3, NVR_r, 0)!NVR_pr because counting starts from 0
                    call MPI_RECV(GP_scatt, 1, mat2, dest, 1, MPI_COMM_WORLD, status, ierr)
                    call MPI_TYPE_FREE(mat2, ierr)

                    XP(1:3, NVR_pr + 1:NVR_pr + NVR_r) = XP_scatt(1:3, 1:NVR_r)
                    QP(1:neqpm + 1, NVR_pr + 1:NVR_pr + NVR_r) = QP_scatt(1:neqpm + 1, 1:NVR_r)
                    UP(1:3, NVR_pr + 1:NVR_pr + NVR_r) = UP_scatt(1:3, 1:NVR_r)
                    GP(1:3, NVR_pr + 1:NVR_pr + NVR_r) = GP_scatt(1:3, 1:NVR_r)

                    NVR_pr = NVR_pr + NVR_r
                end do
            end if
        else
            if (NVR_p .gt. 0) then
                call mpimat2_pm(mat2, 3, NVR_p, 3, NVR_p, 0)
                call MPI_SEND(XP_scatt, 1, mat2, 0, 1, MPI_COMM_WORLD, ierr)
                call MPI_TYPE_FREE(mat2, ierr)

                call mpimat2_pm(mat2, neqpm + 1, NVR_p, neqpm + 1, NVR_p, 0)
                call MPI_SEND(QP_scatt, 1, mat2, 0, 1, MPI_COMM_WORLD, ierr)
                call MPI_TYPE_FREE(mat2, ierr)

                call mpimat2_pm(mat2, 3, NVR_p, 3, NVR_p, 0)
                call MPI_SEND(UP_scatt, 1, mat2, 0, 1, MPI_COMM_WORLD, ierr)
                call MPI_TYPE_FREE(mat2, ierr)

                call mpimat2_pm(mat2, 3, NVR_p, 3, NVR_p, 0)
                call MPI_SEND(GP_scatt, 1, mat2, 0, 1, MPI_COMM_WORLD, ierr)
                call MPI_TYPE_FREE(mat2, ierr)
            end if
        end if
    end subroutine particles_gath

    subroutine proj_gath
        use vpm_vars, only: neqpm
        use vpm_size, only: fine_grid
        use pmgrid, only: RHS_pm
        ! use parvar
        use mpi_matrices, only: mpimat4
        use MPI

        implicit None
        integer                 :: my_rank, np, ierr, i
        integer                 :: dest, source, mat4
        integer                 :: status(MPI_STATUS_SIZE)
        real(dp), allocatable   :: RHS_pm_tmp(:, :, :, :)
        integer                 :: NN(3)

        NN = fine_grid%NN
        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        if (my_rank .eq. 0) then
            allocate (RHS_pm_tmp(neqpm, NN(1), NN(2), NN(3)))
            do i = 2, np
                source = i - 1
                call mpimat4(mat4, neqpm, NN(1), NN(2), NN(3))
                call MPI_RECV(RHS_pm_tmp, 1, mat4, source, 1, MPI_COMM_WORLD, status, ierr)
                call MPI_TYPE_FREE(mat4, ierr)
                RHS_pm = RHS_pm + RHS_pm_tmp
            end do
            deallocate (RHS_pm_tmp)
        else
            dest = 0
            call mpimat4(mat4, neqpm, NN(1), NN(2), NN(3))
            call MPI_SEND(RHS_pm, 1, mat4, dest, 1, MPI_COMM_WORLD, ierr)
            call MPI_TYPE_FREE(mat4, ierr)
        end if
    end subroutine proj_gath

    subroutine proj_gath_new
        use vpm_vars, only: interf_iproj, neqpm, neqpm, ND
        use vpm_size, only: fine_grid
        use pmgrid, only: RHS_pm
        use parvar, only: XP_scatt
        use mpi_matrices, only: mpimat4
        use MPI

        implicit none
        integer               :: my_rank, np, ierr, NN_proj(6), nn1, nn2, nn3, NN_tmp(6)
        integer               :: source, mat4
        integer               :: status(MPI_STATUS_SIZE)
        integer               :: imax, imin, jmax, jmin, kmax, kmin
        real(dp), allocatable :: RHS_pm_tmp(:, :, :, :)
        real(dp)              :: xpmax, xpmin, ypmax, ypmin, zpmax, zpmin
        integer               :: NN(3)
        real(dp)              :: DXpm, DYpm, DZpm, XMIN_pm, YMIN_pm, ZMIN_pm

        NN = fine_grid%NN

        DXpm = fine_grid%Dpm(1)
        DYpm = fine_grid%Dpm(2)
        DZpm = fine_grid%Dpm(3)

        XMIN_pm = fine_grid%Xbound(1)
        YMIN_pm = fine_grid%Xbound(2)
        ZMIN_pm = fine_grid%Xbound(3)

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
        if (ND .eq. 2) then
            xpmax = maxval(XP_scatt(1, :))
            xpmin = minval(XP_scatt(1, :))
            ypmax = maxval(XP_scatt(2, :))
            ypmin = minval(XP_scatt(2, :))

            imax = int((xpmax - XMIN_pm)/DXpm) + 1
            imin = int((xpmin - XMIN_pm)/DXpm) + 1
            jmax = int((ypmax - YMIN_pm)/DYpm) + 1
            jmin = int((ypmin - YMIN_pm)/DYpm) + 1

            NN_proj(1) = max(imin - interf_iproj, 1)
            NN_proj(2) = max(jmin - interf_iproj, 1)
            NN_proj(3) = 1
            NN_proj(4) = min(imax + interf_iproj, NN(1))
            NN_proj(5) = min(jmax + interf_iproj, NN(2))
            NN_proj(6) = 1
        else
            xpmax = maxval(XP_scatt(1, :))
            xpmin = minval(XP_scatt(1, :))
            ypmax = maxval(XP_scatt(2, :))
            ypmin = minval(XP_scatt(2, :))
            zpmax = maxval(XP_scatt(3, :))
            zpmin = minval(XP_scatt(3, :))

            imax = int((xpmax - XMIN_pm)/DXpm) + 1
            imin = int((xpmin - XMIN_pm)/DXpm) + 1
            jmax = int((ypmax - YMIN_pm)/DYpm) + 1
            jmin = int((ypmin - YMIN_pm)/DYpm) + 1
            kmax = int((zpmax - ZMIN_pm)/DZpm) + 1
            kmin = int((zpmin - ZMIN_pm)/DZpm) + 1

            NN_proj(1) = max(imin - interf_iproj, 1)
            NN_proj(2) = max(jmin - interf_iproj, 1)
            NN_proj(3) = max(kmin - interf_iproj, 1)
            NN_proj(4) = min(imax + interf_iproj, NN(1))
            NN_proj(5) = min(jmax + interf_iproj, NN(2))
            NN_proj(6) = min(kmax + interf_iproj, NN(3))
        end if

        if (my_rank .ne. 0) then
            call MPI_SEND(NN_proj, 6, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ierr)
            nn1 = NN_proj(4) - NN_proj(1) + 1
            nn2 = NN_proj(5) - NN_proj(2) + 1
            nn3 = NN_proj(6) - NN_proj(3) + 1
            allocate (RHS_pm_tmp(neqpm, nn1, nn2, nn3))
            RHS_pm_tmp(1:neqpm, 1:nn1, 1:nn2, 1:nn3) = RHS_pm(1:neqpm, NN_proj(1):NN_proj(4), &
                                                              NN_proj(2):NN_proj(5), &
                                                              NN_proj(3):NN_proj(6))
            call mpimat4(mat4, neqpm, nn1, nn2, nn3)
            call MPI_SEND(RHS_pm_tmp, 1, mat4, 0, 1, MPI_COMM_WORLD, ierr)
            call MPI_TYPE_FREE(mat4, ierr)
            deallocate (RHS_pm_tmp)
        else
            do source = 1, np - 1
                call MPI_RECV(NN_tmp, 6, MPI_INTEGER, source, 1, MPI_COMM_WORLD, status, ierr)
                nn1 = NN_tmp(4) - NN_tmp(1) + 1
                nn2 = NN_tmp(5) - NN_tmp(2) + 1
                nn3 = NN_tmp(6) - NN_tmp(3) + 1
                allocate (RHS_pm_tmp(neqpm, nn1, nn2, nn3))
                call mpimat4(mat4, neqpm, nn1, nn2, nn3)
                call MPI_RECV(RHS_pm_tmp, 1, mat4, source, 1, MPI_COMM_WORLD, status, ierr)
                RHS_pm(1:neqpm, NN_tmp(1):NN_tmp(4), NN_tmp(2):NN_tmp(5), NN_tmp(3):NN_tmp(6)) = &
                    RHS_pm(1:neqpm, NN_tmp(1):NN_tmp(4), &
                           NN_tmp(2):NN_tmp(5), &
                           NN_tmp(3):NN_tmp(6)) + &
                    RHS_pm_tmp(1:neqpm, 1:nn1, 1:nn2, 1:nn3)
                deallocate (RHS_pm_tmp)
                call MPI_TYPE_FREE(mat4, ierr)
            end do
        end if

    end subroutine proj_gath_new
end module vpm_mpi

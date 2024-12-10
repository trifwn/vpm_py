module mpi_matrices
    use MPI
    implicit none

    public :: mpimat2, mpimat2_pm, mpimat3_pm, mpimat4, mpimat4int, mpimat5, mpimat5d
contains
    subroutine mpimat2(mat2, nsize1, nsize2)
        use MPI

        integer ierr
        integer::typelist(2)
        integer ::imat(2), mat(2), start(2)
        ! integer ::istart
        integer :: nsize1, nsize2, mat2 !, orig1, orig2,
        !allocate(struct%AS_ij(nsize,nsize))

        mat(1) = nsize1
        mat(2) = nsize2
        imat = mat
        start(1) = 0
        start(2) = 0

        typelist(1) = MPI_DOUBLE_PRECISION
        typelist(2) = MPI_DOUBLE_PRECISION

        call MPI_TYPE_CREATE_SUBARRAY(2, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat2, ierr)
        call MPI_TYPE_COMMIT(mat2, ierr)
    end subroutine mpimat2

    subroutine mpimat2_pm(mat2, orig1, orig2, nsize1, nsize2, istart)
        use MPI

        implicit none
        integer ierr
        integer::typelist(2)
        integer ::imat(2), mat(2), start(2)
        integer ::istart
        integer ::orig1, orig2, nsize1, nsize2, mat2, my_rank, ieer
        !allocate(struct%AS_ij(nsize,nsize))
        imat(1) = orig1
        imat(2) = orig2

        mat(1) = nsize1
        mat(2) = nsize2

        start(1) = 0
        start(2) = istart

        typelist(1) = MPI_DOUBLE_PRECISION
        typelist(2) = MPI_DOUBLE_PRECISION

        call MPI_TYPE_CREATE_SUBARRAY(2, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat2, ierr)
        if (ierr /= 0) then
            print *, 'Error in MPI_TYPE_CREATE_SUBARRAY. Got error:', ierr
            call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ieer)
            print *, 'my_rank:', my_rank
        end if
        call MPI_TYPE_COMMIT(mat2, ierr)
        if (ierr /= 0) then
            print *, 'Error in MPI_TYPE_COMMIT. Got error:', ierr
            call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ieer)
            print *, 'my_rank:', my_rank
        end if
    end subroutine mpimat2_pm

    subroutine mpimat3_pm(mat3, nsize1, nsize2, nsize3)
        use MPI

        implicit none
        integer, intent(in)::nsize1, nsize2, nsize3
        integer ierr
        integer ::imat(3), mat(3), start(3)
        integer ::mat3
        !allocate(struct%AS_ij(nsize,nsize))
        imat(1) = nsize1
        imat(2) = nsize2
        imat(3) = nsize3
        mat(1) = nsize1
        mat(2) = nsize2
        mat(3) = nsize3
        start(1) = 0
        start(2) = 0
        start(3) = 0

        call MPI_TYPE_CREATE_SUBARRAY(3, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat3, ierr)
        call MPI_TYPE_COMMIT(mat3, ierr)
    end subroutine mpimat3_pm

    subroutine mpimat4(mat4, nsize1, nsize2, nsize3, nsize4)
        use MPI

        implicit none
        integer, intent(in)::nsize1, nsize2, nsize3, nsize4
        integer ierr
        integer ::imat(4), mat(4), start(4)
        integer ::mat4
        !allocate(struct%AS_ij(nsize,nsize))
        imat(1) = nsize1
        imat(2) = nsize2
        imat(3) = nsize3
        imat(4) = nsize4
        mat(1) = nsize1
        mat(2) = nsize2
        mat(3) = nsize3
        mat(4) = nsize4
        start(1) = 0
        start(2) = 0
        start(3) = 0
        start(4) = 0

        call MPI_TYPE_CREATE_SUBARRAY(4, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat4, ierr)
        call MPI_TYPE_COMMIT(mat4, ierr)
    end subroutine mpimat4

    subroutine mpimat4int(mat4int, nsize1, nsize2, nsize3, nsize4)
        use MPI

        implicit none
        integer, intent(in)::nsize1, nsize2, nsize3, nsize4
        integer ierr
        integer ::imat(4), mat(4), start(4)
        integer ::mat4int
        !allocate(struct%AS_ij(nsize,nsize))
        imat(1) = nsize1
        imat(2) = nsize2
        imat(3) = nsize3
        imat(4) = nsize4
        mat(1) = nsize1
        mat(2) = nsize2
        mat(3) = nsize3
        mat(4) = nsize4
        start(1) = 0
        start(2) = 0
        start(3) = 0
        start(4) = 0

        call MPI_TYPE_CREATE_SUBARRAY(4, imat, mat, start, MPI_ORDER_FORTRAN, MPI_INTEGER, mat4int, ierr)
        call MPI_TYPE_COMMIT(mat4int, ierr)
    end subroutine mpimat4int

    subroutine mpimat5_old(mat5, nsize1, nsize2, nsize3, nsize4, nsize5)
        use MPI

        implicit none
        integer, intent(in)::nsize1, nsize2, nsize3, nsize4, nsize5
        integer ierr
        integer ::imat(5), mat(5), start(5)
        integer ::mat5
        !allocate(struct%AS_ij(nsize,nsize))
        imat(1) = nsize1
        imat(2) = nsize2
        imat(3) = nsize3
        imat(4) = nsize4
        imat(5) = nsize5
        mat(1) = nsize1
        mat(2) = nsize2
        mat(3) = nsize3
        mat(4) = nsize4
        mat(5) = nsize5
        start(1) = 0
        start(2) = 0
        start(3) = 0
        start(4) = 0
        start(5) = 0

        call MPI_TYPE_CREATE_SUBARRAY(5, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat5, ierr)
        call MPI_TYPE_COMMIT(mat5, ierr)
    end subroutine mpimat5_old

    subroutine mpimat5(mat5, nsize1, nsize2, nsize3, nsize4, nsize5, isize5, nbstart)
        use MPI

        implicit none
        integer, intent(in)::nsize1, nsize2, nsize3, nsize4, nsize5, isize5, nbstart
        integer ierr
        integer ::imat(5), mat(5), start(5)
        integer ::mat5
        !allocate(struct%AS_ij(nsize,nsize))
        imat(1) = nsize1
        imat(2) = nsize2
        imat(3) = nsize3
        imat(4) = nsize4
        imat(5) = nsize5
        mat(1) = nsize1
        mat(2) = nsize2
        mat(3) = nsize3
        mat(4) = nsize4
        mat(5) = isize5
        start(1) = 0
        start(2) = 0
        start(3) = 0
        start(4) = 0
        start(5) = nbstart

        call MPI_TYPE_CREATE_SUBARRAY(5, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat5, ierr)
        call MPI_TYPE_COMMIT(mat5, ierr)
    end subroutine mpimat5

    subroutine mpimat5d(mat5, original_size, send_size, start_size)
        use MPI

        integer, intent(in)::original_size(5), send_size(5), start_size(5)
        integer ierr
        integer ::imat(5), mat(5), start(5)
        integer ::mat5
        mat5 = 0
        imat = original_size
        mat = send_size
        start = start_size

        call MPI_TYPE_CREATE_SUBARRAY(5, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat5, ierr)
        call MPI_TYPE_COMMIT(mat5, ierr)
    end subroutine mpimat5d

End module mpi_matrices

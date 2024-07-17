module mpi_matrices
   use MPI
   implicit none

   public :: mpimat2, mpimat2_pm, mpimat3_pm, mpimat4, mpimat4int, mpimat5, mpimat5d

contains
   Subroutine mpimat2(mat2, nsize1, nsize2)
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

      !write (*,*) nsize
      call MPI_TYPE_CREATE_SUBARRAY(2, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat2, ierr)
      call MPI_TYPE_COMMIT(mat2, ierr)
   end subroutine mpimat2

   Subroutine mpimat2_pm(mat2, orig1, orig2, nsize1, nsize2, istart)
      use MPI

      Implicit None
      integer ierr
      integer::typelist(2)
      integer ::imat(2), mat(2), start(2)
      integer ::istart
      integer ::orig1, orig2, nsize1, nsize2, mat2
      !allocate(struct%AS_ij(nsize,nsize))
      imat(1) = orig1
      imat(2) = orig2

      mat(1) = nsize1
      mat(2) = nsize2

      start(1) = 0
      start(2) = istart

      typelist(1) = MPI_DOUBLE_PRECISION
      typelist(2) = MPI_DOUBLE_PRECISION

      !write (*,*) nsize
      call MPI_TYPE_CREATE_SUBARRAY(2, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat2, ierr)
      call MPI_TYPE_COMMIT(mat2, ierr)
   End Subroutine mpimat2_pm

   Subroutine mpimat3_pm(mat3, nsize1, nsize2, nsize3)
      use MPI

      Implicit None
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
   End subroutine mpimat3_pm

   Subroutine mpimat4(mat4, nsize1, nsize2, nsize3, nsize4)
      use MPI

      Implicit None
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

      !write (*,*) nsize
      call MPI_TYPE_CREATE_SUBARRAY(4, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat4, ierr)
      call MPI_TYPE_COMMIT(mat4, ierr)
   End subroutine mpimat4

   Subroutine mpimat4int(mat4int, nsize1, nsize2, nsize3, nsize4)
      use MPI

      Implicit None
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
   End subroutine mpimat4int

   Subroutine mpimat5_old(mat5, nsize1, nsize2, nsize3, nsize4, nsize5)
      use MPI

      Implicit None
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

      !write (*,*) nsize
      call MPI_TYPE_CREATE_SUBARRAY(5, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat5, ierr)
      call MPI_TYPE_COMMIT(mat5, ierr)
   End subroutine mpimat5_old

   Subroutine mpimat5(mat5, nsize1, nsize2, nsize3, nsize4, nsize5, isize5, nbstart)
      use MPI

      Implicit None
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
      !print *, imat
      !print *,mat
      !print *,start
      !write (*,*) nsize
      call MPI_TYPE_CREATE_SUBARRAY(5, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat5, ierr)
      call MPI_TYPE_COMMIT(mat5, ierr)
   End subroutine mpimat5

   Subroutine mpimat5d(mat5, original_size, send_size, start_size)
      use MPI

      integer, intent(in)::original_size(5), send_size(5), start_size(5)
      integer ierr
      integer ::imat(5), mat(5), start(5)
      integer ::mat5
      mat5 = 0
      imat = original_size
      mat = send_size
      start = start_size
      !print *, imat
      !print *,mat
      !print *,start
      !write (*,*) nsize
      call MPI_TYPE_CREATE_SUBARRAY(5, imat, mat, start, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mat5, ierr)
      call MPI_TYPE_COMMIT(mat5, ierr)
   End subroutine mpimat5d

End module mpi_matrices

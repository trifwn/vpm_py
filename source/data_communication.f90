module data_communication
    use iso_fortran_env, only: real64
    use MPI
    implicit none

    private
    public :: commit_array, free_vars, distribute, collect

    ! MPI-related variables
    integer :: rank, size_mpi, ierr
    integer :: cart_comm

    ! DATA-related variables
    integer :: ndim
    integer, allocatable :: dims(:), coords(:)
    integer, allocatable :: global_sizes(:), local_sizes(:)
    integer, allocatable :: padded_counts(:), padded_starts(:, :), padded_sizes(:, :)
    integer, allocatable :: true_counts(:), true_starts(:, :), true_sizes(:, :)
    integer :: arr_padding = 1

    interface collect
        module subroutine collect_4D(global_data, local_data)
            real(real64), allocatable, intent(inout) :: global_data(:, :, :, :)
            real(real64), allocatable, intent(in) :: local_data(:, :, :, :)
        end subroutine collect_4D

        module subroutine collect_3D(global_data, local_data)
            real(real64), allocatable, intent(inout) :: global_data(:, :, :)
            real(real64), allocatable, intent(in) :: local_data(:, :, :)
        end subroutine collect_3D

        module subroutine collect_2D(global_data, local_data)
            real(real64), allocatable, intent(inout) :: global_data(:, :)
            real(real64), allocatable, intent(in) :: local_data(:, :)
        end subroutine collect_2D
    end interface collect

    interface distribute
        module subroutine distribute_4D(global_data, local_data)
            real(real64), target, intent(in) :: global_data(:, :, :, :)
            real(real64), allocatable, intent(out) :: local_data(:, :, :, :)
        end subroutine distribute_4D

        module subroutine distribute_3D(global_data, local_data)
            real(real64), target, intent(in) :: global_data(:, :, :)
            real(real64), allocatable, intent(out) :: local_data(:, :, :)
        end subroutine distribute_3D

        module subroutine distribute_2D(global_data, local_data)
            real(real64), target, intent(in) :: global_data(:, :)
            real(real64), allocatable, intent(out) :: local_data(:, :)
        end subroutine distribute_2D
    end interface distribute
contains

    subroutine free_vars()
        if (allocated(padded_counts)) deallocate (padded_counts)
        if (allocated(padded_starts)) deallocate (padded_starts)
        if (allocated(padded_sizes)) deallocate (padded_sizes)
        if (allocated(true_counts)) deallocate (true_counts)
        if (allocated(true_starts)) deallocate (true_starts)
        if (allocated(true_sizes)) deallocate (true_sizes)
        if (allocated(dims)) deallocate (dims)
        if (allocated(coords)) deallocate (coords)
        if (allocated(global_sizes)) deallocate (global_sizes)
        if (allocated(local_sizes)) deallocate (local_sizes)
        call MPI_Comm_free(cart_comm, ierr)
        ndim = 0
        arr_padding = 1
    end subroutine free_vars

    subroutine commit_array(array, dimensions, padding)
        implicit none
        real(real64), intent(in) :: array(..)
        integer, intent(in) :: dimensions(:)
        integer, intent(in) :: padding
        logical, allocatable :: periods(:)
        integer :: neq, global_nx, global_ny, global_nz, local_nx, local_ny, local_nz
        integer :: total_size, i

        ! Initialize MPI
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, size_mpi, ierr)

        arr_padding = padding

        ! Get dimensions of global data that are passed as an intrinsic array
        ndim = size(shape(array))
        if (allocated(global_sizes)) deallocate (global_sizes)
        if (allocated(local_sizes)) deallocate (local_sizes)
        allocate (global_sizes(ndim))
        allocate (local_sizes(ndim))
        allocate (periods(ndim))

        if (allocated(dims)) deallocate (dims)
        if (allocated(coords)) deallocate (coords)
        allocate (dims(ndim))
        allocate (coords(ndim))
        if (rank == 0) then
            do i = 1, ndim
                global_sizes(i) = size(array, i)
            end do
        end if
        ! Broadcast global sizes to all processes
        call MPI_Bcast(global_sizes, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        dims = dimensions

        periods = .False.
        call MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, .true., cart_comm, ierr)
        call MPI_Cart_coords(cart_comm, rank, ndim, coords, ierr)
        call calculate_access()
    end subroutine commit_array

    subroutine calculate_access()
        integer :: local_coords(ndim), local_size
        integer :: proc, i

        if (allocated(padded_counts)) deallocate (padded_counts)
        if (allocated(padded_starts)) deallocate (padded_starts)
        if (allocated(padded_sizes)) deallocate (padded_sizes)
        if (allocated(true_counts)) deallocate (true_counts)
        if (allocated(true_starts)) deallocate (true_starts)
        if (allocated(true_sizes)) deallocate (true_sizes)

        allocate (padded_counts(size_mpi))
        allocate (true_counts(size_mpi))

        allocate (padded_starts(size_mpi, ndim))
        allocate (true_starts(size_mpi, ndim))

        allocate (padded_sizes(size_mpi, ndim))
        allocate (true_sizes(size_mpi, ndim))

        padded_counts = 0
        padded_starts = 0
        do proc = 0, size_mpi - 1
            ! Get the coordinates of the current rank
            call MPI_Cart_coords(cart_comm, proc, ndim, local_coords, ierr)
            do i = 1, ndim
                local_size = (global_sizes(i))/dims(i)
                if (local_coords(i) == dims(i) - 1) then
                    local_size = global_sizes(i) - local_size*(dims(i) - 1)
                    true_starts(proc + 1, i) = global_sizes(i) - local_size
                else
                    true_starts(proc + 1, i) = local_coords(i)*local_size
                end if
                true_sizes(proc + 1, i) = local_size

                ! Check for neighbors and adjust local sizes by including ghost cells
                if (local_coords(i) > 0) then
                    padded_starts(proc + 1, i) = true_starts(proc + 1, i) - arr_padding
                    padded_sizes(proc + 1, i) = true_sizes(proc + 1, i) + arr_padding
                else
                    padded_starts(proc + 1, i) = true_starts(proc + 1, i)
                    padded_sizes(proc + 1, i) = true_sizes(proc + 1, i)
                end if

                if (local_coords(i) < dims(i) - 1) then
                    padded_sizes(proc + 1, i) = padded_sizes(proc + 1, i) + arr_padding
                end if
            end do
            true_counts(proc + 1) = product(true_sizes(proc + 1, :))
            padded_counts(proc + 1) = product(padded_sizes(proc + 1, :))
            if (rank == proc) then
                local_sizes = true_sizes(proc + 1, :)
            end if

        end do
        call MPI_Barrier(cart_comm, ierr)
    end subroutine calculate_access

    module subroutine distribute_4D(global_data, local_data)
        real(real64), target, intent(in) :: global_data(:, :, :, :)
        real(real64), allocatable, intent(out) :: local_data(:, :, :, :)
        integer :: source
        ! MPI

        ! Allocate local data array
        if (allocated(local_data)) deallocate (local_data)
        allocate (local_data(padded_sizes(rank + 1, 1), &
                             padded_sizes(rank + 1, 2), &
                             padded_sizes(rank + 1, 3), &
                             padded_sizes(rank + 1, 4)))
        source = 0
        if (rank == source) then
            ! Root process sends data to all other processes

            call send_to_all(global_data)
            ! Root process copies its own portion of data
            local_data = global_data( &
                         padded_starts(source + 1, 1) + 1:padded_starts(source + 1, 1) + local_sizes(1), &
                         padded_starts(source + 1, 2) + 1:padded_starts(source + 1, 2) + local_sizes(2), &
                         padded_starts(source + 1, 3) + 1:padded_starts(source + 1, 3) + local_sizes(3), &
                         padded_starts(source + 1, 4) + 1:padded_starts(source + 1, 4) + local_sizes(4))
        else
            call recv_from_source(local_data, source)
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine distribute_4D

    module subroutine distribute_3D(global_data, local_data)
        real(real64), target, intent(in) :: global_data(:, :, :)
        real(real64), allocatable, intent(out) :: local_data(:, :, :)
        integer :: source

        ! Allocate local data array
        if (allocated(local_data)) deallocate (local_data)
        allocate (local_data(padded_sizes(rank + 1, 1), &
                             padded_sizes(rank + 1, 2), &
                             padded_sizes(rank + 1, 3)))
        source = 0
        if (rank == source) then
            ! Root process sends data to all other processes

            call send_to_all(global_data)
            ! Root process copies its own portion of data
            local_data = global_data( &
                         padded_starts(source + 1, 1) + 1:padded_starts(source + 1, 1) + local_sizes(1), &
                         padded_starts(source + 1, 2) + 1:padded_starts(source + 1, 2) + local_sizes(2), &
                         padded_starts(source + 1, 3) + 1:padded_starts(source + 1, 3) + local_sizes(3))
        else
            call recv_from_source(local_data, source)
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine distribute_3D

    module subroutine distribute_2D(global_data, local_data)
        real(real64), target, intent(in) :: global_data(:, :)
        real(real64), allocatable, intent(out) :: local_data(:, :)
        integer :: source
        ! Allocate local data array
        if (allocated(local_data)) deallocate (local_data)
        allocate (local_data(padded_sizes(rank + 1, 1), padded_sizes(rank + 1, 2)))
        source = 0
        if (rank == source) then
            ! Root process sends data to all other processes
            call send_to_all(global_data)
            ! Root process copies its own portion of data
            local_data = global_data(padded_starts(1, 1) + 1:padded_starts(1, 1) + local_sizes(1), &
                                     padded_starts(1, 2) + 1:padded_starts(1, 2) + local_sizes(2))
        else
            call recv_from_source(local_data, source)
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine distribute_2D

    module subroutine collect_4D(global_data, local_data)
        real(real64), allocatable, intent(inout) :: global_data(:, :, :, :)
        real(real64), allocatable, intent(in) :: local_data(:, :, :, :)
        integer :: source
        ! MPI

        if (allocated(global_data)) deallocate (global_data)
        if (rank == 0) then
            allocate (global_data(global_sizes(1), global_sizes(2), global_sizes(3), global_sizes(4)))
        end if

        source = 0
        if (rank == source) then
            call recv_from_all(global_data)
            ! Root process copies its own portion of data
            global_data( &
                padded_starts(source + 1, 1) + 1:padded_starts(source + 1, 1) + local_sizes(1), &
                padded_starts(source + 1, 2) + 1:padded_starts(source + 1, 2) + local_sizes(2), &
                padded_starts(source + 1, 3) + 1:padded_starts(source + 1, 3) + local_sizes(3), &
                padded_starts(source + 1, 4) + 1:padded_starts(source + 1, 4) + local_sizes(4)) = local_data
        else
            call send_to_source(local_data, source)
        end if
    end subroutine collect_4D

    module subroutine collect_3D(global_data, local_data)
        real(real64), allocatable, intent(inout) :: global_data(:, :, :)
        real(real64), allocatable, intent(in) :: local_data(:, :, :)
        integer :: source
        if (allocated(global_data)) deallocate (global_data)
        if (rank == 0) then
            allocate (global_data(global_sizes(1), global_sizes(2), global_sizes(3)))
        end if

        source = 0
        if (rank == source) then
            call recv_from_all(global_data)
            ! Root process copies its own portion of data
            global_data( &
                padded_starts(source + 1, 1) + 1:padded_starts(source + 1, 1) + local_sizes(1), &
                padded_starts(source + 1, 2) + 1:padded_starts(source + 1, 2) + local_sizes(2), &
                padded_starts(source + 1, 3) + 1:padded_starts(source + 1, 3) + local_sizes(3)) = local_data
        else
            call send_to_source(local_data, source)
        end if
    end subroutine collect_3D

    module subroutine collect_2D(global_data, local_data)
        real(real64), allocatable, intent(inout) :: global_data(:, :)
        real(real64), allocatable, intent(in) :: local_data(:, :)
        integer :: source
        if (allocated(global_data)) deallocate (global_data)
        if (rank == 0) then
            allocate (global_data(global_sizes(1), global_sizes(2)))
        end if

        source = 0
        if (rank == source) then
            call recv_from_all(global_data)
            ! Root process copies its own portion of data
            global_data( &
                padded_starts(source + 1, 1) + 1:padded_starts(source + 1, 1) + local_sizes(1), &
                padded_starts(source + 1, 2) + 1:padded_starts(source + 1, 2) + local_sizes(2)) = local_data
        else
            call send_to_source(local_data, source)
        end if
    end subroutine collect_2D

    subroutine recv_from_source(local_data, source)
        real(real64), intent(out) :: local_data(..)
        integer :: source, my_starts(ndim), my_size(ndim)
        integer :: status(MPI_STATUS_SIZE)
        integer :: recv_subarray
        my_starts = 0
        my_size = padded_sizes(rank + 1, :)
        call MPI_Type_create_subarray(ndim, my_size, my_size, my_starts, &
                                      MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, recv_subarray, ierr)
        call MPI_Type_commit(recv_subarray, ierr)
        ! Non-root processes receive data
        call MPI_Recv(local_data, 1, recv_subarray, source, 0, MPI_COMM_WORLD, status, ierr)
        call MPI_Type_free(recv_subarray, ierr)
    end subroutine recv_from_source

    subroutine send_to_all(global_data)
        real(real64), intent(in) :: global_data(..)
        integer :: source, dest_starts(ndim), dest_size(ndim)
        integer :: status(MPI_STATUS_SIZE)
        integer :: send_subarray, i

        do i = 0, size_mpi - 1
            if (i == rank) cycle
            dest_starts = padded_starts(i + 1, :)
            dest_size = padded_sizes(i + 1, :)
            ! Create subarray type for sending data
            call MPI_Type_create_subarray(ndim, global_sizes, dest_size, dest_starts, &
                                          MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, send_subarray, ierr)
            call MPI_Type_commit(send_subarray, ierr)
            ! Root process sends data to all other processes
            call MPI_Send(global_data, 1, send_subarray, i, 0, MPI_COMM_WORLD, ierr)
            call MPI_Type_free(send_subarray, ierr)
        end do
    end subroutine send_to_all

    subroutine send_to_source(local_data, dest)
        real(real64), intent(in) :: local_data(..)
        integer, intent(in) :: dest
        integer :: send_subarray
        integer :: my_starts(ndim), my_sizes(ndim), my_true_sizes(ndim)

        ! Non-root processes send their data to root
        my_starts = true_starts(rank + 1, :) - padded_starts(rank + 1, :)
        my_sizes = padded_sizes(rank + 1, :)
        my_true_sizes = true_sizes(rank + 1, :)
        call MPI_Type_create_subarray(ndim, my_sizes, my_true_sizes, my_starts, &
                                      MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, send_subarray, ierr)
        call MPI_Type_commit(send_subarray, ierr)
        call MPI_Send(local_data, 1, send_subarray, dest, 0, cart_comm, ierr)
        call MPI_Type_free(send_subarray, ierr)
    end subroutine send_to_source

    subroutine recv_from_all(global_data)
        real(real64), intent(out) :: global_data(..)
        integer :: dest, true_proc_starts(ndim), true_proc_sizes(ndim)
        integer :: status(MPI_STATUS_SIZE)
        integer :: recv_subarray, local_coords(ndim)

        ! Root process receives data from all other processes
        do dest = 0, size_mpi - 1
            if (dest == rank) cycle
            true_proc_sizes = true_sizes(dest + 1, :)
            true_proc_starts = true_starts(dest + 1, :)
            call MPI_Type_create_subarray(ndim, global_sizes, true_proc_sizes, true_proc_starts, &
                                          MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, recv_subarray, ierr)
            call MPI_Type_commit(recv_subarray, ierr)
            call MPI_Recv(global_data, 1, recv_subarray, dest, 0, cart_comm, status, ierr)
            call MPI_Type_free(recv_subarray, ierr)
        end do
    end subroutine recv_from_all

end module data_communication

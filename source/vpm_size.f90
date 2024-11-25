Module vpm_size
    use base_types, only: dp
    use pmgrid, only: cartesian_grid

    ! Fine grid
    type(cartesian_grid)                 :: fine_grid
    type(cartesian_grid), allocatable    :: block_grids(:)
    type(cartesian_grid)                 :: coarse_grid
    integer                              :: my_block_idx

    ! Block grid
    integer, save               :: nb_i, nb_j, nb_k, NBlocks, ndumcell_coarse, ndumcell_bl
    real(dp)                    :: st, et

    public :: cartesian_grid
    public :: fine_grid, block_grids, coarse_grid, my_block_idx
    public :: print_vpm_size_info, define_sizes, get_domain_bounds_from_particles, get_fine_NN, get_fine_NNbl, &
              get_fine_Xbound
contains
    !> Defines Coarse and Fine GRID
    subroutine define_sizes
        use vpm_vars, only: interf_iproj, idefine, NTIME_PM, ilevmax, NBI, NBJ, NBK, ND
        use pmgrid, only: XMIN_pm, YMIN_pm, ZMIN_pm, &
                          XMAX_pm, YMAX_pm, ZMAX_pm, &
                          NXs_fine_bl, NXf_fine_bl, NYs_fine_bl, &
                          NYf_fine_bl, NZs_fine_bl, NZf_fine_bl, &
                          NXpm_fine, NYpm_fine, NZpm_fine, DVPM, &
                          DXpm, DYpm, DZpm, ndumcell, ncoarse
        use pmlib, only: definepm
        use parvar, only: XP
        use console_io, only: dummy_string, vpm_print, nocolor, blue, yellow, red, tab_level
        use MPI

        Implicit None
        integer    :: nsiz(3), nsiz_bl(3)
        integer    :: i, j, k, np, my_rank, ierr, nb, istep, lev
        integer    :: redifine_pm
        integer    :: NN_bl_tmp(6), NXbl, NYbl, NZbl, NXB, NYB, NZB

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        if ((NTIME_pm .eq. 0) .or. idefine .eq. 0) then
            redifine_pm = 1
        else
            redifine_pm = 0
        end if

        if (my_rank .eq. 0) then
            st = MPI_WTIME()
            write (dummy_string, *) 'Defining Sizes'
            call vpm_print(dummy_string, blue, 2)
            write (dummy_string, '(A, I5, A, I1)') achar(9)//'NTIME_PM=', NTIME_pm, &
                achar(9)//'Redefine=', redifine_pm
            call vpm_print(dummy_string, nocolor, 2)
        end if

        tab_level = tab_level + 1

        if (my_rank .eq. 0) then
            if (redifine_pm .eq. 1) then
                call get_domain_bounds_from_particles
                write (dummy_string, '(A)') 'The computational domain bounds are recalculated from the particle positions'
                call vpm_print(dummy_string, red, 2)
                ! Write the min and max values of the particle positions
                write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') &
                    achar(9)//'Particle XMIN=', minval(XP(1, :)), &
                    achar(9)//'Particle XMAX=', maxval(XP(1, :))
                call vpm_print(dummy_string, yellow, 2)
                write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') &
                    achar(9)//'Particle YMIN=', minval(XP(2, :)), &
                    achar(9)//'Particle YMAX=', maxval(XP(2, :))
                call vpm_print(dummy_string, yellow, 2)
                write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') &
                    achar(9)//'Particle ZMIN=', minval(XP(3, :)), &
                    achar(9)//'Particle ZMAX=', maxval(XP(3, :))
                call vpm_print(dummy_string, yellow, 2)
            end if

            write (dummy_string, '(A)') 'The computational domain bounds are:'
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, '(A, F10.5, A, F8.5, A, F10.5)') &
                achar(9)//'XMIN=', XMIN_pm, &
                achar(9)//'YMIN=', YMIN_pm, &
                achar(9)//'ZMIN=', ZMIN_pm
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') &
                achar(9)//'XMAX=', XMAX_pm, &
                achar(9)//'YMAX=', YMAX_pm, &
                achar(9)//'ZMAX=', ZMAX_pm
            call vpm_print(dummy_string, nocolor, 2)
        end if

        ! BRADCAST THE MIN AND MAX OF THE COMPUTATIONAL DOMAIN
        call MPI_BCAST(XMIN_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(YMIN_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ZMIN_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        call MPI_BCAST(XMAX_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(YMAX_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ZMAX_pm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        ! CREATE THE BOUNDING BOX
        fine_grid%Xbound(1) = XMIN_pm
        fine_grid%Xbound(2) = YMIN_pm
        fine_grid%Xbound(3) = ZMIN_pm
        fine_grid%Xbound(4) = XMAX_pm
        fine_grid%Xbound(5) = YMAX_pm
        fine_grid%Xbound(6) = ZMAX_pm

        fine_grid%Dpm(1) = DXpm
        fine_grid%Dpm(2) = DYpm
        fine_grid%Dpm(3) = DZpm

        nsiz(1) = NBI*ncoarse
        nsiz(2) = NBJ*ncoarse
        nsiz(3) = NBK*ncoarse

        if (my_rank .eq. 0) then
            write (dummy_string, '(A)') 'The extended fine domain is defined (no dummy cells)'
            call vpm_print(dummy_string, red, 2)
            write (dummy_string, '(A)') achar(9)//'The number of cells (nodes-1) must be divisible by the processor subdivision:'
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, '(A, I5, A, I5, A, I5)') &
                achar(9)//achar(9)//'X-dir: ', nsiz(1), &
                achar(9)//achar(9)//'Y-dir: ', nsiz(2), &
                achar(9)//achar(9)//'Z-dir: ', nsiz(3)
            call vpm_print(dummy_string, nocolor, 2)
        end if

        ! First Change Dpm so that the numbers of cells divides (that would be case 2)
        ! Here we change the number of cells
        ! by nsize i.e with NBI,NBJ,NBK, ncoarse,levmax depending on the criterion
        ! thats why ndumcell=0
        ndumcell = 0
        if (redifine_pm .eq. 1) then
            call definepm(3, fine_grid%Xbound, fine_grid%Dpm, ND, ndumcell, nsiz, fine_grid%NN, fine_grid%NN_bl)
        end if

        ! THE new extended domain is defined
        XMIN_pm = fine_grid%Xbound(1)
        YMIN_pm = fine_grid%Xbound(2)
        ZMIN_pm = fine_grid%Xbound(3)
        XMAX_pm = fine_grid%Xbound(4)
        YMAX_PM = fine_grid%Xbound(5)
        ZMAX_pm = fine_grid%Xbound(6)

        NXpm_fine = fine_grid%NN(1)
        NYpm_fine = fine_grid%NN(2)
        NZpm_fine = fine_grid%NN(3)

        NXs_fine_bl = fine_grid%NN_bl(1)
        NYs_fine_bl = fine_grid%NN_bl(2)
        NZs_fine_bl = fine_grid%NN_bl(3)
        NXf_fine_bl = fine_grid%NN_bl(4)
        NYf_fine_bl = fine_grid%NN_bl(5)
        NZf_fine_bl = fine_grid%NN_bl(6)

        DXpm = fine_grid%Dpm(1)
        DYpm = fine_grid%Dpm(2)
        DZpm = fine_grid%Dpm(3)

        DVpm = DXpm*DYpm
        if (ND .eq. 3) then
            DVpm = DVpm*DZpm
        end if

        if (my_rank .eq. 0) then
            write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') &
                achar(9)//'XMIN='//achar(9), XMIN_pm, &
                achar(9)//'YMIN='//achar(9), YMIN_pm, &
                achar(9)//'ZMIN='//achar(9), ZMIN_pm
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, '(A, F10.5, A, F10.5, A, F8.5)') &
                achar(9)//'XMAX='//achar(9), XMAX_pm, &
                achar(9)//'YMAX='//achar(9), YMAX_pm, &
                achar(9)//'ZMAX='//achar(9), ZMAX_pm
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, '(A, F10.5, A, F10.5, A, F10.5)') &
                achar(9)//'DX='//achar(9), DXpm, &
                achar(9)//'DY='//achar(9), DYpm, &
                achar(9)//'DZ='//achar(9), DZpm
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, '(A, I5, A, I5, A, I5)') &
                achar(9)//'Nodes X= ', NXpm_fine, &
                achar(9)//achar(9)//'Nodes Y= ', NYpm_fine, &
                achar(9)//achar(9)//'Nodes Z= ', NZpm_fine
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, "(A)") "The indexes of the coarse grid that do not include the dummy cells are:"
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, '(A, I5, A, I5, A, I5)') &
                achar(9)//'NXs_coarse=', NXs_fine_bl, &
                achar(9)//'NYs_coarse=', NYs_fine_bl, &
                achar(9)//'NZs_coarse=', NZs_fine_bl
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, '(A, I5, A, I5, A, I5)') &
                achar(9)//'NXf_coarse=', NXf_fine_bl, &
                achar(9)//'NYf_coarse=', NYf_fine_bl, &
                achar(9)//'NZf_coarse=', NZf_fine_bl
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, '(A, F8.5)') achar(9)//'DV=', DVpm
            call vpm_print(dummy_string, nocolor, 2)
        end if

        ! define block grids
        ! so they are divided by ncoarse and ilevmax
        ! so as to have the coarse information at the boundaries exactly.
        NBlocks = np
        if (.not. allocated(block_grids)) then
            ! allocate (Xbound_block(6, BLOCKS), NN_block(3, BLOCKS), NN_bl_block(6, BLOCKS))
            allocate (block_grids(NBlocks))
        end if

        ! DEVIDE THE DOMAIN INTO EQUAL BLOCKS
        ! NUMBER OF CELLS IN A BLOCK IN X DIRECTION
        NXB = int(nint(((fine_grid%Xbound(4) - fine_grid%Xbound(1))/fine_grid%Dpm(1))))
        ! NUMBER OF CELLS IN A BLOCK IN Y DIRECTION
        NYB = int(nint(((fine_grid%Xbound(5) - fine_grid%Xbound(2))/fine_grid%Dpm(2))))
        ! NUMBER OF CELLS IN A BLOCK IN Z DIRECTION
        NZB = int(nint(((fine_grid%Xbound(6) - fine_grid%Xbound(3))/fine_grid%Dpm(3))))
        NXbl = NXB/NBI
        NYbl = NYB/NBJ
        NZbl = NZB/NBK

        ndumcell_bl = ncoarse
        nsiz_bl = ncoarse!ndumcell_coarse!*2*2**ilevmax

        if (my_rank .eq. 0) then
            write (dummy_string, '(A)') 'The fine block domains are defined for each processor'
            call vpm_print(dummy_string, red, 2)
            write (dummy_string, '(A)') achar(9)//'The number of cells in each direction for the block grid'
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, '(A, I5, A, I5, A, I5)') achar(9)// &
                achar(9)//'Total Num Cells:'// &
                achar(9)//'X='//achar(9), NXB, &
                achar(9)//'Y='//achar(9), NYB, &
                achar(9)//'Z='//achar(9), NZB
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, '(A, I5, A, I5, A, I5)') achar(9)// &
                achar(9)//'Cells per Block:'// &
                achar(9)//'X='//achar(9), NXbl, &
                achar(9)//'Y='//achar(9), NYbl, &
                achar(9)//'Z='//achar(9), NZbl
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, '(A, I5)') achar(9)//achar(9)//'Block dummy cells:', ndumcell_bl
            call vpm_print(dummy_string, nocolor, 2)
        end if

        nb_i = -1
        nb_j = -1
        nb_k = -1
        do k = 1, NBK
            do j = 1, NBJ
                do i = 1, NBI
                    nb = (k - 1)*NBJ*NBI + (j - 1)*NBI + i
                    block_grids(nb)%Xbound(1) = fine_grid%Xbound(1) + (i - 1)*(NXbl)*fine_grid%Dpm(1)
                    block_grids(nb)%Xbound(4) = fine_grid%Xbound(1) + (i)*(NXbl)*fine_grid%Dpm(1)

                    block_grids(nb)%Xbound(2) = fine_grid%Xbound(2) + (j - 1)*(NYbl)*fine_grid%Dpm(2)
                    block_grids(nb)%Xbound(5) = fine_grid%Xbound(2) + (j)*(NYbl)*fine_grid%Dpm(2)

                    block_grids(nb)%Xbound(3) = fine_grid%Xbound(3) + (k - 1)*(NZbl)*fine_grid%Dpm(3)
                    block_grids(nb)%Xbound(6) = fine_grid%Xbound(3) + (k)*(NZbl)*fine_grid%Dpm(3)

                    block_grids(nb)%Dpm(1:3) = fine_grid%Dpm(1:3)

                    call definepm(1, block_grids(nb)%Xbound, block_grids(nb)%Dpm, ND, ndumcell_bl, nsiz_bl, &
                                  block_grids(nb)%NN, block_grids(nb)%NN_bl)

                    if (my_rank .eq. 0) then
                        write (dummy_string, '(A,I3,A,3I3,A)') achar(9)//'Block ', nb, " = (", i, j, k, ")"
                        call vpm_print(dummy_string, blue, 2)
                        call print_grid_info(block_grids(nb))
                    end if

                    if (nb .eq. my_rank + 1) then
                        nb_i = i
                        nb_j = j
                        nb_k = k
                    end if
                end do
            end do
        end do
        if ((nb_i .eq. -1) .or. (nb_j .eq. -1) .or. (nb_k .eq. -1)) then
            print *, "Processor with rank ", my_rank, " does not have a block"
            stop
        end if

        !define coarse grid must cover block grids
        fine_grid%Xbound(1) = XMIN_pm!minval(Xbound_bl(1,:))
        fine_grid%Xbound(2) = YMIN_pm!minval(Xbound_bl(2,:))
        fine_grid%Xbound(3) = ZMIN_pm!minval(Xbound_bl(3,:))
        fine_grid%Xbound(4) = XMAX_pm!maxval(Xbound_bl(4,:))
        fine_grid%Xbound(5) = YMAX_pm!maxval(Xbound_bl(5,:))
        fine_grid%Xbound(6) = ZMAX_pm!maxval(Xbound_bl(6,:))

        coarse_grid%Xbound(1:6) = fine_grid%Xbound(1:6)
        coarse_grid%Dpm(1:3) = ncoarse*fine_grid%Dpm(1:3)
        ndumcell_coarse = 4!2**ilevmax
        nsiz_bl = 2**ilevmax
        call definepm(1, coarse_grid%Xbound, coarse_grid%Dpm, ND, ndumcell_coarse, nsiz_bl, coarse_grid%NN, coarse_grid%NN_bl)

        if (my_rank .eq. 0) then
            write (dummy_string, '(A)') 'The extended coarse domain is redefined (with dummy cells)'
            call vpm_print(dummy_string, red, 2)
            call print_grid_info(coarse_grid)
        end if

        ! CALCULATE THE NUMBER OF LEVELS
        if (my_rank .eq. 0) then
            do lev = 0, ilevmax
                istep = 2**lev
            !!!!!!!!if not divided exactly dummy cell
                if (ND .eq. 2) then
                    if (int((block_grids(1)%NN_bl(4) - block_grids(1)%NN_bl(1))/istep) .eq. 0 .or. &
                        int((block_grids(1)%NN_bl(5) - block_grids(1)%NN_bl(2))/istep) .eq. 0) then
                        ilevmax = lev - 1
                        print *, 'Changing number of levels', ilevmax
                        exit
                    end if
                else
                    if (int((block_grids(1)%NN_bl(4) - block_grids(1)%NN_bl(1))/istep) .eq. 0 .or. &
                        int((block_grids(1)%NN_bl(5) - block_grids(1)%NN_bl(2))/istep) .eq. 0 .or. &
                        int((block_grids(1)%NN_bl(6) - block_grids(1)%NN_bl(3))/istep) .eq. 0) then
                        ilevmax = lev - 1
                        print *, 'Changing number of levels', ilevmax
                        exit
                    end if
                end if
            end do
        end if
        call MPI_BCAST(ilevmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        tab_level = tab_level - 1
        if (my_rank .eq. 0) then
            et = MPI_WTIME()
            write (dummy_string, *) 'Defining Sizes finished in:', int((et - st)/60), 'm', mod(et - st, 60.d0), 's'
            call vpm_print(dummy_string, yellow, 2)
        end if
    end subroutine define_sizes

    subroutine get_domain_bounds_from_particles
        use parvar, only: NVR, XP
        use vpm_vars, only: interf_iproj
        use pmgrid, only: XMIN_pm, YMIN_pm, ZMIN_pm, DXpm, DYpm, DZpm, XMAX_pm, &
                          YMAX_pm, ZMAX_pm
        real(dp) :: XMIN_pm_old, YMIN_pm_old, ZMIN_pm_old, XMAX_pm_old, YMAX_pm_old, ZMAX_pm_old, &
                    X_mean, Y_mean, Z_mean
        logical  :: bounds_changed

        XMIN_pm_old = XMIN_pm
        YMIN_pm_old = YMIN_pm
        ZMIN_pm_old = ZMIN_pm
        XMAX_pm_old = XMAX_pm
        YMAX_pm_old = YMAX_pm
        ZMAX_pm_old = ZMAX_pm

        ! Get mean X, Y, Z
        X_mean = sum(XP(1, 1:NVR))/NVR
        Y_mean = sum(XP(2, 1:NVR))/NVR
        Z_mean = sum(XP(3, 1:NVR))/NVR

        XMIN_pm = minval(XP(1, 1:NVR)) - interf_iproj*DXpm
        YMIN_pm = minval(XP(2, 1:NVR)) - interf_iproj*DYpm
        ZMIN_pm = minval(XP(3, 1:NVR)) - interf_iproj*DZpm

        XMAX_pm = maxval(XP(1, 1:NVR)) + interf_iproj*DXpm
        YMAX_pm = maxval(XP(2, 1:NVR)) + interf_iproj*DYpm
        ZMAX_pm = maxval(XP(3, 1:NVR)) + interf_iproj*DZpm

        ! Normalize the domain bounds to be multiple of DX, DY, DZ
        XMIN_pm = X_mean - ceiling((X_mean - XMIN_pm)/DXpm)*DXpm
        YMIN_pm = Y_mean - ceiling((Y_mean - YMIN_pm)/DYpm)*DYpm
        ZMIN_pm = Z_mean - ceiling((Z_mean - ZMIN_pm)/DZpm)*DZpm

        XMAX_pm = X_mean + ceiling((XMAX_pm - X_mean)/DXpm)*DXpm
        YMAX_pm = Y_mean + ceiling((YMAX_pm - Y_mean)/DYpm)*DYpm
        ZMAX_pm = Z_mean + ceiling((ZMAX_pm - Z_mean)/DZpm)*DZpm

        ! Check if the domain bounds have changed
        if (XMIN_pm .ne. XMIN_pm_old .or. YMIN_pm .ne. YMIN_pm_old .or. ZMIN_pm .ne. ZMIN_pm_old .or. &
            XMAX_pm .ne. XMAX_pm_old .or. YMAX_pm .ne. YMAX_pm_old .or. ZMAX_pm .ne. ZMAX_pm_old) then
            bounds_changed = .TRUE.
        else
            bounds_changed = .FALSE.
        end if
    end subroutine get_domain_bounds_from_particles

    subroutine get_fine_NN(NN_out) bind(C, name='get_NN')
        use iso_c_binding
        implicit none
        integer(c_int), dimension(3) :: NN_out

        NN_out = fine_grid%NN
    end subroutine get_fine_NN

    subroutine get_fine_NNbl(NN_bl_out) bind(C, name='get_NN_bl')
        use iso_c_binding
        implicit none
        integer(c_int), dimension(6) :: NN_bl_out
        NN_bl_out = fine_grid%NN_bl
    end subroutine get_fine_NNbl

    subroutine get_fine_Xbound(Xbound_out) bind(C, name='get_fine_bounds')
        use iso_c_binding
        implicit none
        real(c_double), dimension(6) :: Xbound_out
        Xbound_out = fine_grid%Xbound
    end subroutine get_fine_Xbound

    subroutine print_vpm_size_info()
        use console_io, only: dp_1d_array_info, i_1d_array_info, dp_2d_alloc_info, i_2d_alloc_info, &
                              dp_1d_alloc_info
        print *, "VPM_SIZE INFO"
        print *, "============"
        print *, ""
        print *, achar(9)//"nb_i", " = ", nb_i
        print *, achar(9)//"nb_j", " = ", nb_j
        print *, achar(9)//"nb_k", " = ", nb_k
        print *, achar(9)//"BLOCKS", " = ", NBlocks
        print *, achar(9)//"ndumcell_coarse", " = ", ndumcell_coarse
        print *, achar(9)//"ndumcell_bl", " = ", ndumcell_bl

    end subroutine print_vpm_size_info

    subroutine print_grid_info(grid_in)
        use console_io, only: dummy_string, vpm_print, nocolor, blue, yellow, red, tab_level
        use vpm_vars, only: ND
        implicit none
        type(cartesian_grid) :: grid_in
        real(dp)    :: DV
        write (dummy_string, '(A, F8.5, A, F10.5, A, F10.5)') &
            achar(9)//'XMIN='//achar(9), grid_in%Xbound(1), &
            achar(9)//'YMIN='//achar(9), grid_in%Xbound(2), &
            achar(9)//'ZMIN='//achar(9), grid_in%Xbound(3)
        call vpm_print(dummy_string, nocolor, 2)
        write (dummy_string, '(A, F8.5, A, F10.5, A, F10.5)') &
            achar(9)//'XMAX='//achar(9), grid_in%Xbound(4), &
            achar(9)//'YMAX='//achar(9), grid_in%Xbound(5), &
            achar(9)//'ZMAX='//achar(9), grid_in%Xbound(6)
        call vpm_print(dummy_string, nocolor, 2)
        write (dummy_string, '(A, F8.5, A, F10.5, A, F10.5)') &
            achar(9)//'DX='//achar(9), grid_in%Dpm(1), &
            achar(9)//'DY='//achar(9), grid_in%Dpm(2), &
            achar(9)//'DZ='//achar(9), grid_in%Dpm(3)
        call vpm_print(dummy_string, nocolor, 2)
        write (dummy_string, '(A, I5, A, I5, A, I5)') &
            achar(9)//'Nodes X= ', grid_in%NN(1), &
            achar(9)//achar(9)//'Nodes Y= ', grid_in%NN(2), &
            achar(9)//achar(9)//'Nodes Z= ', grid_in%NN(3)
        call vpm_print(dummy_string, nocolor, 2)
        write (dummy_string, "(A)") "The indexes of the coarse grid that do not include the dummy cells are:"
        call vpm_print(dummy_string, nocolor, 2)
        write (dummy_string, '(A, I5, A, I5, A, I5)') &
            achar(9)//'NXs=', grid_in%NN_bl(1), &
            achar(9)//'NYs=', grid_in%NN_bl(2), &
            achar(9)//'NZs=', grid_in%NN_bl(3)
        call vpm_print(dummy_string, nocolor, 2)
        write (dummy_string, '(A, I5, A, I5, A, I5)') &
            achar(9)//'NXf=', grid_in%NN_bl(4), &
            achar(9)//'NYf=', grid_in%NN_bl(5), &
            achar(9)//'NZf=', grid_in%NN_bl(6)
        call vpm_print(dummy_string, nocolor, 2)
        DV = grid_in%Dpm(1)*grid_in%Dpm(2)
        if (ND .eq. 3) then
            DV = DV*grid_in%Dpm(3)
        end if
        write (dummy_string, '(A, F8.5)') achar(9)//'DV=', DV
        call vpm_print(dummy_string, nocolor, 2)
    end subroutine print_grid_info

End Module vpm_size

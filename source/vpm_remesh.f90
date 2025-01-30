module vpm_remesh
    use vpm_types, only: dp, cartesian_grid
    use console_io, only: vpm_print, nocolor, tab_level, yellow, dummy_string, red, blue

    implicit none

contains

    subroutine interpolate_and_remesh_particles(npar_per_cell, XPR, QPR, UPR, GPR, NVR_in, NVR_size_in, cutoff_value)
        use vpm_functions, only: project_particles_parallel
        use parvar, only: associate_particles
        use vpm_size, only: fine_grid
        use pmgrid, only: RHS_pm
        use MPI

        integer, intent(in)                                             :: npar_per_cell
        real(dp), intent(inout), allocatable, target, dimension(:, :)   :: XPR, QPR, GPR, UPR
        integer, intent(inout)                                          :: NVR_in, NVR_size_in
        real(dp), intent(in), optional                                  :: cutoff_value
        ! LOCAL VARIABLES
        integer                          :: my_rank, ierr, np

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        ! Project particles on PM grid to get RHS
        call associate_particles(NVR_in, NVR_size_in, XPR, QPR, UPR, GPR)
        call project_particles_parallel

        ! Remesh particles
        if (present(cutoff_value)) then
            call remesh_particles_3d(RHS_pm, fine_grid, npar_per_cell, XPR, QPR, UPR, GPR, NVR_in, cutoff_value)
        else
            call remesh_particles_3d(RHS_pm, fine_grid, npar_per_cell, XPR, QPR, UPR, GPR, NVR_in)
        end if
    end subroutine interpolate_and_remesh_particles

    !--------------------------------------------------------------------------------
    !subroutine remesh_particles
    !This subroutine remeshes particles  on the pm grid (4 particles per cell)
    !---------------------------------------------------------------------------------------------------
    !->PM grid is orthogonal (Volume of PM cell
    !--------------------------------------------------------------------------------!
    !-->The loops starts from 2 because we need cells that DO NOT contain particles  !
    !-->Total Number of Cells                                                        !
    !--------------------------------------------------------------------------------!
    subroutine remesh_particles_3d(RHS_pm, grid, npar_per_cell, XPR, QPR, UPR, GPR, NVRR, cutoff_value)
        use vpm_vars, only: mrem, neqpm, interf_iproj, V_ref
        use vpm_interpolate, only: interpolate_particle_Q
        use parvar, only: NVR, XP, QP, GP, UP, NVR_size
        use MPI

        implicit None

        ! PARAMETERS
        type(cartesian_grid), intent(in)                              :: grid
        real(dp), intent(in)                                          :: RHS_pm(neqpm, grid%NN(1), grid%NN(2), grid%NN(3))
        integer, intent(in)                                           :: npar_per_cell
        real(dp), intent(out), allocatable, target, dimension(:, :)   :: XPR, QPR, GPR, UPR
        integer, intent(inout)                                        :: NVRR
        real(dp), intent(in), optional                                :: cutoff_value

        ! LOCAL VARIABLES
        real(dp), dimension(8)           :: X, Y, Z
        real(dp), allocatable            :: XC(:), YC(:), ZC(:), QC(:)
        real(dp), allocatable, target    :: XP_tmp(:, :), QP_tmp(:, :)
        integer                          :: i, j, k, ncell, npar
        integer                          :: nxstart, nxfin, nystart, nyfin, nzstart, nzfin
        integer                          :: nc
        integer                          :: NVR_old
        integer                          :: my_rank, ierr, np, NN(3), NN_bl(6)
        real(dp)                         :: Xbound(6), Dpm(3), wmag, cutoff, XMIN_pm, YMIN_pm, ZMIN_pm, DVpm
        real(dp)                         :: st, et

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

        if (present(cutoff_value)) then
            cutoff = cutoff_value
        else
            cutoff = 1e-09
        end if

        if (my_rank .eq. 0) then
            st = MPI_WTIME()
            write (dummy_string, "(A, E8.3)") 'Remeshing with cutoff value: ', cutoff
            call vpm_print(dummy_string, blue, 0)
            write (dummy_string, "(A)") 'RHS_PM will be used to remesh the particles'
            call vpm_print(dummy_string, nocolor, 2)
            ! Print the Dimensions of RHS_PM and the Maximum and Minimum values
            write (dummy_string, "(A, I5, A, I5, A, I5)") achar(9)//'Dimensions of RHS_PM:', grid%NN(1), 'x', grid%NN(2), 'x', grid%NN(3)
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, "(A, F8.2)") achar(9)//'Maximal value of RHS_PM:', maxval(RHS_pm)
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, "(A, F8.2)") achar(9)//'Minimal value of RHS_PM:', minval(RHS_pm)
            call vpm_print(dummy_string, nocolor, 2)

        end if
        tab_level = tab_level + 1

        NVR_old = NVRR

        ! Get PM grid parameters
        Dpm = grid%Dpm
        Xbound = grid%Xbound
        NN = grid%NN
        NN_bl = grid%NN_bl

        NN = NN*mrem
        NN_bl = NN_bl*mrem
        Dpm(1) = (Xbound(4) - Xbound(1))/(NN(1) - 1)
        Dpm(2) = (Xbound(5) - Xbound(2))/(NN(2) - 1)
        Dpm(3) = (Xbound(6) - Xbound(3))/(NN(3) - 1)
        DVpm = Dpm(1)*Dpm(2)*Dpm(3)
        XMIN_pm = Xbound(1) 
        YMIN_pm = Xbound(2)
        ZMIN_pm = Xbound(3)

        ! Particle remeshing
        if (my_rank .eq. 0) then
            ncell = npar_per_cell
            allocate (XC(ncell), YC(ncell), ZC(ncell), QC(ncell))

            if (ncell .eq. 1) then
                nxstart = NN_bl(1) + interf_iproj/2
                nystart = NN_bl(2) + interf_iproj/2
                nzstart = NN_bl(3) + interf_iproj/2
                nxfin = NN_bl(4) - interf_iproj/2
                nyfin = NN_bl(5) - interf_iproj/2
                nzfin = NN_bl(6) - interf_iproj/2
            else
                nxstart = NN_bl(1)
                nystart = NN_bl(2)
                nzstart = NN_bl(3)
                nxfin = NN_bl(4) - 1
                nyfin = NN_bl(5) - 1
                nzfin = NN_bl(6) - 1
            end if

            write (dummy_string, "(A)") 'The sampling for the particles happens on the intervals:'
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, "(A, I5, A, I5)") '-X: start=', nxstart, " finish=", nxfin
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, "(A, I5, A, I5)") '-Y: start=', nystart, " finish=", nyfin
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, "(A, I5, A, I5)") '-Z: start=', nzstart, " finish=", nzfin
            call vpm_print(dummy_string, nocolor, 2)

            ! Total number of particles
            NVR = (nxfin - nxstart + 1)
            NVR = (nyfin - nystart + 1)*NVR
            NVR = (nzfin - nzstart + 1)*NVR
            NVR = NVR*ncell

            write (dummy_string, "(A,I5)") 'Allocating XP_tmp and QP_tmp with NVR:', NVR
            call vpm_print(dummy_string, nocolor, 2)

            allocate (XP_tmp(3, NVR), QP_tmp(neqpm + 1, NVR))
            XP_tmp = 0
            QP_tmp = 0
            npar = 0
            V_ref = 1.d0/float(ncell)*DVpm

            !!$omp parallel private(i,j,k,npar,X,Y,Z) num_threads(OMPTHREADS)
            !!$omp do
            do k = nzstart, nzfin
                do j = nystart, nyfin
                    do i = nxstart, nxfin

                        X(1) = XMIN_pm + Dpm(1)*(i - 1) !
                        Y(1) = YMIN_pm + Dpm(2)*(j - 1) !
                        Z(1) = ZMIN_pm + Dpm(3)*(k - 1) !
                        !-> Get PM cell nodes orthogonal structured grid
                        if (ncell .gt. 1) then
                            X(2) = XMIN_pm + Dpm(1)*(i)
                            X(3) = XMIN_pm + Dpm(1)*(i)
                            X(4) = XMIN_pm + Dpm(1)*(i - 1)
                            X(5) = XMIN_pm + Dpm(1)*(i - 1)
                            X(6) = XMIN_pm + Dpm(1)*(i)
                            X(7) = XMIN_pm + Dpm(1)*(i)
                            X(8) = XMIN_pm + Dpm(1)*(i - 1)

                            Y(2) = YMIN_pm + Dpm(2)*(j - 1)
                            Y(3) = YMIN_pm + Dpm(2)*(j)
                            Y(4) = YMIN_pm + Dpm(2)*(j)
                            Y(5) = YMIN_pm + Dpm(2)*(j - 1)
                            Y(6) = YMIN_pm + Dpm(2)*(j - 1)
                            Y(7) = YMIN_pm + Dpm(2)*(j)
                            Y(8) = YMIN_pm + Dpm(2)*(j)

                            Z(2) = ZMIN_pm + Dpm(3)*(k - 1)
                            Z(3) = ZMIN_pm + Dpm(3)*(k - 1)
                            Z(4) = ZMIN_pm + Dpm(3)*(k - 1)
                            Z(5) = ZMIN_pm + Dpm(3)*(k)
                            Z(6) = ZMIN_pm + Dpm(3)*(k)
                            Z(7) = ZMIN_pm + Dpm(3)*(k)
                            Z(8) = ZMIN_pm + Dpm(3)*(k)

                            YC = cell3d_interp_euler(Y, ncell, 2)
                            ZC = cell3d_interp_euler(Z, ncell, 2)
                            XC = cell3d_interp_euler(X, ncell, 2)
                            do nc = 1, ncell
                                npar = npar + 1
                                XP_tmp(1, npar) = XC(nc)
                                XP_tmp(2, npar) = YC(nc)
                                XP_tmp(3, npar) = ZC(nc)
                                QP_tmp(neqpm + 1, npar) = DVpm/float(ncell)
                            end do
                        else
                            wmag = sqrt(RHS_pm(1, i, j, k)**2 + RHS_pm(2, i, j, k)**2 + RHS_pm(3, i, j, k)**2)
                            if (wmag .lt. cutoff) cycle
                            npar = npar + 1
                            XP_tmp(1, npar) = X(1)
                            XP_tmp(2, npar) = Y(1)
                            XP_tmp(3, npar) = Z(1)

                            QP_tmp(neqpm + 1, npar) = DVpm
                            QP_tmp(1:neqpm, npar) = RHS_pm(1:neqpm, i, j, k)*DVpm
                        end if
                    end do
                end do
            end do
            !!$omp enddo
            !!$omp endparallel

            NVR = npar
        end if

        ! BCAST NEW NVR
        call MPI_BCAST(NVR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        NVRR = NVR
        NVR_size = NVR
        if (allocated(XPR)) deallocate (XPR)
        if (allocated(QPR)) deallocate (QPR)
        allocate (XPR(3, NVR), QPR(neqpm + 1, NVR))
        
        if (allocated(GPR)) deallocate (GPR)
        if (allocated(UPR)) deallocate (UPR)
        allocate (GPR(3, NVR), UPR(3, NVR))

        if (allocated(XP_tmp).and.allocated(QP_tmp)) then
            XPR = XP_tmp(:, 1:NVR)
            QPR = QP_tmp(:, 1:NVR)
        else
            XPR = 0
            QPR = 0
        end if

        GPR = 0
        UPR = 0

        nullify (XP, QP, GP, UP)
        XP => XPR
        QP => QPR
        UP => UPR
        GP => GPR

        if (my_rank .eq. 0) then
            deallocate (XP_tmp, QP_tmp)
            deallocate (XC, YC, ZC, QC)

            if (ncell .gt. 1) call interpolate_particle_Q(RHS_pm, XP, QP(1:neqpm, :), NVR, 4, NVR_size)
            write (dummy_string, *) 'After remesh'
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9), 'Number of particles before', NVR_old
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9), 'Number of particles after', NVR
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9), 'Volume of a cell', DVpm
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9), 'Number of cells', grid%NN(1), grid%NN(2), grid%NN(3) 
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9), 'Maximal value of QPR', maxval(QP(:, :))
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9), 'Minimal value of QPR', minval(QP(:, :))
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9), 'Maximal value of XP', maxval(XP(:, :))
            call vpm_print(dummy_string, nocolor, 2)
            write (dummy_string, *) achar(9), 'Minimal value of XP', minval(XP(:, :))
            call vpm_print(dummy_string, nocolor, 2)
        end if

        if (my_rank .eq. 0) then
            et = MPI_WTIME()
            write (dummy_string, "(A,I5,A,F8.2,A)") &
                'Remeshing finished in:'//achar(9), int((et - st)/60), ' m', mod(et - st, 60.d0), ' s'
            call vpm_print(dummy_string, yellow, 1)
        end if
        tab_level = tab_level - 1
    end subroutine remesh_particles_3d

    !--------------------------------------------------------------------------------
    !>@function
    ! subroutine    cell3d_interp_euler
    !>
    !>@author Papis
    !>
    !>@brief
    !>subroutine cell3d_interp_euler creates 4 or more particles per cell using ksi ita
    !>coordinates
    !REVISION HISTORY
    !> 17/7/2013 - Initial Version
    !> TODO_dd
    !>
    !>@param [in]  F is the value at the global coordinates
    !>@param [out] FC is the value at global coordinates of the interpolated value
    !--------------------------------------------------------------------------------
    module function cell3d_interp_euler(F, N, M) result(FC)
        use iso_fortran_env
        implicit none

        integer, parameter :: dp = real64

        real(dp), dimension(8), intent(in) :: F
        integer, intent(in) :: N, M
        real(dp), dimension(N) :: FC

        ! LOCAL VARIABLES
        real(dp) :: KSIC(8), HTAC(8), ZETAC(8)
        real(dp), dimension(N) :: KSI, HTA, ZETA
        integer :: i, j
        real(dp) :: addit
        integer :: NTEMP

        NTEMP = N
        !-->Define KSI,HTA corners
        KSIC = [-1.d0, 1.d0, 1.d0, -1.d0, -1.d0, 1.d0, 1.d0, -1.d0]
        HTAC = [-1.d0, -1.d0, 1.d0, 1.d0, -1.d0, -1.d0, 1.d0, 1.d0]
        ZETAC = [-1.d0, -1.d0, -1.d0, -1.d0, 1.d0, 1.d0, 1.d0, 1.d0]

        FC = 0.0_dp
        write (*, *) 'NTEMP', NTEMP, N, M
        call get_ksi_ita_pos_3d(N, M, KSIC, HTAC, ZETAC, KSI, HTA, ZETA)

        do i = 1, 8 !cell nodes
            do j = 1, NTEMP
                addit = F(i)
                addit = addit*(1.0_dp + KSI(j)*KSIC(i))
                addit = addit*(1.0_dp + HTA(j)*HTAC(i))
                addit = addit*(1.0_dp + ZETA(j)*ZETAC(i))
                FC(j) = FC(j) + addit
            end do
        end do

        FC = 0.125d0*FC ! 1/8
    end function cell3d_interp_euler

    !--------------------------------------------------------------------------------
    !>@function
    ! subroutine    get_ksi_ita_pos
    !>
    !>@author Papis
    !>
    !>@brief
    !>subroutine get_ksi_ita_pos  depending on the defined number of particles(must be perfect square)
    !REVISION HISTORY
    !> 22/7/2013 - Initial Version
    !> TODO_dd
    !>
    !>@param [in]  N is the number of particles
    !>@param [in]  KSIC(4),HTAC(4) is the corner coordinates in the KSI,HTA
    !>@param [out] KSI(2*N),HTA(2*N) local position
    !--------------------------------------------------------------------------------
    subroutine get_ksi_ita_pos_3d(N, M, KSIC, HTAC, ZETAC, KSI, HTA, ZETA)
        use vpm_types, only: dp
        implicit none

        integer, intent(in)           :: M
        integer, intent(in)           :: N
        ! N is the number of particles to remesh
        ! M is 2
        real(dp), intent(in)  :: KSIC(8), HTAC(8), ZETAC(8)
        real(dp), intent(out) :: KSI(N), HTA(N), ZETA(N)
        real(dp)              :: DKSI, DHTA, DZETA
        integer                       :: i, j, k, nod

        !--> find position minus ksi minus ita quadrant and then by symmerty * 4
        KSI = 0.d0
        HTA = 0.d0

        DKSI = dabs(2.d0/float(M))
        DHTA = dabs(2.d0/float(M))
        DZETA = dabs(2.d0/float(M))
        write (*, *) 'SIZE of KSI ', size(KSI)
        do k = 1, M
            do j = 1, M
                do i = 1, M
                    nod = (k - 1)*M*M + (j - 1)*M + i
                    write (*, *) 'nod', nod
                    KSI(nod) = KSIC(1) + (i - 1./2.)*DKSI
                    HTA(nod) = HTAC(1) + (j - 1./2.)*DHTA
                    ZETA(nod) = ZETAC(1) + (k - 1./2.)*DZETA
                end do
            end do
        end do

    end subroutine get_ksi_ita_pos_3d
end module vpm_remesh

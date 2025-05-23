!--------------------------------------------------------------------------------
!>@file
!>@brief Poisson solver library.This source code file(pmlib.f90) contains the library
!!for solving poisson problem on a structured grid with constant DX,DY,DZ
!--------------------------------------------------------------------------------

!>@brief This module defines the variables that will be used internally in the library
!>       All variables are private
module pmlib
    use vpm_types, only: dp
    use constants, only: pi, pi2, pi4
    use console_io, only: vpm_print, red, blue, green, nocolor, yellow, dummy_string, tab_level, VERBOCITY
    implicit none
    ! integer, parameter            :: FISHPACK  = 0,   &
    !                                  MUDPACK   = 1,   &
    !                                  INTEL_MKL = 2
    ! integer                       :: SOLVER = MUDPACK
    
    real(dp), save                :: XMIN_pm, XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
    real(dp), save                :: DXpm, DYpm, DZpm

    integer, save                 :: NVR, NXpm, NYpm, NZpm, ND
    integer, save                 :: NXs_bl, NYs_bl, NXf_bl, NYf_bl, NZs_bl, NZf_bl

    integer, save                 :: nbound, levmax
    integer, save                 :: itree, ibctyp

    !Here pointers are defined which will be assigned in the external data to save up space
    real(dp), pointer             :: SOL_pm(:, :, :, :), RHS_pm(:, :, :, :), QP(:, :), XP(:, :)

    real(dp), allocatable         :: SOL_0_pm(:, :, :, :), source_bound(:, :), x_s(:, :), y_s(:, :), &
                                     z_s(:, :), d_s(:)
    real(dp), allocatable, save   :: source_bound_lev(:, :, :), xs_lev(:, :), ys_lev(:, :), &
                                     zs_lev(:, :), ds_lev(:, :)
    integer, allocatable, save   :: nbound_lev(:), ilev_t(:, :, :)

    private :: XMIN_pm, XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm, DXpm, DYpm, DZpm
    private :: NVR, NXpm, NYpm, NZPm, ND
    private :: NXs_bl, NYs_bl, NXf_bl, NYf_bl, NZs_bl, NZf_bl
    private :: itree, ibctyp
    private :: SOL_pm, RHS_pm, QP, XP
    private :: SOL_0_pm, source_bound, x_s, y_s, z_s, d_s
    private :: source_bound_lev, xs_lev, ys_lev, zs_lev, ds_lev
    private :: nbound, ilev_t

! pinfdomain.f90
    interface infdomain
        module subroutine infdomain(neqs, neqf)
            implicit none
            integer, intent(in) :: neqs, neqf
        end subroutine infdomain
    end interface
    interface build_level_nbound
        module subroutine build_level_nbound(NXs, NXf, NYs, NYf, neqs, neqf)
            implicit none
            integer, intent(in)     :: NXs, NXf, NYs, NYf, neqs, neqf
        end subroutine build_level_nbound
    end interface
    interface infdomain_3D
        module subroutine infdomain_3D(neqs, neqf)
            implicit none
            integer, intent(in) :: neqs, neqf
        end subroutine infdomain_3D
    end interface

    interface build_level_nbound_3d
        module subroutine build_level_nbound_3d(NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
            implicit none
            integer, intent(in)     :: NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf
        end subroutine build_level_nbound_3d
    end interface

! pmsolve.f90
    interface solve_eq
        module subroutine solve_eq(NXs, NXf, NYs, NYf, neq)
            implicit none
            integer, intent(in)  :: NXs, NXf, NYs, NYf, neq
        end subroutine solve_eq
    end interface

    interface solve_eq_0
        module subroutine solve_eq_0(NXs, NXf, NYs, NYf, neq)
            implicit none
            integer, intent(in)   :: NXs, NXf, NYs, NYf, neq
        end subroutine solve_eq_0
    end interface

    interface solve_eq_3D
        module subroutine solve_eq_3D(NXs, NXf, NYs, NYf, NZs, NZf, neq)
            implicit none
            integer, intent(in)  :: NXs, NXf, NYs, NYf, NZs, NZf, neq
        end subroutine solve_eq_3D
    end interface

    interface solve_eq_0_3D
        module subroutine solve_eq_0_3d(NXs, NXf, NYs, NYf, NZs, NZf, neq)
            implicit none
            integer, intent(in)  :: NXs, NXf, NYs, NYf, NZs, NZf, neq
        end subroutine solve_eq_0_3D
    end interface

! pmbound.f90
    interface Bounds2D
        module subroutine Bounds2d(itype, NXs, NXf, NYs, NYf, neqs, neqf)
            implicit none
            integer, intent(in)  :: itype, NXs, NXf, NYs, NYf, neqs, neqf
        end subroutine Bounds2D
    end interface

    interface Bounds2d_lev
        module subroutine Bounds2d_lev(itype, NXs, NXf, NYs, NYf, neqs, neqf)
            implicit none
            integer, intent(in):: itype, NXs, NXf, NYs, NYf, neqs, neqf
        end subroutine Bounds2d_lev
    end interface

    interface Bounds3d
        module subroutine Bounds3d(itype, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
            implicit none
            integer, intent(in):: itype, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf
        end subroutine Bounds3d
    end interface

    interface Bounds3d_lev
        module subroutine Bounds3d_lev(itype, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
            implicit none
            integer, intent(in):: itype, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf
        end subroutine Bounds3d_lev
    end interface

contains
    !--------------------------------------------------------------------------------
    !>@function
    ! subroutine  pmesh
    !>
    !>@author Papis
    !>
    !>@brief
    !>This is the main subroutine of the Poisson solver library
    !>The input of the subroutine is DSOL_pm,DRHS_pm,DQP,DXP and d velocities.These variables are assigned the
    !>specified pointers.What is also needed is
    !> Dpm,NN,NN_bl,,Xbound which define the grid
    !> Dpm(3)   is DX,DY,DZ
    !> NN(3)    is the size of extend domain which coincides with the size of DSOL_pm,DRHS_pm
    !> NN_bl(6) is the size of the original grid (smaller domain) in which the poisson problem will be solved
    !>          (NN_bl(1:3),the starting nodes of the grid(X,Y,Z) with respect to the extended domain)
    !>          (NN_bl(4:6),the last nodes of the grid(X,Y,Z)  with respect to the extended domain)
    !>Xbound(6) Xmin,Ymin,Zmin,Xmax,Ymax,Zmax of the extended domain
    !>ibctyp (1 for bio savart law) (2 for infinite domain boundary conditions)
    !>neqs,neqf is an option if we want to solve more than one equation.neqs,neqf shoud coresspond to
    !>DSOL_pm(:,:,:,neqs,neqf)...
    !>iynbc is in case we want to add externally bc's for our equation.0 means normal solve (which means the
    !>solution at boundaries is from ibctyp.1 means keeps the solution at the boundaries at what is defined
    !> externally)
    !>IMPORTANT NOTE: The following library assumes that there are two domains.
    !>                -->the original domain and
    !>                -->an extended domain(that's why NN differ from NN_bl)
    !>                The solution is found on the extended domain.
    !>                IF you don't want an extended domain then set NN_bl(1)=1,NN_bl(4)=NN(1)..NN_bl(2)..
    !>                EXTERNALLY
    !REVISION HISTORY
    !> TODO_dd
    !>
    !>@param [in]
    !>@param [out]
    !>--------------------------------------------------------------------------------
    subroutine pmesh(DSOL_pm, DRHS_pm, DQP, DXP, Xbound, Dpm, NN, NN_bl, ND_in, ibctyp_in, &
                     neqs, neqf, iynbc, NVR_in, itree_in, levmax_in)
        ! use parvar, only : XP, QP
        use MPI
        implicit none
        real(dp), intent(inout), target  :: DSOL_pm(:, :, :, :), DRHS_pm(:, :, :, :), DQP(:, :), DXP(:, :)
        integer, intent(in)              :: ibctyp_in, itree_in, NVR_in, levmax_in, iynbc, neqs, neqf
        real(dp), intent(in)             :: Xbound(6), Dpm(3)
        integer, intent(in)              :: NN_bl(6), NN(3), ND_in
        ! Local variables
        integer                          ::  nb, NXs, NYs, NXf, NYf, NZs, NZf, neq

        ibctyp = ibctyp_in
        itree = itree_in
        NVR = NVR_in
        levmax = levmax_in
        ND = ND_in

        ! integer                        :: ierr, my_rank, np, rank
        ! real(dp)                       :: XPM, YPM, velx, vely
        ! real(dp)                       :: xi, yi, ksi1, ksi2, th1, th2, w1, w2
        ! real(dp)                       :: R, DX, DY, GreenF, nv

        !-->Pass the external data to the namespace of the library
        XMIN_pm = Xbound(1)
        YMIN_pm = Xbound(2)
        ZMIN_pm = Xbound(3)
        XMAX_pm = Xbound(4)
        YMAX_pm = Xbound(5)
        ZMAX_pm = Xbound(6)

        DXpm = Dpm(1)
        DYpm = Dpm(2)
        DZpm = Dpm(3)
        !NXs_bl,NXf_bl refer to the starting ending node of the domain we want to solve
        NXpm = NN(1)
        NYpm = NN(2)
        NZpm = NN(3)

        NXs_bl = NN_bl(1)
        NYs_bl = NN_bl(2)
        NZs_bl = NN_bl(3)
        NXf_bl = NN_bl(4)
        NYf_bl = NN_bl(5)
        NZf_bl = NN_bl(6)
        !Assign the pointers to the external data
        SOL_pm => DSOL_pm
        RHS_pm => DRHS_pm
        QP => DQP; 
        XP => DXP

        !iynbc 1 normal poisson solve.(Calculation of bc's is done here)
        if (iynbc .eq. 1) then
            !NXs,NXf is the points we want to calculate boundary conditions
            if (ND .eq. 2) then
                if (ibctyp .eq. 1) then
                    !Biot Savart boundary conditions(set on the ext. domain because that's what will be solved)
                    NXs = 1              ! NXs_bl
                    NXf = NXpm           !  NXf_bl
                    NYs = 1              !NYs_bl
                    NYf = NYpm           !NYf_bl
                    call Bounds2D(ibctyp, NXs, NXf, NYs, NYf, neqs, neqf)
                else
                    !Infinite domain boundary conditions(asume zero bc' at the boundary
                    call infdomain(neqs, neqf)
                end if
            else if (ND .eq. 3) then
                if (ibctyp .eq. 1) then
                    !Biot Savart boundary conditions(set on the ext. domain because that's what will be solved)
                    NXs = 1           !NXs_bl
                    NXf = NXpm        !NXf_bl
                    NYs = 1           !NYs_bl
                    NYf = NYpm        !NYf_bl
                    NZs = 1           !NZs_bl
                    NZf = NZpm        !NZf_bl
                    call Bounds3D(ibctyp, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
                else
                    !Infinite domain boundary conditions(asume zero bc' at the boundary
                    call infdomain_3D(neqs, neqf)
                end if
            end if

            !Boundary calculations have finished now we have the boundary conditions at NXs,NXf...
            if (ND .eq. 2) then
                nb = 1
                NXs = 1           !NXs_bl(1)
                NXf = NXpm        !NXf_bl
                NYs = 1           !NYs_bl
                NYf = NYpm        !NYf_bl
                !Solve poisson problem with the boundary conditions set at Nxs,NXf.....
                do neq = neqs, neqf
                    call solve_eq(NXs, NXf, NYs, NYf, neq)
                end do
            else
                nb = 1
                NXs = 1         !NXs_bl
                NXf = NXpm      !NXf_bl
                NYs = 1         !NYs_bl
                NYf = NYpm      !NYf_bl
                NZs = 1         !NZs_bl
                NZf = NZpm
                !Solve poisson problem with the boundary conditions set at Nxs,NXf.....
                do neq = neqs, neqf
                    call solve_eq_3D(NXs, NXf, NYs, NYf, NZs, NZf, neq)
                end do
            end if
            !In case iynbc = 0 then the boundary calculations are asummed ok from the external data and thus
            !going fo the final poisson solve.
            !IMPORTANT:: because the original idea was to have a solution at the smaller domain (extension of
            !the domain used for solving the small).The poisson solver is now solved for the small domain,since
            !boundary conditions are assumed ok.
        else
            if (ND .eq. 2) then
                nb = 1
                NXs = NXs_bl
                NXf = NXf_bl
                NYs = NYs_bl
                NYf = NYf_bl
                !Solve poisson problem with the boundary conditions set at Nxs,NXf.....
                do neq = neqs, neqf
                    call solve_eq(NXs, NXf, NYs, NYf, neq)
                end do
            else
                nb = 1
                NXs = NXs_bl
                NXf = NXf_bl
                NYs = NYs_bl
                NYf = NYf_bl
                NZs = NZs_bl
                NZf = NZf_bl
                !Solve poisson problem with the boundary conditions set at Nxs,NXf.....
                do neq = neqs, neqf
                    call solve_eq_3D(NXs, NXf, NYs, NYf, NZs, NZf, neq)
                end do
            end if
        end if

        nullify (SOL_pm, RHS_pm, QP, XP)
    end subroutine pmesh

    !> @brief Defines the characteristics of the Poisson solver grid.
    !>
    !> This subroutine extends the computational domain and adjusts its parameters
    !> based on the specified extension type.
    !>
    !> @param[in]    itype                 Type of domain extension (1, 2, or 3)
    !> @param[inout] Xbound                Domain boundaries [Xmin, Ymin, Zmin, Xmax, Ymax, Zmax]
    !> @param[inout] Dpm                   Grid spacings [DX, DY, DZ]
    !> @param[in]    num_dimensions        Number of dimensions (2 or 3)
    !> @param[in]    ndum                  Number of cells to extend the domain by
    !> @param[in]    nsize                 Desired divisibility of domain size [NX, NY, NZ]
    !> @param[out]   num_nodes             Number of nodes in extended domain [NX, NY, NZ]
    !> @param[out]   original_domain_idxs  Indexs of start and end nodes of original domain [IXs, JYs, KZs, IXe, JYe, KZe]
    !>
    !> @details
    !> Extension types:
    !> - Type 1: Extends the domain by ndum cells
    !> - Type 2: Extends the domain and adjusts Dpm to make it divisible by nsize
    !> - Type 3: Extends the domain by adding nodes to make it divisible by nsize, keeping original Dpm
    subroutine definepm(itype, Xbound, Dpm, num_dimensions, n_dum, nsize, num_nodes, original_domain_idxs)
        implicit None
        integer, intent(in)     :: itype, num_dimensions, nsize(3), n_dum
        real(dp), intent(inout) :: Xbound(6), Dpm(3)
        integer, intent(out)    :: num_nodes(3), original_domain_idxs(6)
        integer                 :: ndum_new(3), nn1, nn2, i, remainder

        ! Extend domain boundaries and calculate initial number of nodes
        do i = 1, num_dimensions
            Xbound(i) = Xbound(i) - (n_dum*Dpm(i))
            Xbound(i + 3) = Xbound(i + 3) + (n_dum*Dpm(i))
            num_nodes(i) = int(nint(abs(Xbound(i + 3) - Xbound(i))/Dpm(i))) + 1
        end do

        select case (itype)
            !Itype 1:
            ! extends the domain by ndum_new cells
        case (1)
            ndum_new = 0
            ! Already done
            !Itype 2:
            !extends the domain by ndum_new cells
            !Makes the domain divisible by nsize by changing the DPM
        case (2)
            ndum_new = 0
            do i = 1, num_dimensions
                num_nodes(i) = num_nodes(i) + nsize(i) - mod(num_nodes(i) - 1, nsize(i))
                Dpm(i) = (abs(Xbound(i + 3) - Xbound(i))/(num_nodes(i) - 1))
            end do
            !Itype 3:
            ! extends the domain by ndum_new cells
            ! Makes the domain divisible by nsize by adding new dummy nodes and keeps the original DPM
        case (3)
            do i = 1, num_dimensions
                ! The new number of nodes is increased so that it is divided by nsize
                remainder = mod(num_nodes(i) - 1, nsize(i))
                if ((remainder .eq. 0) .or. (remainder .eq. nsize(i))) then
                    ndum_new(i) = 0
                    cycle
                end if
                ! Add remainder to make nodes divisible by nsize
                ndum_new(i) = nsize(i) - remainder
                ! Here we add the new nodes symmetrically or asymmetrically to the domain
                nn1 = mod(ndum_new(i), 2)  ! Asymmetric addition part
                nn2 = int(ndum_new(i)/2)  ! Symmetric addition part
                Xbound(i) = Xbound(i) - ((nn1 + nn2)*Dpm(i))
                Xbound(i + 3) = Xbound(i + 3) + (nn2*Dpm(i))
                num_nodes(i) = int(nint(abs(Xbound(i + 3) - Xbound(i))/Dpm(i))) + 1
            end do
        end select

        !Define the nodes which the original domain starts (corresponding to Xbound_in)
        do i = 1, num_dimensions
            nn1 = mod(ndum_new(i), 2)
            if (nn1 .eq. 2) nn1 = 0
            nn2 = int(ndum_new(i)/2)
            original_domain_idxs(i) = n_dum + 1 ! + nn1 + nn2
            original_domain_idxs(i + 3) = num_nodes(i) - n_dum !- nn2
        end do

        ! In the 2D case, zero the 3d dimension
        if (num_dimensions .eq. 2) then
            num_nodes(3) = 1
            original_domain_idxs(3) = 1
            original_domain_idxs(6) = 1
            Xbound(6) = 0
            Xbound(3) = 0
        end if
    end subroutine definepm

end module pmlib

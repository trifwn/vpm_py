
submodule(pmlib) pmsolve_mudpack
   use mud3sp
   implicit none
contains
   !-----------------------------------------------------------------------!
   !-> subroutine solve_phiz                                                !
   !   This subroutines calls the fft library to solve for Phi poisson     !
   !   in all the points of Particle mesh.Dirichlet Boundary Cond. are used!
   !-----------------------------------------------------------------------!
   subroutine solve_eq_3d(NXs, NXf, NYs, NYf, NZs, NZf, neq)
        implicit none
        integer, intent(in)     :: NXs, NXf, NYs, NYf, NZs, NZf, neq
        logical                 :: has_bc

        has_bc = .true.
        call solve_mudpack_3D(NXs, NXf, NYs, NYf, NZs, NZf, neq, has_bc)
   end subroutine solve_eq_3d

    subroutine solve_eq_0_3d(NXs, NXf, NYs, NYf, NZs, NZf, neq)
        implicit none
        integer, intent(in)     :: NXs, NXf, NYs, NYf, NZs, NZf, neq
        logical                 :: has_bc

        has_bc = .false.
        call solve_mudpack_3D(NXs, NXf, NYs, NYf, NZs, NZf, neq, has_bc)

    end subroutine solve_eq_0_3d

    module subroutine solve_eq(NXs, NXf, NYs, NYf, neq)
        implicit none
        integer, intent(in)   :: NXs, NXf, NYs, NYf, neq

        write(*,*) 'solve_eq not implemented'
        stop
    end subroutine solve_eq

   module subroutine solve_eq_0(NXs, NXf, NYs, NYf, neq)
        implicit none
        integer, intent(in)           :: NXs, NXf, NYs, NYf, neq

        ! Raise NotImplementedError
        write(*,*) 'solve_eq_0 not implemented'
        stop
   end subroutine solve_eq_0

    subroutine solve_mudpack_3D(NXs, NXf, NYs, NYf, NZs, NZf, neq, has_bc)
        use console_io, only: print_stats_rank4, print_stats_rank3
        implicit none
        !--------------------------------------------------------------------
        ! Arguments:
        !   SOL_pm(3, nx, ny, nz)  : solution array (to be updated)
        !   RHS_pm(3, nx, ny, nz)  : right-hand side array
        !--------------------------------------------------------------------
        
        ! Get domain dimensions from fine_grid
        integer, intent(in)  :: NXs, NXf, NYs, NYf, NZs, NZf, neq
        logical, intent(in)  :: has_bc

        ! Local arrays for solution and RHS 
        real(dp), allocatable   :: sol_tmp(:,:,:)
        real(dp), allocatable   :: rhs(:,:,:)
        real(dp), allocatable   :: work(:)
        integer :: llwork
        ! Parameter arrays for Mudpack
        real(dp), dimension(8)  :: fprm
        integer, dimension(22)  :: iprm
        integer, dimension(4)   :: mgopt
        ! Local variables
        real(dp)    :: XMinCalc, XmaxCalc, YMinCalc, YmaxCalc, ZminCalc, ZmaxCalc, CX, CY, CZ
        integer     :: nx, ny, nz
        integer     :: ierror, system
        Integer     :: i, j, k
        integer     :: ixp, jyq, kzr, iex, jey, kez
        integer     :: temp, e

        XminCalc = XMIN_pm + (NXs - 1)*DXpm
        XmaxCalc = XMIN_pm + (NXf - 1)*DXpm

        YminCalc = YMIN_pm + (NYs - 1)*DYpm
        YmaxCalc = YMIN_pm + (NYf - 1)*DYpm

        ZminCalc = ZMIN_pm + (NZs - 1)*DZpm
        ZmaxCalc = ZMIN_pm + (NZf - 1)*DZpm

        NX = NXf - NXs + 1 - 2
        NY = NYf - NYs + 1 - 2
        NZ = NZf - NZs + 1 - 2

        ! Compute required workspace size (based on tmud3sp.f example)
        llwork = 3 * (7 * (nx+2) * (ny+2) * (nz+2)) / 2
        
        ! Initialize mgopt to zeros
        mgopt(1) = 2
        mgopt(2) = 2
        mgopt(3) = 1
        mgopt(4) = 3
        
        ! Set domain boundaries from fine_grid%Xbound
        ! Assumed ordering: [xa, xb, yc, yd, ze, zf]
        fprm(1)   = XMinCalc 
        fprm(2)   = XmaxCalc

        fprm(3)   = YMinCalc
        fprm(4)   = YmaxCalc
        
        fprm(5)   = ZminCalc
        fprm(6)   = ZmaxCalc

        ! Set up fprm parameters
        fprm(7)   = 1e-8     ! Tolerance
        fprm(8)   = 0.0      ! Output 
        
        ! Set up iprm parameters
        ! iprm(1) is used to control the call (0: initialization, 1: solve)
        iprm(2) = 1   ! Lower X boundary: Dirichlet
        iprm(3) = 1   ! Upper X boundary: Dirichlet
        iprm(4) = 1   ! Lower Y boundary: Dirichlet
        iprm(5) = 1   ! Upper Y boundary: Dirichlet
        iprm(6) = 1   ! Lower Z boundary: Dirichlet
        iprm(7) = 1   ! Upper Z boundary: Dirichlet
        
        ! Decompose nx: find odd part (ixp) and exponent+1 (iex)
            temp = nx - 1
            e = 0
            do while (temp /= 0 .and. mod(temp, 2) == 0)
                temp = temp / 2
                e = e + 1
            end do
            ixp = temp
            iex = e + 1

            ! Decompose ny: find odd part (jyq) and exponent+1 (jey)
            temp = ny - 1
            e = 0
            do while (temp /= 0 .and. mod(temp, 2) == 0)
                temp = temp / 2
                e = e + 1
            end do
            jyq = temp
            jey = e + 1

            ! Decompose nz: find odd part (kzr) and exponent+1 (kez)
            temp = nz - 1
            e = 0
            do while (temp /= 0 .and. mod(temp, 2) == 0)
                temp = temp / 2
                e = e + 1
            end do
            kzr = temp
            kez = e + 1

        iprm(8)  = ixp
        iprm(9)  = jyq
        iprm(10) = kzr
        iprm(11) = iex
        iprm(12) = jey
        iprm(13) = kez

        iprm(14) = ixp*(2**(iex-1)) + 1
        iprm(15) = jyq*(2**(jey-1)) + 1 
        iprm(16) = kzr*(2**(kez-1)) + 1
        iprm(17) = 1      ! Use the initial guess (phi initially zero)
        iprm(18) = 1000   ! Maximum number of iterations/cycles
        iprm(19) = 0      ! Relaxation method (e.g., Gauss-Seidel)
        iprm(20) = llwork ! Length of the work array
        iprm(21) = 0      ! (Output: minimum required workspace; set by mud3sp)
        iprm(22) = 0      ! (Output: iteration count; set by mud3sp)
                    
        ! Allocate workspace and local solution arrays
        allocate(sol_tmp(nx, ny, nz))
        allocate(rhs(nx, ny, nz))
        allocate(work(llwork))

        rhs    (1:NX, 1:NY, 1:NZ) = RHS_pm(neq, NXs + 1:NXf - 1, NYs + 1:NYf - 1, NZs + 1:NZf - 1)
        SOL_tmp(1:NX, 1:NY, 1:NZ) = RHS_pm(neq, NXs + 1:NXf - 1, NYs + 1:NYf - 1, NZs + 1:NZf - 1)

        if (has_bc) then
            CX = 1.d0!/DXpm**2
            CY = 1.d0!/DYpm**2
            CZ = 1.d0!/DZpm**2
            ! Set boundary conditions
            !-->Set Boundary Conditions
            !---> XMIN,XMAX
            do k = 1, NZ
                do j = 1, NY
                    SOL_tmp(1, j, k)  = SOL_tmp(1, j, k)  - CX*SOL_pm(neq, NXs, j + NYs - 1 + 1, k + NZs - 1 + 1)
                    SOL_tmp(NX, j, k) = SOL_tmp(NX, j, k) - CX*SOL_pm(neq, NXf, j + NYs - 1 + 1, k + NZs - 1 + 1)
                end do
            end do

            !---> YMIN,YMAX
            do k = 1, NZ
                do i = 1, NX
                    SOL_tmp(i, 1, k)  = SOL_tmp(i, 1, k)  - CY*SOL_pm(neq, i + NXs - 1 + 1, NYs, k + NZs - 1 + 1)
                    SOL_tmp(i, NY, k) = SOL_tmp(i, NY, k) - CY*SOL_pm(neq, i + NXs - 1 + 1, NYf, k + NZs - 1 + 1)
                end do
            end do
    
            do j = 1, NY
                do i = 1, NX
                    SOL_tmp(i, j, 1)  = SOL_tmp(i, j, 1)  - CZ*SOL_pm(neq, i + NXs - 1 + 1, j + NYs - 1 + 1, NZs)
                    SOL_tmp(i, j, NZ) = SOL_tmp(i, j, NZ) - CZ*SOL_pm(neq, i + NXs - 1 + 1, j + NYs - 1 + 1, NZf)
                end do
            end do
        end if

        !----------------------------------------------------------------
        ! First call: initialization/discretization step
        iprm(1) = 0
        call mud3sp(iprm, fprm, work, cfx, cfy, cfz, bndc, rhs, sol_tmp, mgopt, ierror)
        if (ierror /= 0) then
            write(*,*) 'Error during mud3sp initialization for system ', system, &
                        ' ierror = ', ierror
        endif
        
        !----------------------------------------------------------------
        ! Second call: iterative solve
        iprm(1) = 1
        call mud3sp(iprm, fprm, work, cfx, cfy, cfz, bndc, rhs, sol_tmp, mgopt, ierror)
        if (ierror /= 0) then
            write(*,*) 'Error during mud3sp solve for system ', system, ' ierror = ', ierror
        endif

        ! Store the computed solution back to SOL_pm
        if (has_bc) then
            SOL_pm(neq, NXs + 1:NXf - 1, NYs + 1:NYf - 1, NZs + 1:NZf - 1) = SOL_tmp(1:NX, 1:NY, 1:NZ)
        else
            SOL_0_pm(neq, NXs + 1:NXf - 1, NYs + 1:NYf - 1, NZs + 1:NZf - 1) = SOL_tmp(1:NX, 1:NY, 1:NZ)
        end if
        deallocate(work, sol_tmp, rhs)
    end subroutine solve_mudpack_3D

    !------------------------------------------------------------------
    ! Subroutine cfx(x,cxx,cx,cex) which provides the
    !     known real coefficients of the x derivative terms for the pde
    !     at any grid point x.  the name chosen in the calling routine
    !     may be different where the coefficient routine must be declared
    !     external.
    subroutine cfx(x, cxx, cx, cex)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: cxx, cx, cex
        ! For the laplacian operator, the coefficients are constant
        cxx = 1.0
        cx = 0.0
        cex = 0.0
    end subroutine cfx

    subroutine cfy(y, cyy, cy, cey)
        real(dp), intent(in) :: y
        real(dp), intent(out) :: cyy, cy, cey
        ! For the laplacian operator, the coefficients are constant
        cyy = 1.0
        cy = 0.0
        cey = 0.0
    end subroutine cfy

    subroutine cfz(z, czz, cz, cez)
        real(dp), intent(in) :: z
        real(dp), intent(out) :: czz, cz, cez
        ! For the laplacian operator, the coefficients are constant
        czz = 1.0
        cz = 0.0
        cez = 0.0
    end subroutine cfz

    subroutine bndc(kbdy, xory, yorz, alfa, phi_val)
        ! Subroutine bndyc(kbdy,xory,yorz,alfa,gbdy).
        !       which are used to input mixed boundary conditions to mud3sp.
        !       the boundaries are numbered one thru six and the form of
        !       conditions are described below.
        !                                          
        integer, intent(in) :: kbdy, xory, yorz, alfa
        real, intent(inout) :: phi_val
        ! For Dirichlet boundary conditions with fixed value 0, nothing is required.
        phi_val = 0.0_dp
    end subroutine bndc
end submodule pmsolve_mudpack
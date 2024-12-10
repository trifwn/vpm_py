submodule(pmlib) pmsolve_fish
    implicit none
contains
    !-----------------------------------------------------------------------!
    !-> subroutine solve_phiz                                                !
    !   This subroutines calls the fft library to solve for Phi poisson     !
    !   in all the points of Particle mesh.Dirichlet Boundary Cond. are used!
    !-----------------------------------------------------------------------!
    subroutine solve_eq_3d(NXs, NXf, NYs, NYf, NZs, NZf, neq)
        use fishpack
        implicit none
        integer, intent(in)  :: NXs, NXf, NYs, NYf, NZs, NZf, neq
        Integer              :: i, j, k, NWORK, INFO, NX, NY, NZ, nbj, LPEROD, MPEROD, NPEROD, IERROR
        real(dp)             :: XPM, YPM, XMinCalc, XmaxCalc, YMinCalc, YmaxCalc, ZminCalc, ZmaxCalc, CX, CY, CZ
        real(dp), allocatable::SOL_tmp(:, :, :), Apois(:), Bpois(:), Cpois(:)

        !--> Assignment of Boundary Values

        XminCalc = XMIN_pm + (NXs - 1)*DXpm
        XmaxCalc = XMIN_pm + (NXf - 1)*DXpm

        YminCalc = YMIN_pm + (NYs - 1)*DYpm
        YmaxCalc = YMIN_pm + (NYf - 1)*DYpm

        ZminCalc = ZMIN_pm + (NZs - 1)*DZpm
        ZmaxCalc = ZMIN_pm + (NZf - 1)*DZpm

        NX = NXf - NXs + 1 - 2
        NY = NYf - NYs + 1 - 2
        NZ = NZf - NZs + 1 - 2
        allocate (SOL_tmp(NX, NY, NZ))
        allocate (Apois(NZ), Bpois(NZ), Cpois(NZ))

        CX = 1.d0/DXpm**2
        CY = 1.d0/DYpm**2
        CZ = 1.d0/DZpm**2
        Apois(1:NZ) = CZ
        Bpois(1:NZ) = -2.d0*CZ
        Cpois(1:NZ) = CZ
        Apois(1) = 0.d0
        Cpois(NZ) = 0.d0
        !-->Set Right hand Side term (except for boundary conditions)
        SOL_tmp(1:NX, 1:NY, 1:NZ) = RHS_pm(neq, NXs + 1:NXf - 1, NYs + 1:NYf - 1, NZs + 1:NZf - 1)
        !-->Set Boundary Conditions
        !---> XMIN,XMAX

        !    Psiz_pm2(1:5,:)  = 0.d0
        !    Psiz_pm2(NX:NX-5,:) = 0.d0
        !    Psiz_pm2(:,1:5)  = 0.d0
        !    Psiz_pm2(:,NY:NY-5) = 0.d0
        do k = 1, NZ
            do j = 1, NY
                SOL_tmp(1, j, k) = SOL_tmp(1, j, k) - CX*SOL_pm(neq, NXs, j + NYs - 1 + 1, k + NZs - 1 + 1)
                SOL_tmp(NX, j, k) = SOL_tmp(NX, j, k) - CX*SOL_pm(neq, NXf, j + NYs - 1 + 1, k + NZs - 1 + 1)
            end do
        end do
        !---> YMIN,YMAX
        do k = 1, NZ
            do i = 1, NX
                SOL_tmp(i, 1, k) = SOL_tmp(i, 1, k) - CY*SOL_pm(neq, i + NXs - 1 + 1, NYs, k + NZs - 1 + 1)
                SOL_tmp(i, NY, k) = SOL_tmp(i, NY, k) - CY*SOL_pm(neq, i + NXs - 1 + 1, NYf, k + NZs - 1 + 1)
            end do
        end do

        do j = 1, NY
            do i = 1, NX
                SOL_tmp(i, j, 1) = SOL_tmp(i, j, 1) - CZ*SOL_pm(neq, i + NXs - 1 + 1, j + NYs - 1 + 1, NZs)
                SOL_tmp(i, j, NZ) = SOL_tmp(i, j, NZ) - CZ*SOL_pm(neq, i + NXs - 1 + 1, j + NYs - 1 + 1, NZf)
            end do
        end do

        call POIS3D(1, NX, CX, 1, NY, CY, 1, NZ, Apois, Bpois, Cpois, &
                    NX, NY, SOL_tmp, IERROR)
        SOL_pm(neq, NXs + 1:NXf - 1, NYs + 1:NYf - 1, NZs + 1:NZf - 1) = SOL_tmp(1:NX, 1:NY, 1:NZ)
        if (IERROR .ne. 0) then
            write (*, *) 'POISSON SOLVER ERROR', ierror
            STOP
        end if
    end subroutine solve_eq_3d

    subroutine solve_eq_0_3d(NXs, NXf, NYs, NYf, NZs, NZf, neq)
        use fishpack
        implicit none
        integer, intent(in)     :: NXs, NXf, NYs, NYf, NZs, NZf, neq
        Integer                 :: i, j, k, NWORK, INFO, NX, NY, NZ, nbj, LPEROD, MPEROD, NPEROD, IERROR
        real(dp)                :: XPM, YPM, XMinCalc, XmaxCalc, YMinCalc, YmaxCalc, ZminCalc, ZmaxCalc, CX, CY, CZ
        real(dp), allocatable   :: SOL_tmp(:, :, :), Apois(:), Bpois(:), Cpois(:)

        !--> Assignment of Boundary Values

        XminCalc = XMIN_pm + (NXs - 1)*DXpm
        XmaxCalc = XMIN_pm + (NXf - 1)*DXpm

        YminCalc = YMIN_pm + (NYs - 1)*DYpm
        YmaxCalc = YMIN_pm + (NYf - 1)*DYpm

        ZminCalc = ZMIN_pm + (NZs - 1)*DZpm
        ZmaxCalc = ZMIN_pm + (NZf - 1)*DZpm

        NX = NXf - NXs + 1 - 2
        NY = NYf - NYs + 1 - 2
        NZ = NZf - NZs + 1 - 2
        allocate (SOL_tmp(NX, NY, NZ))
        allocate (Apois(NZ), Bpois(NZ), Cpois(NZ))

        CX = 1.d0/DXpm**2
        CY = 1.d0/DYpm**2
        CZ = 1.d0/DZpm**2
        Apois(1:NZ) = CZ
        Bpois(1:NZ) = -2.d0*CZ
        Cpois(1:NZ) = CZ
        Apois(1) = 0.d0
        Cpois(NZ) = 0.d0
        !-->Set Right hand Side term (except for boundary conditions)

        SOL_tmp(1:NX, 1:NY, 1:NZ) = RHS_pm(neq, NXs + 1:NXf - 1, NYs + 1:NYf - 1, NZs + 1:NZf - 1)
        !-->Set Boundary Conditions
        !---> XMIN,XMAX

        call POIS3D(1, NX, CX, 1, NY, CY, 1, NZ, Apois, Bpois, Cpois, &
                    NX, NY, SOL_tmp, IERROR)

        SOL_0_pm(neq, NXs + 1:NXf - 1, NYs + 1:NYf - 1, NZs + 1:NZf - 1) = SOL_tmp(1:NX, 1:NY, 1:NZ)

        if (IERROR .ne. 0) then
            write (*, *) 'POISSON SOLVER ERROR', ierror
            STOP
        end if
    end subroutine solve_eq_0_3d

    module subroutine solve_eq(NXs, NXf, NYs, NYf, neq)
        use fishpack
        implicit none
        integer, intent(in)            :: NXs, NXf, NYs, NYf, neq
        Integer                       :: i, j, INFO, NX, NY
        integer                       :: NN
        double precision              :: XMinCalc, XmaxCalc, YMinCalc, YmaxCalc, pertrb
        double precision, allocatable  :: SOL_pm2(:, :)
        double precision, allocatable  :: bd_ax(:), bd_bx(:), bd_ay(:), bd_by(:)

        !--> Assignment of Boundary Values
        XminCalc = XMIN_pm + (NXs - 1)*DXpm
        XmaxCalc = XMIN_pm + (NXf - 1)*DXpm

        YminCalc = YMIN_pm + (NYs - 1)*DYpm
        YmaxCalc = YMIN_pm + (NYf - 1)*DYpm

        NX = NXf - NXs + 1
        NY = NYf - NYs + 1

        !-->Set Right hand Side term (except for boundary conditions)
        NN = NX*NY
        allocate (Sol_pm2(NX, NY))
        SOL_pm2 = -RHS_pm(neq, NXs + 1:NXf - 1, NYs + 1:NYf - 1, 1)

        !-->Set Boundary Conditions
        NN = NY
        allocate (bd_ax(NN), bd_bx(NN))
        do j = 1, NY
            bd_ax(j) = SOL_pm(neq, NXs, j + NYs - 1, 1)
            bd_bx(j) = SOL_pm(neq, NXf, j + NYs - 1, 1)
        end do
        !---> YMIN,YMAX
        NN = NX
        allocate (bd_ay(NN), bd_by(NN))

        do i = 1, NX
            bd_ay(i) = SOL_pm(neq, i + NXs - 1, NYs, 1)
            bd_by(i) = SOL_pm(neq, i + NXs - 1, NYf, 1)
        end do

        !--SOLVE PHI ON PMESH
        call hwscrt(XminCalc, XmaxCalc, NX - 1, 1, bd_ax, bd_bx, & ! Xmin,Xmax, Number of panels, Type of BC, BC at Xmin, BC at Xmax
                    YminCalc, YmaxCalc, NY - 1, 1, bd_ay, bd_by, & ! Same for Y
                    0.d0, & ! The constant LAMBDA in the Helmholtz equation.
                    SOL_pm2, & ! The right-hand side of the Helmholtz equation. F(I,J) = F(X(I),Y(J)).
                    NX, & ! The row (or first) dimension of the array F
                    pertrb, & ! PERTRB is a constant, calculated and subtracted from F, which ensures that a solution exists.
                    INFO)      ! An error flag that indicates invalid input parameters.

        SOL_pm(neq, NXs:NXf, NYs:NYf, 1) = SOL_pm2(1:NX, 1:NY)
        if (INFO .ne. 0) then
            write (*, *) 'POISSON SOLVER ERROR', INFO
            STOP
        end if
    end subroutine solve_eq

    module subroutine solve_eq_0(NXs, NXf, NYs, NYf, neq)
        use fishpack
        implicit none
        integer, intent(in)            :: NXs, NXf, NYs, NYf, neq
        Integer                       :: i, j, INFO, NX, NY
        integer                       :: NN
        double precision              :: XMinCalc, XmaxCalc, YMinCalc, YmaxCalc, pertrb
        double precision, allocatable  :: SOL_pm2(:, :)
        double precision, allocatable  :: bd_ax(:), bd_bx(:), bd_ay(:), bd_by(:)

        !--> Assignment of Boundary Values
        XminCalc = XMIN_pm + (NXs - 1)*DXpm
        XmaxCalc = XMIN_pm + (NXf - 1)*DXpm

        YminCalc = YMIN_pm + (NYs - 1)*DYpm
        YmaxCalc = YMIN_pm + (NYf - 1)*DYpm

        NX = NXf - NXs + 1
        NY = NYf - NYs + 1
        !-->Set Right hand Side term (except for boundary conditions)

        NN = NX*NY
        allocate (Sol_pm2(NX, NY))
        SOL_pm2 = -RHS_pm(neq, NXs + 1:NXf - 1, NYs + 1:NYf - 1, 1)

        !-->Set Boundary Conditions
        NN = NY
        allocate (bd_ax(NN), bd_bx(NN))
        bd_ax = 0.d0; bd_bx = 0.d0; !  SOL_pm(NXs,j + NYs -1,nb,neq)
        !---> YMIN,YMAX
        NN = NX
        allocate (bd_ay(NN), bd_by(NN))
        bd_ay = 0.d0; bd_by = 0.d0! SOL_pm(i + NXs-1,NYf ,nb,neq)

        ! Boundary conditions types:
        !    = 0  If the solution is periodic in X, i.e., U(I,J) = U(M+I,J).
        !    = 1  If the solution is specified at X = A and X = B.
        !    = 2  If the solution is specified at X = A and the derivative of
        !         the solution with respect to X is specified at X = B.
        !    = 3  If the derivative of the solution with respect to X is
        !         specified at X = A and X = B.
        !    = 4  If the derivative of the solution with respect to X is
        !         specified at X = A and the solution is specified at X = B.
        call hwscrt(XminCalc, XmaxCalc, NX - 1, 1, bd_ax, bd_bx, & ! Xmin,Xmax, Number of panels, Type of BC, BC at Xmin, BC at Xmax
                    YminCalc, YmaxCalc, NY - 1, 1, bd_ay, bd_by, & ! Same for Y
                    0.d0, & ! The constant LAMBDA in the Helmholtz equation.
                    SOL_pm2, & ! The right-hand side of the Helmholtz equation. F(I,J) = F(X(I),Y(J)).
                    NX, & ! The row (or first) dimension of the array F
                    pertrb, & ! PERTRB is a constant, calculated and subtracted from F, which ensures that a solution exists.
                    INFO)      ! An error flag that indicates invalid input parameters.
        if (INFO .ne. 0) then
            write (*, *) 'POISSON SOLVER ERROR', INFO
            STOP
        end if
        SOL_pm(neq, NXs:NXf, NYs:NYf, 1) = SOL_pm2(1:NX, 1:NY)

    end subroutine solve_eq_0
end submodule pmsolve_fish

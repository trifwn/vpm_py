submodule(pmlib) pmbound
    implicit none
contains
    !--------------------------------------------------------------------------------
    !>@file
    !>@brief Calculations of boundary conditions using particles or sources at domain
    !> boundaries
    !------------------------------------------------------------------------!
    !----------------------------------------------------------------------------!
    !-->subroutine Bounds2d                                                      !
    !   This subroutine calculates the boundary conditions for the Solver on PM  !
    !   The boundary conditions change for the 2d case                           !
    !   The equation solved div(grad)F = P  needs the exact values of F in PM    !
    !   boundaries.In 2d:                                                        !
    !   For one particle : F = P * (-lnr / (2pi)  )                              !
    !   The boundary values is the sum for all particles at each i,j             !
    !----------------------------------------------------------------------------!
    module subroutine Bounds2d(itype, NXs, NXf, NYs, NYf, neqs, neqf)
        implicit none
        integer, intent(in)  :: itype, NXs, NXf, NYs, NYf, neqs, neqf
        integer              :: iconst, jconst, iplane, i, j
        real(dp)             :: X, Y
        !-->Calculate boundary conditions for each boundary (XMIN,XMAX,YMIN,YMAX)

        !-->In case of infinite domain bc's(sources are used),In case
        !-->XMIN
        iconst = Nxs
        iplane = 2
        if (itype .eq. 1) then
            call calc_bound2d(iplane, iconst, NYs, NYf, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_2d(iplane, iconst, NYs, NYf, neqs, neqf)
        else if (itype .eq. 3) then
            call calc_boundinf_2d_s(iplane, iconst, NYs, NYf, neqs, neqf)
        end if
        !-->XMAX
        iconst = NXf
        iplane = 2
        if (itype .eq. 1) then
            call calc_bound2d(iplane, iconst, NYs, NYf, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_2d(iplane, iconst, NYs, NYf, neqs, neqf)
        else if (itype .eq. 3) then
            iplane = -iplane
            call calc_boundinf_2d_s(iplane, iconst, NYs, NYf, neqs, neqf)
        end if

        !-->YMIN
        !We use Nxs + 1,NXf - 1 For corners since they already calculated
        jconst = NYs
        iplane = 1
        if (itype .eq. 1) then
            call calc_bound2d(iplane, jconst, NXs + 1, NXf - 1, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_2d(iplane, jconst, NXs + 1, NXf - 1, neqs, neqf)
        else if (itype .eq. 3) then
            call calc_boundinf_2d_s(iplane, jconst, NXs + 1, NXf - 1, neqs, neqf)
        end if

        !-->YMAX
        jconst = NYf
        iplane = 1
        if (itype .eq. 1) then
            call calc_bound2d(iplane, jconst, NXs + 1, NXf - 1, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_2d(iplane, jconst, NXs + 1, NXf - 1, neqs, neqf)
        else if (itype .eq. 3) then
            iplane = -iplane
            call calc_boundinf_2d_s(iplane, jconst, NXs + 1, NXf - 1, neqs, neqf)
        end if

    end subroutine Bounds2d

    !>----------------------------------------------------------------------------!
    !>-->subroutine Bounds3d                                                      !
    !>   This subroutine calculates the boundary conditions for the Solver on PM  !
    !>   The boundary conditions change for the 3d case                           !
    !>   The equation solved div(grad)F = P  needs the exact values of F in PM    !
    !>   boundaries.In 2d:                                                        !
    !>   For one particle : F = P * (-lnr / (2pi)  )                              !
    !>   The boundary values is the sum for all particles at each i,j             !
    !>----------------------------------------------------------------------------!
    module subroutine Bounds3d(itype, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
        implicit none
        integer, intent(in):: itype, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf
        integer           :: iconst, jconst, kconst, iplane
        !-->Calculate boundary conditions for each boundary (XMIN,XMAX,YMIN,YMAX, ZMIN,ZMAX)
        !-->iplane is the plane of calculation of the bc's (i.e. for X=const a Y plane is defined)
        !-->iconst defines the poisition of the plane to be calculated
        !-->N*s,N*f is the nodes on the plane to be calculated
        !-->neqs,neqf is the bc's for more than one equations

        !-->XMIN
        iconst = Nxs
        iplane = 2
        if (itype .eq. 1) then
            call calc_bound3d(iplane, iconst, NYs, NYf, NZs, NZf, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_3d(iplane, iconst, NYs, NYf, NZs, NZf, neqs, neqf)
        else if (itype .eq. 3) then
            call calc_boundinf_3d_s(iplane, iconst, NYs, NYf, NZs, NZf, neqs, neqf)
        end if
        !-->XMAX
        iconst = NXf
        iplane = 2
        if (itype .eq. 1) then
            call calc_bound3d(iplane, iconst, NYs, NYf, NZs, NZf, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_3d(iplane, iconst, NYs, NYf, NZs, NZf, neqs, neqf)
        else if (itype .eq. 3) then
            iplane = -iplane
            call calc_boundinf_3d_s(iplane, iconst, NYs, NYf, NZs, NZf, neqs, neqf)
        end if

        !We use Nxs + 1,NXf - 1 For corners since they already calculated
        !-->YMIN
        jconst = NYs
        iplane = 1
        if (itype .eq. 1) then
            call calc_bound3d(iplane, jconst, NXs + 1, NXf - 1, NZs, NZf, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_3d(iplane, jconst, NXs + 1, NXf - 1, NZs, NZf, neqs, neqf)
        else if (itype .eq. 3) then
            call calc_boundinf_3d_s(iplane, jconst, NXs + 1, NXf - 1, NZs, NZf, neqs, neqf)
        end if
        !-->YMAX
        jconst = NYf
        if (itype .eq. 1) then
            call calc_bound3d(iplane, jconst, NXs + 1, NXf - 1, NZs, NZf, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_3d(iplane, jconst, NXs + 1, NXf - 1, NZs, NZf, neqs, neqf)
        else if (itype .eq. 3) then
            iplane = -iplane
            call calc_boundinf_3d_s(iplane, jconst, NXs + 1, NXf - 1, NZs, NZf, neqs, neqf)
        end if

        !-->ZMIN
        kconst = NZs
        iplane = 3
        if (itype .eq. 1) then
            call calc_bound3d(iplane, kconst, NXs + 1, NXf - 1, NYs + 1, NYf - 1, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_3d(iplane, kconst, NXs + 1, NXf - 1, NYs + 1, NYf - 1, neqs, neqf)
        else if (itype .eq. 3) then
            call calc_boundinf_3d_s(iplane, kconst, NXs + 1, NXf - 1, NYs + 1, NYf - 1, neqs, neqf)
        end if
        !-->ZMAX
        kconst = NZf
        if (itype .eq. 1) then
            call calc_bound3d(iplane, kconst, NXs + 1, NXf - 1, NYs + 1, NYf - 1, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_3d(iplane, kconst, NXs + 1, NXf - 1, NYs + 1, NYf - 1, neqs, neqf)
        else if (itype .eq. 3) then
            iplane = -iplane
            call calc_boundinf_3d_s(iplane, kconst, NXs + 1, NXf - 1, NYs + 1, NYf - 1, neqs, neqf)
        end if
    end subroutine Bounds3d

    !-----------------------------------------------------------------------------!
    !->subroutine calc_bound2d                                                    !
    !  This subroutine calculates Biot-Savart Law for 2D applications             !
    !  Input:                                                                     !
    !        iplane : axis of plane                                               !
    !        iconst : Calculation side in nodal value                             !
    !-----------------------------------------------------------------------------!
    subroutine calc_bound2d(iplane, iconst, Ns, Nf, neqs, neqf)
        implicit none

        integer, intent(in)  :: iplane, iconst, Ns, Nf, neqs, neqf
        real(dp)             :: X, Y, XR, YR, r
        integer              :: i, j, nv

        !calculate bc's of all particles on the specified plane defined at iconst
        !-->Y=constant plane
        if (iplane .eq. 1) then
            Y = YMIN_pm + (iconst - 1)*DYpm
            do nv = 1, NVR
                YR = XP(2, nv) - Y

                do i = Ns, Nf

                    X = XMIN_pm + (i - 1)*DXpm
                    XR = XP(1, nv) - X
                    r = sqrt(XR**2 + YR**2)
                    !-->Remember Bio Savart :  logr/2pr for grad phi = delta
                    SOL_pm(neqs:neqf, i, iconst, 1) = SOL_pm(neqs:neqf, i, iconst, 1) + QP(neqs:neqf, nv)*log(r)/pi2
                end do
            end do
            !-->X=constant plane
        else if (iplane .eq. 2) then
            X = XMIN_pm + (iconst - 1)*DXpm

            do nv = 1, NVR
                XR = XP(1, nv) - X

                do j = Ns, Nf !Because corner is calculated TWICE

                    Y = YMIN_pm + (j - 1)*DYpm
                    YR = XP(2, nv) - Y
                    r = sqrt(XR**2 + YR**2)
                    !-->Remember Bio Savart :  logr/2pi
                    SOL_pm(neqs:neqf, iconst, j, 1) = SOL_pm(neqs:neqf, iconst, j, 1) + QP(neqs:neqf, nv)*log(r)/pi2
                end do
            end do

        end if

    end subroutine calc_bound2d

    !-----------------------------------------------------------------------------!
    !->subroutine calc_bound3d                                                    !
    !  This subroutine calculates Biot-Savart Law for 3D applications             !
    !  Input:                                                                     !
    !        iplane : axis of plane                                               !
    !        iconst : Calculation side in nodal value                             !
    !-----------------------------------------------------------------------------!
    subroutine calc_bound3d(iplane, iconst, Ns, Nf, Ns2, Nf2, neqs, neqf)
        implicit none

        integer, intent(in) :: iplane, iconst, Ns, Nf, Ns2, Nf2, neqs, neqf

        real(dp)            :: X, Y, XR, YR, Z, ZR, r
        integer             :: i, j, k, nv

        !calculate bc's of all particles on the specified plane defined at iconst
        !-->Y=constant plane
        if (iplane .eq. 1) then
            Y = YMIN_pm + (iconst - 1)*DYpm

            do nv = 1, NVR
                YR = XP(2, nv) - Y
                do k = Ns2, Nf2
                    do i = Ns, Nf
                        X = XMIN_pm + (i - 1)*DXpm
                        XR = XP(1, nv) - X
                        Z = ZMIN_pm + (k - 1)*DZpm
                        ZR = XP(3, nv) - Z
                        r = sqrt(XR**2 + YR**2 + ZR**2)

                        !-->Remember Bio Savart : - 1/(PI4*R) for gradphi= delta
                        SOL_pm(neqs:neqf, i, iconst, k) = SOL_pm(neqs:neqf, i, iconst, k) - &
                                                          QP(neqs:neqf, nv)/(PI4*r)
                    end do
                end do
            end do
            !-->X=constant plane
        else if (iplane .eq. 2) then
            X = XMIN_pm + (iconst - 1)*DXpm

            do nv = 1, NVR
                XR = XP(1, nv) - X
                do k = Ns2, Nf2
                    do j = Ns, Nf !Because corner is calculated TWICE

                        Y = YMIN_pm + (j - 1)*DYpm
                        YR = XP(2, nv) - Y
                        Z = ZMIN_pm + (k - 1)*DZpm
                        ZR = XP(3, nv) - Z
                        r = sqrt(XR**2 + YR**2 + ZR**2)

                        SOL_pm(neqs:neqf, iconst, j, k) = SOL_pm(neqs:neqf, iconst, j, k) - &
                                                          QP(neqs:neqf, nv)/(PI4*r)
                    end do
                end do
            end do
            !-->Z=constant plane
        else if (iplane .eq. 3) then
            Z = ZMIN_pm + (iconst - 1)*DZpm

            do nv = 1, NVR
                ZR = XP(3, nv) - Z
                do j = Ns2, Nf2 !Because corner is calculated TWICE
                    do i = Ns, Nf

                        X = XMIN_pm + (i - 1)*DXpm
                        XR = XP(1, nv) - X
                        Y = YMIN_pm + (j - 1)*DYpm
                        YR = XP(2, nv) - Y
                        r = sqrt(XR**2 + YR**2 + ZR**2)

                        SOL_pm(neqs:neqf, i, j, iconst) = SOL_pm(neqs:neqf, i, j, iconst) - QP(neqs:neqf, nv)/(PI4*r)
                    end do
                end do
            end do
        end if

    end subroutine calc_bound3d

    !-------------------------------------------------------------------------------!
    !-> subroutine calc_boundinf                                                    !
    !   This subroutine calculates boundary conditions for the sources              !
    !   Same as particles
    !-------------------------------------------------------------------------------!
    subroutine calc_boundinf_2d(iplane, iconst, Ns, Nf, neqs, neqf)
        implicit none

        integer, intent(in) :: iplane, iconst, Ns, Nf, neqs, neqf

        real(dp)            :: X, Y, XR, YR, r, DS
        integer             :: i, j, nv

        !calculate bc's of all sources on the specified plane defined at iconst
        !-->Y=constant plane
        if (abs(iplane) .eq. 1) then
            Y = YMIN_pm + (iconst - 1)*DYpm
            do nv = 1, nbound
                YR = 0.5*(y_s(1, nv) + y_s(2, nv)) - Y

                do i = Ns, Nf
                    X = XMIN_pm + (i - 1)*DXpm
                    XR = 0.5d0*(x_s(1, nv) + x_s(2, nv)) - X

                    r = sqrt(XR**2 + YR**2)
                    DS = d_s(nv)

                    !-->Remember Bio Savart :  logr/2pr for grad phi = delta
                    SOL_pm(neqs:neqf, i, iconst, 1) = SOL_pm(neqs:neqf, i, iconst, 1) + &
                                                      source_bound(neqs:neqf, nv)*log(r)*DS/pi2

                end do
            end do
            !-->X=constant plane
        else if (abs(iplane) .eq. 2) then
            X = XMIN_pm + (iconst - 1)*DXpm

            do nv = 1, nbound
                XR = 0.5*(x_s(1, nv) + x_s(2, nv)) - X

                do j = Ns, Nf !Because corners are calculated twice
                    Y = YMIN_pm + (j - 1)*DYpm
                    YR = 0.5*(y_s(1, nv) + y_s(2, nv)) - Y

                    r = sqrt(XR**2 + YR**2)
                    DS = d_s(nv)

                    SOL_pm(neqs:neqf, iconst, j, 1) = SOL_pm(neqs:neqf, iconst, j, 1) + &
                                                      source_bound(neqs:neqf, nv)*log(r)*DS/pi2

                end do
            end do

        end if

    end subroutine calc_boundinf_2d

    !-------------------------------------------------------------------------------!
    !-> subroutine calc_boundinf                                                    !
    !   This subroutine calculates boundary conditions for the sources              !
    !   Same as particles
    !-------------------------------------------------------------------------------!
    subroutine calc_boundinf_2d_s(iplane, iconst, Ns, Nf, neqs, neqf)
        implicit none

        integer, intent(in) :: iplane, iconst, Ns, Nf, neqs, neqf

        real(dp)            :: X, Y, XR, YR, greenint, cosb, sinb, DS
        integer             :: i, j, nv

        !calculate bc's of all sources on the specified plane defined at iconst
        !-->Y=constant plane

        if (abs(iplane) .eq. 1) then
            Y = YMIN_pm + (iconst - 1)*DYpm
            cosb = -1.d0*sign(1, iplane)
            sinb = 0.d0
            do nv = 1, nbound
                XR = x_s(1, nv)
                YR = y_s(1, nv)
                do i = Ns, Nf
                    X = XMIN_pm + (i - 1)*DXpm

                    DS = d_s(nv)

                    call PHIELS(X, Y, XR, YR, DS, cosb, sinb, greenint)
                    SOL_pm(neqs:neqf, i, iconst, 1) = SOL_pm(neqs:neqf, i, iconst, 1) + &
                                                      source_bound(neqs:neqf, nv)*greenint
                end do
            end do
            !-->X=constant plane
        else if (abs(iplane) .eq. 2) then
            X = XMIN_pm + (iconst - 1)*DXpm
            cosb = 0.d0
            sinb = 1.d0*sign(1, iplane)

            do nv = 1, nbound
                XR = x_s(1, nv)
                YR = y_s(1, nv)
                do j = Ns, Nf !Because corners are calculated twice
                    Y = YMIN_pm + (j - 1)*DYpm
                    DS = d_s(nv)

                    call PHIELS(X, Y, XR, YR, DS, cosb, sinb, greenint)
                    SOL_pm(neqs:neqf, iconst, j, 1) = SOL_pm(neqs:neqf, iconst, j, 1) + &
                                                      source_bound(neqs:neqf, nv)*greenint

                end do
            end do

        end if

    end subroutine calc_boundinf_2d_s

    !-------------------------------------------------------------------------------!
    !-> subroutine calc_boundinf                                                    !
    !   This subroutine calculates boundary conditions for the sources              !
    !-------------------------------------------------------------------------------!
    subroutine calc_boundinf_3d(iplane, iconst, Ns, Nf, Ns2, Nf2, neqs, neqf)
        implicit none

        integer, intent(in) :: iplane, iconst, Ns, Nf, Ns2, Nf2, neqs, neqf

        real(dp)            :: X, Y, XR, YR, Z, ZR, r, a, b, ra, rb, greenint, racos, rasin, DS
        integer             :: i, j, k, nv
        !-->Y=constant plane
        if (abs(iplane) .eq. 1) then
            Y = YMIN_pm + (iconst - 1)*DYpm

            !calculate bc's of all sources on the specified plane defined at iconst
            do nv = 1, nbound
                YR = 0.25*(y_s(1, nv) + y_s(2, nv) + y_s(3, nv) + y_s(4, nv)) - Y
                do k = Ns2, Nf2
                    do i = Ns, Nf

                        Z = ZMIN_pm + (k - 1)*DZpm
                        ZR = 0.25d0*(z_s(1, nv) + z_s(2, nv) + z_s(3, nv) + z_s(4, nv)) - Z

                        X = XMIN_pm + (i - 1)*DXpm
                        XR = 0.25d0*(x_s(1, nv) + x_s(2, nv) + x_s(3, nv) + x_s(4, nv)) - X

                        r = sqrt(XR**2 + YR**2 + ZR**2)
                        DS = d_s(nv)
                        !Green function -1/(4PIR)
                        if (r .gt. 1d-05) then
                            SOL_pm(neqs:neqf, i, iconst, k) = SOL_pm(neqs:neqf, i, iconst, k) - &
                                                              source_bound(neqs:neqf, nv)*DS/(PI4*r)
                        else

                            SOL_pm(neqs:neqf, i, iconst, k) = SOL_pm(neqs:neqf, i, iconst, k) + &
                                                              source_bound(neqs:neqf, nv)*(DXpm + DZpm)*log((sqrt(2.d0) - 1)/ &
                                                                                                            (sqrt(2.d0) + 1))/PI4
                        end if

                    end do
                end do
            end do
            !-->X=constant plane
        else if (abs(iplane) .eq. 2) then
            X = XMIN_pm + (iconst - 1)*DXpm

            do nv = 1, nbound
                XR = 0.25d0*(x_s(1, nv) + x_s(2, nv) + x_s(3, nv) + x_s(4, nv)) - X
                !XR = x_s(1,nv) -  X

                do k = Ns2, Nf2
                    do j = Ns, Nf !Because corners are calculated twice
                        Z = ZMIN_pm + (k - 1)*DZpm
                        ZR = 0.25d0*(z_s(1, nv) + z_s(2, nv) + z_s(3, nv) + z_s(4, nv)) - Z
                        !ZR = z_s(1,nv) -  Z
                        Y = YMIN_pm + (j - 1)*DYpm
                        YR = 0.25*(y_s(1, nv) + y_s(2, nv) + y_s(3, nv) + y_s(4, nv)) - Y
                        !YR =y_s(1,nv)- Y

                        r = sqrt(XR**2 + YR**2 + ZR**2)
                        DS = d_s(nv)

                        if (r .gt. 1d-05) then
                            SOL_pm(neqs:neqf, iconst, j, k) = SOL_pm(neqs:neqf, iconst, j, k) - &
                                                              source_bound(neqs:neqf, nv)*DS/(PI4*r)
                        else
                            SOL_pm(neqs:neqf, iconst, j, k) = SOL_pm(neqs:neqf, iconst, j, k) + &
                                                              source_bound(neqs:neqf, nv)*(DYpm + DZpm)*log((sqrt(2.d0) - 1)/ &
                                                                                                            (sqrt(2.d0) + 1))/PI4
                        end if

                    end do
                end do
            end do

        else if (abs(iplane) .eq. 3) then
            Z = ZMIN_pm + (iconst - 1)*DZpm

            do nv = 1, nbound
                ZR = 0.25d0*(z_s(1, nv) + z_s(2, nv) + z_s(3, nv) + z_s(4, nv)) - Z
                do j = Ns2, Nf2 !Because corners are calculated twice
                    do i = Ns, Nf
                        X = XMIN_pm + (i - 1)*DXpm
                        XR = 0.25d0*(x_s(1, nv) + x_s(2, nv) + x_s(3, nv) + x_s(4, nv)) - X
                        Y = YMIN_pm + (j - 1)*DYpm
                        YR = 0.25*(y_s(1, nv) + y_s(2, nv) + y_s(3, nv) + y_s(4, nv)) - Y

                        r = sqrt(XR**2 + YR**2 + ZR**2)
                        DS = d_s(nv)

                        if (r .gt. 1d-05) then
                            SOL_pm(neqs:neqf, i, j, iconst) = SOL_pm(neqs:neqf, i, j, iconst) - &
                                                              source_bound(neqs:neqf, nv)*DS/(PI4*r)
                        else
                            SOL_pm(neqs:neqf, i, j, iconst) = SOL_pm(neqs:neqf, i, j, iconst) + &
                                                              source_bound(neqs:neqf, nv)*(DXpm + DYpm)*log((sqrt(2.d0) - 1)/ &
                                                                                                            (sqrt(2.d0) + 1))/PI4
                        end if

                    end do
                end do
            end do

        end if

    end subroutine calc_boundinf_3d

    !-------------------------------------------------------------------------------!
    !-> subroutine calc_boundinf                                                    !
    !   This subroutine calculates boundary conditions for the sources              !
    !-------------------------------------------------------------------------------!
    subroutine calc_boundinf_3d_s(iplane, iconst, Ns, Nf, Ns2, Nf2, neqs, neqf)
        implicit none

        integer, intent(in) :: iplane, iconst, Ns, Nf, Ns2, Nf2, neqs, neqf
        real(dp)            :: X, Y, XR, YR, Z, ZR, r, a, b, ra, rb, greenint, racos, rasin, DS
        real(dp)            :: XO(3), RG(3), E1(3), E2(3), E3(3), S(4), T(4), SINB(4), COSB(4), D(4), &
                               AREA, DIAG, EPSS, FIS
        integer             :: i, j, k, nv
        integer             :: ISING, NSIDE, si

        ISING = 0
        NSIDE = 0
        EPSS = 1d-14
        !-->Y=constant plane
        if (abs(iplane) .eq. 1) then
            Y = YMIN_pm + (iconst - 1)*DYpm
            si = sign(1, iplane)
            DIAG = sqrt(DXpm**2 + DZpm**2)
            E1 = 0.d0; E1(1) = 1.d0*si
            E2 = 0.d0; E2(3) = 1.d0
            E3 = 0.d0; E3(2) = -1.d0*si
            COSB = 0.d0; COSB(2) = 1.d0; COSB(4) = -1.d0
            SINB = 0.d0; SINB(1) = 1.d0; SINB(3) = -1.d0
            AREA = DXpm*DZpm
            S(1) = -0.5*DXpm; S(2) = -0.5*DXpm; S(3) = 0.5d0*DXpm; S(4) = 0.5d0*DXpm
            T(1) = -0.5*DZpm; T(2) = 0.5*DZpm; T(3) = 0.5d0*DZpm; T(4) = -0.5d0*DZpm
            D(1) = DZpm; D(2) = DXpm; D(3) = DZpm; D(4) = DXpm
            !calculate bc's of all sources on the specified plane defined at iconst
            do nv = 1, nbound
                RG(1) = 0.25d0*(x_s(1, nv) + x_s(2, nv) + x_s(3, nv) + x_s(4, nv))
                RG(2) = 0.25d0*(y_s(1, nv) + y_s(2, nv) + y_s(3, nv) + y_s(4, nv))
                RG(3) = 0.25d0*(z_s(1, nv) + z_s(2, nv) + z_s(3, nv) + z_s(4, nv))
                do k = Ns2, Nf2
                    do i = Ns, Nf
                        X = XMIN_pm + (i - 1)*DXpm
                        Z = ZMIN_pm + (k - 1)*DZpm
                        XO(1) = X; XO(2) = Y; XO(3) = Z
                        call FSOUR_A4(XO, RG, E1, E2, E3, &
                                      S, T, D, SINB, COSB, &
                                      DIAG, AREA, NSIDE, EPSS, ISING, FIS)

                        !Green function -1/(4PIR)
                        SOL_pm(neqs:neqf, i, iconst, k) = SOL_pm(neqs:neqf, i, iconst, k) + &
                                                          FIS*source_bound(neqs:neqf, nv)
                    end do
                end do
            end do
            !-->X=constant plane
        else if (abs(iplane) .eq. 2) then
            X = XMIN_pm + (iconst - 1)*DXpm
            si = sign(1, iplane)
            DIAG = sqrt(DYpm**2 + DZpm**2)
            E1 = 0.d0; E1(3) = 1.d0
            E2 = 0.d0; E2(2) = -1.d0*si
            E3 = 0.d0; E3(1) = -1.d0*si
            COSB = 0.d0; COSB(2) = 1.d0; COSB(4) = -1.d0
            SINB = 0.d0; SINB(1) = 1.d0; SINB(3) = -1.d0
            AREA = DYpm*DZpm
            S(1) = -0.5*DYpm; S(2) = -0.5*DYpm; S(3) = 0.5d0*DYpm; S(4) = 0.5d0*DYpm
            T(1) = -0.5*DZpm; T(2) = 0.5*DZpm; T(3) = 0.5d0*DZpm; T(4) = -0.5d0*DZpm
            D(1) = DZpm; D(2) = DYpm; D(3) = DZpm; D(4) = DYpm
            !calculate bc's of all sources on the specified plane defined at iconst
            do nv = 1, nbound
                RG(1) = 0.25d0*(x_s(1, nv) + x_s(2, nv) + x_s(3, nv) + x_s(4, nv))
                RG(2) = 0.25d0*(y_s(1, nv) + y_s(2, nv) + y_s(3, nv) + y_s(4, nv))
                RG(3) = 0.25d0*(z_s(1, nv) + z_s(2, nv) + z_s(3, nv) + z_s(4, nv))
                do k = Ns2, Nf2
                    do j = Ns, Nf !Because corners are calculated twice
                        Z = ZMIN_pm + (k - 1)*DZpm
                        Y = YMIN_pm + (j - 1)*DYpm
                        XO(1) = X; XO(2) = Y; XO(3) = Z
                        !Green function -1/(4PIR)
                        call FSOUR_A4(XO, RG, E1, E2, E3, &
                                      S, T, D, SINB, COSB, &
                                      DIAG, AREA, NSIDE, EPSS, ISING, FIS)

                        SOL_pm(neqs:neqf, iconst, j, k) = SOL_pm(neqs:neqf, iconst, j, k) + &
                                                          FIS*source_bound(neqs:neqf, nv)

                    end do
                end do
            end do
            !-->Z=constant plane
        else if (abs(iplane) .eq. 3) then
            Z = ZMIN_pm + (iconst - 1)*DZpm
            si = sign(1, iplane)
            DIAG = sqrt(DXpm**2 + DYpm**2)
            E1 = 0.d0; E1(1) = 1.d0
            E2 = 0.d0; E2(2) = -1.d0*si
            E3 = 0.d0; E3(3) = -1.d0*si
            COSB = 0.d0; COSB(2) = 1.d0; COSB(4) = -1.d0
            SINB = 0.d0; SINB(1) = 1.d0; SINB(3) = -1.d0
            AREA = DXpm*DYpm
            S(1) = -0.5*DXpm; S(2) = -0.5*DXpm; S(3) = 0.5d0*DXpm; S(4) = 0.5d0*DXpm
            T(1) = -0.5*DYpm; T(2) = 0.5*DYpm; T(3) = 0.5d0*DYpm; T(4) = -0.5d0*DYpm
            D(1) = DYpm; D(2) = DXpm; D(3) = DYpm; D(4) = DXpm
            !calculate bc's of all sources on the specified plane defined at iconst
            do nv = 1, nbound
                RG(1) = 0.25d0*(x_s(1, nv) + x_s(2, nv) + x_s(3, nv) + x_s(4, nv))
                RG(2) = 0.25d0*(y_s(1, nv) + y_s(2, nv) + y_s(3, nv) + y_s(4, nv))
                RG(3) = 0.25d0*(z_s(1, nv) + z_s(2, nv) + z_s(3, nv) + z_s(4, nv))
                do j = Ns2, Nf2 !Because corners are calculated twice
                    do i = Ns, Nf
                        X = XMIN_pm + (i - 1)*DXpm
                        Y = YMIN_pm + (j - 1)*DYpm
                        XO(1) = X; XO(2) = Y; XO(3) = Z
                        !Green function -1/(4PIR)
                        call FSOUR_A4(XO, RG, E1, E2, E3, &
                                      S, T, D, SINB, COSB, &
                                      DIAG, AREA, NSIDE, EPSS, ISING, FIS)

                        SOL_pm(neqs:neqf, i, j, iconst) = SOL_pm(neqs:neqf, i, j, iconst) + &
                                                          FIS*source_bound(neqs:neqf, nv)

                    end do
                end do
            end do

        end if

    end subroutine calc_boundinf_3d_s

    module subroutine Bounds2d_lev(itype, NXs, NXf, NYs, NYf, neqs, neqf)
        implicit none
        integer, intent(in):: itype, NXs, NXf, NYs, NYf, neqs, neqf
        integer           :: iconst, jconst, iplane, i, j
        real(dp)          :: X, Y
        !-->Calculate boundary conditions for each boundary (XMIN,XMAX,YMIN,YMAX)

        !-->In case of infinite domain bc's(sources are used),In case
        !-->XMIN
        iconst = Nxs
        iplane = 2
        if (itype .eq. 1) then
            call calc_bound2d(iplane, iconst, NYs, NYf, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_2d_lev(iplane, iconst, NYs, NYf, neqs, neqf)
        else if (itype .eq. 3) then
            call calc_boundinf_2d_lev_s(iplane, iconst, NYs, NYf, neqs, neqf)
        end if
        !-->XMAX
        iconst = NXf
        iplane = 2
        if (itype .eq. 1) then
            call calc_bound2d(iplane, iconst, NYs, NYf, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_2d_lev(iplane, iconst, NYs, NYf, neqs, neqf)
            !call calc_boundinf_2d_lev(iplane,iconst,NYs,NYf,neqs,neqf)
        else if (itype .eq. 3) then
            iplane = -iplane
            call calc_boundinf_2d_lev_s(iplane, iconst, NYs, NYf, neqs, neqf)
            !call calc_boundinf_2d_lev(iplane,iconst,NYs,NYf,neqs,neqf)
        end if

        !-->YMIN
        !We use Nxs + 1,NXf - 1 For corners since they already calculated
        jconst = NYs
        iplane = 1
        if (itype .eq. 1) then
            call calc_bound2d(iplane, jconst, NXs + 1, NXf - 1, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_2d_lev(iplane, jconst, NXs + 1, NXf - 1, neqs, neqf)
        else if (itype .eq. 3) then
            call calc_boundinf_2d_lev_s(iplane, jconst, NXs + 1, NXf - 1, neqs, neqf)
        end if

        !-->YMAX
        jconst = NYf
        iplane = 1
        if (itype .eq. 1) then
            call calc_bound2d(iplane, jconst, NXs + 1, NXf - 1, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_2d_lev(iplane, jconst, NXs + 1, NXf - 1, neqs, neqf)
        else if (itype .eq. 3) then
            iplane = -iplane
            call calc_boundinf_2d_lev_s(iplane, jconst, NXs + 1, NXf - 1, neqs, neqf)
        end if

    end subroutine Bounds2d_lev

    !-------------------------------------------------------------------------------!
    !-> subroutine calc_boundinf                                                    !
    !   This subroutine calculates boundary conditions for the sources              !
    !   Same as particles
    !-------------------------------------------------------------------------------!
    subroutine calc_boundinf_2d_lev(iplane, iconst, Ns, Nf, neqs, neqf)
        implicit none

        integer, intent(in) :: iplane, iconst, Ns, Nf, neqs, neqf

        real(dp)            :: X, Y, XR, YR, r, a, b, ra, rb, greenint, racos, rasin, DS, SOURCE(neqf)
        integer             :: i, j, nv
        integer             :: leafstart, leaffin, lev, nlev, nleaf, branch
        !calculate bc's of all sources on the specified plane defined at iconst
        !-->Y=constant plane
        if (abs(iplane) .eq. 1) then
            Y = YMIN_pm + (iconst - 1)*DYpm
            ! do nv = 1, nbound_lev(levmax)
            nv = 1
            do i = Ns, Nf
                X = XMIN_pm + (i - 1)*DXpm
                SOURCE = 0.d0
                leafstart = 0
                branch = 1
                call tree_calc_2d(nv, levmax, leafstart, X, Y, SOURCE, neqs, neqf)
                SOL_pm(neqs:neqf, i, iconst, 1) = SOL_pm(neqs:neqf, i, iconst, 1) + SOURCE(neqs:neqf)/pi2
            end do
            ! enddo
            !-->X=constant plane
        else if (abs(iplane) .eq. 2) then
            X = XMIN_pm + (iconst - 1)*DXpm
            ! do nv = 1, nbound_lev(levmax)
            nv = 1
            do j = Ns, Nf !Because corners are calculated twice
                Y = YMIN_pm + (j - 1)*DYpm
                SOURCE = 0.d0
                leafstart = 0
                branch = 1
                call tree_calc_2d(nv, levmax, leafstart, X, Y, SOURCE, neqs, neqf)
                SOL_pm(neqs:neqf, iconst, j, 1) = SOL_pm(neqs:neqf, iconst, j, 1) + SOURCE(neqs:neqf)/pi2

            end do
            ! enddo

        end if

    end subroutine calc_boundinf_2d_lev

    Recursive subroutine tree_calc_2d(nv, nlev, leafstart, X, Y, SOURCE, neqs, neqf)
        implicit none
        integer, intent(in) :: nlev, nv, neqs, neqf
        integer, intent(inout) :: leafstart
        real(dp), intent(in)  :: X, Y
        real(dp), intent(inout) :: SOURCE(neqf)

        integer                         :: newlev, nleaf, leaffin, leafs, leaff
        integer                         :: listleaf(4), nlf, nn, nlist, nj
        integer                         :: nmax, npre
        real(dp)                        :: YR, XR, r, DS
        integer                         ::ierr, my_rank
        !The loop for all bounds happens here
        !lev4 is the coarsest one.We start searching the coarsest and then move in finer and finer levels
        !The tree is a set of structured grids so we use that information to go deeper and deeper in the tree

        !for the coarsest level do all elements(The loop for all bounds happens here so we have to do for all coarsest
        !level elements
        if (nlev .eq. levmax) then
            leafstart = 1
            leaffin = nbound_lev(nlev)
            nlist = 0

            !in case of dummy element
        else if (ilev_t(nv, nlev + 1, 6) .lt. 0) then
            listleaf(1:1) = ilev_t(nv, nlev + 1, 1:1)
            !find the 4 finer cells.To do that we use the structured grid information
            leafstart = 1; leaffin = 1
            nlist = 1
            !in case of normal element
        else
            listleaf(1:2) = ilev_t(nv, nlev + 1, 1:2)
            leafstart = 1; leaffin = 2
            nlist = 1
        end if

        newlev = nlev - 1

        !if (my_rank.eq.0) then
        !    write(*,*) leafstart,nlev,nlist,nv
        !endif
        !if (my_rank.eq.0) then
        !if (nlev.lt.levmax.and.nlist.eq.1.and.ilev_t(nv,nlev+1,6).gt.0) then
        !  do nlf= leafstart,leaffin
        !     nleaf=listleaf(nlf)
        !   print *,nv,nlev,nleaf
        !     write(17,*) xs_lev(nleaf,nlev),&
        !                 ys_lev(nleaf,nlev)

        !  enddo
        !    write(17,*)  xs_lev(listleaf(1),nlev),&
        !                 ys_lev(listleaf(1),nlev)
        !    write(17,*)
        !    write(17,*)
        !    write(18,*)  xs_lev(nv,nlev+1),&
        !                 ys_lev(nv,nlev+1)

        !else if (nlev.lt.levmax.and.nlist.eq.1.and.ilev_t(nv,nlev+1,6).lt.0) then
        !  do nlf= leafstart,leaffin
        !     nleaf=listleaf(nlf)
        !     print *,nv,nlev,nleaf
        !     write(16,*) xs_lev(nleaf,nlev),&
        !                 ys_lev(nleaf,nlev)

        !  enddo
        !     write(16,*)
        !     write(16,*)
        !endif
        !read(*,*)
        !endif
        ! write(*,*) nlev,leafstart,leaffin
        ! read(*,*)
        do nlf = leafstart, leaffin
            if (nlist .eq. 0) then
                nleaf = nlf
            else
                nleaf = listleaf(nlf)
            end if
            if (ilev_t(nleaf, nlev, 6) .eq. 0) then
                write (*, *) 'Something is wrong with the tree'
                write (*, *) nleaf, nlev
                write (*, *)
                STOP
            end if
            YR = ys_lev(nleaf, nlev) - Y
            XR = xs_lev(nleaf, nlev) - X
            r = sqrt(XR**2 + YR**2)
            DS = ds_lev(nleaf, nlev)
            if ((r .lt. 10*DS .and. nlev .gt. 0) .or. ilev_t(nleaf, nlev, 6) .lt. 0) then
                call tree_calc_2d(nleaf, newlev, leafstart, X, Y, SOURCE, neqs, neqf)
            else
                SOURCE(neqs:neqf) = SOURCE(neqs:neqf) + source_bound_lev(nleaf, neqs:neqf, nlev)*DS*log(r)
            end if
        end do
    end subroutine tree_calc_2d

    !-------------------------------------------------------------------------------!
    !-> subroutine calc_boundinf                                                    !
    !   This subroutine calculates boundary conditions for the sources              !
    !   Same as particles
    !-------------------------------------------------------------------------------!

    subroutine calc_boundinf_2d_lev_s(iplane, iconst, Ns, Nf, neqs, neqf)
        implicit none

        integer, intent(in) :: iplane, iconst, Ns, Nf, neqs, neqf

        real(dp)            :: X, Y, XR, YR, r, a, b, ra, rb, greenint, cosb, sinb, DS, SOURCE(neqf)
        integer             :: i, j, nv
        integer             :: leafstart, leaffin, lev, nlev, nleaf, branch
        !calculate bc's of all sources on the specified plane defined at iconst
        !-->Y=constant plane
        if (abs(iplane) .eq. 1) then
            Y = YMIN_pm + (iconst - 1)*DYpm
            cosb = -1.d0*sign(1, iplane)
            sinb = 0.d0
            nv = 1
            do i = Ns, Nf
                X = XMIN_pm + (i - 1)*DXpm
                SOURCE = 0.d0
                leafstart = 0
                branch = 1
                call tree_calc_2d_s(nv, levmax, leafstart, X, Y, cosb, sinb, SOURCE, neqs, neqf)
                SOL_pm(neqs:neqf, i, iconst, 1) = SOL_pm(neqs:neqf, i, iconst, 1) + SOURCE(neqs:neqf)
            end do
            !-->X=constant plane
        else if (abs(iplane) .eq. 2) then
            X = XMIN_pm + (iconst - 1)*DXpm
            cosb = 0.d0
            sinb = 1.d0*sign(1, iplane)
            nv = 1
            do j = Ns, Nf !Because corners are calculated twice
                Y = YMIN_pm + (j - 1)*DYpm
                SOURCE = 0.d0
                leafstart = 0
                branch = 1
                call tree_calc_2d_s(nv, levmax, leafstart, X, Y, cosb, sinb, SOURCE, neqs, neqf)
                SOL_pm(neqs:neqf, iconst, j, 1) = SOL_pm(neqs:neqf, iconst, j, 1) + SOURCE(neqs:neqf)

            end do

        end if

    end subroutine calc_boundinf_2d_lev_s

    Recursive subroutine tree_calc_2d_s(nv, nlev, leafstart, X, Y, cosb, sinb, SOURCE, neqs, neqf)
        implicit none
        integer, intent(in)     :: nlev, nv, neqs, neqf
        integer, intent(inout)  :: leafstart
        real(dp), intent(in)    :: X, Y, cosb, sinb
        real(dp), intent(inout) :: SOURCE(neqf)

        integer                         :: newlev, nleaf, leaffin, leafs, leaff
        integer                         :: listleaf(4), nlf, nn, nlist, nj
        integer                         :: nmax, npre
        real(dp)                        :: YR, XR, r, DS, greenint

        if (nlev .eq. levmax) then
            leafstart = 1
            leaffin = nbound_lev(nlev)
            nlist = 0

            !in case of dummy element
        else if (ilev_t(nv, nlev + 1, 6) .lt. 0) then
            listleaf(1:1) = ilev_t(nv, nlev + 1, 1:1)
            !find the 4 finer cells.To do that we use the structured grid information
            leafstart = 1; leaffin = 1
            nlist = 1
            !in case of normal element
        else
            listleaf(1:2) = ilev_t(nv, nlev + 1, 1:2)
            leafstart = 1; leaffin = 2
            nlist = 1
        end if
        newlev = nlev - 1

        do nlf = leafstart, leaffin
            if (nlist .eq. 0) then
                nleaf = nlf
            else
                nleaf = listleaf(nlf)
            end if
            if (ilev_t(nleaf, nlev, 6) .eq. 0) then
                write (*, *) 'Something is wrong with the tree'
                write (*, *) nleaf, nlev
                write (*, *)
                STOP
            end if

            DS = ds_lev(nleaf, nlev)
            YR = ys_lev(nleaf, nlev)
            XR = xs_lev(nleaf, nlev)
            r = sqrt((XR - X)**2 + (YR - Y)**2)
            if ((r .lt. 10*sqrt(DS) .and. nlev .gt. 0) .or. ilev_t(nleaf, nlev, 6) .eq. -1) then
                call tree_calc_2d_s(nleaf, newlev, leafstart, X, Y, cosb, sinb, SOURCE, neqs, neqf)
            else
                if (nlev .gt. 0) then
                    SOURCE(neqs:neqf) = SOURCE(neqs:neqf) + source_bound_lev(nleaf, neqs:neqf, nlev) &
                                        *DS*log(r)/PI2
                else
                    call PHIELS(X, Y, XR, YR, DS, cosb, sinb, greenint)
                    SOURCE(neqs:neqf) = SOURCE(neqs:neqf) + source_bound_lev(nleaf, neqs:neqf, nlev) &
                                        *greenint
                end if
            end if
        end do

    end subroutine tree_calc_2d_s

    !----------------------------------------------------------------------------!
    !-->subroutine Bounds3d                                                      !
    !   This subroutine calculates the boundary conditions for the Solver on PM  !
    !   The boundary conditions change for the 3d case                           !
    !   The equation solved div(grad)F = P  needs the exact values of F in PM    !
    !   boundaries.In 2d:                                                        !
    !   For one particle : F = P * (-lnr / (2pi)  )                                 !
    !   The boundary values is the sum for all particles at each i,j             !
    !----------------------------------------------------------------------------!
    module subroutine Bounds3d_lev(itype, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf)
        implicit none
        integer, intent(in):: itype, NXs, NXf, NYs, NYf, NZs, NZf, neqs, neqf
        integer           :: iconst, jconst, kconst, iplane
        !-->Calculate boundary conditions for each boundary (XMIN,XMAX,YMIN,YMAX)
        !-->iplane is the plane of calculation of the bc's (i.e. for X=const a Y plane is defined)
        !-->iconst defines the poisition of the plane to be calculated
        !-->N*s,N*f is the nodes on the plane to be calculated
        !-->neqs,neqf is the bc's for more than one equations

        !-->XMIN
        iconst = Nxs
        iplane = 2
        if (itype .eq. 1) then
            call calc_bound3d(iplane, iconst, NYs, NYf, NZs, NZf, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_3d_lev(iplane, iconst, NYs, NYf, NZs, NZf, neqs, neqf)
        else if (itype .eq. 3) then
            call calc_boundinf_3d_lev_s(iplane, iconst, NYs, NYf, NZs, NZf, neqs, neqf)
        end if
        !-->XMAX
        iconst = NXf
        if (itype .eq. 1) then
            call calc_bound3d(iplane, iconst, NYs, NYf, NZs, NZf, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_3d_lev(iplane, iconst, NYs, NYf, NZs, NZf, neqs, neqf)
        else if (itype .eq. 3) then
            iplane = -iplane
            call calc_boundinf_3d_lev_s(iplane, iconst, NYs, NYf, NZs, NZf, neqs, neqf)
        end if

        !We use Nxs + 1,NXf - 1 For corners since they already calculated
        !-->YMIN
        jconst = NYs
        iplane = 1
        if (itype .eq. 1) then
            call calc_bound3d(iplane, jconst, NXs + 1, NXf - 1, NZs, NZf, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_3d_lev(iplane, jconst, NXs + 1, NXf - 1, NZs, NZf, neqs, neqf)
        else if (itype .eq. 3) then
            call calc_boundinf_3d_lev_s(iplane, jconst, NXs + 1, NXf - 1, NZs, NZf, neqs, neqf)
        end if
        !-->YMAX
        jconst = NYf
        if (itype .eq. 1) then
            call calc_bound3d(iplane, jconst, NXs + 1, NXf - 1, NZs, NZf, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_3d_lev(iplane, jconst, NXs + 1, NXf - 1, NZs, NZf, neqs, neqf)
        else if (itype .eq. 3) then
            iplane = -iplane
            call calc_boundinf_3d_lev_s(iplane, jconst, NXs + 1, NXf - 1, NZs, NZf, neqs, neqf)
        end if

        !-->ZMIN
        kconst = NZs
        iplane = 3
        if (itype .eq. 1) then
            call calc_bound3d(iplane, kconst, NXs + 1, NXf - 1, NYs + 1, NYf - 1, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_3d_lev(iplane, kconst, NXs + 1, NXf - 1, NYs + 1, NYf - 1, neqs, neqf)
        else if (itype .eq. 3) then
            call calc_boundinf_3d_lev_s(iplane, kconst, NXs + 1, NXf - 1, NYs + 1, NYf - 1, neqs, neqf)
        end if
        !-->ZMAX
        kconst = NZf
        if (itype .eq. 1) then
            call calc_bound3d(iplane, kconst, NXs + 1, NXf - 1, NYs + 1, NYf - 1, neqs, neqf)
        else if (itype .eq. 2) then
            call calc_boundinf_3d_lev(iplane, kconst, NXs + 1, NXf - 1, NYs + 1, NYf - 1, neqs, neqf)
        else if (itype .eq. 3) then
            iplane = -iplane
            call calc_boundinf_3d_lev_s(iplane, kconst, NXs + 1, NXf - 1, NYs + 1, NYf - 1, neqs, neqf)
        end if

    end subroutine Bounds3d_lev

    subroutine calc_boundinf_3d_lev(iplane, iconst, Ns, Nf, Ns2, Nf2, neqs, neqf)
        implicit none

        integer, intent(in) :: iplane, iconst, Ns, Nf, Ns2, Nf2, neqs, neqf

        real(dp)            :: X, Y, XR, YR, Z, ZR, r, a, b, ra, rb, greenint, racos, rasin, DS, SOURCE(1:neqf)
        integer             :: i, j, k, nv
        integer             :: leafstart, branch
        !-->Y=constant plane
        if (abs(iplane) .eq. 1) then
            Y = YMIN_pm + (iconst - 1)*DYpm
            !calculate bc's of all sources on the specified plane defined at iconst
            ! do nv = 1, nbound_lev(0)
            nv = 1
            do k = Ns2, Nf2
                do i = Ns, Nf
                    Z = ZMIN_pm + (k - 1)*DZpm
                    X = XMIN_pm + (i - 1)*DXpm

                    SOURCE = 0.d0
                    leafstart = 0
                    branch = 1

                    call tree_calc_3d(nv, levmax, leafstart, X, Y, Z, SOURCE, neqs, neqf)
                    SOL_pm(neqs:neqf, i, iconst, k) = SOL_pm(neqs:neqf, i, iconst, k) + SOURCE(neqs:neqf)/PI4
                end do
            end do
            !  enddo
            !-->X=constant plane
        else if (abs(iplane) .eq. 2) then
            X = XMIN_pm + (iconst - 1)*DXpm

            !  do nv = 1, nbound_lev(0)
            nv = 1
            do k = Ns2, Nf2
                do j = Ns, Nf !Because corners are calculated twice
                    Z = ZMIN_pm + (k - 1)*DZpm
                    Y = YMIN_pm + (j - 1)*DYpm

                    SOURCE = 0.d0
                    leafstart = 0
                    branch = 1

                    call tree_calc_3d(nv, levmax, leafstart, X, Y, Z, SOURCE, neqs, neqf)
                    SOL_pm(neqs:neqf, iconst, j, k) = SOL_pm(neqs:neqf, iconst, j, k) + SOURCE(neqs:neqf)/PI4

                end do
            end do
            !  enddo
            !-->Z=constant plane
        else if (abs(iplane) .eq. 3) then
            Z = ZMIN_pm + (iconst - 1)*DZpm

            !do nv = 1, nbound_lev(0)
            nv = 1
            do j = Ns2, Nf2 !Because corners are calculated twice
                do i = Ns, Nf
                    X = XMIN_pm + (i - 1)*DXpm
                    Y = YMIN_pm + (j - 1)*DYpm

                    SOURCE = 0.d0
                    leafstart = 0
                    branch = 1
                    call tree_calc_3d(nv, levmax, leafstart, X, Y, Z, SOURCE, neqs, neqf)
                    SOL_pm(neqs:neqf, i, j, iconst) = SOL_pm(neqs:neqf, i, j, iconst) + SOURCE(neqs:neqf)/PI4

                end do
            end do
            !enddo

        end if

    end subroutine calc_boundinf_3d_lev

    Recursive subroutine tree_calc_3d(nv, nlev, leafstart, X, Y, Z, SOURCE, neqs, neqf)
        use MPI
        implicit none
        integer, intent(in) :: nlev, nv, neqs, neqf
        integer, intent(inout) :: leafstart
        real(dp), intent(in)  :: X, Y, Z
        real(dp), intent(inout) :: SOURCE(neqf)

        integer                         :: newlev, nleaf, leaffin, leafs, leaff, listleaf(4), nlf, nn, nlist, nj
        integer                         :: nmax, npre
        real(dp)                        :: XR, YR, ZR, r, DS, ss(neqf)
        integer                         :: my_rank, ierr

        call MPI_Comm_Rank(MPI_COMM_WORLD, my_rank, ierr)
        listleaf = 0
        !The loop for all bounds happens here
        !lev4 is the coarsest one.We start searching the coarsest and then move in finer and finer levels
        !The tree is a set of structured grids so we use that information to go deeper and deeper in the tree

        !for the coarsest level do all elements(The loop for all bounds happens here so we have to do for all coarsest
        !level elements
        if (nlev .eq. levmax) then
            leafstart = 1
            leaffin = nbound_lev(nlev)
            nlist = 0

            !in case of dummy element
        else if (ilev_t(nv, nlev + 1, 6) .lt. 0) then
            listleaf(1:2) = ilev_t(nv, nlev + 1, 1:2)
            !find the 4 finer cells.To do that we use the structured grid information
            leafstart = 1; leaffin = 2
            nlist = 1
            !in case of normal element
        else
            listleaf(1:4) = ilev_t(nv, nlev + 1, 1:4)
            leafstart = 1; leaffin = 4
            nlist = 1
        end if

        newlev = nlev - 1

        !if (my_rank.eq.0) then
        !    write(*,*) leafstart,nlev,nlist,nv
        !endif
        !if (my_rank.eq.0) then
        !if (nlev.lt.levmax.and.nlist.eq.1.and.ilev_t(nv,nlev+1,6).gt.0) then
        !  do nlf= leafstart,leaffin
        !     nleaf=listleaf(nlf)
        !   print *,nv,nlev,nleaf
        !     write(17,*) xs_lev(nleaf,nlev),&
        !                 ys_lev(nleaf,nlev),&
        !                 zs_lev(nleaf,nlev)

        !  enddo
        !    write(17,*)  xs_lev(listleaf(1),nlev),&
        !                 ys_lev(listleaf(1),nlev),&
        !                 zs_lev(listleaf(1),nlev)
        !    write(17,*)
        !    write(17,*)
        !    write(18,*)  xs_lev(nv,nlev+1),&
        !                 ys_lev(nv,nlev+1),&
        !                 zs_lev(nv,nlev+1)

        !else if (nlev.lt.levmax.and.nlist.eq.1.and.ilev_t(nv,nlev+1,6).lt.0) then
        !  do nlf= leafstart,leaffin
        !     nleaf=listleaf(nlf)
        !     print *,nv,nlev,nleaf
        !     write(16,*) xs_lev(nleaf,nlev),&
        !                 ys_lev(nleaf,nlev),&
        !                 zs_lev(nleaf,nlev)

        !  enddo
        !     write(16,*)
        !     write(16,*)
        !endif
        !read(*,*)
        !endif

        !search all elements in coarser grid and go to finer if needed.Initialy at coarsest level
        !we search all elements
        !at a finer level we search the 4
        do nlf = leafstart, leaffin
            if (nlist .eq. 0) then
                nleaf = nlf
            else
                nleaf = listleaf(nlf)
            end if
            if (ilev_t(nleaf, nlev, 6) .eq. 0) then
                write (*, *) 'Something is wrong with the tree'
                write (*, *) nleaf, nlev
                write (*, *)
                STOP
            end if
            XR = xs_lev(nleaf, nlev) - X
            YR = ys_lev(nleaf, nlev) - Y
            ZR = zs_lev(nleaf, nlev) - Z
            r = sqrt(XR**2 + YR**2 + ZR**2)
            DS = ds_lev(nleaf, nlev)
            if ((r .lt. 10*sqrt(DS) .and. nlev .gt. 0) .or. ilev_t(nleaf, nlev, 6) .lt. 0) then
                call tree_calc_3d(nleaf, newlev, leafstart, X, Y, Z, SOURCE, neqs, neqf)
            else
                !Green function  -1/(PI4*R)
                if (r .gt. 1d-05) then
                    SOURCE(neqs:neqf) = SOURCE(neqs:neqf) - source_bound_lev(nleaf, neqs:neqf, nlev)*DS/r
                else
                    SOURCE(neqs:neqf) = SOURCE(neqs:neqf) + &
                                        source_bound_lev(nleaf, neqs:neqf, nlev)*2*sqrt(DS)*log((sqrt(2.d0) - 1)/(sqrt(2.d0) + 1))
                end if

            end if
        end do

    end subroutine tree_calc_3d

    !-------------------------------------------------------------------------------!
    !-> subroutine calc_boundinf                                                    !
    !   This subroutine calculates boundary conditions for the sources              !
    !-------------------------------------------------------------------------------!
    subroutine calc_boundinf_3d_lev_s(iplane, iconst, Ns, Nf, Ns2, Nf2, neqs, neqf)
        implicit none

        integer, intent(in) :: iplane, iconst, Ns, Nf, Ns2, Nf2, neqs, neqf

        real(dp)            :: X, Y, XR, YR, Z, ZR, r, a, b, ra, rb, greenint, racos, rasin, DS, SOURCE(1:neqf)
        integer             :: i, j, k, nv
        integer             :: leafstart, branch

        real(dp)            :: XO(3), RG(3), E1(3), E2(3), E3(3), S(4), T(4), SINB(4), COSB(4), D(4), &
                               AREA, DIAG, EPSS, FIS
        integer             :: ISING, NSIDE, si

        !-->Y=constant plane
        ISING = 0
        NSIDE = 0
        EPSS = 1d-14
        !-->Y=constant plane
        if (abs(iplane) .eq. 1) then
            Y = YMIN_pm + (iconst - 1)*DYpm
            si = sign(1, iplane)
            DIAG = sqrt(DXpm**2 + DZpm**2)
            E1 = 0.d0; E1(1) = 1.d0*si
            E2 = 0.d0; E2(3) = 1.d0
            E3 = 0.d0; E3(2) = -1.d0*si
            COSB = 0.d0; COSB(2) = 1.d0; COSB(4) = -1.d0
            SINB = 0.d0; SINB(1) = 1.d0; SINB(3) = -1.d0
            AREA = DXpm*DZpm
            S(1) = -0.5*DXpm; S(2) = -0.5*DXpm; S(3) = 0.5d0*DXpm; S(4) = 0.5d0*DXpm
            T(1) = -0.5*DZpm; T(2) = 0.5*DZpm; T(3) = 0.5d0*DZpm; T(4) = -0.5d0*DZpm
            D(1) = DZpm; D(2) = DXpm; D(3) = DZpm; D(4) = DXpm

            !calculate bc's of all sources on the specified plane defined at iconst
            nv = 1
            do k = Ns2, Nf2
                do i = Ns, Nf
                    Z = ZMIN_pm + (k - 1)*DZpm
                    X = XMIN_pm + (i - 1)*DXpm

                    XO(1) = X; XO(2) = Y; XO(3) = Z

                    SOURCE = 0.d0
                    leafstart = 0
                    branch = 1

                    call tree_calc_3d_s(nv, levmax, leafstart, XO, SOURCE, neqs, neqf, &
                                        DIAG, E1, E2, E3, COSB, SINB, S, T, D, NSIDE, EPSS, ISING)
                    SOL_pm(neqs:neqf, i, iconst, k) = SOL_pm(neqs:neqf, i, iconst, k) + SOURCE(neqs:neqf)
                end do
            end do
            !-->X=constant plane
        else if (abs(iplane) .eq. 2) then
            X = XMIN_pm + (iconst - 1)*DXpm
            si = sign(1, iplane)
            DIAG = sqrt(DYpm**2 + DZpm**2)
            E1 = 0.d0; E1(3) = 1.d0
            E2 = 0.d0; E2(2) = -1.d0*si
            E3 = 0.d0; E3(1) = -1.d0*si
            COSB = 0.d0; COSB(2) = 1.d0; COSB(4) = -1.d0
            SINB = 0.d0; SINB(1) = 1.d0; SINB(3) = -1.d0
            AREA = DYpm*DZpm
            S(1) = -0.5*DYpm; S(2) = -0.5*DYpm; S(3) = 0.5d0*DYpm; S(4) = 0.5d0*DYpm
            T(1) = -0.5*DZpm; T(2) = 0.5*DZpm; T(3) = 0.5d0*DZpm; T(4) = -0.5d0*DZpm
            D(1) = DZpm; D(2) = DYpm; D(3) = DZpm; D(4) = DYpm

            nv = 1
            do k = Ns2, Nf2
                do j = Ns, Nf !Because corners are calculated twice
                    Z = ZMIN_pm + (k - 1)*DZpm
                    Y = YMIN_pm + (j - 1)*DYpm

                    XO(1) = X; XO(2) = Y; XO(3) = Z

                    SOURCE = 0.d0
                    leafstart = 0
                    branch = 1

                    call tree_calc_3d_s(nv, levmax, leafstart, XO, SOURCE, neqs, neqf, &
                                        DIAG, E1, E2, E3, COSB, SINB, S, T, D, NSIDE, EPSS, ISING)
                    SOL_pm(neqs:neqf, iconst, j, k) = SOL_pm(neqs:neqf, iconst, j, k) + SOURCE(neqs:neqf)

                end do
            end do

        else if (abs(iplane) .eq. 3) then
            Z = ZMIN_pm + (iconst - 1)*DZpm
            Z = ZMIN_pm + (iconst - 1)*DZpm
            si = sign(1, iplane)
            DIAG = sqrt(DXpm**2 + DYpm**2)
            E1 = 0.d0; E1(1) = 1.d0
            E2 = 0.d0; E2(2) = -1.d0*si
            E3 = 0.d0; E3(3) = -1.d0*si
            COSB = 0.d0; COSB(2) = 1.d0; COSB(4) = -1.d0
            SINB = 0.d0; SINB(1) = 1.d0; SINB(3) = -1.d0
            AREA = DXpm*DYpm
            S(1) = -0.5*DXpm; S(2) = -0.5*DXpm; S(3) = 0.5d0*DXpm; S(4) = 0.5d0*DXpm
            T(1) = -0.5*DYpm; T(2) = 0.5*DYpm; T(3) = 0.5d0*DYpm; T(4) = -0.5d0*DYpm
            D(1) = DYpm; D(2) = DXpm; D(3) = DYpm; D(4) = DXpm

            nv = 1
            do j = Ns2, Nf2 !Because corners are calculated twice
                do i = Ns, Nf
                    X = XMIN_pm + (i - 1)*DXpm
                    Y = YMIN_pm + (j - 1)*DYpm

                    XO(1) = X; XO(2) = Y; XO(3) = Z

                    SOURCE = 0.d0
                    leafstart = 0
                    branch = 1
                    call tree_calc_3d_s(nv, levmax, leafstart, XO, SOURCE, neqs, neqf, &
                                        DIAG, E1, E2, E3, COSB, SINB, S, T, D, NSIDE, EPSS, ISING)

                    SOL_pm(neqs:neqf, i, j, iconst) = SOL_pm(neqs:neqf, i, j, iconst) + SOURCE(neqs:neqf)

                end do
            end do

        end if

    end subroutine calc_boundinf_3d_lev_s

    Recursive subroutine tree_calc_3d_s(nv, nlev, leafstart, XO, SOURCE, neqs, neqf, &
                                        DIAG, E1, E2, E3, COSB, SINB, S, T, D, NSIDE, EPSS, ISING)
        implicit none
        integer, intent(in)     :: nlev, nv, neqs, neqf
        integer, intent(inout)  :: leafstart
        real(dp), intent(in)    :: DIAG, E1(3), E2(3), E3(3)
        real(dp), intent(in)    :: COSB(4), SINB(4), S(4), T(4), D(4)
        real(dp), intent(in)    :: EPSS
        real(dp), intent(in)    :: XO(3)
        real(dp), intent(inout) :: SOURCE(neqf)
        integer, intent(in)     :: ISING, NSIDE

        integer                 :: newlev, nleaf, leaffin, leafs, leaff
        integer                 :: listleaf(4), nlf, nn, nlist, nj
        integer                 :: nmax, npre
        real(dp)                ::r, DS, RG(3), FIS, RATIO

        listleaf = 0
        !The loop for all bounds happens here
        !lev4 is the coarsest one.We start searching the coarsest and then move in finer and finer levels
        !The tree is a set of structured grids so we use that information to go deeper and deeper in the tree

        !for the coarsest level do all elements(The loop for all bounds happens here so we have to do for all coarsest
        !level elements
        if (nlev .eq. levmax) then
            leafstart = 1
            leaffin = nbound_lev(nlev)
            nlist = 0

            !in case of dummy element
        else if (ilev_t(nv, nlev + 1, 6) .lt. 0) then
            listleaf(1:2) = ilev_t(nv, nlev + 1, 1:2)
            !find the 4 finer cells.To do that we use the structured grid information
            leafstart = 1; leaffin = 2
            nlist = 1
            !in case of normal element
        else
            listleaf(1:4) = ilev_t(nv, nlev + 1, 1:4)
            leafstart = 1; leaffin = 4
            nlist = 1
        end if

        newlev = nlev - 1

        !if (my_rank.eq.0) then
        !    write(*,*) leafstart,nlev,nlist,nv
        !endif
        !if (my_rank.eq.0) then
        !if (nlev.lt.levmax.and.nlist.eq.1.and.ilev_t(nv,nlev+1,6).gt.0) then
        !  do nlf= leafstart,leaffin
        !     nleaf=listleaf(nlf)
        !   print *,nv,nlev,nleaf
        !     write(17,*) xs_lev(nleaf,nlev),&
        !                 ys_lev(nleaf,nlev),&
        !                 zs_lev(nleaf,nlev)

        !  enddo
        !    write(17,*)  xs_lev(listleaf(1),nlev),&
        !                 ys_lev(listleaf(1),nlev),&
        !                 zs_lev(listleaf(1),nlev)
        !    write(17,*)
        !    write(17,*)
        !    write(18,*)  xs_lev(nv,nlev+1),&
        !                 ys_lev(nv,nlev+1),&
        !                 zs_lev(nv,nlev+1)

        !else if (nlev.lt.levmax.and.nlist.eq.1.and.ilev_t(nv,nlev+1,6).lt.0) then
        !  do nlf= leafstart,leaffin
        !     nleaf=listleaf(nlf)
        !     print *,nv,nlev,nleaf
        !     write(16,*) xs_lev(nleaf,nlev),&
        !                 ys_lev(nleaf,nlev),&
        !                 zs_lev(nleaf,nlev)

        !  enddo
        !     write(16,*)
        !     write(16,*)
        !endif
        !read(*,*)
        !endif

        !search all elements in coarser grid and go to finer if needed.Initialy at coarsest level
        !we search all elements
        !at a finer level we search the 4
        do nlf = leafstart, leaffin
            if (nlist .eq. 0) then
                nleaf = nlf
            else
                nleaf = listleaf(nlf)
            end if
            if (ilev_t(nleaf, nlev, 6) .eq. 0) then
                write (*, *) 'Something is wrong with the tree'
                write (*, *) nleaf, nlev
                write (*, *)
                STOP
            end if

            RG(1) = xs_lev(nleaf, nlev)
            RG(2) = ys_lev(nleaf, nlev)
            RG(3) = zs_lev(nleaf, nlev)
            r = sqrt((XO(1) - RG(1))**2 + (XO(2) - RG(2))**2 + (XO(3) - RG(3))**2)
            DS = ds_lev(nleaf, nlev)
            RATIO = r/sqrt(DS)
            if ((RATIO .lt. 10 .and. nlev .gt. 0) .or. ilev_t(nleaf, nlev, 6) .eq. -1) then
                call tree_calc_3d_s(nleaf, newlev, leafstart, XO, SOURCE, neqs, neqf, &
                                    DIAG, E1, E2, E3, COSB, SINB, S, T, D, NSIDE, EPSS, ISING)
            else
                !write(*,*) newlev,ratio,r,sqrt(DS)
                if (nlev .ne. levmax) then
                    SOURCE = SOURCE - source_bound_lev(nleaf, neqs:neqf, nlev)*DS/(PI4*r)
                else
                    call FSOUR_A4(XO, RG, E1, E2, E3, &
                                  S, T, D, SINB, COSB, &
                                  DIAG, DS, NSIDE, EPSS, ISING, FIS)
                    SOURCE = SOURCE + source_bound_lev(nv, neqs:neqf, nlev)*FIS
                end if
            end if

        end do

    end subroutine tree_calc_3d_s

    !subroutine PHIELS calculates the potential induced by constant panels'
    ! XO is the point of calculation X1,Y1 is the first corner of the constant panel
    ! DS is the area of the face cosb,sinb give the direction assumed.
    subroutine PHIELS(X0, Y0, X1, Y1, DS, COSB, SINB, PHILS)
        implicit none
        real(dp), intent(in) :: X0, Y0, X1, Y1, DS, COSB, SINB
        real(dp), intent(out):: PHILS
        real(dp)             :: AKSIL, HTAL, TA1, TA2
        AKSIL = (X0 - X1)*COSB + (Y0 - Y1)*SINB !this is the vector X0-X1 in local coor
        HTAL = -(X0 - X1)*SINB + (Y0 - Y1)*COSB
        TA1 = AKSIL ! vector XO-X1 from the first point
        TA2 = AKSIL - DS!vector XO-X2 (since in l.c. X2=X1+DS then X0-X2 = XO-X1 - DS
        if (dabs(HTAL) .gt. 1.d-08) then
            PHILS = 1./(PI2)*(TA1*DLOG(DSQRT(TA1**2 + HTAL**2)) &
                              - TA2*DLOG(DSQRT(TA2**2 + HTAL**2)) &
                              + HTAL*(DATAN(TA1/HTAL) - DATAN(TA2/HTAL)) - DS)
        else
            if (abs(TA1) .lt. 1d-08) then
                PHILS = 1.d0/(2.d0*PI2)*(-TA2*DLOG(TA2**2) - 2.d0*DS)
            else if (abs(TA2) .lt. 1d-08) then
                PHILS = 1.d0/(2.d0*PI2)*(TA1*DLOG(TA1**2) - 2d0*DS)
            else
                PHILS = 1.d0/(2.d0*PI2)*(TA1*DLOG(TA1**2) - TA2*DLOG(TA2**2) - 2.d0*DS)
            end if
        end if

        return
    end subroutine PHIELS

!-----------------------------------------------------------------------
!     Chapter 1. FLAT CONSTANT SOURCE ELEMENT
!-----------------------------------------------------------------------
!  Subr      :FSOUR_A4
    subroutine FSOUR_A4(XO, RG, E1, E2, E3, &
                        S, T, D, SINB, COSB, &
                        DIAG, AREA, NSIDE, EPSS, ISING, FIS)

        implicit none

        integer, intent(in)   :: NSIDE, ISING
        real(dp), intent(in)  :: DIAG, AREA, EPSS
        real(dp), intent(out) :: FIS
        real(dp), intent(in)  :: XO(3), RG(3), E1(3), E2(3), E3(3), S(4), T(4), &
                                 D(4), SINB(4), COSB(4)
        real(dp)              :: X, Y, Z, RO, RATIO, AREA1, RO3, P, Q, FIL, &
                                 XK, YK, ZK, A1, A2, AK, ZP, AZ, R(4), E(4), H(4), UL(3), TOO(4), TINY

        integer             :: K, K1, K2

        TINY = 1d-10
        X = (XO(1) - RG(1))*E1(1) + (XO(2) - RG(2))*E1(2) + (XO(3) - RG(3))*E1(3)
        Y = (XO(1) - RG(1))*E2(1) + (XO(2) - RG(2))*E2(2) + (XO(3) - RG(3))*E2(3)
        Z = (XO(1) - RG(1))*E3(1) + (XO(2) - RG(2))*E3(2) + (XO(3) - RG(3))*E3(3)

        RO = sqrt(X*X + Y*Y + Z*Z)
        RATIO = RO/DIAG
        FIS = 0
        UL = 0
        if (RATIO .gt. 10) then
            AREA1 = AREA/PI4
            FIS = -AREA1/RO
            return
        end if

        UL(1) = 0.d0
        UL(2) = 0.d0
        FIL = 0.d0

        ZP = 1.d0
        if (Z .lt. 0.0d0) ZP = -1.d0
        AZ = abs(Z)
        AZ = dmax1(EPSS, AZ)
        Z = AZ*ZP
        ZK = Z*Z

        do K = 1, 4
            XK = (X - S(K))*(X - S(K))
            YK = (Y - T(K))*(Y - T(K))
            H(K) = (X - S(K))*(Y - T(K))
            R(K) = sqrt(XK + YK + ZK)
            E(K) = ZK + XK
        end do

        do K = 1, 4
            if (K .ne. NSIDE) then
                K1 = K
                K2 = K + 1
                if (K .eq. 4) K2 = 1
                A2 = R(K1) + R(K2) + D(K)
                A1 = R(K1) + R(K2) - D(K)
                if (abs(A1*A2) .lt. TINY) cycle
                AK = dlog(A1/A2)
                TOO(K) = -(X - S(K1))*SINB(K) + (Y - T(K1))*COSB(K)
                UL(1) = UL(1) + (T(K2) - T(K1))*AK/D(K)
                UL(2) = UL(2) - (S(K2) - S(K1))*AK/D(K)
                FIL = FIL - TOO(K)*AK
            end if
        end do

        UL(3) = 0.5d0*PI4
        if (ISING .ne. 0) then
            go to 4
        end if
        AZ = dmax1(TINY, AZ)
        Z = AZ*ZP

        UL(3) = 0.d0
        do K = 1, 4
        if (K .ne. NSIDE .and. abs(COSB(K)) .ge. TINY) then
            K1 = K
            K2 = K + 1
            if (K1 .eq. 4) K2 = 1
            UL(3) = UL(3) &
                    + atan((SINB(K)*E(K1) - COSB(K)*H(K1))/(Z*COSB(K)*R(K1))) &
                    - atan((SINB(K)*E(K2) - COSB(K)*H(K2))/(Z*COSB(K)*R(K2)))
        end if
        end do

4       FIS = FIS + (FIL + Z*UL(3))/PI4

        return
    end subroutine FSOUR_A4

!---- subroutine :GPANEL ----------------------------------------------
!        Data :  XCORN(3,i)     global co-ordinates of the vertices
!                IYNCP          if =1 then the c.p. is calculated
!       Output:  SK(4),TK(4)    local co-ordinates of the vertices
!                R(3)           g. c. of the center of the local system
!                EX,EY,EZ(3)    the local base
!                DIAG           diameter of the element
!                AREA           area
!                XCP(3)         g. c. of the control point
!                AMOMX,AMOMY    cartesian 1st order moments
!                AMOMXX,AMOMXY  cartesian 2nd order moments
!                       AMOMYY
!                D(4)           lengths of the element's sides
!                CBA(4), SBA(4) cosines and sines of the element's sides
!                               with respect to the local system
!                NCORN          =0 then the element is a quadrilateral
!                                k then the k-th side of the element is
!                                  of zero length
!
!
    subroutine GPANEL(XCORN, R, EX, EY, EZ, &
                      SK, TK, DIAG, AREA, XCP, &
                      D, CBA, SBA, &
                      IYNCP, NCORN, IERR)

        implicit none

        real(8), intent(in)  :: XCORN(3, 4)

        integer, intent(in)  :: IYNCP

        real(8), intent(out) :: SK(4), TK(4), R(3), EX(3), EY(3), EZ(3), &
                                DIAG, AREA, XCP(3), &
                                D(4), CBA(4), SBA(4)

        integer, intent(out)::NCORN, IERR

        real(8) :: XM(3), RR(4), RX(4), RY(4), DD(4), A(4), XCORN1(3, 4), &
                   DEX, DEY, DN, DE1, SO, TO, &
                   SCP, TCP, U, V, Ux, Uy, Vx, Vy, EPS, DET, &
                   DETS, DETT, DSCP, DTCP, TINY, TINYs

        integer :: J, K, NOD, K1, K2, ITER, L
        TINY = 1d-014
        TINYs = 1d-014

        !-- Define the local coordinate base     EX(), EY(), EZ()
        !          the element's center          XM()
        !      and the maximum diagonal          DIAG
        IERR = 0
        EX(1:3) = XCORN(1:3, 3) - XCORN(1:3, 1)
        EY(1:3) = XCORN(1:3, 4) - XCORN(1:3, 2)
        XM(1:3) = 0.25d0*(XCORN(1:3, 1) + XCORN(1:3, 2) + XCORN(1:3, 3) + XCORN(1:3, 4))
        DEX = dsqrt(EX(1)*EX(1) + EX(2)*EX(2) + EX(3)*EX(3))
        DEY = dsqrt(EY(1)*EY(1) + EY(2)*EY(2) + EY(3)*EY(3))
        DIAG = dmax1(DEX, DEY)
        EZ(1) = -EX(2)*EY(3) + EX(3)*EY(2)
        EZ(2) = -EX(3)*EY(1) + EX(1)*EY(3)
        EZ(3) = -EX(1)*EY(2) + EX(2)*EY(1)
        DN = dsqrt(dot_product(EZ, EZ))
        if (DN <= TINYs) then
            write (*, *) 'GPANEL:NCORN=-1, DN=', DN
            write (*, *) XCORN(1, 1), XCORN(2, 1), XCORN(3, 1)
            write (*, *) XCORN(1, 2), XCORN(2, 2), XCORN(3, 2)
            write (*, *) XCORN(1, 3), XCORN(2, 3), XCORN(3, 3)
            write (*, *) XCORN(1, 4), XCORN(2, 4), XCORN(3, 4)
            NCORN = -1
            return
        end if
        EZ(1:3) = EZ(1:3)/DN
        EY(1) = -EX(2)*EZ(3) + EX(3)*EZ(2)
        EY(2) = -EX(3)*EZ(1) + EX(1)*EZ(3)
        EY(3) = -EX(1)*EZ(2) + EX(2)*EZ(1)
        DEY = dsqrt(dot_product(EY, EY))
        EX(1:3) = EX(1:3)/DEX
        EY(1:3) = EY(1:3)/DEY

!-- Define the plane coordinates    XCORN1(3,.)
        DE1 = dot_product(EZ(1:3), XM(1:3) - XCORN(1:3, 1))
        do NOD = 1, 4
            XCORN1(1:3, NOD) = XCORN(1:3, NOD) + EZ(1:3)*((-1.d0)**(NOD - 1))*DE1
        end do ! NOD

!-- Define the local coordinates    SK() , TK()
        do J = 1, 4
            SK(J) = dot_product(EX(1:3), XCORN1(1:3, J) - XM(1:3))
            TK(J) = dot_product(EY(1:3), XCORN1(1:3, J) - XM(1:3))
        end do
        SO = ((TK(1) - TK(2))*SK(4) + (TK(4) - TK(1))*SK(2))/(3.d0*(TK(2) - TK(4)))
        TO = -TK(1)/3.d0

        SK(1:4) = SK(1:4) - SO
        TK(1:4) = TK(1:4) - TO

        R(1:3) = EX(1:3)*SO + EY(1:3)*TO + XM(1:3)

        NCORN = 0
        do K = 1, 4
            L = K + 1
            if (K == 4) L = 1
            D(K) = dsqrt((SK(L) - SK(K))**2 + (TK(L) - TK(K))**2)
            if (D(K) .gt. TINY) then
                SBA(K) = (TK(L) - TK(K))/D(K)
                CBA(K) = (SK(L) - SK(K))/D(K)
            else
                SBA(K) = 0.d0
                CBA(K) = 0.d0
                NCORN = K
            end if
        end do

!-- Calculate the element's control point
        if (IYNCP == 0) then
            XCP(1:3) = (XCORN(1:3, 1) + XCORN(1:3, 2) + XCORN(1:3, 3) + XCORN(1:3, 4))/4.d0
            goto 16
        end if

        SCP = 0.d0
        TCP = 0.d0
        EPS = 0.0001d0
        ITER = 0
10      ITER = ITER + 1
        do K = 1, 4
            RR(K) = dsqrt((SK(K) - SCP)**2 + (TK(K) - TCP)**2)
        end do
        U = 0.d0
        V = 0.d0
        UX = 0.d0
        UY = 0.d0
        VX = 0.d0
        VY = 0.d0

        do K = 1, 4
            L = K + 1
            if (K == 4) L = 1
            if (NCORN == K) cycle
            RX(K) = (SCP - SK(K))/RR(K) + (SCP - SK(L))/RR(L)
            RY(K) = (TCP - TK(K))/RR(K) + (TCP - TK(L))/RR(L)
            DD(K) = ((RR(K) + RR(L))*(RR(K) + RR(L)) - D(K)**2)/2.d0
            A(K) = (RR(K) + RR(L) - D(K))/(RR(K) + RR(L) + D(K))
            A(K) = dlog(A(K))
            U = U + (TK(L) - TK(K))*A(K)/D(K)
            V = V + (SK(K) - SK(L))*A(K)/D(K)
            UX = UX + (TK(L) - TK(K))*RX(K)/DD(K)
            UY = UY + (TK(L) - TK(K))*RY(K)/DD(K)
            VX = VX + (SK(K) - SK(L))*RX(K)/DD(K)
            VY = VY + (SK(K) - SK(L))*RY(K)/DD(K)
        end do

        DET = UX*VY - UY*VX
        DETS = UY*V - U*VY
        DETT = U*VX - UX*V
        DSCP = DETS/DET
        DTCP = DETT/DET
        SCP = SCP + DSCP
        TCP = TCP + DTCP
        if (ITER >= 10) goto 14
        if ((DSCP < EPS) .and. (DTCP < EPS)) goto 13
        goto 10
! the maximum number of iterations has been exceeded
14      IERR = 1
        SCP = 0.d0
        TCP = 0.d0
! convergence accomplished
13      XCP(1:3) = R(1:3) + EX(1:3)*SCP + EY(1:3)*TCP

16      AREA = 0.5d0*(SK(3) - SK(1))*(TK(2) - TK(4))

    end subroutine GPANEL
!--------------------------------------------------------------------
end submodule pmbound
! end module pmlib
module vpm_interpolate
    implicit none
contains
    !---------------------------------------------------------------------------!
    !-> subroutine back_to_particles                                            !
    !   This subroutine interpolates PM grid values back to particles at the    !
    !   positions they ARE.Convections takes place afterwards.                  !
    !   Input :                                                                 !
    !          itype (1,2) defines what value to interpolate to the particles   !
    !---------------------------------------------------------------------------!
    subroutine back_to_particles_3D(XP, QP, UP, GP, velocity_pm, deform_pm, &
                                    NVR, iproj, itype, NVRM)
        ! TODO: MAKE FLAGS FOR CALC U AND CALC G
        use projlib, only: projection_fun
        use vpm_types, only: dp
        use vpm_size, only: fine_grid
        use vpm_vars, only: neqpm
        implicit none
        integer, intent(in)     :: NVR, iproj, NVRM, itype
        real(dp), intent(in)    :: velocity_pm(3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3))
        real(dp), intent(in)    :: deform_pm(3, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3))
        real(dp), intent(inout) :: QP(neqpm + 1, NVRM), XP(3, NVRM), UP(3, NVRM), GP(3, NVRM)

        real(dp)         :: fx, fy, fz, f, x, y, z
        integer          :: inode, jnode, knode, i, j, k, nv, ips, ipf
        real(dp)         :: XBound(6), Dpm(3)

        Xbound = fine_grid%XBound
        Dpm = fine_grid%Dpm
        if (iproj .eq. 2) then
            ips = 0
            ipf = 1
        else if (iproj .eq. 3) then
            ips = 1
            ipf = 2
        else if (iproj .eq. 4) then
            ips = 1
            ipf = 2
        end if

        if (itype == 1) then
            UP(1:3, :) = 0.d0
        end if
        GP(1:3, :) = 0.d0
        do nv = 1, NVR
            !-->Find the cell/node  the  particle belongs for X and Y and Z direction.
            inode = int((XP(1, nv) - XBound(1))/Dpm(1)) + 1
            jnode = int((XP(2, nv) - XBound(2))/Dpm(2)) + 1
            knode = int((XP(3, nv) - XBound(3))/Dpm(3)) + 1
            ! Stop the program if the nodes are outside the domain bounds
            if (                                                        &
                inode - ips < 1 .or. inode + ipf > fine_grid%NN(1) .or. &
                jnode - ips < 1 .or. jnode + ipf > fine_grid%NN(2) .or. &
                knode - ips < 1 .or. knode + ipf > fine_grid%NN(3)      &
            ) then
                write (*, *) 'Particle outside domain bounds'
                write (*, *) 'inode, jnode, knode = ', inode, jnode, knode
                write (*, *) 'XP = ', XP(:, nv)
                stop
            end if

            !--We search the 4 nodes close to the particles
            do k = knode - ips, knode + ipf
                do j = jnode - ips, jnode + ipf
                    do i = inode - ips, inode + ipf
                        x = XP(1, nv) - (XBound(1) + (i - 1)*Dpm(1))
                        x = x / Dpm(1) ! Normalize the distance
                        fx = projection_fun(iproj, x)

                        y = XP(2, nv) - (XBound(2) + (j - 1)*Dpm(2))
                        y = y / Dpm(2) ! Normalize the distance
                        fy = projection_fun(iproj, y)

                        z = XP(3, nv) - (XBound(3) + (k - 1)*Dpm(3))
                        z = z / Dpm(3) ! Normalize the distance
                        fz = projection_fun(iproj, z)

                        f = fx*fy*fz
                        if (itype == 1) then
                            UP(1, nv) = UP(1, nv) + f*(velocity_pm(1, i, j, k))
                            UP(2, nv) = UP(2, nv) + f*(velocity_pm(2, i, j, k))
                            UP(3, nv) = UP(3, nv) + f*(velocity_pm(3, i, j, k))
                        end if
                        GP(2, nv) = GP(2, nv) + f*(deform_pm(2, i, j, k)) 
                        GP(1, nv) = GP(1, nv) + f*(deform_pm(1, i, j, k)) 
                        GP(3, nv) = GP(3, nv) + f*(deform_pm(3, i, j, k)) 
                    end do
                end do
            end do

            ! Normalize GP using the particle "volume"
            GP(1:3, nv) = GP(1:3, nv)*QP(neqpm + 1, nv)
        end do
    end subroutine back_to_particles_3D

    !---------------------------------------------------------------------------!
    !-> subroutine back_to_particles                                            !
    !   This subroutine interpolates a quantity from the PM grid to the         !
    !     particle locations XP                                                 !
    !  Input :                                                                 !
    !        Q_pm : PM grid values                                              !
    !        XP   : Particle positions                                          !
    !        QP   : Particle values                                             !
    !        NVR   : Number of particles                                        !
    !        iproj : Projection type                                            !
    !---------------------------------------------------------------------------!
    subroutine interpolate_particle_Q(Q_pm, XP, QP, NVR, iproj, NVR_size)
        use projlib, only: projection_fun
        use vpm_types, only: dp
        use vpm_size, only: fine_grid!NN, XBound, Dpm
        use vpm_vars, only: neqpm

        implicit None
        integer, intent(in)     :: NVR, iproj, NVR_size
        real(dp), intent(in)    :: Q_pm(neqpm, fine_grid%NN(1), fine_grid%NN(2), fine_grid%NN(3))
        real(dp), intent(in)    :: XP(3, NVR_size)
        real(dp), intent(out)   :: QP(neqpm + 1, NVR_size)

        real(dp)                :: fx, fy, fz, f, x, y, z
        integer                 :: inode, jnode, knode, i, j, k, nv, ips, ipf
        real(dp)                :: Xbound(6), Dpm(3)

        Xbound  = fine_grid%XBound
        Dpm     = fine_grid%Dpm

        if (iproj .eq. 2) then
            ips = 0
            ipf = 1
        else if (iproj .eq. 3) then
            ips = 1
            ipf = 2
        else if (iproj .eq. 4) then
            ips = 1
            ipf = 2
        else 
            ips = 1
            ipf = 2
        end if

        QP(1:neqpm, :) = 0
        do nv = 1, NVR
            !-->Find the cell/node  the  particle belongs for X and Y and Z direction.
            inode = int((XP(1, nv) - XBound(1))/Dpm(1)) + 1
            jnode = int((XP(2, nv) - XBound(2))/Dpm(2)) + 1
            knode = int((XP(3, nv) - XBound(3))/Dpm(3)) + 1
            !--We search the 4 nodes close to the particles
            do k = knode - ips, knode + ipf
                do j = jnode - ips, jnode + ipf
                    do i = inode - ips, inode + ipf
                        x = (XP(1, nv) - XBound(1) - (i - 1)*Dpm(1))/Dpm(1)
                        fx = projection_fun(iproj, x)

                        y = (XP(2, nv) - XBound(2) - (j - 1)*Dpm(2))/Dpm(2)
                        fy = projection_fun(iproj, y)

                        z = (XP(3, nv) - XBound(3) - (k - 1)*Dpm(3))/Dpm(3)
                        fz = projection_fun(iproj, z)

                        f = fx*fy*fz
                        QP(1:neqpm, nv) = QP(1:neqpm, nv) + f*Q_pm(1:neqpm, i, j, k)
                    end do
                end do
            end do
        end do
        do i = 1, neqpm
            QP(i, :) = QP(i, :)*QP(neqpm + 1, :)
        end do
    end subroutine interpolate_particle_Q

end module vpm_interpolate

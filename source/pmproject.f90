module projlib
    use vpm_types, only: dp, cartesian_grid

    private
    real(dp), save       :: XMIN_pm, XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
    real(dp), save       :: DXpm, DYpm, DZpm, DVpm
    integer, save        :: NXpm, NYpm, NZpm, NXs, NXf, NYs, NYf, NZs, NZf
    integer, save        :: IDVPM, ND

    public :: print_projlib_info, projlibinit,                          &
              project_particles_3D, project_particles_2D,               &
              project_vol3d, project_vol2d, project_particles_2D_vol,   &
              projection_fun

contains
    subroutine projlibinit(grid, IDVPM_in, ND_in)
        implicit none
        type(cartesian_grid), intent(in) :: grid
        integer, intent(in)              :: IDVPM_in, ND_in

        IDVPM = IDVPM_in
        ND = ND_in

        XMIN_pm = grid%Xbound(1)
        YMIN_pm = grid%Xbound(2)
        ZMIN_pm = grid%Xbound(3)
        XMAX_pm = grid%Xbound(4)
        YMAX_pm = grid%Xbound(5)
        ZMAX_pm = grid%Xbound(6)

        DXpm = grid%Dpm(1)
        DYpm = grid%Dpm(2)
        DZpm = grid%Dpm(3)

        NXpm = grid%NN(1)
        NYpm = grid%NN(2)
        NZpm = grid%NN(3)

        NXs = grid%NN_bl(1)
        NYs = grid%NN_bl(2)
        NZs = grid%NN_bl(3)
        NXf = grid%NN_bl(4)
        NYf = grid%NN_bl(5)
        NZf = grid%NN_bl(6)

        DVpm = DXpm*DYpm
        if (ND .eq. 3) then
            DVpm = DVpm*DZpm
        end if
    end subroutine projlibinit

    ! --------------------------------------------------------------------------!
    !-->subroutine project_particles_3D                                          !
    !   This subroutine projects particle values on the PM grid                !
    !   The values projected are :                                             !
    !      - Mass  -> becomes Density on the grid                              !
    !      - Volume                                                            !
    !      - Dilatation                                                        !
    !      - Phi                                                               !
    !      - PsiX , PsiY                                                       !
    !--------------------------------------------------------------------------!
    subroutine project_particles_3D( &
        Q_pm, QP, XP, Qprojtype, ipar, isize, ieq, neq, QINF, iparsize &
        )
        ! Qproj -> RHS of the PM grid
        ! Qpar -> Particle values QP_scat
        ! QpX -> Particle positions XP_scat
        ! Qprojtype -> Type of projection function
        ! ipar -> Number of particles NVR_p
        ! isize -> Number of variables
        ! ieq -> Variables to project
        ! neq -> Number of equations
        ! QINF -> Inflow values
        ! iparsize -> Number of particles NVR_p

        use MPI

        implicit none
        integer, intent(in)     :: ipar, isize, ieq(neq), iparsize, neq
        real(dp), intent(out)   :: Q_pm(neq, NXpm, NYpm, NZpm)
        real(dp), intent(in)    :: QP(neq + 1, iparsize), XP(3, iparsize), QINF(neq)
        integer, intent(in)     :: Qprojtype(iparsize)

        real(dp)                :: fx, fy, fz, f, x, y, z
        integer                 :: inode, jnode, knode, i, j, k, nv, itype, ips, ipf
        integer                 :: ierr
        real(dp)                :: fpriv(isize, NXpm, NYpm, NZpm)

        !-->Projection function (TSC)
        Q_pm = 0.d0
        fpriv = 0.d0

        !$omp parallel private(nv, inode, jnode, knode, i, j, k, x, y, z, fx, fy, fz, f, itype, ips, ipf) &
        !$omp shared(QP,fpriv, Q_pm, QpX, Qprojtype, ipar, isize, ieq, neq)
        !$omp do
        do nv = 1, ipar
            itype = Qprojtype(nv)
            !-->Find the cell/node  the  particle belongs for X and Y and Z direction.
            inode = int((XP(1, nv) - XMIN_pm)/DXpm) + 1
            jnode = int((XP(2, nv) - YMIN_pm)/DYpm) + 1
            if (ND .eq. 3) then
                knode = int((XP(3, nv) - ZMIN_pm)/DZpm) + 1
            else
                knode = 1
            end if

            if (itype .eq. 2) then
                ! For type 2: project to the nearest node and the next (ips = 0, ipf = 1)
                ips = 0
                ipf = 1
            else
                ! For other types: project to the previous, nearest, and next two nodes (ips = 1, ipf = 2)
                ips = 1
                ipf = 2
            end if

            ! Check if particle is within the PM grid
            if (inode - ips < 1 .or. inode + ipf > NXpm) then
                print *, 'Particle out of bounds in X direction', inode, ips, ipf, nv, XP(1, nv)
                read *, ierr
                STOP
            end if
            if (jnode - ips < 1 .or. jnode + ipf > NYpm) then
                print *, 'Particle out of bounds in Y direction', jnode, ips, ipf, nv, XP(2, nv)
                read *, ierr
                STOP
            end if
            if (ND == 3 .and. (knode - ips < 1 .or. knode + ipf > NZpm)) then
                print *, 'Particle out of bounds in Z direction', knode, ips, ipf, nv, XP(3, nv)
                read *, ierr
                STOP
            end if
            !--We search the 4 nodes close to the particles
            do k = knode - ips, knode + ipf
                do j = jnode - ips, jnode + ipf
                    do i = inode - ips, inode + ipf
                        x = (XP(1, nv) - XMIN_pm - (i - 1)*DXpm)/DXpm
                        fx = projection_fun(itype, x)

                        y = (XP(2, nv) - YMIN_pm - (j - 1)*DYpm)/DYpm
                        fy = projection_fun(itype, y)

                        if (ND .eq. 3) then
                            z = (XP(3, nv) - ZMIN_pm - (k - 1)*DZpm)/DZpm
                            fz = projection_fun(itype, z)
                            f = fx*fy*fz
                        else
                            f = fx*fy
                        end if

                        ! Qproj(ieq(1:neq-1),i,j,k) = Qproj(ieq(1:neq-1),i,j,k) +&
                        !           f * (QPar(ieq(1:neq-1),nv))!-QINF(1:neq-1)*QPar(ieq(neq),nv))

                        ! For the last equation, we subtract the inflow value
                        ! Qproj(ieq(neq),i,j,k) = Qproj(neq,i,j,k) +&
                        !           f * (QPar(ieq(neq),nv)-QINF(neq))
                        fpriv(ieq(1:neq), i, j, k) = fpriv(ieq(1:neq), i, j, k) &
                                                     + f*(QP(ieq(1:neq), nv) - QINF(1:neq))

                    end do
                end do
            end do
        end do
        !$omp end do

        ! Reduce the thread-private results into Qproj
        !$omp critical
        Q_pm = Q_pm + fpriv
        !$omp end critical
        !$omp end parallel

        !--After the projection all the values will be divided by Vol_pm to take into account the volume
        !  effect in the interpolation.We extend volume and density values to no-particle regions
    end subroutine project_particles_3D

    !--------------------------------------------------------------------------!
    !-->subroutine project_vol3d                                               !
    !   This subroutine projects particle values on the PM grid                !
    !   The values projected are :                                             !
    !      - Mass  -> becomes Density on the grid                              !
    !      - Volume                                                            !
    !      - Dilatation                                                        !
    !      - Phi                                                               !
    !      - PsiX , PsiY                                                       !
    !--------------------------------------------------------------------------!
    subroutine project_vol3d(Qproj, isize, ieq, neq, iflag)
        implicit none
        integer, intent(in)                 :: neq, isize, iflag
        integer, intent(in), dimension(neq) :: ieq
        real(dp), intent(inout)             :: Qproj(isize, NXpm, NYpm, NZpm)
        integer                             :: i, j, k, IDVPMt

        IDVPMt = IDVPM
        if (iflag .eq. 1) IDVPMt = 1

        if (IDVPMt .eq. 0) then
            do k = 1, NZpm
                do j = 1, NYpm
                    do i = 1, NXpm
                        if (                                                      &
                            (i .lt. NXs) .or. (i .gt. NXf) .or. (j .lt. NYs) .or. &
                            (j .ge. NYf) .or. (k .lt. NZs) .or. (k .gt. NZf)      &
                        ) then
                            Qproj(ieq(1:neq), i, j, k) = 0.d0
                            Qproj(ieq(neq), i, j, k) = DVpm
                            cycle
                        end if
                        Qproj(ieq(1:neq - 1), i, j, k) = Qproj(ieq(1:neq - 1), i, j, k)/(Qproj(ieq(neq), i, j, k))
                    end do
                end do
            end do
        else if (IDVPMt .eq. 1) then
            do k = 1, NZpm
                do j = 1, NYpm
                    do i = 1, NXpm
                        ! if (i.lt.NXs.or.i.gt.NXf.or.j.lt.NYs.or.j.ge.NYf.or.k.lt.NZs.or.k.gt.NZf) then
                        !     Qproj (ieq(1:neq-1),i,j,k) = 0.d0
                        !     cycle
                        ! endif
                        Qproj(ieq(1:neq), i, j, k) = Qproj(ieq(1:neq), i, j, k)/DVpm
                        ! Qproj(ieq(neq), i, j, k) = DVpm
                    end do
                end do
            end do
        else
            write (*, *) 'WRONG IDVPM'
            STOP
        end if
    end subroutine project_vol3d
    
    !--------------------------------------------------------------------------!
    !-->subroutine project_particles                                           !
    !   This subroutine projects particle values on the PM grid                !
    !   The values projected are :                                             !
    !      - Mass  -> becomes Density on the grid                              !
    !      - Volume                                                            !
    !      - Dilatation                                                        !
    !      - Phi                                                               !
    !      - PsiX , PsiY                                                       !
    !--------------------------------------------------------------------------!
    subroutine project_particles_2D( &
        Qproj, Qpar, QpX, Qprojtype, ipar, isize, ieq, neq, QINF &
        )
        implicit none
        integer, intent(in)     :: ipar, isize, ieq(neq)
        real(dp), intent(out)   :: Qproj(isize, NXpm, NYpm, NZpm)
        real(dp), intent(in)    :: Qpar(isize, ipar), QpX(3, ipar), QINF(neq)
        integer, intent(in)     :: Qprojtype(ipar)

        real(dp)                :: fx, fy, f, x, y
        integer                 :: inode, jnode, i, j, k, nv, itype, neq, ips, ipf

        !-->Projection function (TSC)
        Qproj = 0.d0
        !!$omp parallel private(nv,Qprojpriv,jnode,inode,itype,ips,ipf,i,j,k,f,fx,fy,x,y)
        !!allocate(Qprojpriv(NXpm,NYpm,NZpm,isize))
        !!Qprojpriv=0.d0
        !!$omp do
        do nv = 1, ipar
            itype = Qprojtype(nv)
            !-->Find the cell/node  the  particle belongs for X and Y and Z direction.
            inode = int((QpX(1, nv) - XMIN_pm)/DXpm) + 1
            jnode = int((QpX(2, nv) - YMIN_pm)/DYpm) + 1
            ! knode = int(ZVR(nv) / DZpm) 3D

            !    enddo

            !--We search the 4 nodes close to the particles
            !    do k = knode - 1, knode + 2 3D
            if (itype .eq. 2) then
                ips = 0
                ipf = 1
            else
                ips = 1
                ipf = 2
            end if
            k = 1 !nbj
            do j = jnode - ips, jnode + ipf
                do i = inode - ips, inode + ipf

                    x = (QpX(1, nv) - XMIN_pm - (i - 1)*DXpm)/DXpm
                    fx = projection_fun(itype, x)

                    y = (QpX(2, nv) - YMIN_pm - (j - 1)*DYpm)/DYpm
                    fy = projection_fun(itype, y)

                    ! z  = (ZVR(nv) - (k-1) * DZpm) / DZpm 3D
                    ! fz = projection_fun(itype,z) * z     3D

                    f = fx*fy !* fz 3D
                    Qproj(ieq(1:neq - 1), i, j, k) = Qproj(ieq(1:neq - 1), i, j, k) + &
                                                     f*(QPar(ieq(1:neq - 1), nv) - QINF(1:neq - 1)*QPar(ieq(neq), nv))
                    Qproj(ieq(neq), i, j, k) = Qproj(ieq(neq), i, j, k) + &
                                               f*(QPar(ieq(neq), nv) - QINF(neq))
                end do
            end do
        end do
        !!$omp end do
        !!$omp critical
        ! !Qproj=Qproj+Qprojpriv
        !!$omp end critical
        !!deallocate(Qprojpriv)
        !!$omp end parallel

        !--After the projection all the values will be divided by Vol_pm to take into account the volume
        !  effect in the interpolation.We extend volume and density values to no-particle regions

    end subroutine project_particles_2D

    subroutine project_particles_2D_vol(Qproj, Qpar, QpX, Qprojtype, ipar, isize, ieq, neq)
        implicit none
        integer, intent(in) :: ipar, isize, ieq(neq)
        ! real(dp), intent(out), dimension(:,:,:) :: Qproj
        real(dp), intent(out):: Qproj(NXpm, NYpm, NZpm)
        !f2py depend(NXpm, NYpm, NZpm) :: Qproj(NXpm, NYpm, NZpm)
        ! real(dp), allocatable:: Qprojpriv(:, :, :, :)
        real(dp), intent(in) :: Qpar(isize, ipar), QpX(3, ipar) !, QINF(neq)
        integer, intent(in) :: Qprojtype(ipar)
        real(dp)   :: fx, fy, f, x, y !, z, fz
        integer            :: i, j, k, nv, inode, jnode, itype, neq, ips, ipf
        ! integer            ::omp_get_max_threads, omp_get_num_threads
        ! integer            :: knode, nbj, nb
        !-->Projection function (TSC)
        Qproj = 0.d0
        !!$omp parallel private(nv,Qprojpriv,jnode,inode,itype,ips,ipf,i,j,k,f,fx,fy,x,y)
        !!allocate(Qprojpriv(NXpm,NYpm,NZpm,isize))
        !!Qprojpriv=0.d0
        !!$omp do
        do nv = 1, ipar
            itype = Qprojtype(nv)
            if (abs(Qpar(1, nv)) .lt. 0.5*DVpm) cycle
            !-->Find the cell/node  the  particle belongs for X and Y and Z direction.
            inode = int((QpX(1, nv) - XMIN_pm)/DXpm) + 1
            jnode = int((QpX(2, nv) - YMIN_pm)/DYpm) + 1
            ! knode = int(ZVR(nv) / DZpm) 3D

            !    enddo

            !--We search the 4 nodes close to the particles
            !    do k = knode - 1, knode + 2 3D
            if (itype .eq. 2) then
                ips = 0
                ipf = 1
            else
                ips = 1
                ipf = 2
            end if
            k = 1 !nbj
            do j = jnode - ips, jnode + ipf
                do i = inode - ips, inode + ipf

                    x = (QpX(1, nv) - XMIN_pm - (i - 1)*DXpm)/DXpm
                    fx = projection_fun(itype, x)

                    y = (QpX(2, nv) - YMIN_pm - (j - 1)*DYpm)/DYpm
                    fy = projection_fun(itype, y)

                    ! z  = (ZVR(nv) - (k-1) * DZpm) / DZpm 3D
                    ! fz = projection_fun(itype,z) * z     3D

                    f = fx*fy !* fz 3D
                    Qproj(i, j, k) = Qproj(i, j, k) + &
                                     f*(QPar(ieq(neq), nv))
                end do
            end do
        end do
        !!$omp end do
        !!$omp critical
        ! !Qproj=Qproj+Qprojpriv
        !!$omp end critical
        !!deallocate(Qprojpriv)
        !!$omp end parallel

        !--After the projection all the values will be divided by Vol_pm to take into account the volume
        !  effect in the interpolation.We extend volume and density values to no-particle regions

    end subroutine project_particles_2D_vol

    !--------------------------------------------------------------------------!
    !-->subroutine project_vol2d                                               !
    !   This subroutine projects particle values on the PM grid                !
    !   The values projected are :                                             !
    !      - Mass  -> becomes Density on the grid                              !
    !      - Volume                                                            !
    !      - Dilatation                                                        !
    !      - Phi                                                               !
    !      - PsiX , PsiY                                                       !
    !--------------------------------------------------------------------------!
    subroutine project_vol2d(Qproj, isize, ieq, neq, iflag)
        implicit none
        integer, intent(in) :: isize, ieq(neq), iflag
        ! real(dp) , intent(out), dimension(:,:, :, :) :: Qproj
        real(dp)   :: Qproj(isize, NXpm, NYpm, NZpm)
        !f2py depend(isize, NXpm, NYpm, NZpm) :: Qproj(isize, NXpm, NYpm, NZpm)
        ! real(dp)   :: fx, fy, fz, f, x, y, z,
        integer            ::  i, j, k, neq, IDVPMt
        ! integer            :: inode, jnode, knode, nv, itype, nbj, nb
        IDVPMt = IDVPM
        if (iflag .eq. 1) IDVPMt = 1
        if (IDVPMt .eq. 0) then
            !$omp parallel private(i,j,k)
            k = 1 !nbj
            !$omp do
            do j = 1, NYpm
                do i = 1, NXpm
                    if (i .lt. NXs .or. i .gt. NXf .or. j .lt. NYs .or. j .ge. NYf) then
                        Qproj(ieq(1:neq - 1), i, j, k) = 0.d0
                        cycle
                    end if
                    Qproj(ieq(1:neq - 1), i, j, k) = Qproj(ieq(1:neq - 1), i, j, k)/(Qproj(ieq(neq), i, j, k))
                    !else
                    !    Qproj (i,j,k,ieq(1:neq-1)) = 0.d0
                    !    Qproj (i,j,k,ieq(neq)) =DVpm
                    !endif
                end do
            end do
            !$omp end do
            !$omp end parallel
        else if (IDVPMt .eq. 1) then
            !$omp parallel private(i,j,k)
            k = 1 !nbj
            !$omp do
            do j = 1, NYpm
                do i = 1, NXpm
                    if (i .lt. NXs .or. i .gt. NXf .or. j .lt. NYs .or. j .ge. NYf) then
                        Qproj(ieq(neq), i, j, k) = DVpm
                        Qproj(ieq(1:neq - 1), i, j, k) = 0.d0
                        cycle
                    end if
                    Qproj(ieq(1:neq - 1), i, j, k) = Qproj(ieq(1:neq - 1), i, j, k)/DVpm
                    Qproj(ieq(neq), i, j, k) = DVpm
                    !endif
                end do
            end do
            !$omp end do
            !$omp end parallel
        else
            write (*, *) 'WRONG IDVPM'
            STOP
        end if
    end subroutine project_vol2d

    !-----------------------------------------------------------------------------!
    !-->Function projection_fun                                                   !
    !   This subroutine defines the projection functions and finds the projected  !
    !   value of the Particle Value on the PMgrid                                 !
    !   Input :                                                                   !
    !      itype       : type of projection function                              !
    !      x           : position of projection
    !-----------------------------------------------------------------------------!
    real(dp) function projection_fun(itype, x) 
        implicit none
        real(dp), intent(in) :: x
        integer, intent(in)          :: itype

        real(dp)             :: xabs

        xabs = abs(x)
        if (itype .eq. 2) then
            ! Rectified Linear Unit
            projection_fun = max(0.d0, 1.d0 - xabs)
        else if (itype .eq. 3) then
            !--Triangular-Shaped Cloud function
            if (xabs .lt. 0.5d0) projection_fun = 0.5d0*(xabs + 3.d0/2.d0)**2 - 3.d0/2.d0*(xabs + 0.5d0)**2
            if (xabs .ge. 0.5d0 .and. xabs .le. 3.d0/2.d0) projection_fun = 0.5d0*(-xabs + 3.d0/2.d0)**2
            if (xabs .gt. 3.d0/2.d0) projection_fun = 0.d0
        else if (itype .eq. 4) then
            !--Quadratic function
            if (xabs .lt. 1.d0) projection_fun = 1 - 2.5d0*xabs**2 + 1.5d0*xabs**3
            if (xabs .ge. 1.d0 .and. xabs .le. 2.d0) projection_fun = 0.5d0*(2.d0 - xabs)**2*(1.d0 - xabs)
            if (xabs .gt. 2.d0) projection_fun = 0.d0
        else if (itype .eq. 5) then
            !--Quadratic function
            if (xabs .le. 0.5d0) projection_fun = 1.d0 - xabs**2
            if (xabs .gt. 0.5d0 .and. xabs .le. 1.5d0) projection_fun = 0.5d0*(1.d0 - xabs)*(2.d0 - xabs)
            if (xabs .gt. 1.5d0) projection_fun = 0.d0
        else if (itype .eq. 6) then
            !--Quadratic function
            if (xabs .le. 1.d0) projection_fun = 1.d0 - xabs**2
            if (xabs .gt. 1.d0 .and. xabs .le. 2.d0) projection_fun = 0.5d0*(1.d0 - xabs)*(2.d0 - xabs)
            if (xabs .gt. 2.d0) projection_fun = 0.d0
        else
            write (*, *) 'No such projection function', itype
            STOP
        end if

    End Function projection_fun

   !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Printers
   !!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine print_projlib_info()
        print *, "PROJLIB INFO"
        print *, "============"
        print *, ""
        print *, achar(9), 'XMIN_pm', XMIN_pm
        print *, achar(9), 'XMAX_pm', XMAX_pm
        print *, achar(9), 'YMIN_pm', YMIN_pm
        print *, achar(9), 'YMAX_pm', YMAX_pm
        print *, achar(9), 'ZMIN_pm', ZMIN_pm
        print *, achar(9), 'ZMAX_pm', ZMAX_pm
        print *, achar(9), 'DXpm', DXpm
        print *, achar(9), 'DYpm', DYpm
        print *, achar(9), 'DZpm', DZpm
        print *, achar(9), 'DVpm', DVpm
        print *, achar(9), 'NXpm', NXpm
        print *, achar(9), 'NYpm', NYpm
        print *, achar(9), 'NZpm', NZpm
        print *, achar(9), 'NXs', NXs
        print *, achar(9), 'NXf', NXf
        print *, achar(9), 'NYs', NYs
        print *, achar(9), 'NYf', NYf
        print *, achar(9), 'NZs', NZs
        print *, achar(9), 'NZf', NZf
        print *, achar(9), 'IDVPM', IDVPM
        print *, achar(9), 'ND', ND
    end subroutine print_projlib_info
End module projlib

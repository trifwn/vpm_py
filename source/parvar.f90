module parvar
    use vpm_types, only: dp
    integer, save              :: neq_par
    ! PARTICLE VARIABLES
    integer                    :: NVR, NVR_size
    real(dp), pointer, save    :: XP(:, :), QP(:, :), UP(:, :), GP(:, :)
    real(dp), pointer, save    :: VP(:, :)

    ! SCATTERED QUANTITIES
    integer                    :: NVR_p
    real(dp), allocatable      :: XP_scatt(:, :), QP_scatt(:, :), UP_scatt(:, :), GP_scatt(:, :)
    integer, allocatable       :: NVR_projtype_scatt(:)

    !   XP: Particle positions

    !   QP: Particle quantities typically:
    !   -->1 Vorticity X
    !   -->2 Vorticity Y
    !   -->3 Vorticity Z
    !   -->4 Dilatation
    !   -->5 Pseudopressure
    !   -->6 Mass
    !   -->7 Volume


    !   -->neq+1: WAS The volume fraction of the particle
    !   VP: Particle Volumes (old QP(neq+1, nv))

    !   UP: Particle velocites                -> U = grad x psi =
    !   GP: Particle deformation of vorticity -> deformation = - omega*\grad u

    ! Printers
    public :: print_parvar_info

contains

    subroutine associate_particles(NVR_in, NVR_size_in, neq_in, XP_in, QP_in, UP_in, GP_in)
        use MPI
        implicit none
        integer, intent(in)                    :: NVR_in, NVR_size_in, neq_in
        real(dp), intent(in), target           :: XP_in(:, :), QP_in(:, :)
        real(dp), intent(in), target, optional :: UP_in(:, :), GP_in(:, :)
        integer                          :: ierr, my_rank
        call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

        if (associated(XP)) nullify (XP)
        if (associated(QP)) nullify (QP)
        if (associated(UP)) nullify (UP)
        if (associated(GP)) nullify (GP)

        if (my_rank .eq. 0) then
            NVR = NVR_in
            NVR_size = NVR_size_in
            neq_par = neq_in
            XP => XP_in
            QP => QP_in
            if (present(UP_in)) then
                UP => UP_in
            end if
            if (present(GP_in)) then
                GP => GP_in
            end if
        end if

        ! BCAST NVR
        call MPI_BCAST(NVR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(NVR_size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(neq_par, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (NVR .eq. 0) return
    end subroutine associate_particles

   !!!!!!!!!!!!!!!!!!!!!!!
    ! Printers
   !!!!!!!!!!!!!!!!!!!!!!!
    subroutine print_parvar_info()
        use console_io, only: dp_2d_ptr_info, i_1d_alloc_info, dp_2d_alloc_info
        print *, "PARVAR INFO"
        print *, "============"
        print *, ""
        print *, achar(9)//"NVR = ", NVR
        call dp_2d_ptr_info("XP", XP)
        call dp_2d_ptr_info("QP", QP)
        call dp_2d_ptr_info("UP", UP)
        call dp_2d_ptr_info("GP", GP)
        print *, "Scattered Quantities"
        print *, ""
        call i_1d_alloc_info("NVR_projtype_scatt", NVR_projtype_scatt)
        print *, achar(9)//"NVR_p = ", NVR_p
        call dp_2d_alloc_info("XP_scatt", XP_scatt)
        call dp_2d_alloc_info("QP_scatt", QP_scatt)
        call dp_2d_alloc_info("UP_scatt", UP_scatt)
        call dp_2d_alloc_info("GP_scatt", GP_scatt)
    end subroutine print_parvar_info

    subroutine print_particle_info() bind(C, name='print_particles')
        integer :: i
        do i = 1, NVR
            write (*, "(A, 1I5, A, 3F15.5)") "Particle X", i, "->", XP(1, i), XP(2, i), XP(3, i)
            write (*, "(A, 1I5, A, 6F15.5)") "Particle U", i, "->", UP(1, i), UP(2, i), UP(3, i)
            write (*, "(A, 1I5, A, 6F15.5)") "Particle G", i, "->", GP(1, i), GP(2, i), GP(3, i)
            write (*, "(A, 1I5, A, 6F15.5)") "Particle Q", i, "->", QP(:, i)
            write (*, *) ""
        end do
    end subroutine print_particle_info

    subroutine print_particle_positions() bind(C, name='print_particle_positions')
        integer :: i
        do i = 1, NVR
            write (*, "(A, 1I5, A, 3F15.5)") "Particle ", i, " X ->", XP(1, i), XP(2, i), XP(3, i)
        end do
    end subroutine print_particle_positions

   !!!!!!!!!!!!!!!!!!!!!!!
    ! Setters and Getters
   !!!!!!!!!!!!!!!!!!!!!!!

    subroutine set_neq(neq_in)
        implicit none
        integer, intent(in):: neq_in
        neq_par = neq_in
    end subroutine set_neq

    subroutine set_particle_positions(XP_in) bind(C, name='set_particle_positions')
        use iso_c_binding
        implicit none
        real(c_double), dimension(3, NVR_size) :: XP_in
        XP = XP_in(1:3, 1:NVR)
    end subroutine set_particle_positions

    subroutine set_particle_strengths(QP_in) bind(C, name='set_particle_strengths')
        use iso_c_binding
        implicit none
        real(c_double), dimension(neq_par + 1, NVR_size) :: QP_in
        QP = QP_in(1:neq_par + 1, 1:NVR)
    end subroutine set_particle_strengths

    subroutine set_particle_velocities(UP_in) bind(C, name='set_particle_velocities')
        use iso_c_binding
        implicit none
        real(c_double), dimension(3, NVR_size) :: UP_in
        UP = UP_in(1:3, 1:NVR)
    end subroutine set_particle_velocities

    subroutine set_particle_deformation(GP_in) bind(C, name='set_particle_deformation')
        use iso_c_binding
        implicit none
        real(c_double), dimension(3, NVR_size) :: GP_in
        GP = GP_in(1:3, 1:NVR)
    end subroutine set_particle_deformation

    subroutine get_nvr(NVR_out) bind(C, name='get_num_particles')
        use iso_c_binding
        implicit none
        integer(c_int) :: NVR_out
        NVR_out = NVR
    end subroutine get_nvr

    subroutine get_nvr_size(NVR_out) bind(C, name='get_nvr_size')
        use iso_c_binding
        implicit none
        integer(c_int) :: NVR_out
        NVR_out = NVR
    end subroutine get_nvr_size
end module parvar

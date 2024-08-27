module parvar
   use base_types, only: dp
   integer                            :: NVR
   real(dp), pointer, save    :: XP(:, :), QP(:, :)
   real(dp), pointer, save    :: UP(:, :), GP(:, :)
   !     QP
   !   -->1 Vorticity X
   !   -->2 Vorticity Y
   !   -->3 Vorticity Z
   !   -->4 Dilatation
   !   -->5 Pseudopressure
   !   -->6 Mass
   !   -->7 Volume
   integer, allocatable, save         :: NVR_projtype(:)
   integer, save                      :: neq

   ! Printers
   public :: print_parvar_info

contains

   !!!!!!!!!!!!!!!!!!!!!!!
   ! Printers
   !!!!!!!!!!!!!!!!!!!!!!!
   subroutine print_parvar_info()
      use io, only: dp_2d_ptr_info, i_1d_alloc_info
      print *, "PARVAR INFO"
      print *, "============"
      print *, ""
      print *, achar(9)//"NVR = ", NVR
      call dp_2d_ptr_info("XP", XP)
      call dp_2d_ptr_info("QP", QP)
      call dp_2d_ptr_info("UP", UP)
      call dp_2d_ptr_info("GP", GP)
      call i_1d_alloc_info("NVR_projtype", NVR_projtype)

   end subroutine print_parvar_info

   subroutine print_particle_info() bind(C, name='print_particles')
      integer :: i
      do i = 1, NVR
         write(*, "(A, 1I5, A, 3F15.5)") "Particle X", i, "->", XP(1, i), XP(2, i), XP(3, i)
         write(*, "(A, 1I5, A, 6F15.5)") "Particle U", i, "->", UP(1, i), UP(2, i), UP(3, i)
         write(*, "(A, 1I5, A, 6F15.5)") "Particle G", i, "->", GP(1, i), GP(2, i), GP(3, i)
         write(*, "(A, 1I5, A, 6F15.5)") "Particle Q", i, "->", QP(1, i), QP(2, i), QP(3, i), QP(4, i)
         write (*, *) ""
      end do
   end subroutine print_particle_info

   subroutine print_particle_positions() bind(C, name='print_particle_positions')
      integer :: i
      do i = 1, NVR
         write(*, "(A, 1I5, A, 3F15.5)") "Particle ", i, " X ->", XP(1, i), XP(2, i), XP(3, i)
      end do
   end subroutine print_particle_positions

   !!!!!!!!!!!!!!!!!!!!!!!
   ! Setters and Getters
   !!!!!!!!!!!!!!!!!!!!!!!

   Subroutine set_neq(neq_in)
      implicit none
      integer, intent(in):: neq_in
      neq = neq_in
   end subroutine set_neq

   subroutine get_particle_positions(XP_out) bind(C, name='get_particle_positions')
      use ND_Arrays
      implicit none
      type(ND_Array), intent(out) :: XP_out
      XP_out = from_intrinsic(XP, shape(XP))
   end subroutine get_particle_positions

   subroutine get_particle_strengths(QP_out) bind(C, name='get_particle_strengths')
      use iso_c_binding
      use ND_Arrays
      implicit none
      type(ND_Array), intent(out) :: QP_out
      QP_out = from_intrinsic(QP, shape(QP))
   end subroutine get_particle_strengths

   subroutine get_particle_deformation(GP_out) bind(C, name='get_particle_deformation')
      use iso_c_binding
      use ND_Arrays
      implicit none
      type(ND_Array), intent(out) :: GP_out
      GP_out = from_intrinsic(GP, shape(GP))
   end subroutine get_particle_deformation

   subroutine get_particle_velocities(UP_out) bind(C, name='get_particle_velocities')
      use iso_c_binding
      use ND_Arrays
      implicit none
      type(ND_Array), intent(out) :: UP_out
      UP_out = from_intrinsic(UP, shape(UP))
   end subroutine get_particle_velocities

   subroutine set_particle_positions(XP_in) bind(C, name='set_particle_positions')
      use iso_c_binding
      implicit none
      real(c_double), dimension(3, NVR) :: XP_in
      XP = XP_in
   end subroutine set_particle_positions

   subroutine set_particle_strengths(QP_in) bind(C, name='set_particle_strengths')
      use iso_c_binding
      implicit none
      real(c_double), dimension(neq + 1, NVR) :: QP_in
      QP = QP_in
   end subroutine set_particle_strengths

   subroutine set_particle_velocities(UP_in) bind(C, name='set_particle_velocities')
      use iso_c_binding
      implicit none
      real(c_double), dimension(3, NVR) :: UP_in
      UP = UP_in
   end subroutine set_particle_velocities

   subroutine set_particle_deformation(GP_in) bind(C, name='set_particle_deformation')
      use iso_c_binding
      implicit none
      real(c_double), dimension(3, NVR) :: GP_in
      GP = GP_in
   end subroutine set_particle_deformation

   subroutine get_nvr(NVR_out) bind(C, name='get_num_particles')
      use iso_c_binding
      implicit none
      integer(c_int) :: NVR_out

      NVR_out = NVR
   end subroutine get_nvr
end module parvar

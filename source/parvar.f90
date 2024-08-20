module parvar
   integer                            :: NVR
   double precision, pointer, save    :: XP(:, :), QP(:, :)
   double precision, pointer, save    :: UP(:, :), GP(:, :)
   !     QP
   !   -->1 Vorticity X
   !   -->2 Vorticity Y
   !   -->3 Vorticity Z
   !   -->4 Dilatation
   !   -->5 Pseudopressure
   !   -->6 Mass
   !   -->7 Volume
   integer, allocatable, save         :: NVR_projtype(:)

   ! double precision, allocatable        :: XP_CFD_Sa(:, :), UT_CFD_Sa(:), UN_CFD_Sa(:), &
      ! DS_CFD_Sa(:), Vo_CFD_Sa(:)

   ! Getters
   public :: get_NVR, get_XP, get_QP, get_UP, get_GP, get_NVR_projtype
   ! Setters
   public :: set_NVR, set_XP, set_QP, set_UP, set_GP, set_NVR_projtype
   ! Printers
   public :: print_parvar_info

contains

   !!!!!!!!!!!!!!!!!!!!!!!
   ! Getters
   !!!!!!!!!!!!!!!!!!!!!!!
   function get_NVR() result(NVR_)
      integer :: NVR_
      NVR_ = NVR
   end function get_NVR

   function get_XP() result(XP_)
      double precision, pointer :: XP_(:, :)
      XP_ => XP
   end function get_XP

   function get_QP() result(QP_)
      double precision, pointer :: QP_(:, :)
      QP_ => QP
   end function get_QP

   function get_UP() result(UP_)
      double precision, pointer :: UP_(:, :)
      UP_ => UP
   end function get_UP

   function get_GP() result(GP_)
      double precision, pointer :: GP_(:, :)
      GP_ => GP
   end function get_GP

   function get_NVR_projtype() result(NVR_projtype_)
      integer, allocatable :: NVR_projtype_(:)
      NVR_projtype_ = NVR_projtype
   end function get_NVR_projtype

   !!!!!!!!!!!!!!!!!!!!!!!
   ! Setters
   !!!!!!!!!!!!!!!!!!!!!!!
   subroutine set_NVR(NVR_)
      integer, intent(in) :: NVR_
      NVR = NVR_
   end subroutine set_NVR

   subroutine set_XP(XP_)
      double precision, pointer, intent(in) :: XP_(:, :)
      XP => XP_
   end subroutine set_XP

   subroutine set_QP(QP_)
      double precision, pointer, intent(in) :: QP_(:, :)
      QP => QP_
   end subroutine set_QP

   subroutine set_UP(UP_)
      double precision, pointer, intent(in) :: UP_(:, :)
      UP => UP_
   end subroutine set_UP

   subroutine set_GP(GP_)
      double precision, pointer, intent(in) :: GP_(:, :)
      GP => GP_
   end subroutine set_GP

   subroutine set_NVR_projtype(NVR_projtype_)
      integer, intent(in) :: NVR_projtype_(:)
      NVR_projtype = NVR_projtype_
   end subroutine set_NVR_projtype

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



end module parvar

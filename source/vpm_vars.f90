
Module vpm_vars
   use base_types, only: dp
   real(dp), allocatable       :: XP_scatt(:, :), QP_scatt(:, :), UP_scatt(:, :), GP_scatt(:, :)
   real(dp)                    :: DT_c, V_ref, NI
   integer, allocatable                :: NVR_projscatt(:)
   integer                             :: interf_iproj, ncell_rem

   integer                             :: ncoarse, nparcell1d
   integer                             :: neqpm, NVR_p, NVR_size, NTIME_pm

   integer, save                       :: IPMWRITE, idefine, iynslice
   integer, save                       :: mrem = 1
   integer, save                       :: IPMWSTART(10), IPMWSTEPS(10)

   ! Printer
   public :: print_vpm_vars_info

contains

   Subroutine print_vpm_vars_info
      print *, "VPM_VARS INFO"
      print *, "============"
      print *, ""
      if (allocated(XP_scatt)) then
         print *, achar(9)//"XP_scatt", " (2D): Size = (", size(XP_scatt,1), ",", size(XP_scatt,2), ")"
         print *, achar(9)//"Sample values: ", XP_scatt(1,1:min(4,size(XP_scatt,2)))
      else
         print '(A)', achar(9)//"XP_scatt Not allocated"
      end if

      if (allocated(QP_scatt)) then
         print *, achar(9)//"QP_scatt", " (2D): Size = (", size(QP_scatt,1), ",", size(QP_scatt,2), ")"
         print *, achar(9)//"Sample values: ", QP_scatt(1,1:min(4,size(QP_scatt,2)))
      else
         print '(A)', achar(9)//"QP_scatt Not allocated"
      end if

      if (allocated(UP_scatt)) then
         print *, achar(9)//"UP_scatt", " (2D): Size = (", size(UP_scatt,1), ",", size(UP_scatt,2), ")"
         print *, achar(9)//"Sample values: ", UP_scatt(1,1:min(4,size(UP_scatt,2)))
      else
         print '(A)', achar(9)//"UP_scatt Not allocated"
      end if

      if (allocated(GP_scatt)) then
         print *, achar(9)//"GP_scatt", " (2D): Size = (", size(GP_scatt,1), ",", size(GP_scatt,2), ")"
         print *, achar(9)//"Sample values: ", GP_scatt(1,1:min(4,size(GP_scatt,2)))
      else
         print '(A)', achar(9)//"GP_scatt Not allocated"
      end if

      print *, achar(9), 'DT_c = ', DT_c
      print *, achar(9), 'V_ref = ', V_ref
      print *, achar(9), 'NI = ', NI

      if (allocated(NVR_projscatt)) then
         print *, achar(9)//"NVR_projscatt", " (1D): Size = (", size(NVR_projscatt), ")"
         print '(A,4I12)', achar(9)//"Sample values: ", NVR_projscatt(1:min(size(NVR_projscatt), 4))
      else
         print '(A)', achar(9)//"NVR_projscatt Not allocated"
      end if

      print *, achar(9), 'interf_iproj = ', interf_iproj
      print *, achar(9), 'ncell_rem = ', ncell_rem
      print *, achar(9), 'ncoarse = ', ncoarse
      print *, achar(9), 'nparcell1d = ', nparcell1d
      print *, achar(9), 'neqpm = ', neqpm
      print *, achar(9), 'NVR_p = ', NVR_p
      print *, achar(9), 'NVR_size = ', NVR_size
      print *, achar(9), 'NTIME_pm = ', NTIME_pm
      print *, achar(9), 'IPMWRITE = ', IPMWRITE
      print *, achar(9), 'mrem = ', mrem
      print *, achar(9), 'idefine = ', idefine
      print *, achar(9), 'iynslice = ', iynslice
   
      print  *, achar(9)//"IPMWSTART", " (1D): Size = (", size(IPMWSTART), ")"
      print '(A,4I12)', achar(9)//"Sample values: ", IPMWSTART(1:min(size(IPMWSTART), 4))
      print '(A)', ""

      print  *, achar(9)//"IPMWSTEPS", " (1D): Size = (", size(IPMWSTEPS), ")"
      print '(A,4I12)', achar(9)//"Sample values: ", IPMWSTEPS(1:min(size(IPMWSTEPS), 4))
      print '(A)', ""


   End Subroutine print_vpm_vars_info

End Module vpm_vars
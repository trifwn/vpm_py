Module vpm_vars
   use base_types, only: dp
   real(dp)                   :: V_ref, NI
   integer, save              :: interf_iproj, idefine
   integer, save              :: mrem = 1
   integer                    :: neqpm, NTIME_pm
   integer, save              :: IPMWRITE
   integer, save              :: IPMWSTART(10), IPMWSTEPS(10)
   integer                    :: OMPTHREADS
   public :: print_vpm_vars_info

contains
   subroutine print_vpm_vars_info
      print *, "VPM_VARS INFO"
      print *, "============"
      print *, ""
      print *, achar(9), 'V_ref = ', V_ref
      print *, achar(9), 'NI = ', NI

      print *, achar(9), 'interf_iproj = ', interf_iproj
      print *, achar(9), 'neqpm = ', neqpm
      print *, achar(9), 'NTIME_pm = ', NTIME_pm
      print *, achar(9), 'IPMWRITE = ', IPMWRITE
      print *, achar(9), 'mrem = ', mrem
      print *, achar(9), 'idefine = ', idefine
   
      print  *, achar(9)//"IPMWSTART", " (1D): Size = (", size(IPMWSTART), ")"
      print '(A,4I12)', achar(9)//"Sample values: ", IPMWSTART(1:min(size(IPMWSTART), 4))
      print '(A)', ""

      print  *, achar(9)//"IPMWSTEPS", " (1D): Size = (", size(IPMWSTEPS), ")"
      print '(A,4I12)', achar(9)//"Sample values: ", IPMWSTEPS(1:min(size(IPMWSTEPS), 4))
      print '(A)', ""
   End subroutine print_vpm_vars_info
End Module vpm_vars
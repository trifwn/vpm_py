Module vpm_vars
   use base_types, only: dp
   integer, save              :: ND = 3
   integer                    :: NTIME_pm
   integer                    :: neqpm
   real(dp)                   :: V_ref, NI
   integer, save              :: idefine
   integer, save              :: interf_iproj
   integer, save              :: mrem = 1
   integer, save              :: nremesh
   integer, save              :: iynbc 
   integer, save              :: ibctyp
   integer, save              :: iyntree, ilevmax, itree
   integer, save              :: NBI, NBJ, NBK 
   integer                    :: OMPTHREADS
   integer, save              :: iret

   integer, save              :: IPMWRITE
   integer, save              :: IPMWSTART(10), IPMWSTEPS(10)
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
      print  *, achar(9)//"iynbc", " = ", iynbc
      print  *, achar(9)//"iret", " = ", iret
      print  *, achar(9)//"NBI", " = ", NBI
      print  *, achar(9)//"NBJ", " = ", NBJ
      print  *, achar(9)//"NBK", " = ", NBK
      print  *, achar(9)//"NREMESH", " = ", nremesh
      print  *, achar(9)//"iyntree", " = ", iyntree
      print  *, achar(9)//"ilevmax", " = ", ilevmax
      print  *, achar(9)//"itree", " = ", itree
      print  *, achar(9)//"ibctyp", " = ", ibctyp
   
      print  *, achar(9)//"IPMWSTART", " (1D): Size = (", size(IPMWSTART), ")"
      print '(A,4I12)', achar(9)//"Sample values: ", IPMWSTART(1:min(size(IPMWSTART), 4))
      print '(A)', ""

      print  *, achar(9)//"IPMWSTEPS", " (1D): Size = (", size(IPMWSTEPS), ")"
      print '(A,4I12)', achar(9)//"Sample values: ", IPMWSTEPS(1:min(size(IPMWSTEPS), 4))
      print '(A)', ""
   End subroutine print_vpm_vars_info
End Module vpm_vars
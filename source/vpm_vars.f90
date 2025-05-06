module vpm_vars
    use vpm_types, only: dp, timestepInformation, solveInformation

    integer, save              :: ND = 3
    integer                    :: NTIME_pm
    integer                    :: neqpm
    real(dp)                   :: V_ref
    integer, save              :: idefine
    integer, save              :: interf_iproj
    integer, save              :: mrem = 1
    integer, save              :: iynbc
    integer, save              :: ibctyp
    integer, save              :: iyntree, ilevmax, itree
    integer, save              :: NBI, NBJ, NBK
    integer                    :: OMPTHREADS
    integer, save              :: iret

    integer, parameter         :: SOLVER_SERIAL_PMESH = 0,  &
                                  SOLVER_YAPS = 1,          &
                                  SOLVER_MUDPACK = 2
    integer                    :: SOLVER = 2 
    type(timestepInformation)  :: timestep_info
    type(solveInformation)     :: solve_info

    public :: print_vpm_vars_info

contains
    subroutine print_vpm_vars_info
        print *, "VPM_VARS INFO"
        print *, "============"
        print *, ""
        print *, achar(9), 'V_ref = ', V_ref

        print *, achar(9), 'interf_iproj = ', interf_iproj
        print *, achar(9), 'neqpm = ', neqpm
        print *, achar(9), 'NTIME_pm = ', NTIME_pm
        print *, achar(9), 'mrem = ', mrem
        print *, achar(9), 'idefine = ', idefine
        print *, achar(9)//"iynbc", " = ", iynbc
        print *, achar(9)//"iret", " = ", iret
        print *, achar(9)//"NBI", " = ", NBI
        print *, achar(9)//"NBJ", " = ", NBJ
        print *, achar(9)//"NBK", " = ", NBK
        print *, achar(9)//"iyntree", " = ", iyntree
        print *, achar(9)//"ilevmax", " = ", ilevmax
        print *, achar(9)//"itree", " = ", itree
        print *, achar(9)//"ibctyp", " = ", ibctyp
    end subroutine print_vpm_vars_info
End module vpm_vars

# Function to find compilers
set(CMAKE_Fortran_COMPILER_INIT "mpifort")
function(find_compilers)
    # Message in RED
    message(STATUS "${Green}Finding compilers${ColourReset}")
    
    # Intel Fortran Compiler
    find_program(IFORT ifx)
    if (IFORT)
        message("\tFound Intel Fortran Compiler: ${IFORT}")
        set(USE_INTEL_COMPILER TRUE PARENT_SCOPE)

        # Get the ifx version
        execute_process(COMMAND ${IFORT} --version OUTPUT_VARIABLE IFORT_VERSION)
        string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" IFORT_VERSION_MATCH ${IFORT_VERSION})
        message("\tIntel Fortran Compiler version: ${IFORT_VERSION_MATCH}")

        # Extract the major version number
        string(REGEX MATCH "^([0-9]+)" IFORT_VERSION_MAJOR ${IFORT_VERSION_MATCH})

        # find_package(MKL)
        set(MKL_FOUND TRUE)
        if (MKL_FOUND)
            message("\tFound MKL: ${MKL_INCLUDE_DIRS}")
            set(USE_MKL TRUE PARENT_SCOPE)
            # Set MKL_FLAG based on the version
            if (IFORT_VERSION_MAJOR GREATER 2022)
                set(MKL_FLAG "-qmkl" PARENT_SCOPE)
                message("\tSetting MKL_FLAG to '-qmkl' for Intel Fortran Compiler version ${IFORT_VERSION_MATCH}")
            else()
                set(MKL_FLAG "-mkl" PARENT_SCOPE)
                message("\tSetting MKL_FLAG to '-mkl' for Intel Fortran Compiler version ${IFORT_VERSION_MATCH}")
            endif()
        else()
            message("\tMKL not found")
            set(USE_MKL FALSE PARENT_SCOPE)
        endif()

        # GET MKL_ROOT
        if (DEFINED ENV{MKLROOT})
            set(MKLROOT $ENV{MKLROOT})
            message("\tMKLROOT: ${MKLROOT}")
        else()
            message("\tMKLROOT not found")
        endif()
        
    else()
        message("\tIntel Fortran Compiler not found")
        set(USE_INTEL_COMPILER FALSE PARENT_SCOPE)
        set(USE_MKL FALSE PARENT_SCOPE)
    endif()
endfunction()
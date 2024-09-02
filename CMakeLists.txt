cmake_minimum_required(VERSION 3.22)

# Pretty print messages
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

# Define options
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Setting directories
set(SRC_HP ${CMAKE_SOURCE_DIR}/source)
set(SRC_TEST ${CMAKE_SOURCE_DIR}/test)
set(LIB_DIR ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

set(EXENAMEIN vpm)

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

# Grab Python, 3.8 or newer
# find_package(Python 3.11 REQUIRED
#   COMPONENTS Interpreter Development.Module NumPy)

# Grab the variables from a local Python installation
# F2PY headers
execute_process(
COMMAND "${Python_EXECUTABLE}"
-c "import numpy.f2py; print(numpy.f2py.get_include())"
OUTPUT_VARIABLE F2PY_INCLUDE_DIR
OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Print out the discovered paths
include(CMakePrintHelpers)
message(STATUS "${Green}Python Configuration:${ColourReset}")
cmake_print_variables(Python_INCLUDE_DIRS)
cmake_print_variables(F2PY_INCLUDE_DIR)
cmake_print_variables(Python_NumPy_INCLUDE_DIRS)

# Function to find compilers
function(find_compilers)
    # Message in RED
    message(STATUS "${Green}Finding compilers${ColourReset}")
    
    # Intel Fortran Compiler
    find_program(IFORT ifx)
    if (IFORT)
        message("\tFound Intel Fortran Compiler: ${IFORT}")
        set(HAS_INTEL TRUE PARENT_SCOPE)

        # Get the ifx version
        execute_process(COMMAND ${IFORT} --version OUTPUT_VARIABLE IFORT_VERSION)
        string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" IFORT_VERSION_MATCH ${IFORT_VERSION})
        message("\tIntel Fortran Compiler version: ${IFORT_VERSION_MATCH}")

        # Extract the major version number
        string(REGEX MATCH "^([0-9]+)" IFORT_VERSION_MAJOR ${IFORT_VERSION_MATCH})

        # Set MKL_FLAG based on the version
        if (IFORT_VERSION_MAJOR GREATER 2022)
            set(MKL_FLAG "-qmkl" PARENT_SCOPE)
            message("\tSetting MKL_FLAG to '-qmkl' for Intel Fortran Compiler version ${IFORT_VERSION_MATCH}")
        else()
            set(MKL_FLAG "-mkl" PARENT_SCOPE)
            message("\tSetting MKL_FLAG to '-mkl' for Intel Fortran Compiler version ${IFORT_VERSION_MATCH}")
        endif()
    else()
        message("\tIntel Fortran Compiler not found")
        set(HAS_INTEL FALSE PARENT_SCOPE)
    endif()

    # If Intel is found and CMP is not 1 then set CMAKE_Fortran_COMPILER to mpif90 and INTEL to TRUE
    if (IFORT)
        set(CMAKE_Fortran_COMPILER_ID "Intel" PARENT_SCOPE)
        set(USE_INTEL TRUE PARENT_SCOPE)
    else()
        # Display error message
        set(USE_INTEL FALSE PARENT_SCOPE)
        set(CMAKE_Fortran_COMPILER_ID "GNU" PARENT_SCOPE)
    endif()
endfunction()

# Call the function
find_compilers()
message(STATUS "${Green}System Configuration Detected:${ColourReset}")
message("\tFortran compiler: ${CMAKE_Fortran_COMPILER}")
message("\tHas intel ifx: ${HAS_INTEL}")
if (HAS_INTEL)
    message("\tIFX location: ${IFORT}")
endif()


# Determine the compiler name (e.g., "Intel", "GNU")
string(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" COMPILER_NAME)
string(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)
set(EXENAME "${EXENAMEIN}_${BUILD_TYPE}_${COMPILER_NAME}")
set(PATHOBJHP "${CMAKE_BINARY_DIR}/bin_${BUILD_TYPE}_${COMPILER_NAME}")
set(MODULE_PATH "${CMAKE_BINARY_DIR}/mod_${BUILD_TYPE}_${COMPILER_NAME}")

file(MAKE_DIRECTORY ${PATHOBJHP})
file(MAKE_DIRECTORY ${MODULE_PATH})

# Set output directories
set(CMAKE_Fortran_MODULE_DIRECTORY ${MODULE_PATH})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PATHOBJHP})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PATHOBJHP})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PATHOBJHP})

# Print configuration
message(STATUS "${Green}Folder Configuration:${ColourReset}")
message("\tCMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message("\tSRCHP: ${SRC_HP}")
message("\tBASE EXENAMEIN: ${EXENAMEIN}")
message("\tCMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message("\tBuild type: ${CMAKE_BUILD_TYPE}")
message("\tBinary directory: ${CMAKE_BINARY_DIR}")
message("\tCMAKE ARCHIVE OUTPUT DIRECTORY: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
message("\tCMAKE LIBRARY OUTPUT DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message("\tCMAKE RUNTIME OUTPUT DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Setting compiler and flags before project()
message(STATUS "${Green}Compiler Configuration:${ColourReset}")
if (USE_INTEL)
    message("\tConfiguring for Intel compiler")
    set(CMAKE_Fortran_COMPILER mpifort)
    set(FLAGS "-DASCII=1 DHAS_INTEL=1")
    set(MOD_FLAG "-module ${MODULE_PATH}")

    set(compile_MKL_headers, TRUE)
    
    # Extend the debug flags
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 ${MKL_FLAG} -qopenmp -g -traceback -fpe0 -fPIC")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all,nouninit")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -init=snan -init=arrays")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ftrapuv")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -debug full")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -diag-enable=all")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -implicitnone")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ffast-math -march=native")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fstack-protector-all")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -warn all")
    
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ${MKL_FLAG} -qopenmp -ffast-math -march=native -fPIC")
    
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-O0 ${MKL_FLAG}  -qopenmp -g -traceback -fpe0 -check all,nouninit -fPIC") 
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-O3 ${MKL_FLAG} -qopenmp -ffast-math -march=native -fPIC")
else()
    message("\tConfiguring for GNU compiler")
    set(CMAKE_Fortran_COMPILER mpifort)
    set(FLAGS "-DASCII=1")
    set(MOD_FLAG "-J${MODULE_PATH}")
    
    set(compile_MKL_headers, FALSE)
    # Extend the debug flags
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbacktrace -finit-real=nan -finit-integer=nan -fPIC")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all ")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ffpe-trap=invalid,overflow,underflow")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fimplicit-none ")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ffast-math -march=native ")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fstack-protector-all")
    # set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall")
    # set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wextra") 

    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fPIC -ffast-math -march=native")

    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-O0 -g -fbacktrace -fcheck=all -finit-real=nan -finit-integer=nan -fPIC")

    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-O3 -fPIC")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
endif()

# Print configuration for compiler
message("\tCMAKE_Fortran_COMPILER: ${CMAKE_Fortran_COMPILER}")
message("\tMPI_FORTRAN_COMPILER: ${MPI_FORTRAN_COMPILER}")
message("\tFLAGS: ${FLAGS}")
message("\tCMAKE_Fortran_FLAGS_DEBUG: ${CMAKE_Fortran_FLAGS_DEBUG}")
message("\tCMAKE_Fortran_FLAGS_RELEASE: ${CMAKE_Fortran_FLAGS_RELEASE}")
message("\tCMAKE_EXE_LINKER_FLAGS_DEBUG: ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
message("\tCMAKE_EXE_LINKER_FLAGS_RELEASE: ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")

# Project information
project(vpm_py VERSION 1.0 LANGUAGES Fortran)


set(VPM_LIB_INCL 
    ${SRC_HP}/vpm_remesh.f90
    ${SRC_HP}/vpm_interpolate.f90
    ${SRC_HP}/vpm_gcalc.f90
    ${SRC_HP}/vpm_mpi.f90
)

set(PM_LIB_FILES_MKL
    ${SRC_HP}/pmlib.f90
    ${SRC_HP}/pmsolve.f90
    ${SRC_HP}/pmbound.f90
    ${SRC_HP}/pinfdomain.f90
)

set(PM_LIB_FILES_FISHPACK
    ${SRC_HP}/pmlib.f90
    ${SRC_HP}/pmbound.f90
    ${SRC_HP}/pinfdomain.f90
    ${SRC_HP}/pmsolve_fish.f90
)

set(YAPSLIB_FILES
    ${SRC_HP}/yaps.f90
    ${SRC_HP}/yaps2d.f90
    ${SRC_HP}/yaps3d.f90
)


set(TEST_EXE_SRC
    ${SRC_TEST}/test.f90
    ${SRC_TEST}/test_problems.f90
)

add_library(constants OBJECT ${SRC_HP}/base_types.f90)
add_library(types OBJECT ${SRC_HP}/constants.f90)
target_link_libraries(types PRIVATE constants)

add_library(arrays SHARED ${SRC_HP}/arrays.f90)
target_link_libraries(arrays PRIVATE types)

add_library(io OBJECT ${SRC_HP}/io.f90)

# MKL
if (compile_MKL_headers)
    add_library(mkl_dfti OBJECT ${SRC_HP}/mkl_dfti.f90)
    add_library(mkl_poisson OBJECT ${SRC_HP}/mkl_poisson.f90)
    target_link_libraries(mkl_poisson PRIVATE mkl_dfti)
    
    # Build the PM library with MKL
    add_library(pmlib OBJECT ${PM_LIB_FILES_MKL})
    target_link_libraries(pmlib PRIVATE mkl_poisson) 
else()
    # set(FISHPACK_SRC_DIR ${CMAKE_SOURCE_DIR}/third_party/fishpack)
    set(FISHPACK_SRC_DIR ${CMAKE_SOURCE_DIR}/../fishpack)
    # We need to find the fishpack library and link it
    find_library(FISHPACK_LIB fishpack HINTS ${FISHPACK_SRC_DIR})
    if (FISHPACK_LIB)
        message(STATUS "Found fishpack library: ${FISHPACK_LIB}")
    else()
        # Search for the directory of the fishpack source code and compile it
        message(STATUS "${Green}Fishpack library not found. Compiling from source${ColourReset}")
        # Define unique output directories for fishpack
        set(FISHPACK_BUILD_DIR ${CMAKE_BINARY_DIR}/fishpack_build)

        # Add the fishpack subdirectory
        add_subdirectory(${FISHPACK_SRC_DIR} fishpack_build)

        # Build the PM library with Fishpack
        add_library(pmlib OBJECT ${PM_LIB_FILES_FISHPACK})
        target_link_libraries(pmlib PRIVATE fishpack)
        # Add the fish_lib to the include directories
        target_include_directories(pmlib PRIVATE ${CMAKE_BINARY_DIR}/fish_lib)
    endif()
endif()
add_library(mpi_matrices OBJECT  ${SRC_HP}/mpi_matrices.f90)

add_library(parvar OBJECT  ${SRC_HP}/parvar.f90)
target_link_libraries(parvar PRIVATE io arrays)

add_library(pmgrid OBJECT  ${SRC_HP}/pmgrid.f90)
target_link_libraries(pmgrid PRIVATE io)

add_library(pmeshpar OBJECT  ${SRC_HP}/pmeshpar.f90)
target_link_libraries(pmeshpar PRIVATE io)

add_library(pmproject OBJECT ${SRC_HP}/pmproject.f90)

add_library(yaps OBJECT ${YAPSLIB_FILES})
target_link_libraries(yaps PRIVATE mpi_matrices pmlib)

add_library(vpm_vars OBJECT ${SRC_HP}/vpm_vars.f90)
add_library(vpm_size OBJECT ${SRC_HP}/vpm_size.f90)

# Add all object files as sources to a dummy target
add_custom_target(dummy ALL
    DEPENDS yaps pmlib pmproject parvar pmeshpar pmgrid vpm_vars vpm_size
)

add_library(vpm OBJECT ${SRC_HP}/vpm.f90 )
add_dependencies(vpm dummy)
target_include_directories(vpm PRIVATE ${SRC_HP}) # ${VPM_LIB_INCL})
target_link_libraries(vpm PRIVATE vpm_vars vpm_size io )

# -------------------------------------------------------------------------------------------------
#                                            Library
# -------------------------------------------------------------------------------------------------
add_library(${EXENAME} SHARED ${SRC_HP}/api.f90)
add_dependencies(${EXENAME} dummy vpm) # Ensure proper dependency resolution
target_link_libraries(${EXENAME} 
    PRIVATE mpi_matrices parvar pmgrid pmeshpar pmproject pmlib yaps vpm vpm_vars vpm_size io arrays types constants
)
target_include_directories(${EXENAME} PUBLIC ${PATHOBJHP}) # Set include directories using PUBLIC

# -------------------------------------------------------------------------------------------------
#                                            VPM TEST Executable
# -------------------------------------------------------------------------------------------------
# Add executable for testing purposes (not necessary)
add_executable(${EXENAME}_exe ${TEST_EXE_SRC})
add_dependencies(${EXENAME}_exe dummy vpm)
target_link_libraries(${EXENAME}_exe 
    PRIVATE mpi_matrices parvar pmgrid pmeshpar pmproject pmlib yaps vpm vpm_vars vpm_size io arrays types constants
)
target_include_directories(${EXENAME}_exe PUBLIC ${PATHOBJHP}) 

# # -------------------------------------------------------------------------------------------------
#                                           #  Unit tests
# # -------------------------------------------------------------------------------------------------
add_executable(test_arrays_exe ${SRC_TEST}/test_arrays.f90)
add_dependencies(test_arrays_exe arrays)
target_link_libraries(test_arrays_exe PRIVATE arrays)


###########################################################
#    Custom targets 
############################################################

# Add a custom target to clean the build
add_custom_target(clean_build
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
    COMMENT "Cleaning build"
)

# # Add custom target to build Release 
# add_custom_target(build_release
#     COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${EXENAME}
#     COMMENT "Building Release"
# )

# # Add custom target to build Debug
# add_custom_target(build_debug
#     COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${EXENAME}
#     COMMENT "Building Debug"
# )

# # Add target to build Release and Debug
# add_custom_target(build_all
#     COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${EXENAME}
#     COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${EXENAME}
#     COMMENT "Building Release and Debug"
# )
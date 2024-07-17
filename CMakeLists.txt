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
set(EXENAMEIN vpm)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

# Grab Python, 3.8 or newer
find_package(Python 3.8 REQUIRED
  COMPONENTS Interpreter Development.Module NumPy)

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
        message(FATAL_ERROR "No Fortran compiler found")
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

message(STATUS "${Green}Compiler Configuration:${ColourReset}")
# Print compiler information
# Setting compiler and flags before project()

if (USE_INTEL)
    set(CMAKE_Fortran_COMPILER mpifort)
    message("\tConfiguring for Intel compiler")
    set(FLAGS "-DASCII=1")
    set(MOD_FLAG "-module")
    
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -qmkl -qopenmp -g -traceback")
    # Extend the debug flags
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fpe0 ")
    # set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all ")
    # set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check bounds -check format -check output_conversion -check pointers -check stack")
    # set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check uninit")
    # set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fstack-security-check")
    # set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ftrapuv")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -warn all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -qmkl -qopenmp")
    
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-module ${PATHOBJHP} -O0 -qmkl -qopenmp")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-module ${PATHOBJHP} -O3 -qmkl -qopenmp")

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(EXENAME "${EXENAMEIN}_debug")
        set(PATHOBJHP "${CMAKE_BINARY_DIR}/bin_debug_ifo")
    else()
        set(EXENAME "${EXENAMEIN}")
        set(PATHOBJHP "${CMAKE_BINARY_DIR}/bin_release_ifo")
    endif()
else()
    message(FATAL_ERROR "No Fortran compiler found")
endif()

# Print configuration for compiler
message("\tCMAKE_Fortran_COMPILER: ${CMAKE_Fortran_COMPILER}")
message("\tMPI_FORTRAN_COMPILER: ${MPI_FORTRAN_COMPILER}")
message("\tFLAGS: ${FLAGS}")
message("\tCMAKE_Fortran_FLAGS_DEBUG: ${CMAKE_Fortran_FLAGS_DEBUG}")
message("\tCMAKE_Fortran_FLAGS_RELEASE: ${CMAKE_Fortran_FLAGS_RELEASE}")
message("\tCMAKE_EXE_LINKER_FLAGS_DEBUG: ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
message("\tCMAKE_EXE_LINKER_FLAGS_RELEASE: ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
message("\tMod flag: ${MOD_FLAG}")
message("\tEXENAME: ${EXENAME}")
message("\tPATHOBJHP: ${PATHOBJHP}")

# Project information
project(vpm_py VERSION 1.0 LANGUAGES Fortran)
# Create necessary directories
file(MAKE_DIRECTORY ${PATHOBJHP})

# Set output directories
set(CMAKE_Fortran_MODULE_DIRECTORY ${PATHOBJHP})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PATHOBJHP})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PATHOBJHP})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PATHOBJHP})

set(VPM_LIB_INCL 
    ${SRC_HP}/vpm_mpi.f90
    ${SRC_HP}/vpm_remesh.f90
    ${SRC_HP}/vpm_time.f90
    ${SRC_HP}/vpm_gcalc.f90
)

set(PM_LIB_INCL
    ${SRC_HP}/pmsolve.f90
    ${SRC_HP}/pmbound.f90
    ${SRC_HP}/pinfdomain.f90
)

set(YAPS_INCL
    ${SRC_HP}/yaps2d.f90
    ${SRC_HP}/yaps3d.f90
)

set(TEST_EXE_SRC
    ${SRC_TEST}/test.f90
    ${SRC_TEST}/test_applic.f90
)

add_library(mkl_dfti SHARED ${SRC_HP}/mkl_dfti.f90)
add_library(mpi_matrices SHARED  ${SRC_HP}/mpi_matrices.f90)
add_library(main_pm SHARED  ${SRC_HP}/main_pm.f90)
add_library(pmproject SHARED ${SRC_HP}/pmproject.f90)
target_link_libraries(pmproject PRIVATE mkl_dfti main_pm)

add_library(pmlib SHARED ${SRC_HP}/pmlib.f90)
target_include_directories(pmlib PUBLIC ${PM_LIB_INCL})
target_link_libraries(pmlib PRIVATE mkl_dfti)

add_library(yaps SHARED ${SRC_HP}/yaps.f90)
target_include_directories(yaps PUBLIC ${YAPS_INCL})
target_link_libraries(yaps PRIVATE mkl_dfti mpi_matrices)

add_library(vpm SHARED ${SRC_HP}/vpm.f90)
target_include_directories(vpm PUBLIC ${VPM_LIB_INCL})

# Add all object files as sources to a dummy target
add_custom_target(dummy ALL
    DEPENDS mkl_dfti mpi_matrices yaps pmlib pmproject vpm 
)

# -------------------------------------------------------------------------------------------------
#                                            Library
# -------------------------------------------------------------------------------------------------
add_library(${EXENAME} SHARED ${SRC_HP}/api.f90)
add_dependencies(${EXENAME} dummy) # Ensure proper dependency resolution
target_link_libraries(${EXENAME} PRIVATE mkl_dfti mpi_matrices  main_pm pmproject pmlib yaps vpm) #PUBLIC ${LINK_FLAGS} 
target_include_directories(${EXENAME} PUBLIC ${PATHOBJHP}) # Set include directories using PUBLIC
target_compile_options(${EXENAME} PUBLIC ${FLAGS} ${EXTRA_FLAGS}) # Set compile options using PUBLIC

# -------------------------------------------------------------------------------------------------
#                                            TEST_EXE
# -------------------------------------------------------------------------------------------------
# Add executable for testing purposes (not necessary)
add_executable(${EXENAME}_exe ${TEST_EXE_SRC})
add_dependencies(${EXENAME}_exe dummy)
target_link_libraries(${EXENAME}_exe PRIVATE mkl_dfti mpi_matrices  main_pm pmproject pmlib yaps vpm) #PUBLIC ${LINK_FLAGS}
target_include_directories(${EXENAME}_exe PUBLIC ${PATHOBJHP}) 
target_compile_options(${EXENAME}_exe PUBLIC ${FLAGS} ${EXTRA_FLAGS})

# # -------------------------------------------------------------------------------------------------
#                                           #  F2PY
# # -------------------------------------------------------------------------------------------------
# message(STATUS "${Green}Creating Python module ${ColourReset}")
# # Name of the F2PY module
# set(f2py_module_name "vpm_f2py")
# set(f2py_module_c "${f2py_module_name}_module.c")
# # List of Fortran source files
# set(FORTRAN_SOURCE_FILES
#     ${SRC_HP}/mkl_dfti.f90
#     ${SRC_HP}/main_pm.f90
#     ${SRC_HP}/pmlib.f90
#     ${SRC_HP}/yaps.f90
    
#     ${SRC_HP}/pmproject.f90
#     ${SRC_HP}/vpm_mpi.f90
#     ${SRC_HP}/vpm_time.f90
#     ${SRC_HP}/vpm_gcalc.f90
#     ${SRC_HP}/vpm_remesh.f90
#     ${SRC_HP}/vpm.f90
#     ${SRC_HP}/api.f90
# )
# # Generate sources using f2py
# add_custom_target(
#   genpyf
#   DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
# )
# add_custom_command(
#   OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
#   COMMAND ${Python_EXECUTABLE} -m "numpy.f2py"
#           ${FORTRAN_SOURCE_FILES}
#           -m "${f2py_module_name}"
#           --include-paths "${F2PY_INCLUDE_DIR} ${SRC_HP}/yaps2d.f90 ${SRC_HP}/yaps3d.f90"
#           # --quiet
#           --lower
#   DEPENDS ${FORTRAN_SOURCE_FILES}  # Fortran source files
# )
# # Create a Python module from the generated sources
# Python_add_library(${f2py_module_name} MODULE WITH_SOABI
#   "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"  # Generated C file
#   "${F2PY_INCLUDE_DIR}/fortranobject.c"  # Required from NumPy
#   ${FORTRAN_SOURCE_FILES}  # Fortran source files
# )

# # Link with NumPy and add dependencies
# target_link_libraries(${f2py_module_name} PRIVATE Python::NumPy)
# add_dependencies(${f2py_module_name} genpyf)
# target_include_directories(${f2py_module_name} PRIVATE "${F2PY_INCLUDE_DIR}")

###########################################################
#    Custom targets 
############################################################
# # Installation target
# install(TARGETS ${EXENAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# # Custom targets
# add_custom_target(prepare
#     COMMAND ${CMAKE_COMMAND} -E make_directory ${PATHOBJHP}
#     COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/bin
# )

# add_custom_target(cleanall
#     COMMAND ${CMAKE_COMMAND} -E remove_directory ${SRCHP}/debug
#     COMMAND ${CMAKE_COMMAND} -E remove_directory ${SRCHP}/release
#     COMMAND ${CMAKE_COMMAND} -E remove ${EXENAMEIN}*
# )

# add_custom_target(remake
#     COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
#     COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
# )

# add_custom_target(complete ALL
#     DEPENDS prepare ${EXENAME}
# )

# add_custom_target(uninstall
#     COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/${EXENAMEIN}*
# )

# # Create symbolic links
# set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
# message("The following symbolic link will be created:")
# message("\tCMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
# message("\tCMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
# add_custom_command(TARGET ${EXENAME} POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXENAME} ${CMAKE_INSTALL_PREFIX}/bin/${EXENAME}
# )

# add_dependencies(${EXENAME} prepare)

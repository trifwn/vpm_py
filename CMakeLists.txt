cmake_minimum_required(VERSION 3.22)
message(NOTICE "Configuring VPM")

find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
# Define options
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Add option for static build
option(BUILD_STATIC "Build static executable" ON)
option(CMAKE_INSTALL_LOCAL_ONLY "Install only the local project" ON)
SET(CMAKE_INSTALL_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/vpm_py/shared_libs")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Include custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(FindCompilers)
find_compilers()

project(vpm VERSION 1.0)
      
get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
if (NOT "Fortran" IN_LIST languages)
    enable_language(Fortran)
endif()
if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "Intel|GNU")
  message(WARNING "Unexpected Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

include(VpmSetupSolvers)
# Call the function to find the compilers
if(USE_MKL)
    message(STATUS "Using MKL")
    setup_mkl_solver()
else()
    message(STATUS "Using Fishpack")
    setup_fishpack()
endif()

# Define the installation directories
# If the Install directory is not set, then set it to the current directory
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
endif()
if (CMAKE_INSTALL_LOCAL_ONLY)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
endif()
if (CMAKE_INSTALL_PREFIX STREQUAL "/usr/local")
    set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
endif()
# set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/vpm_py/shared_libs")
set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/")
message(NOTICE "LIB_INSTALL_DIR: ${LIB_INSTALL_DIR}")
message(NOTICE "INCLUDE_INSTALL_DIR: ${INCLUDE_INSTALL_DIR}")

# Setup HDF5
include(SetupHDF5)
setup_hdf5()

# Setup compiler flags
include(VpmCompilerFlags)
setup_vpm_compiler_flags()

# Setting directories
set(SRC_VPM ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(SRC_TEST ${CMAKE_CURRENT_SOURCE_DIR}/test)

# Determine the compiler name (e.g., "Intel", "GNU")
string(TOLOWER "${CMAKE_Fortran_COMPILER_ID}" COMPILER_NAME)
string(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)

set(LIB_DIRECTORY "${CMAKE_BINARY_DIR}/${BUILD_TYPE}/${COMPILER_NAME}/lib/")
set(MOD_DIRECTORY "${CMAKE_BINARY_DIR}/${BUILD_TYPE}/${COMPILER_NAME}/mod/")
set(BIN_DIRECTORY "${CMAKE_BINARY_DIR}/${BUILD_TYPE}/${COMPILER_NAME}/bin/")

make_directory(${LIB_DIRECTORY})
make_directory(${MOD_DIRECTORY})
make_directory(${BIN_DIRECTORY})

# Set output directories
set(CMAKE_Fortran_MODULE_DIRECTORY ${MOD_DIRECTORY})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIRECTORY})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIRECTORY})


set(FISHPACK_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../fishpack)

include(VpmLibraries)
define_vpm_targets()

# -------------------------------------------------------------------------------------------------
#                                            Dependencies
# -------------------------------------------------------------------------------------------------

# Print the configuration summary 
message("-------------------------------------------------------------------------------------------------------------------")
message(NOTICE "VPM Configuration summary:")
message(NOTICE "\tBUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(NOTICE "\tFortran_COMPILER: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")
message(NOTICE "\tFortran_COMPILER_VERSION: ${CMAKE_Fortran_COMPILER_VERSION}")
message(NOTICE "\tFortran_COMPILER_TARGET: ${CMAKE_Fortran_COMPILER_TARGET}")
message(NOTICE "\tFortran_COMPILER_EXTERNAL_TOOLCHAIN: ${CMAKE_Fortran_COMPILER_EXTERNAL_TOOLCHAIN}")
# Print the Fortran flags
message("Compiler flags:")
message(NOTICE "\tCMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
message(NOTICE "\tCMAKE_Fortran_FLAGS_DEBUG: ${Fortran_FLAGS_DEBUG}")
message(NOTICE "\tCMAKE_Fortran_FLAGS_RELEASE: ${Fortran_FLAGS_RELEASE}")
# Print the linker flags
message("Linker flags:")
message(NOTICE "\tCMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
message(NOTICE "\tCMAKE_EXE_LINKER_FLAGS_DEBUG: ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
message(NOTICE "\tCMAKE_EXE_LINKER_FLAGS_RELEASE: ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
# Print the libraries
message("Output directories:")
message(NOTICE "\tCMAKE_Fortran_MODULE_DIRECTORY: ${CMAKE_Fortran_MODULE_DIRECTORY}")
message(NOTICE "\tCMAKE_ARCHIVE_OUTPUT_DIRECTORY: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
message(NOTICE "\tCMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message(NOTICE "\tCMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
# Scalapack flags
message("SOLVER Implementation:")
if(USE_MKL)
    message("\tUsing MKL Implementation")
    message(NOTICE "\tMKL_LINK_FLAGS: ${MKL_LINK_FLAGS}")
else()
    message("\tCompiling and linking with Netlib Fishpack Implementation")
endif()
message("-------------------------------------------------------------------------------------------------------------------")

# Install the libraries to the folder VPM_SOURCE_DIR/../vpm_py/shared_libs
# install(TARGETS operators_api DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/vpm_py/shared_libs)
# install(TARGETS vpm DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/vpm_py/shared_libs)
# install(TARGETS arrays DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/vpm_py/shared_libs)
# install(TARGETS vpm_py_api DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/vpm_py/shared_libs)


# Install the operators_api library
install(TARGETS operators_api
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
    RUNTIME DESTINATION ${LIB_INSTALL_DIR}
    INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}
)

# Install the vpm library
install(TARGETS vpm
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
    RUNTIME DESTINATION ${LIB_INSTALL_DIR}
    INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}
)

# Install the arrays library
install(TARGETS arrays
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
    RUNTIME DESTINATION ${LIB_INSTALL_DIR}
    INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}
)

# Install the vpm_py_api library
install(TARGETS vpm_py_api
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
    RUNTIME DESTINATION ${LIB_INSTALL_DIR}
    INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}
)

install(TARGETS h5fortran
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
    RUNTIME DESTINATION ${LIB_INSTALL_DIR}
    INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}
)
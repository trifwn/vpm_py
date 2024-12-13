function(define_vpm_targets)
    set(VPM_LIB_FILES
        ${SRC_VPM}/vpm.f90
    )

    set(PM_LIB_FILES_MKL
        ${SRC_VPM}/pmlib.f90
        ${SRC_VPM}/pmsolve.f90
        ${SRC_VPM}/pmbound.f90
        ${SRC_VPM}/pinfdomain.f90
    )

    set(PM_LIB_FILES_FISHPACK
        ${SRC_VPM}/pmlib.f90
        ${SRC_VPM}/pmbound.f90
        ${SRC_VPM}/pinfdomain.f90
        ${SRC_VPM}/pmsolve_fish.f90
    )

    set(YAPSLIB_FILES
        ${SRC_VPM}/yaps.f90
        ${SRC_VPM}/yaps2d.f90
        ${SRC_VPM}/yaps3d.f90
    )

    set(TEST_EXE_SRC
        ${SRC_TEST}/test.f90
        ${SRC_TEST}/test_problems.f90
    )

    set(ALL_VPM_OBJECTS_SRC
        # VPM
        ${SRC_VPM}/vpm_vars.f90
        ${SRC_VPM}/vpm_size.f90
        ${SRC_VPM}/vpm_interpolate.f90
        ${SRC_VPM}/vpm_remesh.f90
        # VPM_LIB
        ${SRC_VPM}/vpm.f90
        ${SRC_VPM}/vpm_interpolate.f90
        ${SRC_VPM}/vpm_gcalc.f90
        ${SRC_VPM}/vpm_mpi.f90
        ${SRC_VPM}/vpm_functions.f90
        # PM_LIB
        ${SRC_VPM}/pmlib.f90
        ${SRC_VPM}/pmbound.f90
        ${SRC_VPM}/pinfdomain.f90
        $<$<BOOL:${USE_MKL}>:${SRC_VPM}/pmsolve.f90>
        $<$<NOT:$<BOOL:${USE_MKL}>>:${SRC_VPM}/pmsolve_fish.f90>
        
        #  MKL Headers
        $<$<BOOL:${USE_MKL}>:${SRC_VPM}/mkl_poisson.f90>
        $<$<BOOL:${USE_MKL}>:${SRC_VPM}/mkl_dfti.f90>

        # YAPSLIB
        ${SRC_VPM}/yaps.f90
        ${SRC_VPM}/yaps2d.f90
        ${SRC_VPM}/yaps3d.f90

        # Other VPM files
        ${SRC_VPM}/vpm_types.f90
        ${SRC_VPM}/constants.f90
        ${SRC_VPM}/data_communication.f90
        ${SRC_VPM}/operators_serial.f90
        ${SRC_VPM}/arrays.f90
        ${SRC_VPM}/console_io.f90
        ${SRC_VPM}/file_io.f90
        ${SRC_VPM}/mpi_matrices.f90
        ${SRC_VPM}/parvar.f90
        ${SRC_VPM}/pmgrid.f90
        ${SRC_VPM}/pmproject.f90
    )
    # -------------------------------------------------------------------------------------------------
    #                                            Arrays Library
    # -------------------------------------------------------------------------------------------------
    add_library(arrays SHARED ${SRC_VPM}/arrays.f90)

    # -------------------------------------------------------------------------------------------------
    #                                            VPM Library
    # -------------------------------------------------------------------------------------------------
    add_library(types OBJECT ${SRC_VPM}/vpm_types.f90)
    add_library(constants OBJECT ${SRC_VPM}/constants.f90)
    target_link_libraries(constants PRIVATE types)
    target_link_libraries(arrays PRIVATE types)
    
    add_library(console_io OBJECT ${SRC_VPM}/console_io.f90)
    target_link_libraries(console_io PRIVATE types)

    add_library(pmgrid OBJECT  ${SRC_VPM}/pmgrid.f90)
    target_link_libraries(pmgrid PRIVATE console_io)

    add_library(file_io OBJECT ${SRC_VPM}/file_io.f90)
    target_link_libraries(
        file_io PRIVATE types console_io pmgrid
        h5fortran::h5fortran
    )

    message(STATUS "Using MKl implementation: ${USE_MKL}")
    if (USE_MKL)
        message("\tCompiling with MKL headers")
        add_library(mkl_dfti OBJECT ${SRC_VPM}/mkl_dfti.f90)
        
        add_library(mkl_poisson OBJECT ${SRC_VPM}/mkl_poisson.f90 )
        target_link_libraries(mkl_poisson PRIVATE mkl_dfti)

        add_library(pmlib OBJECT ${PM_LIB_FILES_MKL})
    else()
        message("\tCompiling with Fishpack headers")
        add_library(pmlib OBJECT ${PM_LIB_FILES_FISHPACK})
    endif()

    target_link_libraries(pmlib PRIVATE 
        console_io types constants
        $<$<BOOL:${USE_MKL}>:mkl_poisson>               # Link with MKL if USE_MKL is true
        $<$<NOT:$<BOOL:${USE_MKL}>>:fishpack>           # Link with Fishpack if USE_MKL is false
    )
    target_compile_options(pmlib PRIVATE 
        $<$<BOOL:${USE_MKL}>:${MKL_LINK_FLAGS}>         # Link with MKL
    )

    add_library(mpi_matrices OBJECT ${SRC_VPM}/mpi_matrices.f90)

    add_library(operators_serial OBJECT ${SRC_VPM}/operators_serial.f90)
    target_link_libraries(operators_serial PRIVATE types)
    add_library(data_com OBJECT ${SRC_VPM}/data_communication.f90)

    add_library(parvar OBJECT  ${SRC_VPM}/parvar.f90)
    target_link_libraries(parvar PRIVATE console_io)

    add_library(pmproject OBJECT ${SRC_VPM}/pmproject.f90)

    add_library(yaps OBJECT ${YAPSLIB_FILES})
    target_link_libraries(yaps PRIVATE mpi_matrices pmlib pmproject types console_io file_io)

    add_library(vpm_vars OBJECT ${SRC_VPM}/vpm_vars.f90)
    target_link_libraries(vpm_vars PRIVATE types console_io)

    add_library(vpm_size OBJECT ${SRC_VPM}/vpm_size.f90)
    target_link_libraries(vpm_size PRIVATE types console_io pmlib)

    add_library(vpm_interpolate OBJECT ${SRC_VPM}/vpm_interpolate.f90)
    target_link_libraries(vpm_interpolate PRIVATE types console_io pmproject vpm_size vpm_vars pmgrid parvar)

    add_library(vpm_gcalc OBJECT ${SRC_VPM}/vpm_gcalc.f90)
    target_link_libraries(vpm_gcalc PRIVATE types constants console_io vpm_size vpm_vars)

    add_library(vpm_mpi OBJECT ${SRC_VPM}/vpm_mpi.f90)
    target_link_libraries(vpm_mpi PRIVATE types console_io vpm_size vpm_vars pmgrid)

    add_library(vpm_functions OBJECT ${SRC_VPM}/vpm_functions.f90)
    target_link_libraries(vpm_functions PRIVATE 
        types constants console_io 
        vpm_vars vpm_size vpm_mpi vpm_interpolate  
        pmgrid pmproject parvar pmlib yaps operators_serial 
    )

    add_library(vpm_remesh OBJECT ${SRC_VPM}/vpm_remesh.f90)
    target_link_libraries(vpm_remesh PRIVATE types console_io parvar pmgrid vpm_interpolate vpm_vars vpm_functions)


    add_library(vpm_lib OBJECT ${VPM_LIB_FILES}) 
    target_link_libraries(vpm_lib 
    PUBLIC 
        types constants                     
        console_io file_io 
        operators_serial 

        vpm_size vpm_vars vpm_interpolate vpm_gcalc 
        vpm_mpi vpm_remesh vpm_functions

        parvar pmgrid  
        yaps pmlib pmproject 
        # mpi_matrices 
        # $<$<BOOL:${USE_MKL}>:mkl_poisson>                       # Link with MKL if USE_MKL is true
        # $<$<NOT:$<BOOL:${USE_MKL}>>:fishpack>                   # Link with Fishpack
        # h5fortran::h5fortran
    )

    if (BUILD_STATIC)
        add_library(vpm STATIC) 
    else()
        add_library(vpm SHARED)
    endif()
    target_link_libraries(vpm 
    PUBLIC 
        vpm_lib arrays yaps pmlib pmproject parvar pmgrid  
        vpm_size vpm_vars vpm_interpolate vpm_gcalc vpm_mpi vpm_remesh vpm_functions
        console_io file_io types constants 
        mpi_matrices operators_serial 
    PRIVATE 
        $<$<BOOL:${USE_MKL}>:mkl_poisson>                       # Link with MKL if USE_MKL is true
        $<$<NOT:$<BOOL:${USE_MKL}>>:fishpack>                   # Link with Fishpack
        h5fortran::h5fortran
    )
    target_include_directories(vpm PUBLIC ${LIB_DIRECTORY}) # Set include directories using PUBLIC
    target_compile_options(vpm PRIVATE 
        $<$<AND:$<BOOL:${BUILD_STATIC}>,$<BOOL:${USE_INTEL_COMPILER}>>:
                -static -static-intel
        >
        $<$<AND:$<BOOL:${BUILD_STATIC}>,$<NOT:$<BOOL:${USE_INTEL_COMPILER}>>>:
                -static -static-libgfortran -static-libgcc -static-libstdc++
        >
    )

    # -------------------------------------------------------------------------------------------------
    #                                            API Serial Vector Field Operators
    # -------------------------------------------------------------------------------------------------
    add_library(operators_api SHARED ${SRC_VPM}/api_operators_serial.f90)
    add_dependencies(operators_serial arrays types) # Ensure proper dependency resolution
    target_link_libraries(operators_api PRIVATE arrays types operators_serial)

    # -------------------------------------------------------------------------------------------------
    #                                            API Shared Library
    # -------------------------------------------------------------------------------------------------
    # add_library(${EXE} SHARED ${SRC_VPM}/api.f90)
    # add_dependencies(${EXE} vpm arrays operators_api) # Ensure proper dependency resolution
    # target_link_libraries(${EXE} PRIVATE vpm)

    python_add_library(vpm_py_api ${SRC_VPM}/api.f90)
    add_dependencies(vpm_py_api vpm arrays operators_api) # Ensure proper dependency resolution
    target_link_libraries(vpm_py_api PRIVATE vpm)
    set_target_properties(vpm_py_api PROPERTIES
        OUTPUT_NAME "vpm_py_api"
        # LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vpm_py_api
    )
    # -------------------------------------------------------------------------------------------------
    #                                            VPM Executable
    # -------------------------------------------------------------------------------------------------
    # We need to add the link options for the MKL library
    add_executable(${EXE}_exe ${TEST_EXE_SRC})
    target_link_libraries(${EXE}_exe PRIVATE 
        vpm
        h5fortran::h5fortran
    )
    target_link_options(${EXE}_exe PRIVATE 
        # $<$<BOOL:${BUILD_STATIC}>:-static>                      # Link statically
        $<$<BOOL:${USE_MKL}>:${MKL_LINK_FLAGS}>                   # Link MKL
    )

    # # -------------------------------------------------------------------------------------------------
    #                                           Arrays Test Executable 
    # # -------------------------------------------------------------------------------------------------
    add_executable(test_arrays_exe ${SRC_TEST}/test_arrays.f90)
    add_dependencies(test_arrays_exe arrays)
    target_link_libraries(test_arrays_exe PRIVATE arrays)

    # -------------------------------------------------------------------------------------------------
    #                                           Operator Test Executable
    # -------------------------------------------------------------------------------------------------

    # add_executable(test_operators ${SRC_TEST}/test_operators.f90)
    # target_link_libraries(test_operators PRIVATE data_com operators_serial)

    # -------------------------------------------------------------------------------------------------
    #                                          Compiler Flags
    # -------------------------------------------------------------------------------------------------

    # Set Fortran compiler flags for all targets
    set_compiler_flags(arrays)
    set_compiler_flags(types)
    set_compiler_flags(constants)
    set_compiler_flags(console_io)
    set_compiler_flags(pmgrid)
    set_compiler_flags(file_io)
    if(USE_MKL)
        set_compiler_flags(mkl_poisson)
    endif()
    set_compiler_flags(pmlib)
    set_compiler_flags(mpi_matrices)
    set_compiler_flags(operators_serial)
    # set_compiler_flags(data_com)
    set_compiler_flags(parvar)
    set_compiler_flags(pmproject)
    set_compiler_flags(yaps)
    set_compiler_flags(vpm_vars)
    set_compiler_flags(vpm_size)
    set_compiler_flags(vpm_interpolate)
    set_compiler_flags(vpm_gcalc)
    set_compiler_flags(vpm_mpi)
    set_compiler_flags(vpm_functions)
    set_compiler_flags(vpm_remesh)
    set_compiler_flags(vpm_lib)
    set_compiler_flags(vpm)
    set_compiler_flags(vpm_py_api)
    set_compiler_flags(${EXE}_exe)
    set_compiler_flags(test_arrays_exe)
    # set_compiler_flags(test_operators)
endfunction()
function(define_vpm_targets)
    set(VPM_LIB_FILES
        ${SRC_VPM}/vpm.f90
        ${SRC_VPM}/vpm_remesh.f90
        ${SRC_VPM}/vpm_interpolate.f90
        ${SRC_VPM}/vpm_gcalc.f90
        ${SRC_VPM}/vpm_mpi.f90
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
        # VPM_LIB
        ${SRC_VPM}/vpm.f90
        ${SRC_VPM}/vpm_remesh.f90
        ${SRC_VPM}/vpm_interpolate.f90
        ${SRC_VPM}/vpm_gcalc.f90
        ${SRC_VPM}/vpm_mpi.f90
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
        ${SRC_VPM}/base_types.f90
        ${SRC_VPM}/constants.f90
        ${SRC_VPM}/arrays.f90
        ${SRC_VPM}/io.f90
        ${SRC_VPM}/mpi_matrices.f90
        ${SRC_VPM}/parvar.f90
        ${SRC_VPM}/pmgrid.f90
        ${SRC_VPM}/pmeshpar.f90
        ${SRC_VPM}/pmproject.f90
    )
    # -------------------------------------------------------------------------------------------------
    #                                            Arrays Library
    # -------------------------------------------------------------------------------------------------
    add_library(arrays SHARED ${SRC_VPM}/arrays.f90)

    # -------------------------------------------------------------------------------------------------
    #                                            VPM Library
    # -------------------------------------------------------------------------------------------------
    add_library(constants OBJECT ${SRC_VPM}/base_types.f90)
    add_library(types OBJECT ${SRC_VPM}/constants.f90)
    target_link_libraries(types PRIVATE constants)
    target_link_libraries(arrays PRIVATE types)
    
    add_library(io OBJECT ${SRC_VPM}/io.f90)
    target_link_libraries(io PRIVATE types)

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
        io types constants
        $<$<BOOL:${USE_MKL}>:mkl_poisson>               # Link with MKL if USE_MKL is true
        $<$<NOT:$<BOOL:${USE_MKL}>>:fishpack>           # Link with Fishpack if USE_MKL is false
    )
    target_compile_options(pmlib PRIVATE 
        $<$<BOOL:${USE_MKL}>:${MKL_LINK_FLAGS}>         # Link with MKL
    )

    add_library(mpi_matrices OBJECT ${SRC_VPM}/mpi_matrices.f90)

    add_library(parvar OBJECT  ${SRC_VPM}/parvar.f90)
    target_link_libraries(parvar PRIVATE io)

    add_library(pmgrid OBJECT  ${SRC_VPM}/pmgrid.f90)
    target_link_libraries(pmgrid PRIVATE io)

    add_library(pmeshpar OBJECT  ${SRC_VPM}/pmeshpar.f90)
    target_link_libraries(pmeshpar PRIVATE io constants)

    add_library(pmproject OBJECT ${SRC_VPM}/pmproject.f90)

    add_library(yaps OBJECT ${YAPSLIB_FILES})
    target_link_libraries(yaps PRIVATE mpi_matrices pmlib pmproject types io)

    add_library(vpm_vars OBJECT ${SRC_VPM}/vpm_vars.f90)
    target_link_libraries(vpm_vars PRIVATE types io)

    add_library(vpm_size OBJECT ${SRC_VPM}/vpm_size.f90)
    target_link_libraries(vpm_size PRIVATE types io)

    add_library(vpm_lib OBJECT ${VPM_LIB_FILES}) 
    target_link_libraries(vpm_lib PRIVATE 
        yaps pmlib pmproject parvar pmeshpar pmgrid  
        vpm_size vpm_vars io types constants mpi_matrices
        $<$<BOOL:${USE_MKL}>:mkl_poisson>                       # Link with MKL if USE_MKL is true
        $<$<NOT:$<BOOL:${USE_MKL}>>:fishpack>                   # Link with Fishpack
        h5fortran::h5fortran
    )

    if (BUILD_STATIC)
        add_library(vpm STATIC ${ALL_VPM_OBJECTS_SRC}) 
    else()
        add_library(vpm SHARED ${ALL_VPM_OBJECTS_SRC})
    endif()
    
    target_link_libraries(vpm PRIVATE 
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
    #                                            API Shared Library
    # -------------------------------------------------------------------------------------------------
    add_library(${EXE} SHARED ${SRC_VPM}/api.f90)
    add_dependencies(${EXE} vpm arrays) # Ensure proper dependency resolution
    target_link_libraries(${EXE} PRIVATE vpm)

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
        $<$<BOOL:${USE_MKL}>:${MKL_LINK_FLAGS}>                 # Link MKL
    )

    # # -------------------------------------------------------------------------------------------------
    #                                           Arrays Test Executable 
    # # -------------------------------------------------------------------------------------------------
    add_executable(test_arrays_exe ${SRC_TEST}/test_arrays.f90)
    add_dependencies(test_arrays_exe arrays)
    set_compiler_flags(test_arrays_exe)
    target_link_libraries(test_arrays_exe PRIVATE arrays)

    # Set Fortran compiler flags for all targets
    set_compiler_flags(arrays)
    set_compiler_flags(constants)
    set_compiler_flags(types)
    set_compiler_flags(io)
    set_compiler_flags(pmlib)
    set_compiler_flags(mpi_matrices)
    set_compiler_flags(parvar)
    set_compiler_flags(pmgrid)
    set_compiler_flags(pmeshpar)
    set_compiler_flags(pmproject)
    set_compiler_flags(yaps)
    set_compiler_flags(vpm_vars)
    set_compiler_flags(vpm_size)
    set_compiler_flags(vpm)
    set_compiler_flags(vpm_lib)
    set_compiler_flags(${EXE})
    set_compiler_flags(${EXE}_exe)
    set_compiler_flags(test_arrays_exe)
    if(USE_MKL)
        set_compiler_flags(mkl_poisson)
    endif()
endfunction()
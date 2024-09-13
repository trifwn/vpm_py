function(setup_mkl_solver)
    if(NOT DEFINED ENV{MKLROOT})
        message(FATAL_ERROR "MKLROOT is not defined. Please set the MKLROOT environment variable.")
    endif()
    
    # find_package(MKL REQUIRED)
    set(MKLROOT $ENV{MKLROOT})
    if (BUILD_STATIC)
        set(MKL_LINK_FLAGS
            -Wl,--start-group 
            ${MKLROOT}/lib/libmkl_intel_lp64.a 
            ${MKLROOT}/lib/libmkl_intel_thread.a 
            ${MKLROOT}/lib/libmkl_core.a 
            -Wl,--end-group 
            -lpthread -lm -ldl -liomp5
            PARENT_SCOPE
        )
    else()
        set(MKL_LINK_FLAGS
            -L${MKLROOT}/lib 
            -lmkl_intel_lp64 
            -lmkl_intel_thread 
            -lmkl_core 
            -liomp5 -lpthread -lm -ldl
            PARENT_SCOPE
        )
    endif()
endfunction()

function(setup_fishpack)
    # set(FISHPACK_SRC_DIR ${CMAKE_SOURCE_DIR}/third_party/fishpack)
    set(FISHPACK_SRC_DIR ${CMAKE_SOURCE_DIR}/../fishpack)
    # We need to find the fishpack library and link it
    find_library(FISHPACK_LIB libfishpack HINTS ${FISHPACK_SRC_DIR})
    if (FISHPACK_LIB)
        message(STATUS "Found fishpack library: ${FISHPACK_LIB}")
    else()
        # Search for the directory of the fishpack source code and compile it
        message(STATUS "${Green}Fishpack library not found. Compiling from source${ColourReset}")
        # Define unique output directories for fishpack
        set(FISHPACK_BUILD_DIR ${CMAKE_BINARY_DIR}/fishpack_build)

        # Add the fishpack subdirectory
        add_subdirectory(${FISHPACK_SRC_DIR} fishpack_build)
    endif()
endfunction()
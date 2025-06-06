function(setup_vpm_compiler_flags)
    if(USE_INTEL_COMPILER) 
        set(COMPILE_FLAGS  -DASCII=1 -DHAVE_OMP -fPIC)
        set(LINK_FLAGS "-qopenmp -qmkl -fPIC")

        set(Fortran_FLAGS_DEBUG ${COMPILE_FLAGS} 
                    -g -traceback -fpe0 
                    -O0 ${MKL_FLAG} 
                    -I${MKLROOT}/include 
                    -L${MKLROOT}/lib/intel64
                    -check all,nouninit
                    -init=snan -init=arrays
                    -ftrapuv
                    -debug full
                    -diag-enable=all
                    -implicitnone
                    -fstack-protector-all
                    -warn all
                    -lunwind
                    PARENT_SCOPE
        )
        set(Fortran_FLAGS_RELEASE ${COMPILE_FLAGS} 
                    -O3 ${MKL_FLAG} 
                    -I${MKLROOT}/include 
                    -L${MKLROOT}/lib/intel64
                    -ffast-math 
                    -march=native
                    PARENT_SCOPE
        )

        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${LINK_FLAGS} ${MKL_FLAG} -O0 -g -traceback -fpe0 -check all,nouninit" PARENT_SCOPE)        
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${LINK_FLAGS} ${MKL_FLAG} -O3 -march=native -flto" PARENT_SCOPE)        
    else()
        set(COMPILE_FLAGS -ffree-line-length-512 -DASCII=1 -llapack -lblas -lmpi -fPIC -lhdf5)
        set(LINK_FLAGS "-fopenmp -fPIC")
        
        set(Fortran_FLAGS_DEBUG ${COMPILE_FLAGS} 
                    -O0 -g -fbacktrace 
                    -finit-real=nan 
                    -finit-integer=nan -fcheck=all 
                    -ffpe-trap=invalid,overflow,underflow
                    -fimplicit-none 
                    -fstack-protector-all
                    -Wall -Wextra -Wconversion -Warray-bounds -Wuninitialized
                    PARENT_SCOPE
        )
        set(Fortran_FLAGS_RELEASE ${COMPILE_FLAGS} -O3 -ffast-math -march=native PARENT_SCOPE)
        
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${LINK_FLAGS} -O0 -g -fbacktrace -fcheck=all -finit-real=nan -finit-integer=nan -fPIC" PARENT_SCOPE)        
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${LINK_FLAGS} -O3 -flto -march=native" PARENT_SCOPE)        
    endif()
    message(STATUS "Set Fortran compiler flags")
endfunction()

# Function to set compiler flags for a target
function(set_compiler_flags target)
    message(STATUS "Setting compiler flags for ${target}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(${target} PRIVATE ${Fortran_FLAGS_DEBUG})
    else()
        target_compile_options(${target} PRIVATE ${Fortran_FLAGS_RELEASE})
    endif()
endfunction()

function(set_preprocessor_flag target)
    # Set the -cpp or -fpp flag for the target
    if (USE_INTEL_COMPILER)
        target_compile_options(${target} PRIVATE -fpp)
    else()
        target_compile_options(${target} PRIVATE -cpp)
    endif() 
endfunction()
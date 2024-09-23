function(setup_hdf5)
    cmake_minimum_required(VERSION 3.14)

    project(H5fortran_example LANGUAGES C Fortran)

    set(CMAKE_POSITION_INDEPENDENT_CODE ON) # -fPIC
    enable_testing()
    # find_package(h5fortran)

        include(FetchContent)

        option(h5fortran_BUILD_TESTING "h5fortran internal tests")
        # Add the -fPIC flag to the Fortran compiler

        set(FETCHCONTENT_UPDATES_DISCONNECTED true)

        FetchContent_Declare(h5fortran
            GIT_REPOSITORY https://github.com/geospace-code/h5fortran.git
            GIT_TAG v4.10.2
            TLS_VERIFY true
            FIND_PACKAGE_ARGS
        )
        FetchContent_MakeAvailable(H5FORTRAN)

        # Get the include directories for h5fortran
        get_target_property(H5FORTRAN_INCLUDES h5fortran INTERFACE_INCLUDE_DIRECTORIES)

        # Extract the BUILD_INTERFACE directory from the H5FORTRAN_INCLUDES
        string(REPLACE ";" "\n" H5FORTRAN_INCLUDES_LIST "${H5FORTRAN_INCLUDES}")
        foreach(INCLUDE_PATH ${H5FORTRAN_INCLUDES_LIST})
            if(INCLUDE_PATH MATCHES "\\$<BUILD_INTERFACE:(.*)>")
                set(BUILD_INTERFACE_DIR ${CMAKE_MATCH_1})
                # remove '<' and '>' from the path
                string(REGEX REPLACE "<|>" "" BUILD_INTERFACE_DIR ${BUILD_INTERFACE_DIR})
                # Search for $ is found remove everything after it including $
                string(REGEX REPLACE "\\$.*" "" BUILD_INTERFACE_DIR ${BUILD_INTERFACE_DIR})

                # Create the build interface directory
                file(MAKE_DIRECTORY ${BUILD_INTERFACE_DIR})
                message(STATUS "Created build interface directory: ${BUILD_INTERFACE_DIR}")
            endif()
        endforeach() 
       
endfunction()
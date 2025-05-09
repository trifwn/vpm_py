function(setup_mudpack_sp)
    cmake_minimum_required(VERSION 3.14)
    # find_package(mudpack_sp)

    include(FetchContent) 
    set(FETCHCONTENT_UPDATES_DISCONNECTED true)

    message(STATUS "Fetching mudpack_sp")
    # Declare the dependency
    FetchContent_Declare(mudpack_sp
        GIT_REPOSITORY https://github.com/trifwn/mudpack_sp.git
        GIT_TAG main  # or specify a version tag/commit
        TLS_VERIFY true
    )
    # FetchContent_MakeAvailable(mudpack_sp)
    FetchContent_GetProperties(mudpack_sp)
    if(NOT mudpack_sp_POPULATED)
        FetchContent_Populate(mudpack_sp)
        add_subdirectory(${mudpack_sp_SOURCE_DIR} ${mudpack_sp_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()

    # Get the include directories for h5fortran
    get_target_property(MUDPACK_SP_INCLUDES mudpack_sp INTERFACE_INCLUDE_DIRECTORIES)


    # Print the targets in the mudpack_sp library
    get_property(mudpack_sp_targets DIRECTORY ${mudpack_sp_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
    # Print the targets
endfunction()
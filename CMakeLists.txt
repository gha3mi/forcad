cmake_minimum_required(VERSION 3.18)

# Metadata
file(READ "${CMAKE_SOURCE_DIR}/VERSION" PROJECT_VERSION_STRING)
string(STRIP "${PROJECT_VERSION_STRING}" PROJECT_VERSION)
set(PROJECT_VERSION ${PROJECT_VERSION})
set(PROJECT_AUTHOR "Seyed Ali Ghasemi")
set(PROJECT_MAINTAINER "info@gha3mi.com")
set(PROJECT_LICENSE "BSD-3-Clause")
project(forcad VERSION ${PROJECT_VERSION} LANGUAGES Fortran)

# Fortran settings
enable_language(Fortran)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Options
option(FORCAD_BUILD_EXAMPLES "Build forcad examples" OFF)
option(FORCAD_BUILD_TESTS "Build forcad tests" OFF)

# OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_Fortran_FOUND)
  message(STATUS "OpenMP found")
endif()

include(FetchContent)

# forIGES dependency
set(BUILD_TEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  forIGES
  GIT_REPOSITORY https://github.com/rweed/forIGES.git
  GIT_TAG main
)
FetchContent_MakeAvailable(forIGES)
set(FORIGES_MOD_DIR ${CMAKE_Fortran_MODULE_DIRECTORY}/forIGES)
set_target_properties(forIGES PROPERTIES
  Fortran_MODULE_DIRECTORY ${FORIGES_MOD_DIR}
)

# forunittest (only used in tests and examples)
if(FORCAD_BUILD_TESTS OR FORCAD_BUILD_EXAMPLES)
  FetchContent_Declare(
    forunittest
    GIT_REPOSITORY https://github.com/gha3mi/forunittest.git
    GIT_TAG main
  )
  FetchContent_MakeAvailable(forunittest)
  set(FORUNITTEST_MOD_DIR ${CMAKE_Fortran_MODULE_DIRECTORY}/forunittest)
  set_target_properties(forunittest PROPERTIES
    Fortran_MODULE_DIRECTORY ${FORUNITTEST_MOD_DIR}
  )
endif()

# Source files
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS
  src/*.f90
  src/*.F90
)

add_library(forcad ${SOURCES})
set_target_properties(forcad PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
add_dependencies(forcad forIGES)
target_link_libraries(forcad PUBLIC forIGES)

target_include_directories(forcad PUBLIC
  ${CMAKE_Fortran_MODULE_DIRECTORY}
  ${FORIGES_MOD_DIR}
  ${FORUNITTEST_MOD_DIR}
)

# Link OpenMP if found
if(OpenMP_Fortran_FOUND)
  target_link_libraries(forcad PUBLIC OpenMP::OpenMP_Fortran)
endif()

# Examples (optional)
if(FORCAD_BUILD_EXAMPLES)
  # Example-only dependencies
  FetchContent_Declare(
    forimage
    GIT_REPOSITORY https://github.com/gha3mi/forimage.git
    GIT_TAG main
  )
  FetchContent_MakeAvailable(forimage)

  FetchContent_Declare(
    forcolormap
    GIT_REPOSITORY https://github.com/vmagnin/forcolormap.git
    GIT_TAG main
  )
  FetchContent_MakeAvailable(forcolormap)

  FetchContent_Declare(
    fortime
    GIT_REPOSITORY https://github.com/gha3mi/fortime.git
    GIT_TAG main
  )
  FetchContent_MakeAvailable(fortime)

  file(GLOB EXAMPLES CONFIGURE_DEPENDS example/*.f90 example/*.F90)
  foreach(example_file ${EXAMPLES})
    get_filename_component(example_name ${example_file} NAME_WE)
    add_executable(${example_name} ${example_file})
    target_link_libraries(${example_name} PRIVATE forcad)

    # Link example-specific dependencies
    if(${example_name} MATCHES "example_ppm[1-3]")
      target_link_libraries(${example_name} PRIVATE forimage forcolormap fortime)
    elseif(${example_name} STREQUAL "nearest_point_2d_bench")
      target_link_libraries(${example_name} PRIVATE fortime)
    endif()

    target_include_directories(${example_name} PRIVATE
      ${CMAKE_Fortran_MODULE_DIRECTORY}
      ${FORIGES_MOD_DIR}
    )

    if(OpenMP_Fortran_FOUND)
      target_link_libraries(${example_name} PRIVATE OpenMP::OpenMP_Fortran)
    endif()
  endforeach()
endif()

# Tests (optional)
if(FORCAD_BUILD_TESTS)
  enable_testing()
  file(GLOB TESTS CONFIGURE_DEPENDS test/*.f90 test/*.F90)
  foreach(test_file ${TESTS})
    get_filename_component(test_name ${test_file} NAME_WE)
    add_executable(${test_name}_exec ${test_file})
    target_link_libraries(${test_name}_exec PRIVATE forcad forunittest)

    add_test(NAME ${test_name}
      COMMAND ${test_name}_exec
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

    target_include_directories(${test_name}_exec PRIVATE
      ${CMAKE_Fortran_MODULE_DIRECTORY}
      ${FORIGES_MOD_DIR}
      ${FORUNITTEST_MOD_DIR}
    )

    if(OpenMP_Fortran_FOUND)
      target_link_libraries(${test_name}_exec PRIVATE OpenMP::OpenMP_Fortran)
    endif()
  endforeach()
endif()

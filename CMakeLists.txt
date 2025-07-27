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

# Set compiler- and OS-specific flags
message(STATUS "Detected Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "Detected OS: ${CMAKE_SYSTEM_NAME}")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -cpp -Wall -Wextra -fmax-errors=1 -g -fcheck=bounds -fcheck=array-temps -fbacktrace -fcoarray=single -fPIC -fimplicit-none -Werror=implicit-interface -ffree-form -Wno-maybe-uninitialized")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -cpp -Wall -Wextra -fmax-errors=1 -g -fcheck=bounds -fcheck=array-temps -fbacktrace -fcoarray=single -fPIC -fimplicit-none -Werror=implicit-interface -ffree-form -Wno-maybe-uninitialized")
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -cpp -Wall -Wextra -fmax-errors=1 -g -fcheck=bounds -fcheck=array-temps -fbacktrace -fcoarray=single -fPIC -fimplicit-none -Werror=implicit-interface -ffree-form -Wno-maybe-uninitialized")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -cpp -Wall -Wextra -fmax-errors=1 -g -fcheck=bounds -fcheck=array-temps -fbacktrace -fcoarray=single -fPIC -fimplicit-none -Werror=implicit-interface -ffree-form -Wno-maybe-uninitialized")
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -cpp -Wall -Wextra -fmax-errors=1 -g -fcheck=bounds -fcheck=array-temps -fbacktrace -fcoarray=single -fPIC -fimplicit-none -Werror=implicit-interface -ffree-form -Wno-maybe-uninitialized")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -cpp -Wall -Wextra -fmax-errors=1 -g -fcheck=bounds -fcheck=array-temps -fbacktrace -fcoarray=single -fPIC -fimplicit-none -Werror=implicit-interface -ffree-form -Wno-maybe-uninitialized")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fpp -diag-error 10006 -check all,nouninit -error-limit 1 -O0 -g -assume byterecl -traceback -fPIC -free")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fpp -diag-error 10006 -check all,nouninit -error-limit 1 -O0 -g -assume byterecl -traceback -fPIC -free")
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /fpp /Qdiag-error:10006 /check:all /check:nouninit /error-limit:1 /Od /Z7 /debug:full /assume:byterecl /traceback /free /heap-arrays")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} /fpp /Qdiag-error:10006 /check:all /check:nouninit /error-limit:1 /Od /Z7 /debug:full /assume:byterecl /traceback /free /heap-arrays")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -cpp -Minform=inform -Mbackslash -Mbounds -Mchkptr -Mchkstk -traceback -fPIC -mp -Mfree -stdpar=multicore,gpu -Minfo=stdpar,accel")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -cpp -Minform=inform -Mbackslash -Mbounds -Mchkptr -Mchkstk -traceback -fPIC -mp -Mfree -stdpar=multicore,gpu -Minfo=stdpar,accel")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -cpp -fPIC -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -cpp -fPIC -g")
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -cpp -fPIC -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -cpp -fPIC -g")
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -cpp -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -cpp -g")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LFortran")
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} --cpp -fPIC --verbose -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} --cpp -fPIC --verbose -g")
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} --cpp -fPIC --verbose -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} --cpp -fPIC --verbose -g")
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} --cpp -fPIC --verbose -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} --cpp -fPIC --verbose -g")
  endif()
endif()

# Options
option(FORCAD_BUILD_EXAMPLES "Build forcad examples" OFF)
option(FORCAD_BUILD_TESTS "Build forcad tests" OFF)

# OpenMP
find_package(OpenMP)
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

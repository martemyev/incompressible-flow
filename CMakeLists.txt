cmake_minimum_required(VERSION 2.8)

project(incompressible-flow CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build: Debug | Release | PDebug | PRelease" FORCE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)

#-------------------------------------------------------------------------------
# MFEM
#-------------------------------------------------------------------------------
set(MFEM_DIR "/u/artemyev/projects/mfem" CACHE PATH "Path to the MFEM root directory")
if(NOT MFEM_BUILD)
  if(BUILD_TYPE STREQUAL DEBUG)
    set(MFEM_BUILD "serialdebug" CACHE STRING
        "Name of MFEM build (it's a subdirectory with config file and library")
  elseif(BUILD_TYPE STREQUAL RELEASE)
    set(MFEM_BUILD "serial")
  elseif(BUILD_TYPE STREQUAL PDEBUG)
    set(MFEM_BUILD "paralleldebug")
  elseif(BUILD_TYPE STREQUAL PRELEASE)
    set(MFEM_BUILD "parallel")
  else()
    message(FATAL_ERROR "Unknown BUILD_TYPE")
  endif()
endif()
include_directories("${MFEM_DIR}")
include_directories("${MFEM_DIR}/${MFEM_BUILD}")
link_directories("${MFEM_DIR}/${MFEM_BUILD}")

#-------------------------------------------------------------------------------
# Parallel mode
#-------------------------------------------------------------------------------
if(BUILD_TYPE STREQUAL PDEBUG OR BUILD_TYPE STREQUAL PRELEASE)
  # MPI
  find_package(MPI REQUIRED)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  # Hypre
  set(HYPRE_DIR "/u/artemyev/software/hypre-2.10.1/src" CACHE PATH
      "Path to the Hypre root directory")
  include_directories(${HYPRE_DIR}/hypre/include)
  if(BUILD_TYPE STREQUAL PDEBUG)
    link_directories(${HYPRE_DIR}/debugbuild)
  else()
    link_directories(${HYPRE_DIR}/cmbuild)
  endif()
  # Metis
  set(METIS_DIR "/u/artemyev/software/metis-5.1.0/build/Linux-x86_64" CACHE PATH
      "Path to the Metis root directory")
  link_directories(${METIS_DIR}/lib)
endif()

file(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(GLOB HDR_LIST "${PROJECT_SOURCE_DIR}/src/*.hpp")

#-------------------------------------------------------------------------------
# Config info
#-------------------------------------------------------------------------------
string(TIMESTAMP CONFIG_TIME "%Y-%m-%d at %H:%M:%S Z" UTC)
if(WIN32)
  set(USER_NAME $ENV{USERNAME})
else()
  set(USER_NAME $ENV{USER})
endif()
site_name(HOST_NAME)
# Get the current working Git branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

#-------------------------------------------------------------------------------
# Config file
#-------------------------------------------------------------------------------
configure_file(
  "${PROJECT_SOURCE_DIR}/config.hpp.in"
  "${PROJECT_SOURCE_DIR}/src/config.hpp")

#-------------------------------------------------------------------------------
# Config file
#-------------------------------------------------------------------------------
unset(MFEM_BUILD CACHE)

#-------------------------------------------------------------------------------
# Executable
#-------------------------------------------------------------------------------
add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
target_link_libraries(${PROJECT_NAME} mfem rt)
if(BUILD_TYPE STREQUAL PDEBUG OR BUILD_TYPE STREQUAL PRELEASE)
  target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES} HYPRE metis)
endif()

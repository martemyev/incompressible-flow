cmake_minimum_required(VERSION 2.8)

project(incompressible-flow)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build: Debug | Release | PDebug | PRelease" FORCE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
message(STATUS "BUILD_TYPE ${BUILD_TYPE}")

#-------------------------------------------------------------------------------
# MFEM
#-------------------------------------------------------------------------------
find_path(MFEM_DIR include/mfem.hpp HINTS ENV MFEM_DIR PATHS $ENV{HOME}/projects/mfem DOC "MFEM directory")
set(MFEM_INCLUDE_DIR ${MFEM_DIR}/include)

unset(MFEM_CONFIG CACHE)
find_path(MFEM_CONFIG config.hpp HINTS ENV MFEM_DIR PATHS ${MFEM_INCLUDE_DIR}/config DOC "MFEM config directory")

list(APPEND MFEM_INCLUDE_DIR ${MFEM_CONFIG})
include_directories(${MFEM_INCLUDE_DIR})

message(STATUS "MFEM_DIR ${MFEM_DIR}")
message(STATUS "MFEM_CONFIG ${MFEM_CONFIG}")
message(STATUS "MFEM_INCLUDE_DIR ${MFEM_INCLUDE_DIR}")

unset(MFEM_LIBRARY CACHE)
find_library(MFEM_LIBRARY mfem HINTS ${MFEM_DIR} PATHS ${MFEM_DIR}/lib PATH_SUFFIXES lib DOC "MFEM library")
if(NOT MFEM_LIBRARY)
  message(FATAL_ERROR "MFEM library was NOT found")
else()
  message(STATUS "MFEM lib ${MFEM_LIBRARY}")
endif()

#-------------------------------------------------------------------------------
# Parallel mode
#-------------------------------------------------------------------------------
if(BUILD_TYPE STREQUAL PDEBUG OR BUILD_TYPE STREQUAL PRELEASE)
  # MPI
  find_package(MPI REQUIRED)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  # Hypre
  find_path(HYPRE_DIR include/HYPRE.h HINTS ENV HYPRE_DIR PATHS "$ENV{HOME}/hypre" DOC "Hypre directory")
  if(EXISTS ${HYPRE_DIR}/include/HYPRE.h)
    set(HYPRE_FOUND YES)
    message(STATUS "Hypre found ${HYPRE_DIR}")
    set(HYPRE_INCLUDES ${HYPRE_DIR})
    find_path(HYPRE_INCLUDE_DIR HYPRE.h HINTS "${HYPRE_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
    list(APPEND HYPRE_INCLUDES ${HYPRE_INCLUDE_DIR})
    include_directories(${HYPRE_INCLUDES})
    find_library(HYPRE_LIBRARY NAMES hypre HYPRE HINTS "${HYPRE_DIR}/lib" DOC "Hypre library")
    message(STATUS "Hypre lib ${HYPRE_LIBRARY}")
    get_filename_component(HYPRE_LIB_DIR ${HYPRE_LIBRARY} DIRECTORY)
    message(STATUS "HYPRE_LIB_DIR ${HYPRE_LIB_DIR}")
    link_directories(${HYPRE_LIB_DIR})
  else()
    set(HYPRE_FOUND NO)
    message(FATAL_ERROR "Hypre is NOT found")
  endif()
  # Metis
  find_path(METIS_DIR include/metis.h HINTS ENV METIS_DIR PATHS "$ENV{HOME}/metis" DOC "Metis directory")
  if(EXISTS ${METIS_DIR}/include/metis.h)
    set(METIS_FOUND YES)
    message(STATUS "Metis found ${METIS_DIR}")
    set(METIS_INCLUDES ${METIS_DIR})
    find_path(METIS_INCLUDE_DIR metis.h HINTS "${METIS_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
    list(APPEND METIS_INCLUDES ${METIS_INCLUDE_DIR})
    include_directories(${METIS_INCLUDES})
    find_library(METIS_LIBRARY NAMES metis METIS HINTS "${METIS_DIR}/lib" DOC "Metis library")
    message(STATUS "Metis lib ${METIS_LIBRARY}")
    string(FIND ${METIS_DIR} "metis-5" IS_METIS_5)
    if(IS_METIS_5)
      set(MFEM_USE_METIS_5 ON)
    else()
      set(MFEM_USE_METIS_5 OFF)
    endif()
    get_filename_component(METIS_LIB_DIR ${METIS_LIBRARY} DIRECTORY)
    message(STATUS "METIS_LIB_DIR ${METIS_LIB_DIR}")
    link_directories(${METIS_LIB_DIR})
  else()
    message(FATAL_ERROR "Metis is NOT found")
    set(METIS_FOUND NO)
  endif()
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
#unset(MFEM_BUILD CACHE) # this results in recompiling the whole project every time
#unset(PARALLEL CACHE)

#-------------------------------------------------------------------------------
# Executable
#-------------------------------------------------------------------------------
add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
target_link_libraries(${PROJECT_NAME} ${MFEM_LIBRARY} rt)
if(BUILD_TYPE STREQUAL PDEBUG OR BUILD_TYPE STREQUAL PRELEASE)
  target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES} ${HYPRE_LIBRARY} ${METIS_LIBRARY})
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    target_link_libraries(${PROJECT_NAME} -lmkl -openmp)
  endif()
endif()

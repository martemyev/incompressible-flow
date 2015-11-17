cmake_minimum_required(VERSION 2.8)

project(incompressible-flow CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build: Debug | Release | RelWithDebInfo" FORCE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)

set(MFEM_DIR "/u/artemyev/projects/mfem")
if(BUILD_TYPE MATCHES DEBUG)
  set(MFEM_BUILD "serialdebug")
else()
  set(MFEM_BUILD "serial")
endif()
include_directories("${MFEM_DIR}")
include_directories("${MFEM_DIR}/${MFEM_BUILD}")

link_directories("${MFEM_DIR}/${MFEM_BUILD}")

file(GLOB SRC_LIST "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(GLOB HDR_LIST "${PROJECT_SOURCE_DIR}/src/*.hpp")

configure_file(
  "${PROJECT_SOURCE_DIR}/config.hpp.in"
  "${PROJECT_SOURCE_DIR}/src/config.hpp")

add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
target_link_libraries(${PROJECT_NAME} mfem rt)


cmake_minimum_required(VERSION 3.20)
project(polycube_parameterization VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(PCP_BUILD_VIEWER "Build the optional polyscope viewer" OFF)
option(PCP_BUILD_TESTS "Build tests" ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(CompilerOptions)
include(Dependencies)

add_subdirectory(src/pcp)
add_subdirectory(apps/pcp_cli)
if(PCP_BUILD_VIEWER)
  add_subdirectory(apps/pcp_view)
endif()

if(PCP_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

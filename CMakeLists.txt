
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(ViennaPDE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

add_compile_definitions(CL_TARGET_OPENCL_VERSION=120) # for the openCL version check
include(ViennaCLCommon)
# set(include_dir "${CMAKE_CURRENT_LIST_DIR}/include")
include_directories("viennapde")

# Set release build:
IF(DEFINED CMAKE_BUILD_TYPE)
 SET (CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
ELSE()
 SET (CMAKE_BUILD_TYPE Release)
ENDIF()


add_subdirectory(src)





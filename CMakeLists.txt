
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(ViennaPDE) 

set(CMAKE_CXX_STANDARD 17) # The project has used a lot of novel characteristics of C++17 本库语法使用了很多 C++17 的要素。
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# add_compile_definitions(CL_TARGET_OPENCL_VERSION=120) # It was needed for the openCL version check, but now it seems to be unnecessary.

# LIBRARY 调用库
# ViennaCL library 库, based on the cmake folder /ViennaPDE/cmake  
include(ViennaCLCommon) 
# set(include_dir "${CMAKE_CURRENT_LIST_DIR}/include")
# ViennaPDE library 库, put in the folder /ViennaPDE/viennapde, essentially the core of this project.
include_directories("viennapde")

# DEFAULT RELEASE build mode, 默认为 Release 模式进行编译
IF(DEFINED CMAKE_BUILD_TYPE)
 SET (CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
ELSE()
 SET (CMAKE_BUILD_TYPE Release)
ENDIF()


# GTEST 开启测试
enable_testing()
add_subdirectory(test)





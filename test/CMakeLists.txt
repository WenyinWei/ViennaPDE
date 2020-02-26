## GTest Folder

# find the GTest library, 查找 GTest 库, the user may download the GTest src zip from github and unzip it, e.g., https://github.com/google/googletest/releases/tag/release-1.10.0.
# Following the steps: mkdir build && cd build && cmake .. && sudo make install, after which the gtest will be installed in your linux system.
# For other system that you are concerned, please search relevant online materials to install it. 
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS}) # GTest header file folder, GTest 的头文件

add_subdirectory(viennapde)

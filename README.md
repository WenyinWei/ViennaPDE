## ViennaPDE
ViennaPDE is written by Wenyin during his study of PDE at Tsinghua University. The project has heavy reliance on the ViennaCL library.


# Compile
It is compiled with a NVIDIA CUDA driver and powered by ViennaCL openCL library, otherwise the user may need to configure the cmake command to use other backend (CUDA or openMP).

The command needs to be edited for the above purpose is located at `/home/ununtu/~/ViennaPDE/test/viennapde/core/CMakeLists.txt `
```
set_target_properties(${exec} PROPERTIES COMPILE_FLAGS "-DVIENNACL_WITH_OPENCL") # Guide Viennacl to use opencl mode.
```

ViennaCL is very easy to install, because it is just a hpp-only library and ViennaPDE is the same. To install ViennaCL in your system, as long as you can include its hpp file when needed.

Due to the fact that ViennaPDE is not an official library maintained by University at Vienna but a student in Beijing, the viennacl folder is not contained in viennapde as it does in viennaFEM. The unfamiliar user can consult https://github.com/viennacl/viennacl-dev for install help.

If you want to run the test module, you need GoogleTest. Otherwise comment the 
```
enable_testing()
add_subdirectory(test)
```
lines in CMakeLists.txt in the main directory.

The full compile procedure is 
```
git clone https://github.com/WenyinWei/ViennaPDE.git
mkdir build && cd build
cmake ..
make -j8
make test 
```

# Test
Google test mechanism is adopted in our project and it is recommended that users offer your typical tests.

The install the GTest library, the user may need to download the GTest src zip from github and unzip it, e.g., https://github.com/google/googletest/releases/tag/release-1.10.0.

Following the steps: 

```
mkdir build && cd build
cmake .. 
sudo make install
```  

, after which the gtest would be installed in your linux system.

For other system that you are concerned, please search relevant online materials to install it. 

The following methods have been achieved and you can `make test` to have the result.
- Godunov Method
- LaxWendroff Method 
- WENO + minimod limiter + LaxFreidrichs Method

# Analyze
The ipython code `test/plot_test_data.ipynb` is the python script to analyze the computed result of Burger's equation.

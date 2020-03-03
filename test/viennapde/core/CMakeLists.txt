


foreach(exec mesh scheme convol_mat convol_mesh) # mesh bmesh convol_mat convol_mesh scheme 
  add_executable(${exec} ${exec}.cpp)
  # 设置 ViennaCL library 用的是 OpenCL 后端，你还可以设置它为 CUDA 或者 OpenMP
  set_target_properties(${exec} PROPERTIES COMPILE_FLAGS "-DVIENNACL_WITH_OPENCL") # Guide Viennacl to use opencl mode
  # 链接测试库
  target_link_libraries(${exec} ${OPENCL_LIBRARIES}) # Viennacl Necessary
  target_link_libraries(${exec}
        ${GTEST_BOTH_LIBRARIES} 
        pthread ) # GTest relies on pthread, which you can install by sudo apt-get install glibc-doc && sudo apt-get install manpages-posix manpages-posix-dev
  # 添加到测试
  gtest_discover_tests(${exec})
endforeach()

foreach(scheme Warmup Godunov Godunov_RK3order LaxWendroff WENO_TVDlimiter_LaxFriedrichs)
  string(TOUPPER ${scheme} SCHEME)
  set(DATAPATH_SCHEME1D_${SCHEME} data/scheme/${scheme}/)
  file(MAKE_DIRECTORY ../../../${DATAPATH_SCHEME1D_${SCHEME}})
  target_compile_definitions(scheme PRIVATE
  DATAPATH_SCHEME1D_${SCHEME}="${DATAPATH_SCHEME1D_${SCHEME}}"
  )
endforeach(scheme Warmup Godunov Godunov_RK3order LaxWendroff WENO_TVDlimiter_LaxFriedrichs)


foreach(files pde_test)
  add_executable(${files} ${files}.cpp)
  set_target_properties(${files} PROPERTIES COMPILE_FLAGS "-DVIENNACL_WITH_OPENCL") # Guide Viennacl to use opencl mode
  target_link_libraries(${files} ${OPENCL_LIBRARIES}) # Viennacl Necessary
endforeach()

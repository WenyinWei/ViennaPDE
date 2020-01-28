
foreach(files pde_test)
  add_executable(${files} ${files}.cpp)
  set_target_properties(${files} PROPERTIES COMPILE_FLAGS "-DVIENNACL_WITH_OPENCL")
  target_link_libraries(${files} ${OPENCL_LIBRARIES})
endforeach()

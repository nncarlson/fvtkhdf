if(ENABLE_MPI)
  add_executable(vtkhdf_ug_demo vtkhdf_ug_demo.F90)
  target_link_libraries(vtkhdf_ug_demo PRIVATE fvtkhdf)

  add_executable(vtkhdf_mb_demo vtkhdf_mb_demo.F90)
  target_link_libraries(vtkhdf_mb_demo PRIVATE fvtkhdf)
else()
  add_executable(vtkhdf_ug_demo vtkhdf_ug_demo-serial.F90)
  target_link_libraries(vtkhdf_ug_demo PRIVATE fvtkhdf)

  add_executable(vtkhdf_mb_demo vtkhdf_mb_demo-serial.F90)
  target_link_libraries(vtkhdf_mb_demo PRIVATE fvtkhdf)
endif()

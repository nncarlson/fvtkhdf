if(ENABLE_MPI)
  add_executable(vtkhdf_mb_test vtkhdf_mb_test.F90)
  target_link_libraries(vtkhdf_mb_test PRIVATE fvtkhdf)

  add_test(vtkhdf_mb_file_type ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./vtkhdf_mb_test ${MPIEXEC_POSTFLAGS})
  set_tests_properties(vtkhdf_mb_file_type PROPERTIES PROCESSORS 4)

  add_executable(vtkhdf_ug_test vtkhdf_ug_test.F90)
  target_link_libraries(vtkhdf_ug_test PRIVATE fvtkhdf)

  add_test(vtkhdf_ug_file_type ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./vtkhdf_ug_test ${MPIEXEC_POSTFLAGS})
  set_tests_properties(vtkhdf_ug_file_type PROPERTIES PROCESSORS 4)
else()
  add_executable(vtkhdf_ug_test vtkhdf_ug_test-serial.F90)
  target_link_libraries(vtkhdf_ug_test PRIVATE fvtkhdf)

  add_test(vtkhdf_ug_file_type ./vtkhdf_ug_test)
  set_tests_properties(vtkhdf_ug_file_type PROPERTIES PROCESSORS 1)

  add_executable(vtkhdf_mb_test vtkhdf_mb_test-serial.F90)
  target_link_libraries(vtkhdf_mb_test PRIVATE fvtkhdf)

  add_test(vtkhdf_mb_file_type ./vtkhdf_mb_test)
  set_tests_properties(vtkhdf_mb_file_type PROPERTIES PROCESSORS 1)
endif()

add_executable(vtkhdf_2d_test vtkhdf_2d_test.F90)
target_link_libraries(vtkhdf_2d_test PRIVATE fvtkhdf)
if(ENABLE_MPI)
  add_test(vtkhdf_2d_test ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ./vtkhdf_2d_test ${MPIEXEC_POSTFLAGS})
else()
  add_test(vtkhdf_2d_test ./vtkhdf_2d_test)
endif()
set_tests_properties(vtkhdf_2d_test PROPERTIES PROCESSORS 1)

add_executable(test_leak test_leak.F90)
target_link_libraries(test_leak PRIVATE fvtkhdf)
add_test(NAME test_leak COMMAND test_leak)

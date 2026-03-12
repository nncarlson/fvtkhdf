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

  add_executable(vtkhdf_moving_ug_test vtkhdf_moving_ug_test-serial.F90)
  target_link_libraries(vtkhdf_moving_ug_test PRIVATE fvtkhdf)

  add_test(vtkhdf_moving_ug_file_type ./vtkhdf_moving_ug_test)
  set_tests_properties(vtkhdf_moving_ug_file_type PROPERTIES PROCESSORS 1)

  add_executable(vtkhdf_mb_test vtkhdf_mb_test-serial.F90)
  target_link_libraries(vtkhdf_mb_test PRIVATE fvtkhdf)

  add_test(vtkhdf_mb_file_type ./vtkhdf_mb_test)
  set_tests_properties(vtkhdf_mb_file_type PROPERTIES PROCESSORS 1)
endif()

add_executable(vtkhdf_dynamic_ug_test vtkhdf_dynamic_ug_test.F90)
target_link_libraries(vtkhdf_dynamic_ug_test PRIVATE fvtkhdf)

add_executable(vtkhdf_chunk_test vtkhdf_chunk_test.F90)
target_link_libraries(vtkhdf_chunk_test PRIVATE fvtkhdf)

find_program(H5DUMP_EXECUTABLE NAMES h5dump)

if(ENABLE_MPI)
  add_test(vtkhdf_dynamic_ug_file_type ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./vtkhdf_dynamic_ug_test ${MPIEXEC_POSTFLAGS})
  set_tests_properties(vtkhdf_dynamic_ug_file_type PROPERTIES PROCESSORS 4)

  add_test(vtkhdf_chunk_settings ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./vtkhdf_chunk_test ${MPIEXEC_POSTFLAGS})
  set_tests_properties(vtkhdf_chunk_settings PROPERTIES
    PROCESSORS 4
    FIXTURES_SETUP vtkhdf_chunk_files
  )
else()
  add_test(vtkhdf_dynamic_ug_file_type vtkhdf_dynamic_ug_test)

  add_test(vtkhdf_chunk_settings vtkhdf_chunk_test)
  set_tests_properties(vtkhdf_chunk_settings PROPERTIES FIXTURES_SETUP vtkhdf_chunk_files)
endif()

if(H5DUMP_EXECUTABLE)
  add_test(
    NAME vtkhdf_chunk_default_layout
    COMMAND ${H5DUMP_EXECUTABLE} -pH -d /VTKHDF/Types
            ${CMAKE_CURRENT_BINARY_DIR}/chunk_default.vtkhdf
  )
  set_tests_properties(vtkhdf_chunk_default_layout PROPERTIES
    FIXTURES_REQUIRED vtkhdf_chunk_files
    PASS_REGULAR_EXPRESSION "CHUNKED \\( 131072 \\)"
  )

  add_test(
    NAME vtkhdf_chunk_custom_layout
    COMMAND ${H5DUMP_EXECUTABLE} -pH -d /VTKHDF/Group_0/Types
            ${CMAKE_CURRENT_BINARY_DIR}/chunk_custom.vtkhdf
  )
  set_tests_properties(vtkhdf_chunk_custom_layout PROPERTIES
    FIXTURES_REQUIRED vtkhdf_chunk_files
    PASS_REGULAR_EXPRESSION "CHUNKED \\( 1048576 \\)"
  )
endif()

add_executable(test_leak test_leak.F90)
target_link_libraries(test_leak PRIVATE fvtkhdf)
add_test(NAME test_leak COMMAND test_leak)

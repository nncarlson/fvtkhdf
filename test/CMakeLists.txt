if(ENABLE_MPI)
  add_executable(vtkhdf_mb_test vtkhdf_mb_test.F90)
  target_link_libraries(vtkhdf_mb_test PRIVATE fvtkhdf)

  add_test(vtkhdf_mb_file_type ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
           ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_mb_test ${MPIEXEC_POSTFLAGS})
  set_tests_properties(vtkhdf_mb_file_type PROPERTIES
    PROCESSORS 4
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_executable(vtkhdf_ug_test vtkhdf_ug_test.F90)
  target_link_libraries(vtkhdf_ug_test PRIVATE fvtkhdf)

  add_test(vtkhdf_ug_file_type ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
           ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_ug_test ${MPIEXEC_POSTFLAGS})
  set_tests_properties(vtkhdf_ug_file_type PROPERTIES
    PROCESSORS 4
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
else()
  add_executable(vtkhdf_ug_test vtkhdf_ug_test-serial.F90)
  target_link_libraries(vtkhdf_ug_test PRIVATE fvtkhdf)

  add_test(vtkhdf_ug_file_type ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_ug_test)
  set_tests_properties(vtkhdf_ug_file_type PROPERTIES
    PROCESSORS 1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_executable(vtkhdf_moving_ug_test vtkhdf_moving_ug_test-serial.F90)
  target_link_libraries(vtkhdf_moving_ug_test PRIVATE fvtkhdf)

  add_test(vtkhdf_moving_ug_file_type ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_moving_ug_test)
  set_tests_properties(vtkhdf_moving_ug_file_type PROPERTIES
    PROCESSORS 1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_executable(vtkhdf_mb_test vtkhdf_mb_test-serial.F90)
  target_link_libraries(vtkhdf_mb_test PRIVATE fvtkhdf)

  add_test(vtkhdf_mb_file_type ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_mb_test)
  set_tests_properties(vtkhdf_mb_file_type PROPERTIES
    PROCESSORS 1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

add_executable(vtkhdf_dynamic_ug_test vtkhdf_dynamic_ug_test.F90)
target_link_libraries(vtkhdf_dynamic_ug_test PRIVATE fvtkhdf)

add_executable(vtkhdf_chunk_test vtkhdf_chunk_test.F90)
target_link_libraries(vtkhdf_chunk_test PRIVATE fvtkhdf)

find_program(H5DUMP_EXECUTABLE NAMES h5dump)
set(VTKHDF_READER_SMOKE_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_reader_smoke.py)

execute_process(
  COMMAND ${Python3_EXECUTABLE} ${VTKHDF_READER_SMOKE_SCRIPT} --check-import
  RESULT_VARIABLE VTKHDF_PYTHON_READER_STATUS
  ERROR_VARIABLE VTKHDF_PYTHON_READER_ERROR
  ERROR_STRIP_TRAILING_WHITESPACE
)

if(VTKHDF_PYTHON_READER_STATUS EQUAL 0)
  set(VTKHDF_PYTHON_READER_AVAILABLE TRUE)
  message(STATUS "VTK Python reader available: enabling VTK HDF readback tests")
else()
  set(VTKHDF_PYTHON_READER_AVAILABLE FALSE)
  if(VTKHDF_PYTHON_READER_ERROR)
    message(STATUS "VTK Python reader not available: skipping VTK HDF readback tests (${VTKHDF_PYTHON_READER_ERROR})")
  else()
    message(STATUS "VTK Python reader not available: skipping VTK HDF readback tests")
  endif()
endif()

if(ENABLE_MPI)
  add_test(vtkhdf_dynamic_ug_file_type ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
           ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_dynamic_ug_test ${MPIEXEC_POSTFLAGS})
  set_tests_properties(vtkhdf_dynamic_ug_file_type PROPERTIES
    PROCESSORS 4
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_test(vtkhdf_chunk_settings ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
           ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_chunk_test ${MPIEXEC_POSTFLAGS})
  set_tests_properties(vtkhdf_chunk_settings PROPERTIES
    PROCESSORS 4
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    FIXTURES_SETUP vtkhdf_chunk_files
  )
else()
  add_test(vtkhdf_dynamic_ug_file_type ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_dynamic_ug_test)
  set_tests_properties(vtkhdf_dynamic_ug_file_type PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_test(vtkhdf_chunk_settings ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_chunk_test)
  set_tests_properties(vtkhdf_chunk_settings PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    FIXTURES_SETUP vtkhdf_chunk_files
  )
endif()

if(H5DUMP_EXECUTABLE)
  if(ENABLE_MPI)
    add_test(
      NAME vtkhdf_ug_static_cell_attribute
      COMMAND ${H5DUMP_EXECUTABLE} -H -d /VTKHDF/CellData/static-cell-scalar
              ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
    )
    set_tests_properties(vtkhdf_ug_static_cell_attribute PROPERTIES
      DEPENDS vtkhdf_ug_file_type
      PASS_REGULAR_EXPRESSION "ATTRIBUTE \"Name\""
    )

    add_test(
      NAME vtkhdf_ug_static_point_attribute
      COMMAND ${H5DUMP_EXECUTABLE} -H -d /VTKHDF/PointData/static-point-vector
              ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
    )
    set_tests_properties(vtkhdf_ug_static_point_attribute PROPERTIES
      DEPENDS vtkhdf_ug_file_type
      PASS_REGULAR_EXPRESSION "ATTRIBUTE \"Name\""
    )

    add_test(
      NAME vtkhdf_ug_temporal_cell_attribute
      COMMAND ${H5DUMP_EXECUTABLE} -H -d /VTKHDF/CellData/cell-scalar
              ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
    )
    set_tests_properties(vtkhdf_ug_temporal_cell_attribute PROPERTIES
      DEPENDS vtkhdf_ug_file_type
      PASS_REGULAR_EXPRESSION "ATTRIBUTE \"Name\""
    )

    add_test(
      NAME vtkhdf_ug_temporal_point_attribute
      COMMAND ${H5DUMP_EXECUTABLE} -H -d /VTKHDF/PointData/point-vector
              ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
    )
    set_tests_properties(vtkhdf_ug_temporal_point_attribute PROPERTIES
      DEPENDS vtkhdf_ug_file_type
      PASS_REGULAR_EXPRESSION "ATTRIBUTE \"Name\""
    )

    add_test(
      NAME vtkhdf_mb_static_cell_attribute
      COMMAND ${H5DUMP_EXECUTABLE} -H -d /VTKHDF/Group_0/CellData/static-cell-scalar
              ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
    )
    set_tests_properties(vtkhdf_mb_static_cell_attribute PROPERTIES
      DEPENDS vtkhdf_mb_file_type
      PASS_REGULAR_EXPRESSION "ATTRIBUTE \"Name\""
    )

    add_test(
      NAME vtkhdf_mb_temporal_point_attribute
      COMMAND ${H5DUMP_EXECUTABLE} -H -d /VTKHDF/Group_0/PointData/point-vector
              ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
    )
    set_tests_properties(vtkhdf_mb_temporal_point_attribute PROPERTIES
      DEPENDS vtkhdf_mb_file_type
      PASS_REGULAR_EXPRESSION "ATTRIBUTE \"Name\""
    )
  else()
    add_test(
      NAME vtkhdf_ug_static_cell_attribute
      COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/CellData/static-cell-scalar/Attribute
              ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
    )
    set_tests_properties(vtkhdf_ug_static_cell_attribute PROPERTIES
      DEPENDS vtkhdf_ug_file_type
      PASS_REGULAR_EXPRESSION "\"Scalars\""
    )

    add_test(
      NAME vtkhdf_ug_static_point_attribute
      COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/PointData/static-point-vector/Attribute
              ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
    )
    set_tests_properties(vtkhdf_ug_static_point_attribute PROPERTIES
      DEPENDS vtkhdf_ug_file_type
      PASS_REGULAR_EXPRESSION "\"Vectors\""
    )

    add_test(
      NAME vtkhdf_ug_temporal_cell_attribute
      COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/CellData/cell-radius/Attribute
              ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
    )
    set_tests_properties(vtkhdf_ug_temporal_cell_attribute PROPERTIES
      DEPENDS vtkhdf_ug_file_type
      PASS_REGULAR_EXPRESSION "\"Scalars\""
    )

    add_test(
      NAME vtkhdf_ug_temporal_point_attribute
      COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/PointData/point-velocity/Attribute
              ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
    )
    set_tests_properties(vtkhdf_ug_temporal_point_attribute PROPERTIES
      DEPENDS vtkhdf_ug_file_type
      PASS_REGULAR_EXPRESSION "\"Vectors\""
    )

    add_test(
      NAME vtkhdf_mb_static_cell_attribute
      COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/Group_0/CellData/static-cell-scalar/Attribute
              ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
    )
    set_tests_properties(vtkhdf_mb_static_cell_attribute PROPERTIES
      DEPENDS vtkhdf_mb_file_type
      PASS_REGULAR_EXPRESSION "\"Scalars\""
    )

    add_test(
      NAME vtkhdf_mb_temporal_point_attribute
      COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/Group_1/PointData/point-velocity/Attribute
              ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
    )
    set_tests_properties(vtkhdf_mb_temporal_point_attribute PROPERTIES
      DEPENDS vtkhdf_mb_file_type
      PASS_REGULAR_EXPRESSION "\"Vectors\""
    )
  endif()

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

add_test(
  NAME vtkhdf_ug_reader_smoke
  COMMAND ${Python3_EXECUTABLE} ${VTKHDF_READER_SMOKE_SCRIPT}
          ug ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
)
set_tests_properties(vtkhdf_ug_reader_smoke PROPERTIES
  DEPENDS vtkhdf_ug_file_type
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  SKIP_RETURN_CODE 77
)

add_test(
  NAME vtkhdf_mb_reader_smoke
  COMMAND ${Python3_EXECUTABLE} ${VTKHDF_READER_SMOKE_SCRIPT}
          mb ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
)
set_tests_properties(vtkhdf_mb_reader_smoke PROPERTIES
  DEPENDS vtkhdf_mb_file_type
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  SKIP_RETURN_CODE 77
)

add_executable(test_leak test_leak.F90)
target_link_libraries(test_leak PRIVATE fvtkhdf)
add_test(NAME test_leak COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_leak)
set_tests_properties(test_leak PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

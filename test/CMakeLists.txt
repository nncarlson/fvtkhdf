add_library(vtkhdf_test_support STATIC
  support/vtkhdf_test_env.F90
  support/vtkhdf_test_data.F90
)
target_link_libraries(vtkhdf_test_support PUBLIC fvtkhdf)

function(vtkhdf_add_test_target target source)
  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE vtkhdf_test_support)
endfunction()

function(vtkhdf_add_runtime_test test_name target)
  if(ENABLE_MPI)
    add_test(${test_name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
             ${CMAKE_CURRENT_BINARY_DIR}/${target} ${MPIEXEC_POSTFLAGS})
    set_tests_properties(${test_name} PROPERTIES
      PROCESSORS 4
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  else()
    add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${target})
    set_tests_properties(${test_name} PROPERTIES
      PROCESSORS 1
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  endif()
endfunction()

vtkhdf_add_test_target(vtkhdf_mb_test vtkhdf_mb_test.F90)
vtkhdf_add_test_target(vtkhdf_mb_basic_blocks_test vtkhdf_mb_basic_blocks_test.F90)
vtkhdf_add_test_target(vtkhdf_ug_test vtkhdf_ug_test.F90)
vtkhdf_add_test_target(vtkhdf_ug_static_test vtkhdf_ug_static_test.F90)
vtkhdf_add_test_target(vtkhdf_ug_fixed_file_test vtkhdf_ug_fixed_test.F90)
vtkhdf_add_test_target(vtkhdf_ug_fail_late_register_test vtkhdf_ug_fail_late_register_test.F90)
vtkhdf_add_test_target(vtkhdf_ug_fail_wrong_cell_count_test vtkhdf_ug_fail_wrong_cell_count_test.F90)
vtkhdf_add_test_target(vtkhdf_ug_fail_data_before_mesh_test vtkhdf_ug_fail_data_before_mesh_test.F90)
vtkhdf_add_test_target(vtkhdf_mb_fail_wrong_block_handle_test vtkhdf_mb_fail_wrong_block_handle_test.F90)
vtkhdf_add_test_target(vtkhdf_mb_fail_late_block_test vtkhdf_mb_fail_late_block_test.F90)
vtkhdf_add_test_target(vtkhdf_moving_ug_test vtkhdf_moving_ug_test.F90)
vtkhdf_add_test_target(vtkhdf_dynamic_ug_test vtkhdf_dynamic_ug_test.F90)
vtkhdf_add_test_target(vtkhdf_chunk_test vtkhdf_chunk_test.F90)

vtkhdf_add_runtime_test(vtkhdf_mb_fixture_file vtkhdf_mb_test)
vtkhdf_add_runtime_test(vtkhdf_mb_basic_blocks_file_type vtkhdf_mb_basic_blocks_test)
vtkhdf_add_runtime_test(vtkhdf_ug_fixture_file vtkhdf_ug_test)
vtkhdf_add_runtime_test(vtkhdf_ug_static_file_type vtkhdf_ug_static_test)
vtkhdf_add_runtime_test(vtkhdf_ug_fixed_file_type vtkhdf_ug_fixed_file_test)
vtkhdf_add_runtime_test(vtkhdf_ug_fail_late_register vtkhdf_ug_fail_late_register_test)
vtkhdf_add_runtime_test(vtkhdf_ug_fail_wrong_cell_count vtkhdf_ug_fail_wrong_cell_count_test)
vtkhdf_add_runtime_test(vtkhdf_ug_fail_data_before_mesh vtkhdf_ug_fail_data_before_mesh_test)
vtkhdf_add_runtime_test(vtkhdf_mb_fail_wrong_block_handle vtkhdf_mb_fail_wrong_block_handle_test)
vtkhdf_add_runtime_test(vtkhdf_mb_fail_late_block vtkhdf_mb_fail_late_block_test)
vtkhdf_add_runtime_test(vtkhdf_moving_ug_file_type vtkhdf_moving_ug_test)
vtkhdf_add_runtime_test(vtkhdf_dynamic_ug_file_type vtkhdf_dynamic_ug_test)
vtkhdf_add_runtime_test(vtkhdf_chunk_settings vtkhdf_chunk_test)

set_tests_properties(
  vtkhdf_ug_fail_late_register
  vtkhdf_ug_fail_wrong_cell_count
  vtkhdf_ug_fail_data_before_mesh
  vtkhdf_mb_fail_wrong_block_handle
  vtkhdf_mb_fail_late_block
  PROPERTIES
  WILL_FAIL TRUE
)

set_tests_properties(vtkhdf_chunk_settings PROPERTIES
  FIXTURES_SETUP vtkhdf_chunk_files
)

if(ENABLE_MPI)
  vtkhdf_add_test_target(vtkhdf_mb_zero_part_test vtkhdf_mb_zero_part_test.F90)
  vtkhdf_add_runtime_test(vtkhdf_mb_zero_part_file_type vtkhdf_mb_zero_part_test)
endif()

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
  set(VTKHDF_DYNAMIC_UG_NUMBER_OF_POINTS_REGEX "\\(0\\): 4, 4, 4, 4, 5, 5, 5, 5")
  set(VTKHDF_DYNAMIC_UG_NUMBER_OF_CELLS_REGEX "\\(0\\): 1, 1, 1, 1, 4, 4, 4, 4")
  set(VTKHDF_DYNAMIC_UG_POINT_OFFSETS_REGEX "\\(0\\): 0, 0, 16, 36, 56")
  set(VTKHDF_DYNAMIC_UG_CELL_OFFSETS_REGEX "\\(0\\): 0, 0, 4, 4, 4")
  set(VTKHDF_DYNAMIC_UG_CONNECTIVITY_ID_OFFSETS_REGEX "\\(0\\): 0, 0, 16, 16, 16")
else()
  set(VTKHDF_DYNAMIC_UG_NUMBER_OF_POINTS_REGEX "\\(0\\): 4, 5")
  set(VTKHDF_DYNAMIC_UG_NUMBER_OF_CELLS_REGEX "\\(0\\): 1, 4")
  set(VTKHDF_DYNAMIC_UG_POINT_OFFSETS_REGEX "\\(0\\): 0, 0, 4, 9, 14")
  set(VTKHDF_DYNAMIC_UG_CELL_OFFSETS_REGEX "\\(0\\): 0, 0, 1, 1, 1")
  set(VTKHDF_DYNAMIC_UG_CONNECTIVITY_ID_OFFSETS_REGEX "\\(0\\): 0, 0, 4, 4, 4")
endif()

if(H5DUMP_EXECUTABLE)

  add_test(
    NAME vtkhdf_ug_static_cell_attribute
    COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/CellData/static-cell-scalar/Attribute
            ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
  )
  set_tests_properties(vtkhdf_ug_static_cell_attribute PROPERTIES
    DEPENDS vtkhdf_ug_fixture_file
    PASS_REGULAR_EXPRESSION "\"Scalars\""
  )

  add_test(
    NAME vtkhdf_ug_static_point_attribute
    COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/PointData/static-point-vector/Attribute
            ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
  )
  set_tests_properties(vtkhdf_ug_static_point_attribute PROPERTIES
    DEPENDS vtkhdf_ug_fixture_file
    PASS_REGULAR_EXPRESSION "\"Vectors\""
  )

  add_test(
    NAME vtkhdf_ug_temporal_cell_attribute
    COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/CellData/cell-radius/Attribute
            ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
  )
  set_tests_properties(vtkhdf_ug_temporal_cell_attribute PROPERTIES
    DEPENDS vtkhdf_ug_fixture_file
    PASS_REGULAR_EXPRESSION "\"Scalars\""
  )

  add_test(
    NAME vtkhdf_ug_temporal_point_attribute
    COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/PointData/point-velocity/Attribute
            ${CMAKE_CURRENT_BINARY_DIR}/ug_test.vtkhdf
  )
  set_tests_properties(vtkhdf_ug_temporal_point_attribute PROPERTIES
    DEPENDS vtkhdf_ug_fixture_file
    PASS_REGULAR_EXPRESSION "\"Vectors\""
  )

  add_test(
    NAME vtkhdf_mb_static_cell_attribute
    COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/Group_0/CellData/static-cell-scalar/Attribute
            ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
  )
  set_tests_properties(vtkhdf_mb_static_cell_attribute PROPERTIES
    DEPENDS vtkhdf_mb_fixture_file
    PASS_REGULAR_EXPRESSION "\"Scalars\""
  )

  add_test(
    NAME vtkhdf_mb_temporal_point_attribute
    COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/Group_1/PointData/point-velocity/Attribute
            ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
  )
  set_tests_properties(vtkhdf_mb_temporal_point_attribute PROPERTIES
    DEPENDS vtkhdf_mb_fixture_file
    PASS_REGULAR_EXPRESSION "\"Vectors\""
  )

  add_test(
    NAME vtkhdf_mb_assembly_block_a_link
    COMMAND ${H5DUMP_EXECUTABLE} -n ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
  )
  set_tests_properties(vtkhdf_mb_assembly_block_a_link PROPERTIES
    DEPENDS vtkhdf_mb_fixture_file
    PASS_REGULAR_EXPRESSION "link +/VTKHDF/Assembly/Block-A -> /VTKHDF/Group_0"
  )

  add_test(
    NAME vtkhdf_mb_assembly_block_b_link
    COMMAND ${H5DUMP_EXECUTABLE} -n ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
  )
  set_tests_properties(vtkhdf_mb_assembly_block_b_link PROPERTIES
    DEPENDS vtkhdf_mb_fixture_file
    PASS_REGULAR_EXPRESSION "link +/VTKHDF/Assembly/Block-B -> /VTKHDF/Group_1"
  )

  add_test(
    NAME vtkhdf_mb_group_0_name
    COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/Group_0/Name
            ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
  )
  set_tests_properties(vtkhdf_mb_group_0_name PROPERTIES
    DEPENDS vtkhdf_mb_fixture_file
    PASS_REGULAR_EXPRESSION "\"Block-A\""
  )

  add_test(
    NAME vtkhdf_mb_group_1_name
    COMMAND ${H5DUMP_EXECUTABLE} -a /VTKHDF/Group_1/Name
            ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
  )
  set_tests_properties(vtkhdf_mb_group_1_name PROPERTIES
    DEPENDS vtkhdf_mb_fixture_file
    PASS_REGULAR_EXPRESSION "\"Block-B\""
  )

  if(ENABLE_MPI)
    add_test(
      NAME vtkhdf_mb_zero_part_sparse_number_of_parts
      COMMAND ${H5DUMP_EXECUTABLE} -d /VTKHDF/Group_1/Steps/NumberOfParts
              ${CMAKE_CURRENT_BINARY_DIR}/mb_zero_part_test.vtkhdf
    )
    set_tests_properties(vtkhdf_mb_zero_part_sparse_number_of_parts PROPERTIES
      DEPENDS vtkhdf_mb_zero_part_file_type
      PASS_REGULAR_EXPRESSION "\\(0\\): 3, 3"
    )

    add_test(
      NAME vtkhdf_mb_zero_part_sparse_number_of_points
      COMMAND ${H5DUMP_EXECUTABLE} -d /VTKHDF/Group_1/NumberOfPoints
              ${CMAKE_CURRENT_BINARY_DIR}/mb_zero_part_test.vtkhdf
    )
    set_tests_properties(vtkhdf_mb_zero_part_sparse_number_of_points PROPERTIES
      DEPENDS vtkhdf_mb_zero_part_file_type
      PASS_REGULAR_EXPRESSION "\\(0\\): 8, 8, 8"
    )

    add_test(
      NAME vtkhdf_mb_zero_part_sparse_point_data_offsets
      COMMAND ${H5DUMP_EXECUTABLE} -d /VTKHDF/Group_1/Steps/PointDataOffsets/point-radius
              ${CMAKE_CURRENT_BINARY_DIR}/mb_zero_part_test.vtkhdf
    )
    set_tests_properties(vtkhdf_mb_zero_part_sparse_point_data_offsets PROPERTIES
      DEPENDS vtkhdf_mb_zero_part_file_type
      PASS_REGULAR_EXPRESSION "\\(0\\): 0, 24"
    )
  endif()

  add_test(
    NAME vtkhdf_dynamic_ug_number_of_points
    COMMAND ${H5DUMP_EXECUTABLE} -d /VTKHDF/NumberOfPoints
            ${CMAKE_CURRENT_BINARY_DIR}/dyn_ug_test.vtkhdf
  )
  set_tests_properties(vtkhdf_dynamic_ug_number_of_points PROPERTIES
    DEPENDS vtkhdf_dynamic_ug_file_type
    PASS_REGULAR_EXPRESSION "${VTKHDF_DYNAMIC_UG_NUMBER_OF_POINTS_REGEX}"
  )

  add_test(
    NAME vtkhdf_dynamic_ug_number_of_cells
    COMMAND ${H5DUMP_EXECUTABLE} -d /VTKHDF/NumberOfCells
            ${CMAKE_CURRENT_BINARY_DIR}/dyn_ug_test.vtkhdf
  )
  set_tests_properties(vtkhdf_dynamic_ug_number_of_cells PROPERTIES
    DEPENDS vtkhdf_dynamic_ug_file_type
    PASS_REGULAR_EXPRESSION "${VTKHDF_DYNAMIC_UG_NUMBER_OF_CELLS_REGEX}"
  )

  add_test(
    NAME vtkhdf_dynamic_ug_point_offsets
    COMMAND ${H5DUMP_EXECUTABLE} -d /VTKHDF/Steps/PointOffsets
            ${CMAKE_CURRENT_BINARY_DIR}/dyn_ug_test.vtkhdf
  )
  set_tests_properties(vtkhdf_dynamic_ug_point_offsets PROPERTIES
    DEPENDS vtkhdf_dynamic_ug_file_type
    PASS_REGULAR_EXPRESSION "${VTKHDF_DYNAMIC_UG_POINT_OFFSETS_REGEX}"
  )

  add_test(
    NAME vtkhdf_dynamic_ug_cell_offsets
    COMMAND ${H5DUMP_EXECUTABLE} -d /VTKHDF/Steps/CellOffsets
            ${CMAKE_CURRENT_BINARY_DIR}/dyn_ug_test.vtkhdf
  )
  set_tests_properties(vtkhdf_dynamic_ug_cell_offsets PROPERTIES
    DEPENDS vtkhdf_dynamic_ug_file_type
    PASS_REGULAR_EXPRESSION "${VTKHDF_DYNAMIC_UG_CELL_OFFSETS_REGEX}"
  )

  add_test(
    NAME vtkhdf_dynamic_ug_connectivity_id_offsets
    COMMAND ${H5DUMP_EXECUTABLE} -d /VTKHDF/Steps/ConnectivityIdOffsets
            ${CMAKE_CURRENT_BINARY_DIR}/dyn_ug_test.vtkhdf
  )
  set_tests_properties(vtkhdf_dynamic_ug_connectivity_id_offsets PROPERTIES
    DEPENDS vtkhdf_dynamic_ug_file_type
    PASS_REGULAR_EXPRESSION "${VTKHDF_DYNAMIC_UG_CONNECTIVITY_ID_OFFSETS_REGEX}"
  )

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
  DEPENDS vtkhdf_ug_fixture_file
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  SKIP_RETURN_CODE 77
)

add_test(
  NAME vtkhdf_mb_reader_smoke
  COMMAND ${Python3_EXECUTABLE} ${VTKHDF_READER_SMOKE_SCRIPT}
          mb ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
)
set_tests_properties(vtkhdf_mb_reader_smoke PROPERTIES
  DEPENDS vtkhdf_mb_fixture_file
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  SKIP_RETURN_CODE 77
)

if(ENABLE_MPI)
  add_test(
    NAME vtkhdf_mb_zero_part_reader_smoke
    COMMAND ${Python3_EXECUTABLE} ${VTKHDF_READER_SMOKE_SCRIPT}
            mb ${CMAKE_CURRENT_BINARY_DIR}/mb_zero_part_test.vtkhdf
  )
  set_tests_properties(vtkhdf_mb_zero_part_reader_smoke PROPERTIES
    DEPENDS vtkhdf_mb_zero_part_file_type
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    SKIP_RETURN_CODE 77
  )
endif()

add_executable(test_leak test_leak.F90)
target_link_libraries(test_leak PRIVATE vtkhdf_test_support)
add_test(NAME test_leak COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_leak)
set_tests_properties(test_leak PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

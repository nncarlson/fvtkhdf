add_library(vtkhdf_test_support STATIC
  support/vtkhdf_test_env.F90
  support/vtkhdf_test_data.F90
)
target_link_libraries(vtkhdf_test_support PUBLIC fvtkhdf)

function(vtkhdf_add_test_target target source)
  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE vtkhdf_test_support)
endfunction()

# Exercise the real fatal handler without linking the application's MPI/HDF5 stack.
add_executable(vtkhdf_assert_returning_mpi_test
  vtkhdf_assert_test.F90
  support/vtkhdf_assert_stubs.F90
  ${PROJECT_SOURCE_DIR}/src/vtkhdf_assert.F90
)
set(ASSERT_TEST_MODULE_DIR ${CMAKE_CURRENT_BINARY_DIR}/assert-returning-mpi-mod)
file(MAKE_DIRECTORY ${ASSERT_TEST_MODULE_DIR})
set_target_properties(vtkhdf_assert_returning_mpi_test PROPERTIES
  Fortran_MODULE_DIRECTORY ${ASSERT_TEST_MODULE_DIR}
)
target_include_directories(vtkhdf_assert_returning_mpi_test PRIVATE
  ${ASSERT_TEST_MODULE_DIR} ${PROJECT_SOURCE_DIR}/src
)
target_compile_definitions(vtkhdf_assert_returning_mpi_test PRIVATE USE_MPI)
add_test(NAME vtkhdf_assert_returning_mpi
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_assert_check.py
          $<TARGET_FILE:vtkhdf_assert_returning_mpi_test>
)
set_tests_properties(vtkhdf_assert_returning_mpi PROPERTIES
  PROCESSORS 1 TIMEOUT 30
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

function(vtkhdf_add_runtime_test test_name target)
  if(ENABLE_MPI)
    add_test(${test_name} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
             ${CMAKE_CURRENT_BINARY_DIR}/${target} ${MPIEXEC_POSTFLAGS} ${ARGN})
    set_tests_properties(${test_name} PROPERTIES
      PROCESSORS 4
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  else()
    add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${target} ${ARGN})
    set_tests_properties(${test_name} PROPERTIES
      PROCESSORS 1
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  endif()
endfunction()

if(ENABLE_MPI)
  set(ASSEMBLY_READER_RANKS 4)
else()
  set(ASSEMBLY_READER_RANKS 1)
endif()

# Both concrete file types exercise the same Assembly implementation.
foreach(format mb pdc)
  vtkhdf_add_test_target(vtkhdf_${format}_assembly_test vtkhdf_assembly_test.F90)
  vtkhdf_add_test_target(vtkhdf_${format}_assembly_fail_test vtkhdf_assembly_fail_test.F90)
  vtkhdf_add_test_target(vtkhdf_${format}_assembly_types_test vtkhdf_mb_collections_test.F90)
  target_compile_definitions(vtkhdf_${format}_assembly_types_test PRIVATE ASSEMBLY_COLLECTION)
  vtkhdf_add_runtime_test(vtkhdf_${format}_assembly_types vtkhdf_${format}_assembly_types_test)
  set_tests_properties(vtkhdf_${format}_assembly_types PROPERTIES FIXTURES_SETUP vtkhdf_${format}_assembly_type_files)
  add_test(NAME vtkhdf_${format}_assembly_types_layout
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_mb_collections_check.py
            layout ${CMAKE_CURRENT_BINARY_DIR} ${format}_assembly_types)
  set_tests_properties(vtkhdf_${format}_assembly_types_layout PROPERTIES
    FIXTURES_REQUIRED vtkhdf_${format}_assembly_type_files SKIP_RETURN_CODE 77)

  if(format STREQUAL "pdc")
    target_compile_definitions(vtkhdf_${format}_assembly_types_test PRIVATE PDC_COLLECTION)
    target_compile_definitions(vtkhdf_${format}_assembly_test PRIVATE PDC_COLLECTION)
    target_compile_definitions(vtkhdf_${format}_assembly_fail_test PRIVATE PDC_COLLECTION)
  endif()
  vtkhdf_add_runtime_test(vtkhdf_${format}_assembly vtkhdf_${format}_assembly_test)
  set_tests_properties(vtkhdf_${format}_assembly PROPERTIES FIXTURES_SETUP vtkhdf_${format}_assembly_files)
  add_test(NAME vtkhdf_${format}_assembly_layout
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_assembly_check.py
            ${CMAKE_CURRENT_BINARY_DIR} ${format})
  set_tests_properties(vtkhdf_${format}_assembly_layout PROPERTIES
    FIXTURES_REQUIRED vtkhdf_${format}_assembly_files SKIP_RETURN_CODE 77)
  add_test(NAME vtkhdf_${format}_assembly_reader_compatibility
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_assembly_reader_probe.py
            ${CMAKE_CURRENT_BINARY_DIR} ${format} ${ASSEMBLY_READER_RANKS})
  set_tests_properties(vtkhdf_${format}_assembly_reader_compatibility PROPERTIES
    FIXTURES_REQUIRED vtkhdf_${format}_assembly_files SKIP_RETURN_CODE 77 LABELS reader-compatibility)

  set(assembly_errors negative_parent future_parent foreign standalone closed_node closed_ref closed_commit
      repeat_commit block_after node_after ref_after duplicate_after)
  if(ENABLE_MPI AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    list(APPEND assembly_errors mismatch_name mismatch_parent mismatch_reference mismatch_custom mismatch_count)
  endif()
  foreach(error IN LISTS assembly_errors)
    set(name vtkhdf_${format}_assembly_fail_${error})
    vtkhdf_add_runtime_test(${name} vtkhdf_${format}_assembly_fail_test ${error})
    set(dir ${CMAKE_CURRENT_BINARY_DIR}/assembly-fail-${format}-${error})
    file(MAKE_DIRECTORY ${dir})
    set_tests_properties(${name} PROPERTIES WILL_FAIL TRUE WORKING_DIRECTORY ${dir} TIMEOUT 30)
  endforeach()
endforeach()

vtkhdf_add_test_target(vtkhdf_mb_test vtkhdf_mb_test.F90)
vtkhdf_add_test_target(vtkhdf_mb_basic_blocks_test vtkhdf_mb_basic_blocks_test.F90)
vtkhdf_add_test_target(vtkhdf_mb_typed_blocks_test vtkhdf_mb_typed_blocks_test.F90)
vtkhdf_add_test_target(vtkhdf_mb_collections_test vtkhdf_mb_collections_test.F90)
vtkhdf_add_test_target(vtkhdf_pdc_collections_test vtkhdf_mb_collections_test.F90)
target_compile_definitions(vtkhdf_pdc_collections_test PRIVATE PDC_COLLECTION)
vtkhdf_add_runtime_test(vtkhdf_pdc_collections vtkhdf_pdc_collections_test)
set_tests_properties(vtkhdf_pdc_collections PROPERTIES FIXTURES_SETUP vtkhdf_pdc_collection_files)
add_test(NAME vtkhdf_pdc_collections_layout
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_mb_collections_check.py
          layout ${CMAKE_CURRENT_BINARY_DIR} pdc)
set_tests_properties(vtkhdf_pdc_collections_layout PROPERTIES
  FIXTURES_REQUIRED vtkhdf_pdc_collection_files SKIP_RETURN_CODE 77)
add_test(NAME vtkhdf_pdc_reader_compatibility
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_pdc_reader_probe.py
          ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(vtkhdf_pdc_reader_compatibility PROPERTIES
  FIXTURES_REQUIRED vtkhdf_pdc_collection_files SKIP_RETURN_CODE 77 LABELS reader-compatibility)


vtkhdf_add_test_target(vtkhdf_mb_duplicate_block_name_test vtkhdf_mb_duplicate_block_name_test.F90)
vtkhdf_add_test_target(vtkhdf_unstructured_test vtkhdf_unstructured_test.F90)
vtkhdf_add_test_target(vtkhdf_unstructured_static_test vtkhdf_unstructured_static_test.F90)
vtkhdf_add_test_target(vtkhdf_unstructured_fixed_file_test vtkhdf_unstructured_fixed_test.F90)
vtkhdf_add_test_target(vtkhdf_unstructured_fail_late_register_test vtkhdf_unstructured_fail_late_register_test.F90)
vtkhdf_add_test_target(vtkhdf_unstructured_fail_wrong_cell_count_test vtkhdf_unstructured_fail_wrong_cell_count_test.F90)
vtkhdf_add_test_target(vtkhdf_unstructured_fail_data_before_mesh_test vtkhdf_unstructured_fail_data_before_mesh_test.F90)
vtkhdf_add_test_target(vtkhdf_mb_fail_wrong_block_handle_test vtkhdf_mb_fail_wrong_block_handle_test.F90)
vtkhdf_add_test_target(vtkhdf_mb_fail_late_block_test vtkhdf_mb_fail_late_block_test.F90)
vtkhdf_add_test_target(vtkhdf_moving_unstructured_test vtkhdf_moving_unstructured_test.F90)
vtkhdf_add_test_target(vtkhdf_dynamic_unstructured_test vtkhdf_dynamic_unstructured_test.F90)
vtkhdf_add_test_target(vtkhdf_chunk_test vtkhdf_chunk_test.F90)
vtkhdf_add_test_target(vtkhdf_fail_foreign_data_handle_test vtkhdf_fail_foreign_data_handle_test.F90)

vtkhdf_add_runtime_test(vtkhdf_mb_fixture_file vtkhdf_mb_test)
vtkhdf_add_runtime_test(vtkhdf_mb_basic_blocks_file_type vtkhdf_mb_basic_blocks_test)
vtkhdf_add_runtime_test(vtkhdf_mb_typed_blocks vtkhdf_mb_typed_blocks_test)
set_tests_properties(vtkhdf_mb_typed_blocks PROPERTIES FIXTURES_SETUP vtkhdf_mb_typed_files)
foreach(check layout reader)
  add_test(NAME vtkhdf_mb_typed_${check}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_mb_typed_blocks_check.py
            ${check} ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(vtkhdf_mb_typed_${check} PROPERTIES
    FIXTURES_REQUIRED vtkhdf_mb_typed_files SKIP_RETURN_CODE 77)
endforeach()
vtkhdf_add_runtime_test(vtkhdf_mb_collections vtkhdf_mb_collections_test)
set_tests_properties(vtkhdf_mb_collections PROPERTIES FIXTURES_SETUP vtkhdf_mb_collection_files)
foreach(check layout reader)
  add_test(NAME vtkhdf_mb_collections_${check}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_mb_collections_check.py
            ${check} ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(vtkhdf_mb_collections_${check} PROPERTIES
    FIXTURES_REQUIRED vtkhdf_mb_collection_files SKIP_RETURN_CODE 77)
endforeach()
foreach(format mb pdc)
  vtkhdf_add_test_target(vtkhdf_${format}_hypertree_test vtkhdf_mb_hypertree_test.F90)
  if(format STREQUAL "pdc")
    target_compile_definitions(vtkhdf_${format}_hypertree_test PRIVATE PDC_COLLECTION)
  endif()
  vtkhdf_add_runtime_test(vtkhdf_${format}_hypertree vtkhdf_${format}_hypertree_test)
  set_tests_properties(vtkhdf_${format}_hypertree PROPERTIES FIXTURES_SETUP vtkhdf_${format}_hypertree_files)
  foreach(check layout reader)
    add_test(NAME vtkhdf_${format}_hypertree_${check}
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_mb_hypertree_check.py
              ${check} ${CMAKE_CURRENT_BINARY_DIR} ${format})
    set_tests_properties(vtkhdf_${format}_hypertree_${check} PROPERTIES
      FIXTURES_REQUIRED vtkhdf_${format}_hypertree_files SKIP_RETURN_CODE 77)
  endforeach()
  if(ENABLE_MPI)
    set(HYPERTREE_TWO_RANK_DIR ${CMAKE_CURRENT_BINARY_DIR}/${format}-hypertree-two-rank)
    file(MAKE_DIRECTORY ${HYPERTREE_TWO_RANK_DIR})
    add_test(NAME vtkhdf_${format}_hypertree_two_rank
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2
              $<TARGET_FILE:vtkhdf_${format}_hypertree_test> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(vtkhdf_${format}_hypertree_two_rank PROPERTIES PROCESSORS 2
      WORKING_DIRECTORY ${HYPERTREE_TWO_RANK_DIR} FIXTURES_SETUP vtkhdf_${format}_hypertree_two_rank_files)
    foreach(check layout reader)
      add_test(NAME vtkhdf_${format}_hypertree_two_rank_${check}
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_mb_hypertree_check.py
                ${check} ${HYPERTREE_TWO_RANK_DIR} ${format})
      set_tests_properties(vtkhdf_${format}_hypertree_two_rank_${check} PROPERTIES
        FIXTURES_REQUIRED vtkhdf_${format}_hypertree_two_rank_files SKIP_RETURN_CODE 77)
    endforeach()
  endif()
endforeach()
vtkhdf_add_runtime_test(vtkhdf_mb_duplicate_block_name_file_type vtkhdf_mb_duplicate_block_name_test)
vtkhdf_add_runtime_test(vtkhdf_unstructured_fixture_file vtkhdf_unstructured_test)
vtkhdf_add_runtime_test(vtkhdf_unstructured_static_file_type vtkhdf_unstructured_static_test)
vtkhdf_add_runtime_test(vtkhdf_unstructured_fixed_file_type vtkhdf_unstructured_fixed_file_test)
vtkhdf_add_runtime_test(vtkhdf_unstructured_fail_late_register vtkhdf_unstructured_fail_late_register_test)
vtkhdf_add_runtime_test(vtkhdf_unstructured_fail_wrong_cell_count vtkhdf_unstructured_fail_wrong_cell_count_test)
vtkhdf_add_runtime_test(vtkhdf_unstructured_fail_data_before_mesh vtkhdf_unstructured_fail_data_before_mesh_test)
vtkhdf_add_runtime_test(vtkhdf_mb_fail_wrong_block_handle vtkhdf_mb_fail_wrong_block_handle_test)
vtkhdf_add_runtime_test(vtkhdf_mb_fail_late_block vtkhdf_mb_fail_late_block_test)
vtkhdf_add_runtime_test(vtkhdf_moving_unstructured_file_type vtkhdf_moving_unstructured_test)
vtkhdf_add_runtime_test(vtkhdf_dynamic_unstructured_file_type vtkhdf_dynamic_unstructured_test)
vtkhdf_add_runtime_test(vtkhdf_chunk_settings vtkhdf_chunk_test)

set_tests_properties(
  vtkhdf_unstructured_fail_late_register
  vtkhdf_unstructured_fail_wrong_cell_count
  vtkhdf_unstructured_fail_data_before_mesh
  vtkhdf_mb_fail_wrong_block_handle
  vtkhdf_mb_fail_late_block
  PROPERTIES
  WILL_FAIL TRUE
)

foreach(format unstructured mb)
  foreach(association cell point field)
    set(test_name vtkhdf_${format}_fail_foreign_${association}_handle)
    vtkhdf_add_runtime_test(${test_name} vtkhdf_fail_foreign_data_handle_test ${format} ${association})
    set_tests_properties(${test_name} PROPERTIES WILL_FAIL TRUE)
  endforeach()
endforeach()

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
  set(VTKHDF_DYNAMIC_UNSTRUCTURED_NUMBER_OF_POINTS_REGEX "\\(0\\): 4, 4, 4, 4, 5, 5, 5, 5")
  set(VTKHDF_DYNAMIC_UNSTRUCTURED_NUMBER_OF_CELLS_REGEX "\\(0\\): 1, 1, 1, 1, 4, 4, 4, 4")
  set(VTKHDF_DYNAMIC_UNSTRUCTURED_POINT_OFFSETS_REGEX "\\(0\\): 0, 0, 16, 36, 56")
  set(VTKHDF_DYNAMIC_UNSTRUCTURED_CELL_OFFSETS_REGEX "\\(0\\): 0, 0, 4, 4, 4")
  set(VTKHDF_DYNAMIC_UNSTRUCTURED_CONNECTIVITY_ID_OFFSETS_REGEX "\\(0\\): 0, 0, 16, 16, 16")
else()
  set(VTKHDF_DYNAMIC_UNSTRUCTURED_NUMBER_OF_POINTS_REGEX "\\(0\\): 4, 5")
  set(VTKHDF_DYNAMIC_UNSTRUCTURED_NUMBER_OF_CELLS_REGEX "\\(0\\): 1, 4")
  set(VTKHDF_DYNAMIC_UNSTRUCTURED_POINT_OFFSETS_REGEX "\\(0\\): 0, 0, 4, 9, 14")
  set(VTKHDF_DYNAMIC_UNSTRUCTURED_CELL_OFFSETS_REGEX "\\(0\\): 0, 0, 1, 1, 1")
  set(VTKHDF_DYNAMIC_UNSTRUCTURED_CONNECTIVITY_ID_OFFSETS_REGEX "\\(0\\): 0, 0, 4, 4, 4")
endif()

# Keep structural checks visible even when the optional HDF5 tools are absent.
function(vtkhdf_add_h5dump_test name)
  if(H5DUMP_EXECUTABLE)
    add_test(NAME ${name} COMMAND "${H5DUMP_EXECUTABLE}" ${ARGN})
  else()
    add_test(NAME ${name}
      COMMAND ${CMAKE_COMMAND} -E echo "h5dump unavailable: skipping structural check")
    set_tests_properties(${name} PROPERTIES
      SKIP_REGULAR_EXPRESSION "^h5dump unavailable:")
  endif()
  set_tests_properties(${name} PROPERTIES LABELS h5dump)
endfunction()

vtkhdf_add_h5dump_test(vtkhdf_unstructured_static_cell_attribute
  -a /VTKHDF/CellData/static-cell-scalar/Attribute
  ${CMAKE_CURRENT_BINARY_DIR}/unstructured_test.vtkhdf
)
set_tests_properties(vtkhdf_unstructured_static_cell_attribute PROPERTIES
  DEPENDS vtkhdf_unstructured_fixture_file
  PASS_REGULAR_EXPRESSION "\"Scalars\""
)

vtkhdf_add_h5dump_test(vtkhdf_unstructured_static_point_attribute
  -a /VTKHDF/PointData/static-point-vector/Attribute
  ${CMAKE_CURRENT_BINARY_DIR}/unstructured_test.vtkhdf
)
set_tests_properties(vtkhdf_unstructured_static_point_attribute PROPERTIES
  DEPENDS vtkhdf_unstructured_fixture_file
  PASS_REGULAR_EXPRESSION "\"Vectors\""
)

vtkhdf_add_h5dump_test(vtkhdf_unstructured_temporal_cell_attribute
  -a /VTKHDF/CellData/cell-radius/Attribute
  ${CMAKE_CURRENT_BINARY_DIR}/unstructured_test.vtkhdf
)
set_tests_properties(vtkhdf_unstructured_temporal_cell_attribute PROPERTIES
  DEPENDS vtkhdf_unstructured_fixture_file
  PASS_REGULAR_EXPRESSION "\"Scalars\""
)

vtkhdf_add_h5dump_test(vtkhdf_unstructured_temporal_point_attribute
  -a /VTKHDF/PointData/point-velocity/Attribute
  ${CMAKE_CURRENT_BINARY_DIR}/unstructured_test.vtkhdf
)
set_tests_properties(vtkhdf_unstructured_temporal_point_attribute PROPERTIES
  DEPENDS vtkhdf_unstructured_fixture_file
  PASS_REGULAR_EXPRESSION "\"Vectors\""
)

vtkhdf_add_h5dump_test(vtkhdf_mb_static_cell_attribute
  -a /VTKHDF/Group_0/CellData/static-cell-scalar/Attribute
  ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
)
set_tests_properties(vtkhdf_mb_static_cell_attribute PROPERTIES
  DEPENDS vtkhdf_mb_fixture_file
  PASS_REGULAR_EXPRESSION "\"Scalars\""
)

vtkhdf_add_h5dump_test(vtkhdf_mb_temporal_point_attribute
  -a /VTKHDF/Group_1/PointData/point-velocity/Attribute
  ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
)
set_tests_properties(vtkhdf_mb_temporal_point_attribute PROPERTIES
  DEPENDS vtkhdf_mb_fixture_file
  PASS_REGULAR_EXPRESSION "\"Vectors\""
)

vtkhdf_add_h5dump_test(vtkhdf_mb_assembly_block_a_link
  -n ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
)
set_tests_properties(vtkhdf_mb_assembly_block_a_link PROPERTIES
  DEPENDS vtkhdf_mb_fixture_file
  PASS_REGULAR_EXPRESSION "link +/VTKHDF/Assembly/Block-A -> /VTKHDF/Group_0"
)

vtkhdf_add_h5dump_test(vtkhdf_mb_assembly_block_b_link
  -n ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
)
set_tests_properties(vtkhdf_mb_assembly_block_b_link PROPERTIES
  DEPENDS vtkhdf_mb_fixture_file
  PASS_REGULAR_EXPRESSION "link +/VTKHDF/Assembly/Block-B -> /VTKHDF/Group_1"
)

vtkhdf_add_h5dump_test(vtkhdf_mb_group_0_name
  -a /VTKHDF/Group_0/Name
  ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
)
set_tests_properties(vtkhdf_mb_group_0_name PROPERTIES
  DEPENDS vtkhdf_mb_fixture_file
  PASS_REGULAR_EXPRESSION "\"Block-A\""
)

vtkhdf_add_h5dump_test(vtkhdf_mb_group_1_name
  -a /VTKHDF/Group_1/Name
  ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf
)
set_tests_properties(vtkhdf_mb_group_1_name PROPERTIES
  DEPENDS vtkhdf_mb_fixture_file
  PASS_REGULAR_EXPRESSION "\"Block-B\""
)

vtkhdf_add_h5dump_test(vtkhdf_mb_duplicate_block_name_first_link
  -n ${CMAKE_CURRENT_BINARY_DIR}/mb_duplicate_block_name_test.vtkhdf
)
set_tests_properties(vtkhdf_mb_duplicate_block_name_first_link PROPERTIES
  DEPENDS vtkhdf_mb_duplicate_block_name_file_type
  PASS_REGULAR_EXPRESSION "link +/VTKHDF/Assembly/Repeated_Block -> /VTKHDF/Group_0"
)

vtkhdf_add_h5dump_test(vtkhdf_mb_duplicate_block_name_second_link
  -n ${CMAKE_CURRENT_BINARY_DIR}/mb_duplicate_block_name_test.vtkhdf
)
set_tests_properties(vtkhdf_mb_duplicate_block_name_second_link PROPERTIES
  DEPENDS vtkhdf_mb_duplicate_block_name_file_type
  PASS_REGULAR_EXPRESSION "link +/VTKHDF/Assembly/Repeated_Block_1 -> /VTKHDF/Group_2"
)

if(ENABLE_MPI)
  vtkhdf_add_h5dump_test(vtkhdf_mb_zero_part_sparse_number_of_parts
    -d /VTKHDF/Group_1/Steps/NumberOfParts
    ${CMAKE_CURRENT_BINARY_DIR}/mb_zero_part_test.vtkhdf
  )
  set_tests_properties(vtkhdf_mb_zero_part_sparse_number_of_parts PROPERTIES
    DEPENDS vtkhdf_mb_zero_part_file_type
    PASS_REGULAR_EXPRESSION "\\(0\\): 3, 3"
  )

  vtkhdf_add_h5dump_test(vtkhdf_mb_zero_part_sparse_number_of_points
    -d /VTKHDF/Group_1/NumberOfPoints
    ${CMAKE_CURRENT_BINARY_DIR}/mb_zero_part_test.vtkhdf
  )
  set_tests_properties(vtkhdf_mb_zero_part_sparse_number_of_points PROPERTIES
    DEPENDS vtkhdf_mb_zero_part_file_type
    PASS_REGULAR_EXPRESSION "\\(0\\): 8, 8, 8"
  )

  vtkhdf_add_h5dump_test(vtkhdf_mb_zero_part_sparse_point_data_offsets
    -d /VTKHDF/Group_1/Steps/PointDataOffsets/point-radius
    ${CMAKE_CURRENT_BINARY_DIR}/mb_zero_part_test.vtkhdf
  )
  set_tests_properties(vtkhdf_mb_zero_part_sparse_point_data_offsets PROPERTIES
    DEPENDS vtkhdf_mb_zero_part_file_type
    PASS_REGULAR_EXPRESSION "\\(0\\): 0, 24"
  )
endif()

vtkhdf_add_h5dump_test(vtkhdf_dynamic_unstructured_number_of_points
  -d /VTKHDF/NumberOfPoints
  ${CMAKE_CURRENT_BINARY_DIR}/dyn_unstructured_test.vtkhdf
)
set_tests_properties(vtkhdf_dynamic_unstructured_number_of_points PROPERTIES
  DEPENDS vtkhdf_dynamic_unstructured_file_type
  PASS_REGULAR_EXPRESSION "${VTKHDF_DYNAMIC_UNSTRUCTURED_NUMBER_OF_POINTS_REGEX}"
)

vtkhdf_add_h5dump_test(vtkhdf_dynamic_unstructured_number_of_cells
  -d /VTKHDF/NumberOfCells
  ${CMAKE_CURRENT_BINARY_DIR}/dyn_unstructured_test.vtkhdf
)
set_tests_properties(vtkhdf_dynamic_unstructured_number_of_cells PROPERTIES
  DEPENDS vtkhdf_dynamic_unstructured_file_type
  PASS_REGULAR_EXPRESSION "${VTKHDF_DYNAMIC_UNSTRUCTURED_NUMBER_OF_CELLS_REGEX}"
)

vtkhdf_add_h5dump_test(vtkhdf_dynamic_unstructured_point_offsets
  -d /VTKHDF/Steps/PointOffsets
  ${CMAKE_CURRENT_BINARY_DIR}/dyn_unstructured_test.vtkhdf
)
set_tests_properties(vtkhdf_dynamic_unstructured_point_offsets PROPERTIES
  DEPENDS vtkhdf_dynamic_unstructured_file_type
  PASS_REGULAR_EXPRESSION "${VTKHDF_DYNAMIC_UNSTRUCTURED_POINT_OFFSETS_REGEX}"
)

vtkhdf_add_h5dump_test(vtkhdf_dynamic_unstructured_cell_offsets
  -d /VTKHDF/Steps/CellOffsets
  ${CMAKE_CURRENT_BINARY_DIR}/dyn_unstructured_test.vtkhdf
)
set_tests_properties(vtkhdf_dynamic_unstructured_cell_offsets PROPERTIES
  DEPENDS vtkhdf_dynamic_unstructured_file_type
  PASS_REGULAR_EXPRESSION "${VTKHDF_DYNAMIC_UNSTRUCTURED_CELL_OFFSETS_REGEX}"
)

vtkhdf_add_h5dump_test(vtkhdf_dynamic_unstructured_connectivity_id_offsets
  -d /VTKHDF/Steps/ConnectivityIdOffsets
  ${CMAKE_CURRENT_BINARY_DIR}/dyn_unstructured_test.vtkhdf
)
set_tests_properties(vtkhdf_dynamic_unstructured_connectivity_id_offsets PROPERTIES
  DEPENDS vtkhdf_dynamic_unstructured_file_type
  PASS_REGULAR_EXPRESSION "${VTKHDF_DYNAMIC_UNSTRUCTURED_CONNECTIVITY_ID_OFFSETS_REGEX}"
)

vtkhdf_add_h5dump_test(vtkhdf_chunk_default_layout
  -pH -d /VTKHDF/Types
  ${CMAKE_CURRENT_BINARY_DIR}/chunk_default.vtkhdf
)
set_tests_properties(vtkhdf_chunk_default_layout PROPERTIES
  FIXTURES_REQUIRED vtkhdf_chunk_files
  PASS_REGULAR_EXPRESSION "CHUNKED \\( 131072 \\)"
)

vtkhdf_add_h5dump_test(vtkhdf_chunk_custom_layout
  -pH -d /VTKHDF/Group_0/Types
  ${CMAKE_CURRENT_BINARY_DIR}/chunk_custom.vtkhdf
)
set_tests_properties(vtkhdf_chunk_custom_layout PROPERTIES
  FIXTURES_REQUIRED vtkhdf_chunk_files
  PASS_REGULAR_EXPRESSION "CHUNKED \\( 1048576 \\)"
)

add_test(
  NAME vtkhdf_unstructured_reader_smoke
  COMMAND ${Python3_EXECUTABLE} ${VTKHDF_READER_SMOKE_SCRIPT}
          unstructured ${CMAKE_CURRENT_BINARY_DIR}/unstructured_test.vtkhdf
)
set_tests_properties(vtkhdf_unstructured_reader_smoke PROPERTIES
  DEPENDS vtkhdf_unstructured_fixture_file
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

add_test(
  NAME vtkhdf_mb_reader_temporal_state
  COMMAND ${Python3_EXECUTABLE} ${VTKHDF_READER_SMOKE_SCRIPT}
          mb ${CMAKE_CURRENT_BINARY_DIR}/mb_test.vtkhdf --check-temporal-state
)
set_tests_properties(vtkhdf_mb_reader_temporal_state PROPERTIES
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

# ImageData fixtures and contract failures.
add_fypp_sources(IMAGE_TEST_SRC SOURCES vtkhdf_image_test.F90.fypp FLAGS --line-numbering)
vtkhdf_add_test_target(vtkhdf_image_test "${IMAGE_TEST_SRC}")
vtkhdf_add_runtime_test(vtkhdf_image_file_type vtkhdf_image_test)
set_tests_properties(vtkhdf_image_file_type PROPERTIES FIXTURES_SETUP vtkhdf_image_files)
vtkhdf_add_test_target(vtkhdf_image_fail_test vtkhdf_image_fail_test.F90)
foreach(case bounds source type components late outside_step missing field_duplicate
             before_mesh register_before_mesh field_outside_step wrong_handle static_type
             static_components mode_moving mode_dynamic mode_invalid)
  vtkhdf_add_runtime_test(vtkhdf_image_fail_${case} vtkhdf_image_fail_test ${case})
  set_tests_properties(vtkhdf_image_fail_${case} PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
endforeach()
# Require the duplicate-write assertion, rather than accepting a setup failure.
set(unstructured_duplicate_command $<TARGET_FILE:vtkhdf_image_fail_test>)
set(unstructured_duplicate_processors 1)
if(ENABLE_MPI)
  list(PREPEND unstructured_duplicate_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4)
  list(APPEND unstructured_duplicate_command ${MPIEXEC_POSTFLAGS})
  set(unstructured_duplicate_processors 4)
endif()
list(APPEND unstructured_duplicate_command unstructured_duplicate)
add_test(NAME vtkhdf_image_fail_unstructured_duplicate
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_expect_assertion.py
          --source ${PROJECT_SOURCE_DIR}/src/vtkhdf_unstructured_writer_type.F90.fypp
          --assertion "no double writing in a single step"
          -- ${unstructured_duplicate_command})
set_tests_properties(vtkhdf_image_fail_unstructured_duplicate PROPERTIES
  PROCESSORS ${unstructured_duplicate_processors}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  TIMEOUT 30)
# Rank-disagreement checks use ASSERT and are absent with NDEBUG.
if(ENABLE_MPI)
  foreach(case geometry_mismatch name_mismatch mode_mismatch)
    set(test_name vtkhdf_image_fail_${case})
    if(CMAKE_CONFIGURATION_TYPES)
      add_test(NAME ${test_name}
        COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_image_fail_test ${MPIEXEC_POSTFLAGS} ${case}
        CONFIGURATIONS Debug)
      set_tests_properties(${test_name} PROPERTIES
        PROCESSORS 4 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
      vtkhdf_add_runtime_test(${test_name} vtkhdf_image_fail_test ${case})
    else()
      continue()
    endif()
    set_tests_properties(${test_name} PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
  endforeach()
endif()

set(VTKHDF_DEBUG_FILL "$<IF:$<CONFIG:Debug>,1,0>")
add_test(NAME vtkhdf_image_layout
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_image_check.py
          layout ${CMAKE_CURRENT_BINARY_DIR} --debug-fill ${VTKHDF_DEBUG_FILL})
set_tests_properties(vtkhdf_image_layout PROPERTIES
  FIXTURES_REQUIRED vtkhdf_image_files SKIP_RETURN_CODE 77)
foreach(case static temporal temporal_fielddata registration_names numeric dim_1 dim_2 dim_3 dim_4 dim_5 dim_6 dim_7)
  add_test(NAME vtkhdf_image_reader_${case}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_image_check.py
            reader ${CMAKE_CURRENT_BINARY_DIR} ${case})
  set_tests_properties(vtkhdf_image_reader_${case} PROPERTIES
    FIXTURES_REQUIRED vtkhdf_image_files SKIP_RETURN_CODE 77)
endforeach()

if(ENABLE_MPI)
  # Two ranks own unequal numbers of blocks, including multiple blocks on
  # one rank. Keep these output files separate from the four-rank fixture.
  set(IMAGE_TWO_RANK_DIR ${CMAKE_CURRENT_BINARY_DIR}/image-two-rank)
  file(MAKE_DIRECTORY ${IMAGE_TWO_RANK_DIR})
  add_test(NAME vtkhdf_image_two_rank
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2
            ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_image_test ${MPIEXEC_POSTFLAGS})
  set_tests_properties(vtkhdf_image_two_rank PROPERTIES PROCESSORS 2
    WORKING_DIRECTORY ${IMAGE_TWO_RANK_DIR} FIXTURES_SETUP vtkhdf_image_two_rank_files)
  add_test(NAME vtkhdf_image_two_rank_layout
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_image_check.py
            layout ${IMAGE_TWO_RANK_DIR} --debug-fill ${VTKHDF_DEBUG_FILL})
  set_tests_properties(vtkhdf_image_two_rank_layout PROPERTIES
    FIXTURES_REQUIRED vtkhdf_image_two_rank_files SKIP_RETURN_CODE 77)
endif()

# Explicit structured geometry: coordinate kinds, distributed regions, and time reuse.
add_fypp_sources(GRID_TEST_SRC SOURCES vtkhdf_grid_test.F90.fypp FLAGS --line-numbering)
vtkhdf_add_test_target(vtkhdf_grid_test "${GRID_TEST_SRC}")
vtkhdf_add_runtime_test(vtkhdf_grid_file_types vtkhdf_grid_test)
set_tests_properties(vtkhdf_grid_file_types PROPERTIES FIXTURES_SETUP vtkhdf_grid_files)
add_test(NAME vtkhdf_grid_layout
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_grid_check.py
          layout ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(vtkhdf_grid_layout PROPERTIES
  FIXTURES_REQUIRED vtkhdf_grid_files SKIP_RETURN_CODE 77)
foreach(mask RANGE 1 7)
  add_test(NAME vtkhdf_grid_reader_dim_${mask}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_grid_check.py
            reader ${CMAKE_CURRENT_BINARY_DIR} ${mask})
  set_tests_properties(vtkhdf_grid_reader_dim_${mask} PROPERTIES
    FIXTURES_REQUIRED vtkhdf_grid_files SKIP_RETURN_CODE 77)
endforeach()
vtkhdf_add_test_target(vtkhdf_grid_fail_test vtkhdf_grid_fail_test.F90)
foreach(case axis_length before_mesh missing_close missing_step components bounds source kind late duplicate
             cell_before_geometry point_before_geometry temporal_before_geometry geometry_after_data)
  vtkhdf_add_runtime_test(vtkhdf_grid_fail_${case} vtkhdf_grid_fail_test ${case})
  set_tests_properties(vtkhdf_grid_fail_${case} PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
endforeach()

vtkhdf_add_test_target(vtkhdf_simple_dataset_test vtkhdf_simple_dataset_test.F90)
vtkhdf_add_runtime_test(vtkhdf_simple_dataset_lifecycle vtkhdf_simple_dataset_test)

add_fypp_sources(UNSTRUCTURED_CONSTANT_DATA_TEST_SRC
  SOURCES vtkhdf_unstructured_constant_data_test.F90.fypp FLAGS --line-numbering)
vtkhdf_add_test_target(vtkhdf_unstructured_constant_data_test "${UNSTRUCTURED_CONSTANT_DATA_TEST_SRC}")
vtkhdf_add_runtime_test(vtkhdf_unstructured_constant_data_files vtkhdf_unstructured_constant_data_test)
set_tests_properties(vtkhdf_unstructured_constant_data_files PROPERTIES FIXTURES_SETUP vtkhdf_unstructured_constant_files)
foreach(check layout reader)
  add_test(NAME vtkhdf_unstructured_constant_data_${check}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_unstructured_constant_data_check.py
            ${check} ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(vtkhdf_unstructured_constant_data_${check} PROPERTIES
    FIXTURES_REQUIRED vtkhdf_unstructured_constant_files SKIP_RETURN_CODE 77)
endforeach()
vtkhdf_add_test_target(vtkhdf_unstructured_constant_data_fail_test vtkhdf_unstructured_constant_data_fail_test.F90)
foreach(case dynamic_point dynamic_cell wrong_point_count wrong_cell_count missing_point missing_cell)
  vtkhdf_add_runtime_test(vtkhdf_unstructured_constant_data_fail_${case}
    vtkhdf_unstructured_constant_data_fail_test ${case})
  set_tests_properties(vtkhdf_unstructured_constant_data_fail_${case} PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
endforeach()

add_fypp_sources(CONSTANT_FIELD_TEST_SRC SOURCES vtkhdf_constant_field_test.F90.fypp FLAGS --line-numbering)
vtkhdf_add_test_target(vtkhdf_constant_field_test "${CONSTANT_FIELD_TEST_SRC}")
vtkhdf_add_runtime_test(vtkhdf_constant_field_files vtkhdf_constant_field_test)
set_tests_properties(vtkhdf_constant_field_files PROPERTIES FIXTURES_SETUP vtkhdf_constant_field_files)
foreach(check layout reader)
  add_test(NAME vtkhdf_constant_field_${check}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_constant_field_check.py
            ${check} ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(vtkhdf_constant_field_${check} PROPERTIES
    FIXTURES_REQUIRED vtkhdf_constant_field_files SKIP_RETURN_CODE 77)
endforeach()

# Moving geometry uses snapshots independently of time steps and data writes.
add_fypp_sources(MOVING_GRID_TEST_SRC SOURCES vtkhdf_moving_grid_test.F90.fypp FLAGS --line-numbering)
vtkhdf_add_test_target(vtkhdf_moving_grid_test "${MOVING_GRID_TEST_SRC}")
vtkhdf_add_runtime_test(vtkhdf_moving_grid_file_types vtkhdf_moving_grid_test)
set_tests_properties(vtkhdf_moving_grid_file_types PROPERTIES FIXTURES_SETUP vtkhdf_moving_grid_files)
add_test(NAME vtkhdf_moving_grid_layout
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_moving_grid_check.py
          layout ${CMAKE_CURRENT_BINARY_DIR} --debug-fill ${VTKHDF_DEBUG_FILL})
set_tests_properties(vtkhdf_moving_grid_layout PROPERTIES
  FIXTURES_REQUIRED vtkhdf_moving_grid_files SKIP_RETURN_CODE 77)
foreach(mask RANGE 1 7)
  add_test(NAME vtkhdf_moving_grid_reader_dim_${mask}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_moving_grid_check.py
            reader ${CMAKE_CURRENT_BINARY_DIR} ${mask} --debug-fill ${VTKHDF_DEBUG_FILL})
  set_tests_properties(vtkhdf_moving_grid_reader_dim_${mask} PROPERTIES
    FIXTURES_REQUIRED vtkhdf_moving_grid_files SKIP_RETURN_CODE 77)
endforeach()
if(ENABLE_MPI)
  set(MOVING_GRID_TWO_RANK_DIR ${CMAKE_CURRENT_BINARY_DIR}/moving-grid-two-rank)
  file(MAKE_DIRECTORY ${MOVING_GRID_TWO_RANK_DIR})
  add_test(NAME vtkhdf_moving_grid_two_rank
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2
            ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_moving_grid_test ${MPIEXEC_POSTFLAGS})
  set_tests_properties(vtkhdf_moving_grid_two_rank PROPERTIES PROCESSORS 2
    WORKING_DIRECTORY ${MOVING_GRID_TWO_RANK_DIR} FIXTURES_SETUP vtkhdf_moving_grid_two_rank_files)
  add_test(NAME vtkhdf_moving_grid_two_rank_layout
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_moving_grid_check.py
            layout ${MOVING_GRID_TWO_RANK_DIR} --debug-fill ${VTKHDF_DEBUG_FILL})
  set_tests_properties(vtkhdf_moving_grid_two_rank_layout PROPERTIES
    FIXTURES_REQUIRED vtkhdf_moving_grid_two_rank_files SKIP_RETURN_CODE 77)
endif()
vtkhdf_add_test_target(vtkhdf_moving_grid_fail_test vtkhdf_moving_grid_fail_test.F90)
foreach(grid rectilinear structured)
  foreach(case before_grid missing data_before_geometry after_data missing_data between_steps kind bounds_change dynamic)
    set(test_name vtkhdf_moving_${grid}_fail_${case})
    vtkhdf_add_runtime_test(${test_name} vtkhdf_moving_grid_fail_test ${grid} ${case})
    set_tests_properties(${test_name} PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
  endforeach()
endforeach()
foreach(case duplicate axis_length)
  set(test_name vtkhdf_moving_rectilinear_fail_${case})
  vtkhdf_add_runtime_test(${test_name} vtkhdf_moving_grid_fail_test rectilinear ${case})
  set_tests_properties(${test_name} PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
endforeach()

# HyperTreeGrid: breadth-first tree descriptors, per-part coordinates, and a mask.
vtkhdf_add_test_target(vtkhdf_hypertree_test vtkhdf_hypertree_test.F90)
vtkhdf_add_runtime_test(vtkhdf_hypertree_file_types vtkhdf_hypertree_test)
set_tests_properties(vtkhdf_hypertree_file_types PROPERTIES FIXTURES_SETUP vtkhdf_hypertree_files)
foreach(check layout reader)
  add_test(NAME vtkhdf_hypertree_${check}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_hypertree_check.py
            ${check} ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(vtkhdf_hypertree_${check} PROPERTIES
    FIXTURES_REQUIRED vtkhdf_hypertree_files SKIP_RETURN_CODE 77)
endforeach()
if(ENABLE_MPI)
  set(HYPERTREE_TWO_RANK_DIR ${CMAKE_CURRENT_BINARY_DIR}/hypertree-two-rank)
  file(MAKE_DIRECTORY ${HYPERTREE_TWO_RANK_DIR})
  add_test(NAME vtkhdf_hypertree_two_rank
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2
            ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_hypertree_test ${MPIEXEC_POSTFLAGS})
  set_tests_properties(vtkhdf_hypertree_two_rank PROPERTIES PROCESSORS 2
    WORKING_DIRECTORY ${HYPERTREE_TWO_RANK_DIR} FIXTURES_SETUP vtkhdf_hypertree_two_rank_files)
  foreach(check layout reader)
    add_test(NAME vtkhdf_hypertree_two_rank_${check}
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_hypertree_check.py
              ${check} ${HYPERTREE_TWO_RANK_DIR})
    set_tests_properties(vtkhdf_hypertree_two_rank_${check} PROPERTIES
      FIXTURES_REQUIRED vtkhdf_hypertree_two_rank_files SKIP_RETURN_CODE 77)
  endforeach()
endif()
vtkhdf_add_test_target(vtkhdf_hypertree_fail_test vtkhdf_hypertree_fail_test.F90)
foreach(case no_grid bad_dims bad_branch_factor grid_twice axis_length tree_id_range missing_cells extra_cells
             masked_refined topology_twice data_before_topology data_before_geometry wrong_cell_count
             geometry_change late_register mesh_after_data static_data missing_data
             invalid_mode cell_outside_step temporal_count missing_initial
             field_outside_step field_late_register)
  vtkhdf_add_runtime_test(vtkhdf_hypertree_fail_${case} vtkhdf_hypertree_fail_test ${case})
  set_tests_properties(vtkhdf_hypertree_fail_${case} PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
endforeach()

# Numeric overloads, cell-name/handle isolation, and HTG field lifecycle wiring.
add_fypp_sources(HYPERTREE_API_TEST_SRC SOURCES vtkhdf_hypertree_api_test.F90.fypp FLAGS --line-numbering)
vtkhdf_add_test_target(vtkhdf_hypertree_api_test "${HYPERTREE_API_TEST_SRC}")
vtkhdf_add_runtime_test(vtkhdf_hypertree_api_files vtkhdf_hypertree_api_test)
set_tests_properties(vtkhdf_hypertree_api_files PROPERTIES FIXTURES_SETUP vtkhdf_hypertree_api_files)
foreach(check layout reader)
  add_test(NAME vtkhdf_hypertree_api_${check}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_hypertree_api_check.py
            ${check} ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(vtkhdf_hypertree_api_${check} PROPERTIES
    FIXTURES_REQUIRED vtkhdf_hypertree_api_files SKIP_RETURN_CODE 77)
endforeach()

# These checks require the intended diagnostic, not merely a failing setup.
foreach(case cell_foreign field_foreign cell_duplicate)
  set(hypertree_fail_command $<TARGET_FILE:vtkhdf_hypertree_fail_test>)
  set(hypertree_fail_processors 1)
  if(ENABLE_MPI)
    list(PREPEND hypertree_fail_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4)
    list(APPEND hypertree_fail_command ${MPIEXEC_POSTFLAGS})
    set(hypertree_fail_processors 4)
  endif()
  list(APPEND hypertree_fail_command ${case})
  set(hypertree_fail_source ${PROJECT_SOURCE_DIR}/src/vtkhdf_hypertree_writer_type.F90.fypp)
  if(case STREQUAL "cell_duplicate")
    set(hypertree_fail_assertion "no double writing in a single step")
  elseif(case STREQUAL "cell_foreign")
    set(hypertree_fail_assertion "INSIST(index >= 1)")
  else()
    set(hypertree_fail_source ${PROJECT_SOURCE_DIR}/src/vtkhdf_field_data_type.F90.fypp)
    set(hypertree_fail_assertion "INSIST(index >= 1 .and. index <= size(this%vars))")
  endif()
  add_test(NAME vtkhdf_hypertree_fail_${case}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_expect_assertion.py
            --source ${hypertree_fail_source} --assertion "${hypertree_fail_assertion}" -- ${hypertree_fail_command})
  set_tests_properties(vtkhdf_hypertree_fail_${case} PROPERTIES
    PROCESSORS ${hypertree_fail_processors} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} TIMEOUT 30)
endforeach()
# Agreement diagnostics are compiled out with NDEBUG.
if(ENABLE_MPI)
  foreach(case mode_mismatch dims_mismatch branch_mismatch)
    set(test_name vtkhdf_hypertree_fail_${case})
    if(CMAKE_CONFIGURATION_TYPES)
      add_test(NAME ${test_name}
        COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                ${CMAKE_CURRENT_BINARY_DIR}/vtkhdf_hypertree_fail_test ${MPIEXEC_POSTFLAGS} ${case}
        CONFIGURATIONS Debug)
      set_tests_properties(${test_name} PROPERTIES
        PROCESSORS 4 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
      vtkhdf_add_runtime_test(${test_name} vtkhdf_hypertree_fail_test ${case})
    else()
      continue()
    endif()
    set_tests_properties(${test_name} PROPERTIES WILL_FAIL TRUE TIMEOUT 30)
  endforeach()
endif()

# Material-interface pairs use the ordinary cell registry with private handles.
add_fypp_sources(HYPERTREE_INTERFACE_TEST_SRC SOURCES vtkhdf_hypertree_interface_test.F90.fypp FLAGS --line-numbering)
vtkhdf_add_test_target(vtkhdf_hypertree_interface_test "${HYPERTREE_INTERFACE_TEST_SRC}")
vtkhdf_add_test_target(vtkhdf_hypertree_interface_fail_test vtkhdf_hypertree_interface_fail_test.F90)
vtkhdf_add_runtime_test(vtkhdf_hypertree_interface_files vtkhdf_hypertree_interface_test)
set_tests_properties(vtkhdf_hypertree_interface_files PROPERTIES FIXTURES_SETUP vtkhdf_hypertree_interface_files)
foreach(check layout reader)
  add_test(NAME vtkhdf_hypertree_interface_${check}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_hypertree_interface_check.py
            ${check} ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(vtkhdf_hypertree_interface_${check} PROPERTIES
    FIXTURES_REQUIRED vtkhdf_hypertree_interface_files SKIP_RETURN_CODE 77)
endforeach()
if(ENABLE_MPI)
  set(interface_two_rank_dir ${CMAKE_CURRENT_BINARY_DIR}/interface-two-rank)
  file(MAKE_DIRECTORY ${interface_two_rank_dir})
  add_test(NAME vtkhdf_hypertree_interface_two_rank
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2
            $<TARGET_FILE:vtkhdf_hypertree_interface_test> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(vtkhdf_hypertree_interface_two_rank PROPERTIES PROCESSORS 2
    WORKING_DIRECTORY ${interface_two_rank_dir} FIXTURES_SETUP vtkhdf_hypertree_interface_two_rank_files)
  foreach(check layout reader)
    add_test(NAME vtkhdf_hypertree_interface_two_rank_${check}
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_hypertree_interface_check.py
              ${check} ${interface_two_rank_dir})
    set_tests_properties(vtkhdf_hypertree_interface_two_rank_${check} PROPERTIES
      FIXTURES_REQUIRED vtkhdf_hypertree_interface_two_rank_files SKIP_RETURN_CODE 77)
  endforeach()
endif()

set(interface_failures
  "closed_write|interface write requires an open writer"
  "closed_register|interface registration requires an open writer"
  "static_register|INSIST(this%is_temporal)"
  "dynamic_static|static interface mode"
  "static_twice|INSIST(this%interface_state /= INTERFACE_STATIC)"
  "register_twice|INSIST(this%interface_state == INTERFACE_UNDEFINED)"
  "register_after_static|INSIST(this%interface_state == INTERFACE_UNDEFINED)"
  "late_register|interface registration precedes time stepping"
  "outside_step|temporal interface requires an open step"
  "between_steps|temporal interface requires an open step"
  "no_grid|INSIST(this%grid_defined .and. this%geometry_defined .and. this%topology_defined)"
  "no_geometry|INSIST(this%grid_defined .and. this%geometry_defined .and. this%topology_defined)"
  "no_topology|INSIST(this%grid_defined .and. this%geometry_defined .and. this%topology_defined)"
  "missing_initial|initial interface pair required, even with no cells"
  "duplicate|interface pair already written"
  "normals_components|INSIST(size(normals,1) == 3)"
  "intercepts_components|INSIST(size(intercepts,1) == 3)"
  "normals_count|INSIST(size(normals,2,kind=int64) == this%ncell)"
  "intercepts_count|INSIST(size(intercepts,2,kind=int64) == this%ncell)"
  "temporal_normals_components|INSIST(size(normals,1) == 3)"
  "temporal_intercepts_count|INSIST(size(intercepts,2,kind=int64) == this%ncell)"
  "kind_64_to_32|INSIST(this%interface_kind =="
  "kind_32_to_64|INSIST(this%interface_kind =="
  "static_barrier|geometry cannot follow cell data"
  "temporal_barrier|topology cannot follow cell data"
  "changed_count|last_entity_count == this%ncell")
if(ENABLE_MPI)
  list(APPEND interface_failures "missing_initial_empty|initial interface pair required, even with no cells")
endif()
set(interface_debug_failures)
if(ENABLE_MPI AND (CMAKE_CONFIGURATION_TYPES OR CMAKE_BUILD_TYPE STREQUAL "Debug"))
  list(APPEND interface_debug_failures
    "kind_agreement_register|interface registration collective agreement"
    "kind_agreement_write|interface write collective agreement"
    "operation_agreement|interface write collective agreement")
endif()
foreach(failure IN LISTS interface_failures interface_debug_failures)
  string(REPLACE "|" ";" fields "${failure}")
  list(GET fields 0 scenario)
  list(GET fields 1 assertion)
  set(command $<TARGET_FILE:vtkhdf_hypertree_interface_fail_test>)
  set(processors 1)
  if(ENABLE_MPI)
    list(PREPEND command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4)
    list(APPEND command ${MPIEXEC_POSTFLAGS})
    set(processors 4)
  endif()
  list(APPEND command ${scenario})
  set(interface_test_configs)
  if(CMAKE_CONFIGURATION_TYPES AND failure IN_LIST interface_debug_failures)
    set(interface_test_configs CONFIGURATIONS Debug)
  endif()
  add_test(NAME vtkhdf_hypertree_interface_fail_${scenario}
    ${interface_test_configs}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/vtkhdf_expect_assertion.py
            --source ${PROJECT_SOURCE_DIR}/src/vtkhdf_hypertree_writer_type.F90.fypp
            --assertion "${assertion}" --diagnostic-glob "interface-assertion-R*.log" -- ${command})
  set(dir ${CMAKE_CURRENT_BINARY_DIR}/interface-fail-${scenario})
  file(MAKE_DIRECTORY ${dir})
  set_tests_properties(vtkhdf_hypertree_interface_fail_${scenario} PROPERTIES
    PROCESSORS ${processors} WORKING_DIRECTORY ${dir} TIMEOUT 30)
endforeach()

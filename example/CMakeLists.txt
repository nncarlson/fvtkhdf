add_executable(unstructured_static unstructured_static.F90)
target_link_libraries(unstructured_static PRIVATE fvtkhdf)

add_executable(unstructured_fixed_mesh unstructured_fixed_mesh.F90)
target_link_libraries(unstructured_fixed_mesh PRIVATE fvtkhdf)

add_executable(unstructured_moving_mesh unstructured_moving_mesh.F90)
target_link_libraries(unstructured_moving_mesh PRIVATE fvtkhdf)

add_executable(unstructured_dynamic_mesh unstructured_dynamic_mesh.F90)
target_link_libraries(unstructured_dynamic_mesh PRIVATE fvtkhdf)

add_executable(mb_example1 mb_example1.F90)
target_link_libraries(mb_example1 PRIVATE fvtkhdf)

add_executable(pdc_example1 pdc_example1.F90)
target_link_libraries(pdc_example1 PRIVATE fvtkhdf)

add_executable(image_fixed_mesh image_fixed_mesh.F90)
target_link_libraries(image_fixed_mesh PRIVATE fvtkhdf)

foreach(grid rectilinear structured)
  add_executable(${grid}_fixed_mesh ${grid}_fixed_mesh.F90)
  target_link_libraries(${grid}_fixed_mesh PRIVATE fvtkhdf)
  add_executable(${grid}_moving_mesh ${grid}_moving_mesh.F90)
  target_link_libraries(${grid}_moving_mesh PRIVATE fvtkhdf)
endforeach()

add_executable(hypertree_dynamic_mesh hypertree_dynamic_mesh.F90)
target_link_libraries(hypertree_dynamic_mesh PRIVATE fvtkhdf)

add_executable(composite_assembly composite_assembly.F90)
target_link_libraries(composite_assembly PRIVATE fvtkhdf)

if(BUILD_TESTING)
  if(ENABLE_MPI)
    set(circle_ranks 2 4)
  else()
    set(circle_ranks 1)
  endif()
  foreach(ranks IN LISTS circle_ranks)
    set(name hypertree_dynamic_mesh_example_${ranks}_rank)
    set(dir ${CMAKE_CURRENT_BINARY_DIR}/circle-${ranks}-rank)
    file(MAKE_DIRECTORY ${dir})
    set(command $<TARGET_FILE:hypertree_dynamic_mesh>)
    if(ENABLE_MPI)
      list(PREPEND command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${ranks})
      list(APPEND command ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(NAME ${name} COMMAND ${command})
    set_tests_properties(${name} PROPERTIES PROCESSORS ${ranks}
      WORKING_DIRECTORY ${dir} FIXTURES_SETUP ${name}_file)
    add_test(NAME ${name}_geometry
      COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/vtkhdf_hypertree_example_check.py
              ${dir}/hypertree_dynamic_mesh.vtkhdf)
    set_tests_properties(${name}_geometry PROPERTIES
      FIXTURES_REQUIRED ${name}_file SKIP_RETURN_CODE 77)
  endforeach()
endif()

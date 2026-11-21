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

add_executable(ug_static ug_static.F90)
target_link_libraries(ug_static PRIVATE fvtkhdf)

add_executable(ug_fixed_mesh ug_fixed_mesh.F90)
target_link_libraries(ug_fixed_mesh PRIVATE fvtkhdf)

add_executable(ug_moving_mesh ug_moving_mesh.F90)
target_link_libraries(ug_moving_mesh PRIVATE fvtkhdf)

add_executable(ug_dynamic_mesh ug_dynamic_mesh.F90)
target_link_libraries(ug_dynamic_mesh PRIVATE fvtkhdf)

add_executable(mb_example1 mb_example1.F90)
target_link_libraries(mb_example1 PRIVATE fvtkhdf)

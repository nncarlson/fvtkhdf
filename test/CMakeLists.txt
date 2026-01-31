add_executable(vtkhdf_demo vtkhdf_demo.F90)
target_link_libraries(vtkhdf_demo PRIVATE fvtkhdf)

add_test(vtkhdf_file_type ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./vtkhdf_demo ${MPIEXEC_POSTFLAGS})
set_tests_properties(vtkhdf_file_type PROPERTIES PROCESSORS 4)

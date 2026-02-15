# CMake generated Testfile for
# Source directory: /app/test
# Build directory: /app/build/test
#
# This file includes the relevant testing commands required for
# testing this directory and lists subdirectories to be tested as well.
add_test([=[vtkhdf_ug_file_type]=] "./vtkhdf_ug_test")
set_tests_properties([=[vtkhdf_ug_file_type]=] PROPERTIES  PROCESSORS "1" _BACKTRACE_TRIPLES "/app/test/CMakeLists.txt;17;add_test;/app/test/CMakeLists.txt;0;")
add_test([=[vtkhdf_mb_file_type]=] "./vtkhdf_mb_test")
set_tests_properties([=[vtkhdf_mb_file_type]=] PROPERTIES  PROCESSORS "1" _BACKTRACE_TRIPLES "/app/test/CMakeLists.txt;23;add_test;/app/test/CMakeLists.txt;0;")
add_test([=[test_leak]=] "/app/build/test/test_leak")
set_tests_properties([=[test_leak]=] PROPERTIES  _BACKTRACE_TRIPLES "/app/test/CMakeLists.txt;29;add_test;/app/test/CMakeLists.txt;0;")

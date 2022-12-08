# CMake generated Testfile for 
# Source directory: /Users/touzani/ofeli/demos/1D
# Build directory: /Users/touzani/ofeli/build/demos/1D
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(1D-1 "elliptic")
set_tests_properties(1D-1 PROPERTIES  _BACKTRACE_TRIPLES "/Users/touzani/ofeli/demos/1D/CMakeLists.txt;39;add_test;/Users/touzani/ofeli/demos/1D/CMakeLists.txt;0;")
add_test(1D-2 "heat" "10" "0.1")
set_tests_properties(1D-2 PROPERTIES  _BACKTRACE_TRIPLES "/Users/touzani/ofeli/demos/1D/CMakeLists.txt;40;add_test;/Users/touzani/ofeli/demos/1D/CMakeLists.txt;0;")
add_test(1D-3 "transport" "10" "0.1")
set_tests_properties(1D-3 PROPERTIES  _BACKTRACE_TRIPLES "/Users/touzani/ofeli/demos/1D/CMakeLists.txt;41;add_test;/Users/touzani/ofeli/demos/1D/CMakeLists.txt;0;")

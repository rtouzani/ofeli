# CMake generated Testfile for 
# Source directory: /Users/touzani/ofeli/util/cfield
# Build directory: /Users/touzani/ofeli/build/util/cfield
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(cfield "cfield" "-f" "gmsh" "-m" "cavity.m" "-i" "cavity.s" "-o" "cavity.pos")
set_tests_properties(cfield PROPERTIES  _BACKTRACE_TRIPLES "/Users/touzani/ofeli/util/cfield/CMakeLists.txt;41;add_test;/Users/touzani/ofeli/util/cfield/CMakeLists.txt;0;")

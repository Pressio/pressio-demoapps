# CMake generated Testfile for 
# Source directory: /Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/sod1d/s7
# Build directory: /Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/buildClang2/tests_cpp/sod1d/s7
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(sod1d_s7 "/usr/local/Cellar/cmake/3.18.4/bin/cmake" "-DMESHDRIVER=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/../meshing_scripts/create_full_mesh.py" "-DOUTDIR=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/buildClang2/tests_cpp/sod1d/s7" "-DEXENAME=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/buildClang2/tests_cpp/sod1d/s7/sod1dexeS7" "-DSTENCILVAL=7" "-P" "/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/sod1d/s7/../test.cmake")
set_tests_properties(sod1d_s7 PROPERTIES  _BACKTRACE_TRIPLES "/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/sod1d/s7/CMakeLists.txt;8;add_test;/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/sod1d/s7/CMakeLists.txt;0;")

# CMake generated Testfile for 
# Source directory: /Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests/riemann2d
# Build directory: /Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/build2/tests/riemann2d
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(riemann2d_s3 "/usr/local/Cellar/cmake/3.18.4/bin/cmake" "-DMESHDRIVER=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests/../meshing_scripts/create_full_mesh.py" "-DOUTDIR=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/build2/tests/riemann2d" "-DEXENAME=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/build2/tests/riemann2d/riemannexe" "-DSTENCILVAL=3" "-P" "/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests/riemann2d/test.cmake")
set_tests_properties(riemann2d_s3 PROPERTIES  _BACKTRACE_TRIPLES "/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests/riemann2d/CMakeLists.txt;11;add_test;/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests/riemann2d/CMakeLists.txt;0;")
add_test(riemann2d_s7 "/usr/local/Cellar/cmake/3.18.4/bin/cmake" "-DMESHDRIVER=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests/../meshing_scripts/create_full_mesh.py" "-DOUTDIR=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/build2/tests/riemann2d" "-DEXENAME=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/build2/tests/riemann2d/riemannexe" "-DSTENCILVAL=7" "-P" "/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests/riemann2d/test.cmake")
set_tests_properties(riemann2d_s7 PROPERTIES  _BACKTRACE_TRIPLES "/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests/riemann2d/CMakeLists.txt;20;add_test;/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests/riemann2d/CMakeLists.txt;0;")

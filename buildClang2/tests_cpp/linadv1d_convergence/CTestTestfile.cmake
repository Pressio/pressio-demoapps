# CMake generated Testfile for 
# Source directory: /Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/linadv1d_convergence
# Build directory: /Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/buildClang2/tests_cpp/linadv1d_convergence
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(linadv1d_periodic "/usr/local/Cellar/cmake/3.18.4/bin/cmake" "-DMESHDRIVER=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/../meshing_scripts/create_full_mesh.py" "-DTESTNAME=linadv1d_periodic" "-DOUTDIR=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/buildClang2/tests_cpp/linadv1d_convergence" "-DEXENAME=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/buildClang2/tests_cpp/linadv1d_convergence/ladvexe" "-P" "/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/linadv1d_convergence/test.cmake")
set_tests_properties(linadv1d_periodic PROPERTIES  _BACKTRACE_TRIPLES "/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/linadv1d_convergence/CMakeLists.txt;6;add_test;/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/linadv1d_convergence/CMakeLists.txt;0;")

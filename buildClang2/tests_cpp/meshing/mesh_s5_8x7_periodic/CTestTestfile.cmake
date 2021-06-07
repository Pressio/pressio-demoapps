# CMake generated Testfile for 
# Source directory: /Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/meshing/mesh_s5_8x7_periodic
# Build directory: /Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/buildClang2/tests_cpp/meshing/mesh_s5_8x7_periodic
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(full_mesh_8x7_s5_periodic "/usr/local/Cellar/cmake/3.18.4/bin/cmake" "-DMESHDRIVER=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/meshing/../../meshing_scripts/create_full_mesh.py" "-DOUTDIR=/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/buildClang2/tests_cpp/meshing/mesh_s5_8x7_periodic" "-Dnx=8" "-Dny=7" "-Dss=5" "-P" "/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/meshing/mesh_s5_8x7_periodic/../test2d.cmake")
set_tests_properties(full_mesh_8x7_s5_periodic PROPERTIES  _BACKTRACE_TRIPLES "/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/meshing/mesh_s5_8x7_periodic/CMakeLists.txt;8;add_test;/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio-demoapps/tests_cpp/meshing/mesh_s5_8x7_periodic/CMakeLists.txt;0;")

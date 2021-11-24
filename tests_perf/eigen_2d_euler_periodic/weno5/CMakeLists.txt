
set(testname perf_test_eigen_2d_euler_periodic_weno5)
set(exename  ${testname}_exe)
add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(${exename} PUBLIC -DUSE_WENO5)

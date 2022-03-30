
set(testname eigen_1d_linadv_default_velo_jacobian_weno5_fd)
set(exename  ${testname}_exe)

configure_file(../plot.py plot.py COPYONLY)
add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(${exename} PUBLIC -DUSE_WENO5)

add_test(NAME ${testname}
  COMMAND ${CMAKE_COMMAND}
  -DMESHDRIVER=${MESHSRC}/create_full_mesh.py
  -DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
  -DEXENAME=$<TARGET_FILE:${exename}>
  -DSTENCILVAL=7
  -P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
  )

set_tests_properties(${testname}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

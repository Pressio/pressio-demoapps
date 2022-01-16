
set(testname eigen_2d_swe_jacobian_fd_with_sample_mesh_firstorder)
set(exename  ${testname}_exe)

configure_file(../sample_mesh_gids.dat sample_mesh_gids.dat COPYONLY)
add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)

add_test(NAME ${testname}
COMMAND ${CMAKE_COMMAND}
-DFMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DSMESHDRIVER=${MESHSRC}/create_sample_mesh.py
-DEXENAME=$<TARGET_FILE:${exename}>
-DSTENCILVAL=3
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)

set_tests_properties(${testname}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

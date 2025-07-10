
set(testname eigen_2d_euler_riemann_ghost_fill_check)
set(exename  ${testname}_exe)

set(A "sample_mesh_s5;sample_mesh_s7")
foreach(id IN LISTS A)
  set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/${id})
  set(TMPDIR ${CMAKE_CURRENT_BINARY_DIR}/${id})
  file(MAKE_DIRECTORY ${TMPDIR})
  configure_file(${SRCDIR}/info.dat         ${TMPDIR}/info.dat COPYONLY)
  configure_file(${SRCDIR}/connectivity.dat ${TMPDIR}/connectivity.dat COPYONLY)
  configure_file(${SRCDIR}/coordinates.dat  ${TMPDIR}/coordinates.dat COPYONLY)
endforeach()

add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
target_compile_definitions(${exename} PUBLIC -DPRESSIODEMOAPPS_ENABLE_TESTS)

add_test(NAME ${testname} COMMAND ${exename})
set_tests_properties(${testname}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

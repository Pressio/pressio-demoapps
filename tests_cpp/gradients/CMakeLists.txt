
set(testname gradients_against_sine_function)
set(exename ${testname}_exe)

set(A "fullmesh_s3;fullmesh_s5;fullmesh_s7;samplemesh_s3;samplemesh_s5;samplemesh_s7")
foreach(id IN LISTS A)
  set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/${id})
  set(TMPDIR ${CMAKE_CURRENT_BINARY_DIR}/${id})
  file(MAKE_DIRECTORY ${TMPDIR})
  configure_file(${SRCDIR}/info.dat         ${TMPDIR}/info.dat COPYONLY)
  configure_file(${SRCDIR}/connectivity.dat ${TMPDIR}/connectivity.dat COPYONLY)
  configure_file(${SRCDIR}/coordinates.dat  ${TMPDIR}/coordinates.dat COPYONLY)
endforeach()

configure_file(plot.py plot.py COPYONLY)

add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
add_test(NAME ${testname} COMMAND ${exename})
set_tests_properties(${testname}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

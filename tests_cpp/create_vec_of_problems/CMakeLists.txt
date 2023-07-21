
set(testname create_vec_of_problems)
set(exename ${testname}_exe)

foreach(id RANGE 2)
  set(TMPDIR ${CMAKE_CURRENT_BINARY_DIR}/domain${id})
  message(STATUS ${TMPDIR})
  file(MAKE_DIRECTORY ${TMPDIR})
  configure_file(info.dat         ${TMPDIR}/info.dat COPYONLY)
  configure_file(connectivity.dat ${TMPDIR}/connectivity.dat COPYONLY)
  configure_file(coordinates.dat  ${TMPDIR}/coordinates.dat COPYONLY)
endforeach()

add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
add_test(NAME ${testname} COMMAND ${exename})
set_tests_properties(${testname}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

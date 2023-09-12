
set(testname gradients)
set(exename ${testname}_exe)

configure_file(info.dat         info.dat COPYONLY)
configure_file(connectivity.dat connectivity.dat COPYONLY)
configure_file(coordinates.dat  coordinates.dat COPYONLY)

add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
add_test(NAME ${testname} COMMAND ${exename})
set_tests_properties(${testname}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

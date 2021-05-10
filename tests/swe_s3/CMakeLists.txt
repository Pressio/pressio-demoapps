
configure_file(connectivity.dat connectivity.dat COPYONLY)
configure_file(coordinates.dat  coordinates.dat  COPYONLY)
configure_file(info.dat	     info.dat	   COPYONLY)

set(test_name swe_s3)
add_executable(${test_name} main.cc)
add_test(NAME ${test_name} COMMAND ${test_name})
set_tests_properties(${test_name}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

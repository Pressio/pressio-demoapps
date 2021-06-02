
add_executable(weno5exe main.cc)
add_test(NAME weno5 COMMAND weno5exe)
set_tests_properties(weno5
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )


add_executable(weno3exe main.cc)
add_test(NAME weno3 COMMAND weno3exe)
set_tests_properties(weno3
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

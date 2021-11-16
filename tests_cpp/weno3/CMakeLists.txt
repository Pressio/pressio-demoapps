
add_executable(weno3exe main.cc)
add_test(NAME weno3 COMMAND weno3exe)
set_tests_properties(weno3
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )
add_executable(weno3jacobiansexe main_jacobians.cc)
add_test(NAME weno3jacobians COMMAND weno3jacobiansexe)
set_tests_properties(weno3jacobians
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

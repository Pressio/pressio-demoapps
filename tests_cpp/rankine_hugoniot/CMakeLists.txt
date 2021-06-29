
add_executable(rankineHugoniotexe main.cc)
add_test(NAME rankineHugoniot COMMAND rankineHugoniotexe)
set_tests_properties(rankineHugoniot
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

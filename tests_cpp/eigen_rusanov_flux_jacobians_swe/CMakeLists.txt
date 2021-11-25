
set(testname eigen_swe_rusanov_flux_jacobian2d)
set(exename  ${testname}_exe)
add_executable(${exename} main2d.cc)
add_test(NAME ${testname} COMMAND ${exename})
set_tests_properties(${testname}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )




set(testname eigen_euler_rusanov_flux_jacobian1d)
set(exename  ${testname}_exe)
add_executable(${exename} main1d.cc)
add_test(NAME ${testname} COMMAND ${exename})
set_tests_properties(${testname}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

set(testname eigen_euler_rusanov_flux_jacobian2d)
set(exename  ${testname}_exe)
add_executable(${exename} main1d.cc)
add_test(NAME ${testname} COMMAND ${exename})
set_tests_properties(${testname}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

set(testname eigen_euler_rusanov_flux_jacobian3d)
set(exename  ${testname}_exe)
add_executable(${exename} main1d.cc)
add_test(NAME ${testname} COMMAND ${exename})
set_tests_properties(${testname}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )
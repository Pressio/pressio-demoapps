
add_executable(euler_jacobians1d_exe testJacobian1d.cc)
add_test(NAME test_euler_jacobians1d COMMAND euler_jacobians1d_exe)
set_tests_properties(test_euler_jacobians1d
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

add_executable(euler_jacobians2d_exe testJacobian2d.cc)
add_test(NAME test_euler_jacobians2d COMMAND euler_jacobians2d_exe)
set_tests_properties(test_euler_jacobians2d
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

add_executable(euler_jacobians3d_exe testJacobian3d.cc)
add_test(NAME test_euler_jacobians3d COMMAND euler_jacobians3d_exe)
set_tests_properties(test_euler_jacobians3d
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )


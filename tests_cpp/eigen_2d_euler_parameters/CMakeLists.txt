
add_executable(eigen_2d_euler_parameters main.cc)
add_test(NAME eigen_2d_euler_parameters COMMAND eigen_2d_euler_parameters)
set_tests_properties(eigen_2d_euler_parameters
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
)

configure_file(info.dat info.dat COPYONLY)
configure_file(coordinates.dat coordinates.dat COPYONLY)
configure_file(connectivity.dat connectivity.dat COPYONLY)
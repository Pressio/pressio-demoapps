
#configure_file(info.dat info.dat COPYONLY)
#configure_file(connectivity.dat connectivity.dat COPYONLY)
#configure_file(coordinates.dat coordinates.dat COPYONLY)

add_executable(eigencrsmatrix main.cc)
add_test(NAME eigencrsmatrix COMMAND eigencrsmatrix)
set_tests_properties(eigencrsmatrix
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )


configure_file(info.dat info.dat COPYONLY)
configure_file(connectivity.dat connectivity.dat COPYONLY)
configure_file(coordinates.dat coordinates.dat COPYONLY)

add_executable(meshobjexe main.cc)
add_test(NAME ccuMeshObj COMMAND meshobjexe)
set_tests_properties(ccuMeshObj
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

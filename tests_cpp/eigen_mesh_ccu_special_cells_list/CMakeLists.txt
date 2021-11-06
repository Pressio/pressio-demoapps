
configure_file(info.dat info.dat COPYONLY)
configure_file(connectivity.dat connectivity.dat COPYONLY)
configure_file(coordinates.dat coordinates.dat COPYONLY)

add_executable(meshobj2exe main.cc)
add_test(NAME eigen_cellCenterUniformMeshObjSpecialCellsList COMMAND meshobj2exe)
set_tests_properties(eigen_cellCenterUniformMeshObjSpecialCellsList
  PROPERTIES
  PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

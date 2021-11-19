
function(make_test reconstructionScheme stencilSize)

  string(TOLOWER ${reconstructionScheme} recSchemeLowerCase)
  set(testname eigen_euler_1d_jac_${recSchemeLowerCase}_fd)
  set(exename  ${testname}_exe)

  configure_file(plot.py plot.py COPYONLY)
  add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
  target_compile_definitions(${exename} PUBLIC -DUSE_${reconstructionScheme})

  add_test(NAME ${testname}
    COMMAND ${CMAKE_COMMAND}
    -DMESHDRIVER=${MESHSRC}/create_full_mesh.py
    -DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
    -DEXENAME=$<TARGET_FILE:${exename}>
    -DSTENCILVAL=${stencilSize}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
    )

  set_tests_properties(${testname}
    PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
    FAIL_REGULAR_EXPRESSION "FAILED"
    )

endfunction()

make_test(FIRSTORDER 3)
#make_test(WENO3 5)
#make_test(WENO5 7)

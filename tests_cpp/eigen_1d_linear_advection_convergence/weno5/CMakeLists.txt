
set(testname eigen_1d_linear_adv_convergence_weno5)
set(exename  ${testname}_exe)

configure_file(conv.py conv.py COPYONLY)
add_executable(${exename} main.cc)

add_test(NAME ${testname}
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DTESTNAME=${testname}
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DSSTEN=7
-DEXENAME=$<TARGET_FILE:${exename}>
-P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
)


set(testname eigen_2d_riemann_firstorder_implicit)
set(exename  ${testname}_exe)

configure_file(../compare.py compare.py COPYONLY)
configure_file(../plot.py plot.py COPYONLY)
configure_file(p_gold.txt p_gold.txt COPYONLY)
configure_file(rho_gold.txt rho_gold.txt COPYONLY)

add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)

add_test(NAME ${testname}
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:${exename}>
-DSTENCILVAL=3
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)

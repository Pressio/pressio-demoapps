
set(testname eigen_2d_cross_shock_weno5_explicit)
set(exename  ${testname}_exe)

configure_file(../compare.py compare.py COPYONLY)
configure_file(rho_gold.txt rho_gold.txt COPYONLY)
configure_file(../plot.py plot.py COPYONLY)

add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(${exename} PUBLIC -DUSE_WENO5)

add_test(NAME ${testname}
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:${exename}>
-DSTENCILVAL=7
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)

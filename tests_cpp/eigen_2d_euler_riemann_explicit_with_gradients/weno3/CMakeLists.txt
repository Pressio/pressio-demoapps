
set(testname eigen_2d_riemann_weno3_explicit_with_gradients)
set(exename  ${testname}_exe)

configure_file(../plot.py plot.py COPYONLY)
configure_file(../compare.py compare.py COPYONLY)
configure_file(grad_gold_init.txt grad_gold_init.txt COPYONLY)
configure_file(grad_gold_final.txt grad_gold_final.txt COPYONLY)

add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(${exename} PUBLIC -DUSE_WENO3 PRESSIODEMOAPPS_ENABLE_GRADIENTS)

add_test(NAME ${testname}
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:${exename}>
-DSTENCILVAL=7
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)

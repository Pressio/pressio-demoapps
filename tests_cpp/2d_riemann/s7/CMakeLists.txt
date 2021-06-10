
add_executable(riemannexeS7 ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(riemannexeS7 PUBLIC -DUSE_WENO5)

configure_file(../compare.py compare.py COPYONLY)
configure_file(p_gold.txt p_gold.txt COPYONLY)
configure_file(rho_gold.txt rho_gold.txt COPYONLY)

add_test(NAME riemann2d_s7
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:riemannexeS7>
-DSTENCILVAL=7
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)



configure_file(compare.py compare.py COPYONLY)
configure_file(p_golds3.txt p_golds3.txt COPYONLY)
configure_file(rho_golds3.txt rho_golds3.txt COPYONLY)
configure_file(p_golds7.txt p_golds7.txt COPYONLY)
configure_file(rho_golds7.txt rho_golds7.txt COPYONLY)

add_executable(riemannexe main.cc)

add_test(NAME riemann2d_s3
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:riemannexe>
-DSTENCILVAL=3
-P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
)

add_test(NAME riemann2d_s7
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:riemannexe>
-DSTENCILVAL=7
-P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
)

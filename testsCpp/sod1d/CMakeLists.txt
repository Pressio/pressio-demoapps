
configure_file(compare.py compare.py COPYONLY)
configure_file(golds3.txt golds3.txt COPYONLY)
configure_file(golds7.txt golds7.txt COPYONLY)

add_executable(sod1dexe main.cc)

add_test(NAME sod1d_s3
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:sod1dexe>
-DSTENCILVAL=3
-P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
)

add_test(NAME sod1d_s7
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:sod1dexe>
-DSTENCILVAL=7
-P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
)

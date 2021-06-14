
add_executable(smooth2dexeS7 ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(smooth2dexeS7 PUBLIC -DUSE_WENO5)

configure_file(../compare.py compare.py COPYONLY)

add_test(NAME smooth2d_s7
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:smooth2dexeS7>
-DSTENCILVAL=7
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)

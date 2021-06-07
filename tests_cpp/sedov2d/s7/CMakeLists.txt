
add_executable(sedovexeS7 ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(sedovexeS7 PUBLIC -DUSE_WENO5)

configure_file(../compare.py compare.py COPYONLY)
configure_file(p_gold.txt p_gold.txt COPYONLY)

add_test(NAME sedov_s7
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:sedovexeS7>
-DSTENCILVAL=7
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)

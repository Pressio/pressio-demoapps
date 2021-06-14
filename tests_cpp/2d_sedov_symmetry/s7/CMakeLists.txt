
add_executable(sedovsymexeS7 ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(sedovsymexeS7 PUBLIC -DUSE_WENO5)

configure_file(../compare.py compare.py COPYONLY)
configure_file(final_state_gold.txt final_state_gold.txt COPYONLY)
configure_file(../plot.py plot.py COPYONLY)

add_test(NAME sedovsym_s7
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:sedovsymexeS7>
-DSTENCILVAL=7
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)

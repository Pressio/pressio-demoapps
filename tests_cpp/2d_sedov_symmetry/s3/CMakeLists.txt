
add_executable(sedovsymexeS3 ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)

configure_file(../compare.py compare.py COPYONLY)
configure_file(final_state_gold.txt final_state_gold.txt COPYONLY)
configure_file(../plot.py plot.py COPYONLY)

add_test(NAME sedovsym_s3
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:sedovsymexeS3>
-DSTENCILVAL=3
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)

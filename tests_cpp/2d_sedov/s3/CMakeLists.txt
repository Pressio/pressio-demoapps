
add_executable(sedovexeS3 ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)

configure_file(../compare.py compare.py COPYONLY)
configure_file(p_gold.txt p_gold.txt COPYONLY)

add_test(NAME sedov_s3
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:sedovexeS3>
-DSTENCILVAL=3
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)


add_executable(sod1dexeS7 ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(sod1dexeS7 PUBLIC -DUSE_WENO5)

configure_file(../compare.py compare.py COPYONLY)
configure_file(gold.txt gold.txt COPYONLY)

add_test(NAME sod1d_s7
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:sod1dexeS7>
-DSTENCILVAL=7
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)

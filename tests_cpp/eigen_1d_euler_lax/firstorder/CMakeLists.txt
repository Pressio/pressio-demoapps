
set(testname eigen_1d_lax_firstorder)
set(exename  ${testname}_exe)

configure_file(../compare.py compare.py COPYONLY)
configure_file(gold.txt gold.txt COPYONLY)
add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)

add_test(NAME ${testname}
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:${exename}>
-DSTENCILVAL=3
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)

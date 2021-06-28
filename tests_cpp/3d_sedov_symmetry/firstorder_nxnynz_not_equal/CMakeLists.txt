
set(testname 3d_sedov_symemtry_firstorder_nxnynz_not_equal)
set(exename  ${testname}_exe)

add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)

configure_file(../compare.py compare.py COPYONLY)
configure_file(../plot.py plot.py COPYONLY)
configure_file(gold_state.txt gold_state.txt COPYONLY)

add_test(NAME ${testname}
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:${exename}>
-DSTENCILVAL=3
-DNX=20
-DNY=25
-DNZ=28
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)


set(testname eigen_3d_sedov_symemtry_weno3_nxnynz_equal)
set(exename  ${testname}_exe)

add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(${exename} PUBLIC -DUSE_WENO3)

configure_file(../compare.py compare.py COPYONLY)
configure_file(../plot.py plot.py COPYONLY)
configure_file(gold_state.txt gold_state.txt COPYONLY)

add_test(NAME ${testname}
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:${exename}>
-DSTENCILVAL=5
-DNX=20
-DNY=20
-DNZ=20
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)

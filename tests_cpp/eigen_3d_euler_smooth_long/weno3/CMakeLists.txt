
set(testname eigen_3d_smooth_weno3)
set(exename  ${testname}_exe)

configure_file(../plot.py plot.py COPYONLY)
configure_file(../compare.py compare.py COPYONLY)

add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)
target_compile_definitions(${exename} PUBLIC -DUSE_WENO3)

add_test(NAME ${testname}
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:${exename}>
-DSTENCILVAL=5
-P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
)

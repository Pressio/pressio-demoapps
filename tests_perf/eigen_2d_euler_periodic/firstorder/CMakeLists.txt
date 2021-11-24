
set(testname perf_test_eigen_2d_euler_periodic_firstorder)
set(exename  ${testname}_exe)
add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/../main.cc)

# add_test(NAME ${testname}
# COMMAND ${CMAKE_COMMAND}
# -DMESHDRIVER=${MESHSRC}/create_full_mesh.py
# -DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
# -DEXENAME=$<TARGET_FILE:${exename}>
# -DSTENCILVAL=3
# -P ${CMAKE_CURRENT_SOURCE_DIR}/../test.cmake
# )

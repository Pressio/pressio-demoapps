
set(testname eigen_2d_swe_slip_wall_parameters)
set(exename  ${testname}_exe)
add_executable(${exename} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)

add_test(NAME ${testname}
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:${exename}>
-DSTENCILVAL=3
-P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
)

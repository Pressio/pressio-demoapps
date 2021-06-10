
configure_file(conv.py conv.py COPYONLY)

add_executable(e1dexe main.cc)

add_test(NAME 1d_euler_convergence
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DTESTNAME=1d_euler_convergence
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:e1dexe>
-P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
)


configure_file(conv.py conv.py COPYONLY)

add_executable(e1dexe main.cc)

add_test(NAME euler1d_periodic
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DTESTNAME=euler1d_periodic
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:e1dexe>
-P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
)

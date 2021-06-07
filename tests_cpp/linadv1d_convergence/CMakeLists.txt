
configure_file(conv.py conv.py COPYONLY)

add_executable(ladvexe main.cc)

add_test(NAME linadv1d_periodic
COMMAND ${CMAKE_COMMAND}
-DMESHDRIVER=${MESHSRC}/create_full_mesh.py
-DTESTNAME=linadv1d_periodic
-DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}
-DEXENAME=$<TARGET_FILE:ladvexe>
-P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
)

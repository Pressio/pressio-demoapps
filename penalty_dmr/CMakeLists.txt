
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include/pressiodemoapps)
set(MESHSRC ${CMAKE_CURRENT_SOURCE_DIR}/../meshing_scripts)

configure_file(plot.py plot.py COPYONLY)
configure_file(svd.py svd.py COPYONLY)

set(fomExeName dmrFomExe)
add_executable(${fomExeName} ${CMAKE_CURRENT_SOURCE_DIR}/main_fom.cc)

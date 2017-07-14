#!/bin/bash

# editted by AT to accomodate MPI

find_package(MPI REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/CToolbox/Modules
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/fi_lib
        ${CMAKE_SOURCE_DIR}/src/rts
		${MPI_INCLUDE_PATH})

link_directories(
  ${CMAKE_SOURCE_DIR/src}
)

add_executable(allzeros allzeros.cpp)
# add_executable(cxsc_mpi_example cxsc_mpi_example.cpp)
add_executable(example example.cpp)
add_executable(example2 example.cpp)
add_executable(inewton inewton.cpp)
add_executable(io io.cpp)
add_executable(lexample lexample.cpp)
add_executable(linewton linewton.cpp)
add_executable(rungekutta rungekutta.cpp)
add_executable(trace trace.cpp)

# add_executable(cxsc_mpicomm cxsc_mpicomm.cpp)
	
target_link_libraries(allzeros cxsc)
# target_link_libraries(cxsc_mpi_example cxsc cxsc_mpicomm)
target_link_libraries(example cxsc)
target_link_libraries(example2 cxsc)
target_link_libraries(inewton cxsc)
target_link_libraries(io cxsc)
target_link_libraries(lexample cxsc)
target_link_libraries(linewton cxsc)
target_link_libraries(rungekutta cxsc)
target_link_libraries(trace cxsc)

# begin add
# if(MPI_COMPILE_FLAGS)
#   set_target_properties(cxsc_mpi_example PROPERTIES
#     COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
# endif()

# if(MPI_LINK_FLAGS)
#   set_target_properties(cxsc_mpi_example PROPERTIES
#     LINK_FLAGS "${MPI_LINK_FLAGS} cxsc_mpicomm.cpp")
# endif()
# end add

file(GLOB example_src "*.cpp")
install(FILES ${example_src} DESTINATION examples)
install(TARGETS allzeros example example2 inewton io lexample linewton rungekutta trace RUNTIME DESTINATION examples)


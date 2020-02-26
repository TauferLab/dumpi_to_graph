cmake_minimum_required(VERSION 3.14)
project(dumpi_to_graph)

list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 17)

add_executable(dumpi_to_graph)

target_sources(dumpi_to_graph PRIVATE src/blocking_collective_callbacks.cpp
src/Channel.cpp src/CommunicatorManager.cpp src/CSMPI_Callstack.cpp
src/Debug.cpp src/EventGraph.cpp src/Glob.cpp
src/matching_function_callbacks.cpp src/Request.cpp src/Trace.cpp
src/Utilities.cpp src/blocking_p2p_callbacks.cpp
src/communicator_management_callbacks.cpp src/Configuration.cpp
src/CSMPI_Trace.cpp src/dumpi_to_graph.cpp src/finalize_callbacks.cpp
src/init_callbacks.cpp src/nonblocking_p2p_callbacks.cpp
src/request_mutating_callbacks.cpp src/undumpi_callbacks.cpp) 

target_include_directories(dumpi_to_graph PRIVATE include include/dumpi_to_graph include/external)

# Boost
find_package(Boost 1.69.0 REQUIRED mpi serialization)
target_include_directories(dumpi_to_graph PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(dumpi_to_graph PRIVATE Boost::mpi)
target_link_libraries(dumpi_to_graph PRIVATE Boost::serialization)
target_link_libraries(dumpi_to_graph PRIVATE Boost::system)

# Variant of sst-dumpi with libdumpi enabled, MPI-3 support, and PAPI support
# not available through spack
target_include_directories(dumpi_to_graph PRIVATE "${DUMPI_ROOT}/include")
target_link_libraries(dumpi_to_graph PRIVATE "${DUMPI_ROOT}/lib/libundumpi.so")

# Igraph
find_package(Igraph)
target_include_directories(dumpi_to_graph PRIVATE ${IGRAPH_INCLUDE_DIR})
target_link_libraries(dumpi_to_graph PRIVATE ${IGRAPH_LIBRARY})

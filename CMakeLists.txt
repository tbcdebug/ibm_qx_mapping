cmake_minimum_required (VERSION 3.0)
project( ibm_qx_mapping )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})

SET(CMAKE_CXX_FLAGS "-Wall -g -O2 -DNDEBUG -DVERSION=\"$(VERSION)\" -std=c++11")

include_directories(src)

FIND_PACKAGE( Boost 1.50 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )


add_executable(ibm_qx_mapping src/main.cpp src/util.cpp src/coupling_graph.cpp src/cost.cpp src/layer_handling.cpp src/circuit_property_handling.cpp src/node_handling.cpp src/mapping.cpp src/qasm/QASMparser.cpp src/qasm/QASMscanner.cpp src/qasm/QASMtoken.cpp)

include_directories(src src/qasm)
target_link_libraries(ibm_qx_mapping boost_program_options)

cmake_minimum_required(VERSION 2.8)
project(FASTQuick)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3  -std=c++0x -Wno-unused-variable -Wno-unused-result -Wno-format-security")

add_subdirectory(libbwa)
add_subdirectory(libmpu)
add_subdirectory(misc)
add_subdirectory(src)
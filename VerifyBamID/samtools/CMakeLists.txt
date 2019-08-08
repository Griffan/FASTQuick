cmake_minimum_required(VERSION 2.8)
project(samtools)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3  -Wno-unused-variable -Wno-unused-result -Wno-format-security")

file(GLOB SOURCE_FILES *.c)

include_directories(${HTS_INCLUDE_DIRS})

add_library(samtools ${SOURCE_FILES})

target_link_libraries(samtools ${HTS_LIBRARIES})

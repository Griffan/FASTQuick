cmake_minimum_required(VERSION 2.8)
project(samtools)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3  -Wno-unused-variable -Wno-unused-result -Wno-format-security -Wno-unused-but-set-variable")

file(GLOB SOURCE_FILES *.c)

add_library(samtools ${SOURCE_FILES})
include_directories(samtools ${HTS_INCLUDE_DIRS})
target_link_libraries(samtools ${HTS_LIBRARIES})

cmake_minimum_required(VERSION 2.8)
file(GLOB CXX_SOURCE_FILES *.cpp)
add_definitions(-DHAVE_PTHREAD)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fpic -pthread")
include_directories(../misc/vcf/ ../misc/general/ ../misc/)

add_executable(FASTQuick ${CXX_SOURCE_FILES})

target_link_libraries(FASTQuick bwa mpu bam tabix vcf_statgen general_statgen misc)

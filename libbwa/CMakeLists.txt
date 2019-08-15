cmake_minimum_required(VERSION 2.8.4)

include_directories(${STATGEN_INCLUDE_DIR})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3  -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-result -Wno-format-security -Wno-unused-function -Wno-write-strings -Wno-sign-compare -Wno-unused-but-set-variable")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -x c++ -pthread -DHAVE_PTHREAD  -mmmx -msse -msse2 -msse3")

set(libbwa_a_SOURCES bamlite.c bntseq.c bwape.c bwase.c bwaseqio.c bwtaln.c bwt.c bwtgap.c bwtio.c bwtmisc.c cs2nt.c is.c kstring.c stdaln.c utils.c bwt_gen.c QSufSort.c)
add_library(bwa  ${libbwa_a_SOURCES})
target_link_libraries(bwa z)

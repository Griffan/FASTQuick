cmake_minimum_required(VERSION 2.8.4)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -x c++ -O3 -pthread  -Wall  -DHAVE_PTHREAD -Wno-write-strings -Wno-unused-result -Wno-sign-compare -Wno-unused-variable -mmmx -msse -msse2 -msse3")
set(libbwa_a_SOURCES bamlite.c bntseq.c bwape.c bwase.c bwaseqio.c bwtaln.c bwt.c bwtgap.c bwtio.c bwtmisc.c cs2nt.c is.c kstring.c stdaln.c utils.c bwt_gen.c QSufSort.c)
add_library(bwa  ${libbwa_a_SOURCES})
target_link_libraries(bwa z)

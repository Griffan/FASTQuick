cmake_minimum_required(VERSION 2.8.4)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fpic")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -fpic -Wno-constant-conversion -Wno-unused-function")


include_directories(${STATGEN_INCLUDE_DIR} ../misc/ ./)

set(libtabix_a_SOURCES bedidx.c mbgzf.c  kstring.c index.c)

set(libmpu_a_SOURCES boolParser.cpp  fVcf.cpp  MathGenMin.cpp  MathGold.cpp  mpuMixupLikelihoods.cpp  mpuPileBases.cpp
        mpuTool.cpp  mpuVerify.cpp  Parameters.cpp  PhredHelper.cpp  Random.cpp)

add_library(tabix ${libtabix_a_SOURCES})
add_library(mpu ${libmpu_a_SOURCES} )



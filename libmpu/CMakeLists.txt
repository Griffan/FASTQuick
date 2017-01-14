cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fpic")

include_directories(../misc/general/ ../misc/ ./)
set(libtabix_a_SOURCES bedidx.c mbgzf.c  kstring.c index.c)

set(libmpu_a_SOURCES boolParser.cpp  fVcf.cpp  MathGenMin.cpp  MathGold.cpp  mpuMixupLikelihoods.cpp  mpuPileBases.cpp
        mpuTool.cpp  mpuVerify.cpp  Parameters.cpp  PhredHelper.cpp  Random.cpp)

add_library(tabix ${libtabix_a_SOURCES})
add_library(mpu ${libmpu_a_SOURCES} )



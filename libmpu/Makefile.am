## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

## Paths to include
AM_CPPFLAGS = -Wall -fpic -O2
#AM_CPPFLAGS = -Wall -fpic -g -O0

#noinst_HEADERS = bam_endian.h boolParser.h kseq.h pException.h tabix.h ksort.h pFile.h bgix.h kstring.h PhredHelper.h bgzf.h fVcf.h khash.h Parameters.h boolExprParser.h knetfile.h wFile.h mpuFile.h Error.h Constant.h mpuLikelihoods.h BaseAsciiMap.h fPed.h MathGold.h MathGenMin.h MathConstant.h MathVector.h MathMatrix.h Random.h Sort.h mpuVerify.h mpuPileBases.h mpuMixupLikelihoods.h mpuSetLikelihoods.h mpuSet.h
noinst_HEADERS = bam_endian.h boolParser.h kseq.h pException.h tabix.h ksort.h pFile.h bgix.h kstring.h PhredHelper.h mbgzf.h fVcf.h khash.h Parameters.h boolExprParser.h knetfile.h wFile.h mpuFile.h Error.h Constant.h BaseAsciiMap.h fPed.h MathGold.h MathGenMin.h MathConstant.h MathVector.h MathMatrix.h Random.h Sort.h mpuVerify.h mpuPileBases.h mpuMixupLikelihoods.h

bin_PROGRAMS =  tabix

noinst_LIBRARIES = libtabix.a libmpu.a

## pEmmax files
#mpuTool_SOURCES = mpuTool.cpp Parameters.cpp Error.cpp PhredHelper.cpp mpuLikelihoods.cpp BaseAsciiMap.cpp MathGold.cpp MathGenMin.cpp MathVector.cpp MathMatrix.cpp Random.cpp Sort.cpp fVcf.cpp boolParser.cpp mpuVerify.cpp mpuPileBases.cpp mpuMixupLikelihoods.cpp mpuSetLikelihoods.cpp
#mpuTool_SOURCES = mpuTool.cpp Parameters.cpp Error.cpp PhredHelper.cpp BaseAsciiMap.cpp MathGold.cpp MathGenMin.cpp MathVector.cpp MathMatrix.cpp Random.cpp Sort.cpp fVcf.cpp boolParser.cpp mpuVerify.cpp mpuPileBases.cpp mpuMixupLikelihoods.cpp
libmpu_a_LIBADD = libtabix.a

#mpuTool_LDADD = libtabix.a

#bgzip_SOURCES = bgzip.c bedidx.c bgzf.c knetfile.c kstring.c

tabix_SOURCES = tabix.c bedidx.c mbgzf.c knetfile.c kstring.c index.c

libtabix_a_SOURCES = bedidx.c mbgzf.c knetfile.c kstring.c index.c

libmpu_a_SOURCES = mpuTool.cpp Parameters.cpp Error.cpp PhredHelper.cpp BaseAsciiMap.cpp MathGold.cpp MathGenMin.cpp MathVector.cpp MathMatrix.cpp Random.cpp Sort.cpp fVcf.cpp boolParser.cpp mpuVerify.cpp mpuPileBases.cpp mpuMixupLikelihoods.cpp

## Define an independent executable script for inclusion in the distribution
## archive. It will not be installed on an end user's system, however.
#dist_noinst_SCRIPTS = autogen.sh

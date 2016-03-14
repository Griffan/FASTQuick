## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = -I  m4 ${ACLOCAL_FLAGS}

## Paths to include
AM_CXXFLAGS = -O3 -Wall -fpic  -I../misc/general/  -I../misc/

noinst_HEADERS = bam_endian.h  boolExprParser.h  Constant.h  fVcf.h   kseq.h   kstring.h     MathGold.h  mpuFile.h              mpuPileBases.h  mpuVerify.h   pFile.h        Random.h  wFile.h bgix.h        boolParser.h      fPed.h      khash.h  ksort.h  MathGenMin.h  mbgzf.h     mpuMixupLikelihoods.h  mpuTool.h       Parameters.h  PhredHelper.h  tabix.h

bin_PROGRAMS =  tabix

noinst_LIBRARIES = libtabix.a libmpu.a

## pEmmax files
#mpuTool_SOURCES = mpuTool.cpp Parameters.cpp Error.cpp PhredHelper.cpp mpuLikelihoods.cpp BaseAsciiMap.cpp MathGold.cpp MathGenMin.cpp MathVector.cpp MathMatrix.cpp Random.cpp Sort.cpp fVcf.cpp boolParser.cpp mpuVerify.cpp mpuPileBases.cpp mpuMixupLikelihoods.cpp mpuSetLikelihoods.cpp
#mpuTool_SOURCES = mpuTool.cpp Parameters.cpp Error.cpp PhredHelper.cpp BaseAsciiMap.cpp MathGold.cpp MathGenMin.cpp MathVector.cpp MathMatrix.cpp Random.cpp Sort.cpp fVcf.cpp boolParser.cpp mpuVerify.cpp mpuPileBases.cpp mpuMixupLikelihoods.cpp
libmpu_a_LIBADD = libtabix.a

#mpuTool_LDADD = libtabix.a

#bgzip_SOURCES = bgzip.c bedidx.c bgzf.c knetfile.c kstring.c

tabix_SOURCES = tabix.c bedidx.c mbgzf.c kstring.c index.c

libtabix_a_SOURCES = bedidx.c mbgzf.c  kstring.c index.c

#libmpu_a_SOURCES = mpuTool.cpp Parameters.cpp  PhredHelper.cpp  MathGold.cpp MathGenMin.cpp  MathMatrix.cpp Random.cpp Sort.cpp fVcf.cpp boolParser.cpp mpuVerify.cpp mpuPileBases.cpp mpuMixupLikelihoods.cpp
libmpu_a_SOURCES = boolParser.cpp  fVcf.cpp  MathGenMin.cpp  MathGold.cpp  mpuMixupLikelihoods.cpp  mpuPileBases.cpp  mpuTool.cpp  mpuVerify.cpp  Parameters.cpp  PhredHelper.cpp  Random.cpp
## Define an independent executable script for inclusion in the distribution
## archive. It will not be installed on an end user's system, however.
#dist_noinst_SCRIPTS = autogen.sh

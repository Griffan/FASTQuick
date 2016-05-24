## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

## Paths to include
#TABIX_VER = 0.2.5
#TABIX = tabix-$(TABIX_VER)
AM_CPPFLAGS = -Wall -fpic -O3 -I../../libmpu/ 

#noinst_LIBRARIES = libtabix.a
#libtabix_a_SOURCES = bedidx.c bgzf.c knetfile.c kstring.c index.c
#noinst_HEADERS = hDist.h \
#sparseHMM.h icGraph.h pcHMM.h \
#AnnotationOutput.h \
#		bam_endian.h boolParser.h Constant.h GeneFormat.h GenomeSequence.h \
#		kseq.h ModelParser.h pException.h Range.h tabix.h \
#		AnnotationResultCollection.h BedReader.h cdsStat.h Error.h Gene.h \
#		GitVersion.h ksort.h ncStat.h pFile.h TabixReader.h AnnotationResult.h \
#		bgix.h Codon.h FreqTable.h genomeLoci.h genomePosition.h IO.h kstring.h OrderedMap.h \
#		PhredHelper.h SequenceUtil.h TypeConversion.h AnnotationString.h bgzf.h \
#		codonHelper.h fVcf.h genomeScore.h khash.h locusStat.h Parameters.h Priority.h \
#		StringTemplate.h Type.h Argument.h boolExprParser.h Common.h GeneAnnotation.h \
#		GenomeScore.h LineBreaker.h knetfile.h LogFile.h pEmmax.h pVcf.h StringUtil.h \
#		wFile.h

bin_PROGRAMS = vcfast

## vcfast files
vcfast_SOURCES = vcfast.cpp Error.cpp PhredHelper.cpp Parameters.cpp codonHelper.cpp genomePosition.cpp
vcfast_LDADD = ../../libmpu/libtabix.a ../../libmpu/libmpu.a


## Define an independent executable script for inclusion in the distribution
## archive. It will not be installed on an end user's system, however.
#dist_noinst_SCRIPTS = autogen.sh

AUTOMAKE_OPTIONS=foreign     

noinst_LIBRARIES=libbam.a

libbam_a_SOURCES= BamIndex.cpp      GenericSamInterface.cpp    PosList.cpp         SamFileHeader.cpp  SamHeaderPG.cpp      SamHeaderSQ.cpp   SamQuerySeqWithRefHelper.cpp  SamRecordPool.cpp     SamTags.cpp BamInterface.cpp  PileupElementBaseQual.cpp  SamCoordOutput.cpp  SamFilter.cpp      SamHeaderRecord.cpp  SamHeaderTag.cpp  SamRecord.cpp                 SamReferenceInfo.cpp  SamValidation.cpp CigarHelper.cpp   PileupElement.cpp          SamFile.cpp         SamHeaderHD.cpp    SamHeaderRG.cpp      SamInterface.cpp  SamRecordHelper.cpp           SamStatistics.cpp

libbam_a_CXXFLAGS= -O3 -I../ -lz -I../general/  -Wall -std=c++0x  -mmmx -msse -msse2 -msse3 -pipe -Werror  -Wno-strict-overflow -I../include -I.   -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS

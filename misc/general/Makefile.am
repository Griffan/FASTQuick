AUTOMAKE_OPTIONS=foreign     

noinst_LIBRARIES=libgeneral.a



libgeneral_a_SOURCES= BaseAsciiMap.cpp   BgzfFileTypeRecovery.cpp  CigarRoller.cpp   FileType.cpp        GzipFileType.cpp  IndexBase.cpp  MathMatrix.cpp      MemoryMap.cpp  StatGenStatus.cpp  StringBasics.cpp  UncompressedFileType.cpp BaseUtilities.cpp  CharBuffer.cpp            Error.cpp         Generic.cpp         GzipHeader.cpp    InputFile.cpp  MathVector.cpp      PhoneHome.cpp  STLUtilities.cpp   StringHash.cpp BgzfFileType.cpp   Cigar.cpp                 ErrorHandler.cpp  GenomeSequence.cpp  Hash.cpp          IntArray.cpp   MemoryMapArray.cpp  Sort.cpp       StringArray.cpp    Tabix.cpp

libgeneral_a_CXXFLAGS= -O3 -I../ -L../ -lmisc -lz  -Wall -std=c++0x   -mmmx -msse -msse2 -msse3 -pipe -Werror  -Wno-strict-overflow  -I.   -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS


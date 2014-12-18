#ifndef __MPU_PILE_BASES__H
#define __MPU_PILE_BASES__H

#include <map>
#include <vector>
#include <string>
#include "mpuFile.h"
#include "Error.h"

class mpuPileBases {
 public:
  // constructor
  mpuPileBases(const char* mpuFile, const char* smID = NULL, int maxDP = 255);

  // read a markers
  int readMarker(const char* chrom, int position, char refB, char altB);

  // parameters (public)
  int minMapQ;
  int maxDepth;
  int minQ;
  int maxQ;

  std::vector<uint8_t> nBases;         // bases
  std::vector<uint8_t> nQuals;         // quals
  std::vector<uint8_t> nMapQs;
  std::vector<uint32_t> nBegins;    // m-th marker is begin <= i < end
  std::vector<uint32_t> nEnds;      // m-th marker is begin <= i < end

 //protected:
  mpuFile inMpu;
  std::string sSM;
};
#endif

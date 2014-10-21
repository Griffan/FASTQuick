#include "mpuPileBases.h"

// construct BamPileBases object from bamFile
mpuPileBases::mpuPileBases(const char* mpuFile, const char* smID, int maxDP) : inMpu(maxDP) 
{
  // open BAM File
  if ( ! inMpu.load( mpuFile )  ) {
    error("Cannot open MPU file %s for reading", mpuFile);
  }
  sSM = smID;
}

int mpuPileBases::readMarker(const char* chrom, int position, char refB, char altB) {
  inMpu.advanceTo(chrom,position, refB, altB);
  nBegins.push_back(nBases.size());
  int totalDepth = 0;
  for(int i=0; i < inMpu.nbase; ++i) {
    if ( ( inMpu.bQs[i] >= minQ ) && ( inMpu.mQs[i] >= minMapQ ) ) {
      if ( inMpu.bQs[i] == maxQ ) inMpu.bQs[i] = maxQ;
      nBases.push_back(inMpu.bases[i]);
      nQuals.push_back(inMpu.bQs[i]);
      nMapQs.push_back(inMpu.mQs[i]);
      ++totalDepth;
    }
  }
  nEnds.push_back(nBases.size());
  return nEnds.back()-nBegins.back();
}

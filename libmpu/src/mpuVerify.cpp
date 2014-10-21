#include "mpuVerify.h"
#include "fVcf.h"

#include <set>

GenMatrixBinary::GenMatrixBinary(const char* vcfFile, bool siteOnly, bool findBest, std::vector<std::string>& subsetInds, double minAF, double minCallRate) {
  // open a VCF file
  fVcf tvcf;
  tvcf.infoAF = ( siteOnly || ( ( subsetInds.size() <= 1 ) && !findBest ) );

  int i = 0, m = 0, M = 0;
  double af, maf;
  int unit = 10000;

  std::set<std::string> idset;
  for(i=0; i < (int)subsetInds.size(); ++i) {
    if ( idset.find(subsetInds[i]) != idset.end() ) {
      error("ERROR: Duplicate individual ID %s",subsetInds[i].c_str());
    }
    idset.insert(subsetInds[i]);
  }
  tvcf.ignoreEmpty = findBest;
  tvcf.load(vcfFile, NULL, "GT", NULL, true, idset);

  indids = tvcf.inds;

  // set bytesPerMarker attribute
  //error("indids.size() = %d",tvcf.nInds);

  if ( siteOnly ) {
    bytesPerMarker = 0;
  }
  else {
    bytesPerMarker = (tvcf.nInds + 3)/4;
  }

  for(M=0; tvcf.readMarkers(unit); ) {
    M += tvcf.nMarkers;
    fprintf(stderr,"Processing %d markers across %d individuals..\n", M, tvcf.nInds);
    for(i=0, m=0; i < tvcf.nMarkers; ++i) { // for each marker
      af = tvcf.alleleFreq(i);
      maf = af > 0.5 ? 1-af : af;
      if ( ( maf >= minAF ) && 
	   ( tvcf.callRate(i) >= minCallRate ) ) {
	addMarker(tvcf.chroms[i].c_str(), tvcf.pos1s[i], tvcf.refs[i][0], tvcf.alts[i][0], af);
	for(int j=0; j < tvcf.nInds; ++j) {
	  setGenotype( tvcf.genos[i*tvcf.nInds + j], j, m );
	}
	++m;
      }
    }
  }
}

int GenMatrixBinary::addMarker(const char* chrom, int position, char refBase, char altBase, double alleleFreq) {
  chroms.push_back(chrom);
  positions.push_back(position);
  refBases.push_back(refBase);
  altBases.push_back(altBase);
  alleleFrequencies.push_back(alleleFreq);

  //notice("%s:%d\t%d",chrom,position,(int)genotypes.size());
  for(int i=0; i < bytesPerMarker; ++i) {
    genotypes.push_back(0);
  }
  //notice("%s:%d\t%d",chrom,position,(int)genotypes.size());

  return (int)chroms.size();
}

int GenMatrixBinary::getGenotype(int indIndex, int markerIndex) {
  int genoIndex = markerIndex * bytesPerMarker + indIndex / 4;
  int offset = (indIndex % 4) * 2;
  return static_cast<int>(( genotypes[genoIndex] >> offset ) & 0x03);
}

void GenMatrixBinary::setGenotype(float geno, int indIndex, int markerIndex) {
  uint8_t ngeno = isnan(geno) ? 0 : ((uint8_t)floor(geno)+1);
  int genoIndex = (chroms.size()-1) * bytesPerMarker + (indIndex / 4);
  int shift = ((indIndex % 4) * 2);
  genotypes[genoIndex] |= (ngeno << shift);
}

void mpuVerify::loadFiles(const char* mpuFile, const char* vcfFile) {
  // create a pile object
  notice("Opening MPU file %s",mpuFile);

  const char* smID = pArgs->sSMID.empty() ? NULL : pArgs->sSMID.c_str();
  pPile = new mpuPileBases(mpuFile, smID, pArgs->maxDepth);
  pPile->minMapQ = pArgs->minMapQ;
  pPile->maxDepth = pArgs->maxDepth;
  pPile->minQ = pArgs->minQ;
  pPile->maxQ = pArgs->maxQ;

  mixOut.init(pArgs->pRefRef,pArgs->pRefHet,pArgs->pRefAlt,pArgs->maxDepth);
  selfOut.init(pArgs->pRefRef,pArgs->pRefHet,pArgs->pRefAlt,0);
  bestOut.init(pArgs->pRefRef,pArgs->pRefHet,pArgs->pRefAlt,0);

  if ( pArgs->bSiteOnly ) {
    if ( subsetInds.size() > 0 ) {
      error("--siteOnly option cannot be combined with --subset option");
    }
  }
  else if ( pArgs->bSelfOnly )  {
    notice("--selfOnly option applied : finding sample ID %s from VCF file", pPile->sSM.c_str());
    if ( ! (pArgs->sSubsetInds.empty()) ) {
      error("--selfOnly option cannot be combined with --subset option");
    }
    subsetInds.push_back(pPile->sSM);
  }
  else if ( pArgs->bFindBest ) { // if --best option is set, load all individuals
  }
  else {
    error("None of --site, --self, --best option was set");
  }

  notice("Opening VCF file %s",vcfFile);
  // create genotype matrix, and load vcfFile
  pGenotypes = new GenMatrixBinary(vcfFile, pArgs->bSiteOnly, pArgs->bFindBest, subsetInds, pArgs->minAF, pArgs->minCallRate);

  notice("Reading MPU file %s",mpuFile);
  // read base information corresponding to each marker
  nMarkers = (int)pGenotypes->chroms.size();
  
  for(int i=0; i < nMarkers; ++i) {
    //if ( ( i > 0 ) && ( i % 10000 == 0 ) ) {
    if ( ( i % 1000 == 0 ) ) {
      notice("Extracting read information in %d markers from BAM file",i);
    }
    //notice("%s:%d",pGenotypes->chroms[i].c_str(),pGenotypes->positions[i]);
    pPile->readMarker(pGenotypes->chroms[i].c_str(), pGenotypes->positions[i],pGenotypes->refBases[i], pGenotypes->altBases[i]);
    //++markerDPs[nb];
  }
  nBases = pPile->nBases.size();
  notice("Finished extracting %d bases in %d markers from BAM file -- Avg Depth = %.3lf", nBases, nMarkers, (double)nBases/nMarkers);
  notice("Finished Reading MPU file %s and VCF file %s\n",mpuFile,vcfFile);
}

// computeMixLLKs() - calculate Pr(Data|fMix) for each SM and RG
//                    given AF and 
double mpuVerify::computeMixLLKs(double fMix) { //{ double* rgLLKs) {
//int nMarkers = (int)(pGenotypes->chroms.size());
//int nRGs = (int)(pPile->vsRGIDs.size());
  double genoFreq[3], baseLKs[9], smMarkerLKs[9];
  double af, baseError, baseMatch;
  uint8_t a1, a2, base;
  double smLLK = 0;
  double tmpf;

  for(int i=0; i < nMarkers; ++i) {
    if ( fVcf::getChromosomeType(pGenotypes->chroms[i].c_str(), pGenotypes->positions[i]) != CT_AUTOSOME ) continue;
    af = pGenotypes->alleleFrequencies[i];
    // force allele-frequency as nonzero to avoid boundary conditions
    if ( af < 0.001 ) af = 0.001;
    if ( af > 0.999 ) af = 0.999;

    // frequency of genotypes under HWE
    genoFreq[0] = (1.-af)*(1.-af);
    genoFreq[1] = 2.*af*(1.-af);
    genoFreq[2] = af*af;
    a1 = BaseAsciiMap::base2int[(int)pGenotypes->refBases[i]];
    a2 = BaseAsciiMap::base2int[(int)pGenotypes->altBases[i]];

    // initialize the marker-level likelihoods
    for(int k=0; k < 9; ++k) {
      smMarkerLKs[k] = (pArgs->bPrecise ? 0. : 1.);
    }

    for(int j=(int)pPile->nBegins[i]; j < (int)pPile->nEnds[i]; ++j) {
      // obtain b (base), (error), and readgroup info
      base = pPile->nBases[j];
      baseError = fPhred2Err[pPile->nQuals[j]];
      baseMatch = 1.-baseError;
	
      for(int k1=0; k1 < 3; ++k1) {
	for(int k2=0; k2 < 3; ++k2) {
	  baseLKs[k1*3+k2] = ( fMix * pSN[0*3+k1] + (1.-fMix) * pSN[0*3+k2] ) * ( base == a1 ? baseMatch : baseError/3.) + ( fMix * pSN[1*3+k1] + (1.-fMix) * pSN[1*3+k2] ) * ( base == a2 ? baseMatch : baseError/3. );
	}
      }
	
      for(int k1=0; k1 < 3; ++k1) {
	for(int k2=0; k2 < 3; ++k2) {
	  if ( pArgs->bPrecise ) {
	    tmpf = log(baseLKs[k1*3+k2]);
	    smMarkerLKs[k1*3+k2] += tmpf;
	  }
	  else {
	    tmpf = baseLKs[k1*3+k2];
	    smMarkerLKs[k1*3+k2] *= tmpf; 
	  }
	}
      }
    }

    //notice("smMarkerLKs[0] = %lf\t%d\t%d",smMarkerLKs[0], pPile->nBegins[i], pPile->nEnds[i]);

    // sample-level per-marker likelihood
    // smLLK = Pr(Data|fMix) = \prod_{markers} Pr(bases|fMix)
    //  = \prod_{markers} \sum_{G1,G2} Pr(bases|G1,G2,fMix) Pr(G1|AF)Pr(G2|AF)
    double perMarkerProb = 0;
    if ( pArgs->bPrecise ) {
      int maxIdx = 0;
      for(int k=1; k < 9; ++k) {
	if ( smMarkerLKs[maxIdx] < smMarkerLKs[k] ) {
	  maxIdx = k;
	}
      }
      for(int k1=0; k1 < 3; ++k1) {
	for(int k2=0; k2 < 3; ++k2) {
	  perMarkerProb += (exp(smMarkerLKs[k1*3+k2] - smMarkerLKs[maxIdx]) * genoFreq[k1] * genoFreq[k2]);
	}
      }
      smLLK += (log(perMarkerProb) + smMarkerLKs[maxIdx]);
    }
    else {
      for(int k1=0; k1 < 3; ++k1) {
	for(int k2=0; k2 < 3; ++k2) {
	  perMarkerProb += (smMarkerLKs[k1*3+k2] * genoFreq[k1] * genoFreq[k2]);
	}
      }
      smLLK += log(perMarkerProb);
    }
  }
  notice("double mpuVerify::computeMixLLK( %.6lf | %.6lf, %.6lf, %.6lf ) = %.6lf",fMix,pSN[0],pSN[1],pSN[2],smLLK);
  return smLLK;
}

// computeIBDLLKs() - calculate Pr(Data|fIBD) for each SM and RG
//                    given AF and call rate threshold
double mpuVerify::computeIBDLLKs(double fIBD, int indIdx) {
  int nMarkers = (int)(pGenotypes->chroms.size());
  //int nRGs = (int)(pPile->vsRGIDs.size());
  double genoFreq[3], genoProb[3], baseLKs[9], smMarkerLKs[9];
  double af, baseError, baseMatch;
  uint8_t a1, a2, base;
  //int rgIdx;
  double smLLK = 0;
  double tmpf;

  for(int i=0; i < nMarkers; ++i) {
    if ( fVcf::getChromosomeType(pGenotypes->chroms[i].c_str(), pGenotypes->positions[i]) != CT_AUTOSOME ) continue;
    af = pGenotypes->alleleFrequencies[i];
    // force allele-frequency as nonzero to avoid boundary conditions
    if ( af < 0.001 ) af = 0.001;
    if ( af > 0.999 ) af = 0.999;

    // frequency of genotypes under HWE
    genoFreq[0] = (1.-af)*(1.-af);
    genoFreq[1] = 2.*af*(1.-af);
    genoFreq[2] = af*af;

    int geno = pGenotypes->getGenotype(indIdx,i);
    switch(geno) {
    case 0: // MISSING
      genoProb[0] = genoFreq[0];
      genoProb[1] = genoFreq[1];
      genoProb[2] = genoFreq[2];
      break;
    case 1: // HOMREF;
      genoProb[0] = 1.-pArgs->genoError;
      genoProb[1] = pArgs->genoError * (GENO_ERR_TO_HET);
      genoProb[2] = pArgs->genoError * (1.-GENO_ERR_TO_HET);
      break;
    case 2: // HET;
      genoProb[0] = pArgs->genoError/2.;
      genoProb[1] = 1.-pArgs->genoError;
      genoProb[2] = pArgs->genoError/2.;
      break;
    case 3: // HET;
      genoProb[0] = pArgs->genoError/2. * (1.-GENO_ERR_TO_HET);
      genoProb[1] = pArgs->genoError/2. * GENO_ERR_TO_HET;
      genoProb[2] = 1.-pArgs->genoError;
      break;
    default:
      error("Unrecognized genotype %d at ind %d, marker %d",indIdx,i);
    }

    a1 = BaseAsciiMap::base2int[(int)pGenotypes->refBases[i]];
    a2 = BaseAsciiMap::base2int[(int)pGenotypes->altBases[i]];

    // initialize the marker-level likelihoods
    for(int k=0; k < 9; ++k) {
      smMarkerLKs[k] = (pArgs->bPrecise ? 0. : 1.);
    }
    //for(int k=0; k < 9 * nRGs; ++k) {
    //  rgMarkerLKs[k] = (pArgs->bPrecise ? 0. : 1.);
    //}

    for(int j=(int)pPile->nBegins[i]; j < (int)pPile->nEnds[i]; ++j) {
      base = pPile->nBases[j];
      baseError = fPhred2Err[pPile->nQuals[j]];
      baseMatch = 1.-baseError;
	
      for(int k1=0; k1 < 3; ++k1) {
	for(int k2=0; k2 < 3; ++k2) {
	  baseLKs[k1*3+k2] = ( fIBD * pSN[0*3+k1] + (1.-fIBD) * pSN[0*3+k2] ) * ( base == a1 ? baseMatch : baseError/3.) + ( fIBD * pSN[1*3+k1] + (1.-fIBD) * pSN[1*3+k2] ) * ( base == a2 ? baseMatch : baseError/3. );
	}
      }
      
      // merge base-level likelihood into marker-level likelihood
      // smMarkerLKs[i*3+j] 
      //   = \prod Pr(b|G1,G2) = \prod (f Pr(b|G1) + (1-f) Pr(b|G2))
      for(int k1=0; k1 < 3; ++k1) {
	for(int k2=0; k2 < 3; ++k2) {
	  if ( pArgs->bPrecise ) {
	    //tmpf = log(fMix * baseLKs[k1] + (1.-fMix) * baseLKs[k2]);
	    tmpf = log(baseLKs[k1*3+k2]);
	    smMarkerLKs[k1*3+k2] += tmpf;
	    //rgMarkerLKs[rgIdx*9+k1*3+k2] += tmpf;
	  }
	  else {
	    //tmpf = (fMix * baseLKs[k1] + (1.-fMix) * baseLKs[k2]);
	    tmpf = baseLKs[k1*3+k2];
	    smMarkerLKs[k1*3+k2] *= tmpf; 
	    //rgMarkerLKs[rgIdx*9+k1*3+k2] *= tmpf;
	  }
	}
      }
    }

    // sample-level per-marker likelihood
    // smLLK = Pr(Data|fMix) = \prod_{markers} Pr(bases|fMix)
    //  = \prod_{markers} \sum_{G1,G2} Pr(bases|G1,G2,fMix) Pr(G1|AF)Pr(G2|AF)
    double perMarkerProb = 0;
    if ( pArgs->bPrecise ) {
      int maxIdx = 0;
      for(int k=1; k < 9; ++k) {
	if ( smMarkerLKs[maxIdx] < smMarkerLKs[k] ) {
	  maxIdx = k;
	}
      }
      for(int k1=0; k1 < 3; ++k1) {
	for(int k2=0; k2 < 3; ++k2) {
	  perMarkerProb += (exp(smMarkerLKs[k1*3+k2] - smMarkerLKs[maxIdx]) * genoProb[k1] * genoFreq[k2]);
	}
      }
      smLLK += (log(perMarkerProb) + smMarkerLKs[maxIdx]);
    }
    else {
      for(int k1=0; k1 < 3; ++k1) {
	for(int k2=0; k2 < 3; ++k2) {
	  perMarkerProb += (smMarkerLKs[k1*3+k2] * genoProb[k1] * genoFreq[k2]);
	}
      }
      smLLK += log(perMarkerProb);
    }
  }
  notice("double mpuVerifyID::computeIBDLLK( %.6lf | %.6lf, %.6lf, %.6lf ) = %.6lf",fIBD,pSN[0],pSN[1],pSN[2],smLLK);
  return smLLK;
}

// computeIBDLLKs() - calculate Pr(Data|fIBD) for each SM and RG
//                    given AF and call rate threshold
double mpuVerify::computePairIBDLLKs(double fIBD, int indIdx1, int indIdx2) {
  int nMarkers = (int)(pGenotypes->chroms.size());
  //int nRGs = (int)(pPile->vsRGIDs.size());
  double genoFreq[3], genoProb1[3], genoProb2[3], baseLKs[9], smMarkerLKs[9];
  double af, baseError, baseMatch;
  uint8_t a1, a2, base;
  double smLLK = 0;
  double tmpf;

  for(int i=0; i < nMarkers; ++i) {
    af = pGenotypes->alleleFrequencies[i];
    if ( af < 0.001 ) af = 0.001;
    if ( af > 0.999 ) af = 0.999;

    // frequency of genotypes under HWE
    genoFreq[0] = (1.-af)*(1.-af);
    genoFreq[1] = 2.*af*(1.-af);
    genoFreq[2] = af*af;

    int geno1 = pGenotypes->getGenotype(indIdx1,i);
    int geno2 = pGenotypes->getGenotype(indIdx2,i);
    switch(geno1) {
    case 0: // MISSING
      genoProb1[0] = genoFreq[0];
      genoProb1[1] = genoFreq[1];
      genoProb1[2] = genoFreq[2];
      break;
    case 1: // HOMREF;
      genoProb1[0] = 1.-pArgs->genoError;
      genoProb1[1] = pArgs->genoError * (GENO_ERR_TO_HET);
      genoProb1[2] = pArgs->genoError * (1.-GENO_ERR_TO_HET);
      break;
    case 2: // HET;
      genoProb1[0] = pArgs->genoError/2.;
      genoProb1[1] = 1.-pArgs->genoError;
      genoProb1[2] = pArgs->genoError/2.;
      break;
    case 3: // HET;
      genoProb1[0] = pArgs->genoError * (1.-GENO_ERR_TO_HET);
      genoProb1[1] = pArgs->genoError * GENO_ERR_TO_HET;
      genoProb1[2] = 1.-pArgs->genoError;
      break;
    default:
      notice("Unrecognized genotype %d at ind %d, marker %d",indIdx1,i);
    }

    if ( indIdx1 == indIdx2 ) geno2 = 0;
    switch(geno2) {
    case 0: // MISSING
      genoProb2[0] = genoFreq[0];
      genoProb2[1] = genoFreq[1];
      genoProb2[2] = genoFreq[2];
      break;
    case 1: // HOMREF;
      genoProb2[0] = 1.-pArgs->genoError;
      genoProb2[1] = pArgs->genoError * (GENO_ERR_TO_HET);
      genoProb2[2] = pArgs->genoError * (1.-GENO_ERR_TO_HET);
      break;
    case 2: // HET;
      genoProb2[0] = pArgs->genoError/2.;
      genoProb2[1] = 1.-pArgs->genoError;
      genoProb2[2] = pArgs->genoError/2.;
      break;
    case 3: // HET;
      genoProb2[0] = pArgs->genoError * (1.-GENO_ERR_TO_HET);
      genoProb2[1] = pArgs->genoError * GENO_ERR_TO_HET;
      genoProb2[2] = 1.-pArgs->genoError;
      break;
    default:
      notice("Unrecognized genotype %d at ind %d, marker %d",indIdx2,i);
    }

    a1 = BaseAsciiMap::base2int[(int)pGenotypes->refBases[i]];
    a2 = BaseAsciiMap::base2int[(int)pGenotypes->altBases[i]];

    // initialize the marker-level likelihoods
    for(int k=0; k < 9; ++k) {
      smMarkerLKs[k] = (pArgs->bPrecise ? 0. : 1.);
    }

    for(int j=(int)pPile->nBegins[i]; j < (int)pPile->nEnds[i]; ++j) {
      // obtain b (base), (error), and readgroup info
      base = pPile->nBases[j];
      baseError = fPhred2Err[pPile->nQuals[j]];
      baseMatch = 1.-baseError;
      //rgIdx = static_cast<int>(pPile->nRGIndices[j]);
      
      for(int k1=0; k1 < 3; ++k1) {
	for(int k2=0; k2 < 3; ++k2) {
	  baseLKs[k1*3+k2] = ( fIBD * pSN[0*3+k1] + (1.-fIBD) * pSN[0*3+k2] ) * ( base == a1 ? baseMatch : baseError/3.) + ( fIBD * pSN[1*3+k1] + (1.-fIBD) * pSN[1*3+k2] ) * ( base == a2 ? baseMatch : baseError/3. );
	}
      }
      
      // merge base-level likelihood into marker-level likelihood
      // smMarkerLKs[i*3+j] 
      //   = \prod Pr(b|G1,G2) = \prod (f Pr(b|G1) + (1-f) Pr(b|G2))
      for(int k1=0; k1 < 3; ++k1) {
	for(int k2=0; k2 < 3; ++k2) {
	  if ( pArgs->bPrecise ) {
	    //tmpf = log(fMix * baseLKs[k1] + (1.-fMix) * baseLKs[k2]);
	    tmpf = log(baseLKs[k1*3+k2]);
	    smMarkerLKs[k1*3+k2] += tmpf;
	    //rgMarkerLKs[rgIdx*9+k1*3+k2] += tmpf;
	  }
	  else {
	    //tmpf = (fMix * baseLKs[k1] + (1.-fMix) * baseLKs[k2]);
	    tmpf = baseLKs[k1*3+k2];
	    smMarkerLKs[k1*3+k2] *= tmpf; 
	    //rgMarkerLKs[rgIdx*9+k1*3+k2] *= tmpf;
	  }
	}
      }
    }

    // sample-level per-marker likelihood
    // smLLK = Pr(Data|fMix) = \prod_{markers} Pr(bases|fMix)
    //  = \prod_{markers} \sum_{G1,G2} Pr(bases|G1,G2,fMix) Pr(G1|AF)Pr(G2|AF)
    double perMarkerProb = 0;
    if ( pArgs->bPrecise ) {
      int maxIdx = 0;
      for(int k=1; k < 9; ++k) {
	if ( smMarkerLKs[maxIdx] < smMarkerLKs[k] ) {
	  maxIdx = k;
	}
      }
      for(int k1=0; k1 < 3; ++k1) {
	for(int k2=0; k2 < 3; ++k2) {
	  //perMarkerProb += (exp(smMarkerLKs[k1*3+k2] - smMarkerLKs[maxIdx]) * genoProb[k1] * genoFreq[k2]);
	  perMarkerProb += (exp(smMarkerLKs[k1*3+k2] - smMarkerLKs[maxIdx]) * genoProb1[k1] * genoProb2[k2]);
	}
      }
      smLLK += (log(perMarkerProb) + smMarkerLKs[maxIdx]);
    }
    else {
      for(int k1=0; k1 < 3; ++k1) {
	for(int k2=0; k2 < 3; ++k2) {
	  //perMarkerProb += (smMarkerLKs[k1*3+k2] * genoProb[k1] * genoFreq[k2]);
	  perMarkerProb += (smMarkerLKs[k1*3+k2] * genoProb1[k1] * genoProb2[k2]);
	}
      }
      smLLK += log(perMarkerProb);
    }
  }
  notice("double mpuVerifyID::computeIBDLLK( %.6lf | %.6lf, %.6lf, %.6lf ) = %.6lf",fIBD,pSN[0],pSN[1],pSN[2],smLLK);
  return smLLK;
}

// computeIBDLLKs() - calculate Pr(Data|fIBD) for each SM and RG
//                    given AF and call rate threshold
void mpuVerify::calculateDepthByGenotype(int indIdx, mpuVerifyOut &vbo) {
  int nMarkers = (int)(pGenotypes->chroms.size());

  for(int i=0; i < nMarkers; ++i) {
    int geno = pGenotypes->getGenotype(indIdx,i);
    int genoIdx = geno;
    ++vbo.numGenos[genoIdx];
    vbo.numReads[genoIdx] += ( pPile->nEnds[i] - pPile->nBegins[i] );
  }
  return;
}

void mpuVerify::calculateDepthDistribution(int maxDepth, mpuVerifyOut &vbo) {
  int dp = 0;
  vbo.numGenos[0] = nMarkers;
  for(int i=0; i < nMarkers; ++i) {
    dp = 0;
    for(int j=(int)pPile->nBegins[i]; j < (int)pPile->nEnds[i]; ++j) {
      ++dp;
      ++vbo.numReads[0];
    }
    ++vbo.depths[dp];
  }
  return;
}

void mpuVerify::printPerMarkerInfo(const char* filename, int indIdx) {
  wFile oFile(filename);
  int nMarkers = (int)(pGenotypes->chroms.size());
  uint8_t base, a1, a2;

  oFile.printf("#CHROM\tPOS\tA1\tA2\tAF\tGENO\t#REF\t#ALT\t#OTHERS\tBASES\tQUALS\tMAPQS\n");
  for(int i=0; i < nMarkers; ++i) {
    int counts[3] = {0,0,0};
    std::vector<char> bases;
    std::vector<char> quals;
    std::vector<char> mqs;

    oFile.printf("%s\t%d\t%c\t%c\t%.4lf\t",pGenotypes->chroms[i].c_str(),pGenotypes->positions[i],pGenotypes->refBases[i],pGenotypes->altBases[i],pGenotypes->alleleFrequencies[i]);
    int geno = pGenotypes->getGenotype(indIdx,i);
    switch(geno) {
    case 0: // MISSING
      oFile.printf("./.");
      break;
    case 1: // HOMREF;
      oFile.printf("0/0");
      break;
    case 2: // HET;
      oFile.printf("0/1");
      break;
    case 3: // HOMALT;
      oFile.printf("1/1");
      break;
    default:
      error("Unrecognized genotype %d at ind %d, marker %d",indIdx,i);
    }

    a1 = BaseAsciiMap::base2int[(int)pGenotypes->refBases[i]];
    a2 = BaseAsciiMap::base2int[(int)pGenotypes->altBases[i]];

    for(int j=(int)pPile->nBegins[i]; j < (int)pPile->nEnds[i]; ++j) {
      base = pPile->nBases[j];
      if ( base == a1 ) {
	++counts[0];
      }
      else if ( base == a2 ) {
	++counts[1];
      }
      else {
	++counts[2];
      }
      bases.push_back(base);
      quals.push_back(pPile->nQuals[j]+33);
      mqs.push_back(((uint8_t)(pPile->nMapQs[j]) > 90) ? '~' : static_cast<char>(pPile->nMapQs[j]+33));
    }
    oFile.printf("\t%d\t%d\t%d\t%.3lf\t",counts[0],counts[1],counts[2],(counts[0]+counts[1] == 0) ? 0.5 : (double)counts[0]/(double)(counts[0]+counts[1]));

    oFile.printf("\t");
    for(int j=0; j < (int)bases.size(); ++j)
      oFile.printf("%c",bases[j]);

    oFile.printf("\t");
    for(int j=0; j < (int)quals.size(); ++j)
      oFile.printf("%c",quals[j]);

    oFile.printf("\t");
    for(int j=0; j < (int)mqs.size(); ++j)
      oFile.printf("%c",mqs[j]);

    oFile.printf("\n");
  }
}

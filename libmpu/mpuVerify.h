#ifndef __MPU_VERIFY_H
#define __MPU_VERIFY_H

#include <vector>
#include <string>
#include <cmath>
#include "MathGold.h"
#include "mpuFile.h"
#include "mpuPileBases.h"

#define MAX_Q 100
#define MAX_DBL 1e99
#define GENO_ERR_TO_HET 0.5

class mpuVerifyArgs {
public:
  std::string sSubsetInds;
  std::string sVcfFile;
  std::string sMpuFile;
  std::string sOutFile;
  std::string sSMID;
  std::string sScanPair;

  bool bSelfOnly;
  bool bSiteOnly;
  bool bFindBest;
  bool bIgnoreSex;

  //bool bIgnoreOverlapPair;
  bool bPrecise;
  bool bVerbose;
  bool bSilent;

  double genoError;
  double minAF;
  double minCallRate;
  double contamThres;

  double pRefRef;
  double pRefHet;
  double pRefAlt;

  bool bFreeNone;
  bool bFreeMixOnly;
  bool bFreeRefBiasOnly;
  bool bFreeFull;

  bool bChipNone;
  bool bChipMixOnly;
  bool bChipRefBiasOnly;
  bool bChipFull;

  int minMapQ;
  int maxDepth;
  int minQ;
  int maxQ;

  double grid;

  mpuVerifyArgs() {
    bSelfOnly = false;
    bSiteOnly = false;
    bFindBest = false;
    bIgnoreSex = false;

    bPrecise = false;

    genoError = 0.001;
    minAF = 0.01;
    minCallRate = 0.50;
    contamThres = 0.02;

    pRefRef = 1;
    pRefHet = 0.5;
    pRefAlt = 0;

    bFreeNone = false;
    bFreeMixOnly = true;
    bFreeRefBiasOnly = false;
    bFreeFull = false;

    bChipNone = false;
    bChipMixOnly = true;
    bChipRefBiasOnly = false;
    bChipFull = false;

    bVerbose = false;

    minMapQ = 10;
    maxDepth = 20;
    minQ = 13;
    maxQ = 40;

    grid = 0.05;
  }
};

class GenMatrixBinary {
public:
  //int nInds;
  std::vector<std::string> chroms;
  std::vector<int> positions;
  std::vector<double> alleleFrequencies;
  std::vector<char> refBases;
  std::vector<char> altBases;

  std::vector<std::string> indids;

  std::vector<uint8_t> genotypes;
  int bytesPerMarker;

  GenMatrixBinary(const char* vcfFile, bool siteOnly, bool findBest, std::vector<std::string>& subsetInds, double minAF, double minCallRate);
  int addMarker(const char* chrom, int position, char refBase, char altBase, double alleleFreq);
  void setGenotype(float genotype, int indIndex, int markerIndex = -1);
  int getGenotype(int indIndex, int markerIndex);
  double computeAlleleFrequency(int markerIndex);
};

class mpuVerifyOut {
 public:
  mpuVerifyOut() {}

  void init(double pRefRef, double pRefHet, double pRefAlt, int maxDepth = -1) {
    numGenos.resize(4,0);
    numReads.resize(4,0);
    depths.resize((maxDepth+1),0);
  }
  
  double llk0;
  double llk1;
  double fMix;
  double refRef;
  double refHet;
  double refAlt;

  std::vector<int> numGenos;
  std::vector<int> numReads;

  std::vector<int> depths;
};

class mpuVerify {
 public:
  ////////////////////////////////////////////
  // core member variables
  ////////////////////////////////////////////
  mpuPileBases* pPile;
  GenMatrixBinary* pGenotypes;   // genotypes containing AFs
  mpuVerifyArgs* pArgs;        // arguments

  mpuVerifyOut mixOut;
  mpuVerifyOut selfOut;
  mpuVerifyOut bestOut;

  int nMarkers;     // # SNPs
  int nBases;

  bool inferProbRefs;
  double pSN[6]; // pSN[i*3+j] == Pr(Sampled=i|True=j,H0); 0(REF), 1(ALT), 2(others)
  std::vector<std::string> subsetInds;
  double fPhred2Err[MAX_Q+1];

  static double logit(double p) { return log(p/(1.-p)); }
  static double invLogit(double x) { double e = exp(x); return e/(1.+e); }

  ///////////////////////////////////////////
  // core member functions
  //////////////////////////////////////////

  mpuVerify(mpuVerifyArgs* p) : pPile(NULL), pGenotypes(NULL), pArgs(p) {
    pSN[0*3+0] = pArgs->pRefRef; 
    pSN[0*3+1] = pArgs->pRefHet;
    pSN[0*3+2] = pArgs->pRefAlt;
    pSN[1*3+0] = 1.-pArgs->pRefRef; 
    pSN[1*3+1] = 1.-pArgs->pRefHet;
    pSN[1*3+2] = 1.-pArgs->pRefAlt;

    for(int i=0; i < MAX_Q+1; ++i) {
      // If Phred>=maxQ, assume that the base quality is
      // overestimated and apply an upper threshold.
      if ( i > static_cast<int>(p->maxQ) ) {
	fPhred2Err[i] = fPhred2Err[p->maxQ]; 
      }
      else {
	fPhred2Err[i] = pow(10.,(0-i)/10.);
      }
    }
    nMarkers = nBases = 0;
    subsetInds.clear();
  }

  ~mpuVerify() {
    //if ( pPile != NULL )
      //delete pPile;
      //if ( pGenotypes != NULL )
      //delete pGenotypes;
  }

  void setRefBiasParams(double pRefBaseRefGeno, double pRefBaseHetGeno, double pRefBaseAltGeno) {
    // Pr(ref||RR) = 1
    pSN[0*3+0] = 1;
    // Pr(ref|RA) = pRefBaseHetGeno
    pSN[0*3+1] = pRefBaseHetGeno;
    // Pr(ref|AA) = pRefBaseAltGeno
    pSN[0*3+2] = pRefBaseAltGeno;
    // Pr(alt|RR) = 0
    pSN[1*3+0] = 0;
    // Pr(alt|RA)
    pSN[1*3+1] = 1.-pRefBaseHetGeno;
    // Pr(alt|AA)
    pSN[1*3+2] = 1.-pRefBaseAltGeno;
    //pSN[2*3+0] = pSN[2*3+1] = pSN[2*3+2] = 0;
  }

  void loadSubsetInds(const char* subsetFile);
  void loadFiles(const char* mpuFile, const char* vcfFile);

  double computeMixLLKs(double fMix);
  double computeIBDLLKs(double fIBD, int indIdx);
  double computePairIBDLLKs(double fIBD, int indIdx1, int indIdx2);
  void calculateDepthByGenotype(int indIdx, mpuVerifyOut &vbo);
  void calculateDepthDistribution(int maxDepth, mpuVerifyOut &vbo);
  void printPerMarkerInfo(const char* outfile, int indIndex);

  class refBiasMixLLKFunc : public VectorFunc {
  public:
    mpuVerify* pVBI;

    double llk0;
    double llk1;
    double pRefHet;
    double pRefAlt;

    refBiasMixLLKFunc(mpuVerify* p) : pVBI(p), llk0(MAX_DBL), llk1(MAX_DBL), pRefHet(0), pRefAlt(0) {
      pVBI->setRefBiasParams(1.,0.5,0);
      llk1 = llk0 = (0-pVBI->computeMixLLKs(0));
      pRefHet = 0.5;
      pRefAlt = 0;
  }

    virtual double Evaluate(Vector& v) {
      if ( v.Length() != 2 )
	error("refBiasLLK(): Input vector must be length of 2");

      double refHet = invLogit(v[0]);
      double refAlt = invLogit(v[1]);

      pVBI->setRefBiasParams(1.,refHet,refAlt);
      double smLLK = 0-pVBI->computeMixLLKs(0);

      if ( smLLK < llk1 ) {
	llk1 = smLLK;
	pRefHet = refHet;
	pRefAlt = refAlt;
      }

      return smLLK;
    }
  };

  class fullMixLLKFunc : public VectorFunc {
  public:
    mpuVerify* pVBI;

    double llk0;
    double llk1;
    double fMix;
    double pRefHet;
    double pRefAlt;

  fullMixLLKFunc(mpuVerify* p) : pVBI(p),llk0(MAX_DBL), llk1(MAX_DBL), fMix(0), pRefHet(0), pRefAlt(0) {
      pVBI->setRefBiasParams(1.,0.5,0);
      llk1 = llk0 = (0-pVBI->computeMixLLKs(0));
      pRefHet = 0.5;
      pRefAlt = 0;
    }

    virtual double Evaluate(Vector& v) {
      if ( v.Length() != 3 ) 
	error("fullMixLLKFunc(): Input vector must be length of 3");

      double mix = invLogit(v[0]);
      double refHet = invLogit(v[1]);
      double refAlt = invLogit(v[2]);

      pVBI->setRefBiasParams(1.,refHet,refAlt);
      double smLLK = 0-pVBI->computeMixLLKs(mix); // rgLLKs);

      if ( smLLK < llk1 ) {
	llk1 = smLLK;
	fMix = mix;
	pRefHet = refHet;
	pRefAlt = refAlt;
      }

      return smLLK;
    }
  };

  class refBiasIbdLLKFunc : public VectorFunc {
  public:
    mpuVerify* pVBI;
    int indIdx;

    double llk0;
    double llk1;
    double pRefHet;
    double pRefAlt;

  refBiasIbdLLKFunc(mpuVerify* p) : pVBI(p), indIdx(0), llk0(MAX_DBL), llk1(MAX_DBL), pRefHet(0), pRefAlt(0) {
      // calculate null likelihood
      pVBI->setRefBiasParams(1.,0.5,0);
      llk1 = llk0 = (0-pVBI->computeIBDLLKs(1., 0));
      pRefHet = 0.5;
      pRefAlt = 0;
  }

    virtual double Evaluate(Vector& v) {
      if ( v.Length() != 2 ) 
	error("refBiasLLK(): Input vector must be length of 2");

      double refHet = invLogit(v[0]);
      double refAlt = invLogit(v[1]);

      pVBI->setRefBiasParams(1.,refHet,refAlt);
      double smLLK = 0-pVBI->computeIBDLLKs(1, indIdx); //rgLLKs, indIdx);

      if ( smLLK < llk1 ) {
	llk1 = smLLK;
	pRefHet = refHet;
	pRefAlt = refAlt;
      }
      return smLLK;
    }
  };

  class fullIbdLLKFunc : public VectorFunc {
  public:
    mpuVerify* pVBI;
    int indIdx;
    double llk0;
    double llk1;
    double fIBD;
    double pRefHet;
    double pRefAlt;
  fullIbdLLKFunc(mpuVerify* p, int ind) : pVBI(p), indIdx(ind), llk0(MAX_DBL), llk1(MAX_DBL), fIBD(1), pRefHet(0), pRefAlt(0) {
      pVBI->setRefBiasParams(1.,0.5,0);
      llk1 = llk0 = (0-pVBI->computeIBDLLKs(1, ind));
      pRefHet = 0.5;
      pRefAlt = 0;
    }

    virtual double Evaluate(Vector& v) {
      if ( v.Length() != 3 ) 
	error("fullIbdLLKFunc(): Input vector must be length of 3");

      double ibd = invLogit(v[0]);
      double refHet = invLogit(v[1]);
      double refAlt = invLogit(v[2]);

      pVBI->setRefBiasParams(1.,refHet,refAlt);
      double smLLK = 0-pVBI->computeIBDLLKs(ibd, indIdx);

      if ( smLLK < llk1 ) {
	llk1 = smLLK;
	fIBD = ibd;
	pRefHet = refHet;
	pRefAlt = refAlt;
      }
      return smLLK;
    }
  };

  class ibdLLK : public ScalarMinimizer {
  public:
    mpuVerify* pVBI;
    int indIdx;
    double llk0;
    double llk1;
    double fIBD;

  ibdLLK(mpuVerify* p) : pVBI(p), indIdx(0), llk0(0), llk1(0), fIBD(0) {
    }

    virtual double f(double fIBD) {
      double smLLK = pVBI->computeIBDLLKs(fIBD, indIdx);
      return 0-smLLK;
    }

    double OptimizeLLK(int ind = 0) {
      indIdx = ind;

      std::vector<double> lks;
      std::vector<double> alphas;
      double minLK = f(1.);
      lks.push_back(minLK);
      alphas.push_back(0.);
      double grid = pVBI->pArgs->grid;

      int minIdx = 0;
      double alpha;
      for(alpha = grid; alpha < 0.999; alpha += grid) {
	lks.push_back(f(1.-alpha));
	alphas.push_back(alpha);
	if ( lks.back() < minLK ) {
	  minLK = lks.back();
	  minIdx = (int)lks.size()-1;
	}
      }

      if ( minIdx == 0 ) {
	a = 1.; fa = lks[0];
	b = 1.-grid/2.; fb = f(b);
	c = 1.-grid; fc = lks[1];
      }
      else if ( minIdx == (int)lks.size()-1 ) {
	a = 0; fa = f(0);
	b = (1.-alphas.back())/2.; fb = f(b);
	c = 1.-alphas.back(); fc = lks.back();
      }
      else {
	a = 1.-alphas[minIdx-1]; fa = lks[minIdx-1];
	b = 1.-alphas[minIdx]; fb = lks[minIdx];
	c = 1.-alphas[minIdx+1]; fc = lks[minIdx+1];
      }

      Brent(0.0001);

      llk0 = lks[0];
      llk1 = fmin;
      fIBD = min;

      return min;
    }
  };

  class pairLLK : public ScalarMinimizer {
  public:
    mpuVerify* pVBI;
    int indIdx1;
    int indIdx2;
    double llk0;
    double llk1;
    double fIBD;

  pairLLK(mpuVerify* p) : pVBI(p), indIdx1(0), indIdx2(0), llk0(0), llk1(0), fIBD(0) {
    }

    virtual double f(double fIBD) {
      double smLLK = pVBI->computePairIBDLLKs(fIBD, indIdx1, indIdx2);
      return 0-smLLK;
    }

    double OptimizeLLK(int ind1 = 0, int ind2 = 0) {
      indIdx1 = ind1;
      indIdx2 = ind2;

      std::vector<double> lks;
      std::vector<double> alphas;
      double minLK = f(1.);
      lks.push_back(minLK);
      alphas.push_back(0.);
      double grid = pVBI->pArgs->grid;

      int minIdx = 0;
      double alpha;
      for(alpha = grid; alpha < 0.999; alpha += grid) {
	lks.push_back(f(1.-alpha));
	alphas.push_back(alpha);
	if ( lks.back() < minLK ) {
	  minLK = lks.back();
	  minIdx = (int)lks.size()-1;
	}
      }

      if ( minIdx == 0 ) {
	a = 1.; fa = lks[0];
	b = 1.-grid/2.; fb = f(b);
	c = 1.-grid; fc = lks[1];
      }
      else if ( minIdx == (int)lks.size()-1 ) {
	a = 0; fa = f(0);
	b = (1.-alphas.back())/2.; fb = f(b);
	c = 1.-alphas.back(); fc = lks.back();
      }
      else {
	a = 1.-alphas[minIdx-1]; fa = lks[minIdx-1];
	b = 1.-alphas[minIdx]; fb = lks[minIdx];
	c = 1.-alphas[minIdx+1]; fc = lks[minIdx+1];
      }

      Brent(0.0001);

      llk0 = lks[0];
      llk1 = fmin;
      fIBD = min;

      return min;
    }
  };

  class mixLLK : public ScalarMinimizer {
  public:
    mpuVerify* pVBI;
    double llk0;
    double llk1;
    double fMix;

  mixLLK(mpuVerify* p) : pVBI(p), llk0(0), llk1(0), fMix(0) {}
    virtual double f(double fMix) {
      double smLLK = pVBI->computeMixLLKs(fMix); 
      return 0-smLLK;
    }
    
    double OptimizeLLK(int rg = -1) {
      std::vector<double> lks;
      std::vector<double> alphas;
      double minLK = f(0);
      lks.push_back(minLK);
      alphas.push_back(0.);
      double grid = pVBI->pArgs->grid;

      int minIdx = 0;
      double alpha;
      for(alpha = grid; alpha < 0.499; alpha += grid) {
	lks.push_back(f(alpha));
	alphas.push_back(alpha);
	if ( lks.back() < minLK ) {
	  minLK = lks.back();
	  minIdx = (int)lks.size()-1;
	}
      }
      if ( minIdx == 0 ) {
	a = 0.; fa = lks[0];
	b = grid/2.; fb = f(b);
	c = grid; fc = lks[1];
      }
      else if ( minIdx == (int)lks.size()-1 ) {
	a = 0.5; fa = f(0.5);
	b = (0.5-alphas.back())/2.; fb = f(b);
	c = alphas.back(); fc = lks.back();
      }
      else {
	a = alphas[minIdx-1]; fa = lks[minIdx-1];
	b = alphas[minIdx]; fb = lks[minIdx];
	c = alphas[minIdx+1]; fc = lks[minIdx+1];
      }

      Brent(0.0001);

      llk0 = lks[0];
      llk1 = fmin;
      fMix = min;

      return min;
    }
  };
};



#endif

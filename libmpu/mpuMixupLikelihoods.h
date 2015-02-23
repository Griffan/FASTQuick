#ifndef __MPU_MIXUP_LIKELIHOODS_H__
#define __MPU_MIXUP_LIKELIHOODS_H__

#include "mpuFile.h"
//#include "MathGold.h"
#include "fPed.h"

#include <math.h>

// Chromosome types
#define CT_AUTOSOME     0
#define CT_CHRX         1
#define CT_CHRY         2
#define CT_MITO         3

class mpuMixupLikelihood { //: public ScalarMinimizer {
 public:
  int        n;
  mpuFile    * mpu;
  int        * pairs;
  //char       * sexes;
  int        chromosomeType;
  bool       hasContam;
  double     inbreedingCoeff;

  double nullLLK, hweAF, hweLLK, hwdAF1, hwdAF2, hwdLLK, Acount, ABcount;
  
  mpuMixupLikelihood(int count, mpuFile * mpuPointers, int* pairPointers, bool contam = false, double fic = 0);
  virtual ~mpuMixupLikelihood();
  
  void SetAlleles(int al1, int al2);
  
  //virtual double f(double freq) { return -Evaluate(freq); }
  
  //virtual double Evaluate(double freq);
  
  void GetPriors(double * priors, double freq, int i);
  void GetMalePriors(double * priors, double freq);
  void GetFemalePriors(double * priors, double freq);
  
  void OptimizeFrequencyEM(bool hwe = true);
  
  static int GenotypeIndex(int base1, int base2)
  {
    return base1 < base2 ? (base1) *(9 - base1) / 2 + (base2 - base1) :
      (base2) *(9 - base2) / 2 + (base1 - base2);
  }
  
 protected:
  int allele1, allele2;
  int geno11, geno12, geno22;
};
#endif

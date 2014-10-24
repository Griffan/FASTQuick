#include "mpuMixupLikelihoods.h"
#include "mpuFile.h"
#include "Constant.h"
#include "Error.h"

mpuMixupLikelihood::mpuMixupLikelihood(int count, mpuFile * mpuPointers, int * pairPointers, bool contam, double fic)
{
  n = count;
  mpu = mpuPointers;
  pairs = pairPointers;
  hasContam = contam;
  inbreedingCoeff = fic;
  
  //sexes = new char [n];

  //for (int i = 0; i < n; i++)
  //sexes[i] = SEX_FEMALE;
  
  chromosomeType = CT_AUTOSOME;
}

mpuMixupLikelihood::~mpuMixupLikelihood()
{
  //if (sexes != NULL)
  //delete [] sexes;
}

void mpuMixupLikelihood::SetAlleles(int al1, int al2)
{
  allele1 = al1;
  allele2 = al2;

  geno11 = GenotypeIndex(allele1, allele1);
  geno12 = GenotypeIndex(allele1, allele2);
  geno22 = GenotypeIndex(allele2, allele2);
}

/*
double mpuMixupLikelihood::Evaluate(double freq)
{
  double prior11 = freq * freq;
  double prior12 = freq * (1.0 - freq) * 2.0;
  double prior22 = (1.0 - freq) * (1.0 - freq);

  double prior1 = freq;
  double prior2 = 1.0 - freq;
  
  double likelihood = 0.0;
  double *lks1, *lks2;

  double fracMale = .5; // assume 50/50 mixture when sex of contaminating genotype is unknown

  if ( inbreedingCoeff > 0 ) {
    double fpq = freq * (1.0-freq) * inbreedingCoeff;
    prior11 += fpq;
    prior12 -= (2*fpq);
    prior22 += fpq;
  }

  switch (chromosomeType) {
  case CT_MITO :     // for mitochondrial genomes, use haploid model ignoring sex information
    prior11 = prior1;
    prior12 = 0.0;
    prior22 = prior2;
  case CT_AUTOSOME : // for autosomes, use diploid model
    for (int i = 0; i < n; i++) {  // calculate likelihood for unmatched pair
      if ( pairs[i] == -1 ) {  
	lks1 = mpu[i].LKs;
	likelihood += log(prior11 * prior11 * lks1[0] +
			  prior11 * prior12 * lks1[1] + 
			  prior11 * prior22 * lks1[2] +
			  prior12 * prior11 * lks1[3] +
			  prior12 * prior12 * lks1[4] +
			  prior12 * prior22 * lks1[5] +
			  prior22 * prior11 * lks1[6] +
			  prior22 * prior12 * lks1[7] +
			  prior22 * prior22 * lks1[8] +
			  1e-30);
      }
      else if ( pairs[i] < i ) { // calculate likelihood for each pair only once
	lks1 = mpu[i].LKs;
	lks2 = mpu[pairs[i]].LKs;
	likelihood += log(prior11 * prior11 * lks1[0] * lks2[0] +
			  prior11 * prior12 * lks1[1] * lks2[3] + 
			  prior11 * prior22 * lks1[2] * lks2[6] +
			  prior12 * prior11 * lks1[3] * lks2[1] +
			  prior12 * prior12 * lks1[4] * lks2[4] +
			  prior12 * prior22 * lks1[5] * lks2[7] +
			  prior22 * prior11 * lks1[6] * lks2[2] +
			  prior22 * prior12 * lks1[7] * lks2[5] +
			  prior22 * prior22 * lks1[8] * lks2[8] +
			  1e-30);
      }
    }
    break;
  case CT_CHRY :   // for Y chromosomes, add up only male genotypes
    for (int i = 0; i < n; i++) {
      if ( pairs[i] == -1 ) {   
	lks1 = mpu[i].LKs;
	if ( mpu[i].sex == SEX_MALE )  // Assume 50/50 mixture  REF/. p REF/REF p2 REF/ALT ALT/. ALT/ALT ALT/ALT
	  likelihood += log(prior1 * ( fracMale * prior1 + (1.-fracMale) ) * lks1[0] +  // REF/REF or REF/.
			    prior1 * ( fracMale * prior2 ) * lks1[2] +       // REF/ALT
			    prior2 * ( fracMale * prior1 ) * lks1[6] +       // ALT/REF
			    prior2 * ( fracMale * prior2 + (1.-fracMale) ) * lks1[8] +   // ALT/ALT or ALT/.
			    1e-30);
      }
      else if ( pairs[i] < i ) {
	lks1 = mpu[i].LKs;
	lks2 = mpu[pairs[i]].LKs;
	if ( mpu[i].sex == SEX_MALE ) {
	  if ( mpu[pairs[i]].sex == SEX_MALE ) {
	    likelihood += log(prior1 * prior1 * lks1[0] * lks2[0] +  // REF/REF
			      prior1 * prior2 * lks1[2] * lks2[6] +  // REF/ALT
			      prior2 * prior1 * lks1[6] * lks2[2] +  // ALT/REF
			      prior2 * prior2 * lks1[8] * lks2[8] +  // ALT/ALT
			      1e-30);
	  }
	  else {   // MALE-FEMALE mixture
	    likelihood += log(prior1 * lks1[0] * lks2[0] +  // REF/.
			      prior2 * lks1[8] * lks2[8] +  // ALT/.
			      1e-30);
	  }
	}
	else {
	  if ( mpu[pairs[i]].sex == SEX_MALE ) { // FEMALE-MALE mixture
	    likelihood += log(prior1 * lks1[0] * lks2[0] +  // REF/.
			      prior2 * lks1[8] * lks2[8] +  // ALT/.
			      1e-30);
	  }
	  //else {} // do nothing for FEMALE-FEMALE mixture
	}
      }
    }
    break;
  case CT_CHRX :  // for non-PAR chromosomes, add male and female genotypes proportionally
    for (int i = 0; i < n; i++) {
      if ( pairs[i] == -1 ) { // assume 50/50 mixture of MALE-FEMALE
	lks1 = mpu[i].LKs;  
	if ( mpu[i].sex == SEX_MALE ) {  // MALE vs UNKNOWN
	  likelihood += log(prior1 * ( fracMale * prior1 + (1.-fracMale) * prior11 ) * lks1[0] +  // REF/REF
			    prior1 * ( fracMale * prior1 + (1.-fracMale) * prior12 ) * lks1[1] +  // REF/HET
			    prior1 * ( fracMale * prior1 + (1.-fracMale) * prior22 ) * lks1[2] +  // REF/ALT
			    prior2 * ( fracMale * prior1 + (1.-fracMale) * prior11 ) * lks1[6] +  // ALT/REF
			    prior2 * ( fracMale * prior1 + (1.-fracMale) * prior12 ) * lks1[7] +  // ALT/HET
			    prior2 * ( fracMale * prior1 + (1.-fracMale) * prior22 ) * lks1[8] +  // ALT/ALT
			    1e-30);
	}
	else { // FEMALE vs UNKNOWN
	  likelihood += log(prior11 * ( fracMale * prior1 + (1.-fracMale) * prior11 ) * lks1[0] +  // REF/REF
			    prior11 * ( fracMale * prior1 + (1.-fracMale) * prior12 ) * lks1[1] +  // REF/HET
			    prior11 * ( fracMale * prior1 + (1.-fracMale) * prior22 ) * lks1[2] +  // REF/ALT
			    prior12 * ( fracMale * prior1 + (1.-fracMale) * prior11 ) * lks1[3] +  // HET/REF
			    prior12 * ( fracMale * prior1 + (1.-fracMale) * prior12 ) * lks1[4] +  // HET/HET
			    prior12 * ( fracMale * prior1 + (1.-fracMale) * prior22 ) * lks1[5] +  // HET/ALT
			    prior22 * ( fracMale * prior1 + (1.-fracMale) * prior11 ) * lks1[6] +  // ALT/REF
			    prior22 * ( fracMale * prior1 + (1.-fracMale) * prior12 ) * lks1[7] +  // ALT/HET
			    prior22 * ( fracMale * prior1 + (1.-fracMale) * prior22 ) * lks1[8] +  // ALT/ALT
			    1e-30);
	}
      }
      else if ( pairs[i] < i ) {
	lks1 = mpu[i].LKs;
	lks2 = mpu[pairs[i]].LKs;
	if ( mpu[i].sex == SEX_MALE ) {
	  if ( mpu[pairs[i]].sex == SEX_MALE ) {
	    likelihood += log(prior1 * prior1 * lks1[0] * lks2[0] +  // REF/REF
			      prior1 * prior2 * lks1[2] * lks2[6] +  // REF/ALT
			      prior2 * prior1 * lks1[6] * lks2[2] +  // ALT/REF
			      prior2 * prior2 * lks1[8] * lks2[8] +  // ALT/ALT
			      1e-30);
	  }
	  else {   // MALE-FEMALE mixture
	    likelihood += log(prior1 * prior11 * lks1[0] * lks2[0] +  // REF/REF
			      prior1 * prior12 * lks1[1] * lks2[3] +  // REF/HET
			      prior1 * prior22 * lks1[2] * lks2[6] +  // REF/ALT
			      prior2 * prior11 * lks1[6] * lks2[2] +  // ALT/REF
			      prior2 * prior12 * lks1[7] * lks2[5] +  // ALT/HET
			      prior2 * prior22 * lks1[8] + lks2[8] +  // ALT/ALT
			      1e-30);
	  }
	}
	else {
	  if ( mpu[pairs[i]].sex == SEX_MALE ) { // FEMALE-MALE mixture
	    likelihood += log(prior11 * prior1 * lks1[0] * lks2[0] +  // REF/REF
			      prior11 * prior2 * lks1[2] * lks2[6] +  // REF/ALT
			      prior12 * prior1 * lks1[3] * lks2[1] +  // HET/REF
			      prior12 * prior2 * lks1[5] * lks2[7] +  // HET/ALT
			      prior22 * prior1 * lks1[6] * lks2[2] +  // ALT/REF
			      prior22 * prior2 * lks1[8] + lks2[8] +  // ALT/ALT
			      1e-30);
	  }
	  else {  // FEMALE-FEMALE mixture is the same to autosomal chromosomes
	    likelihood += log(prior11 * prior11 * lks1[0] * lks2[0] +
			      prior11 * prior12 * lks1[1] * lks2[3] + 
			      prior11 * prior22 * lks1[2] * lks2[6] +
			      prior12 * prior11 * lks1[3] * lks2[1] +
			      prior12 * prior12 * lks1[4] * lks2[4] +
			      prior12 * prior22 * lks1[5] * lks2[7] +
			      prior22 * prior11 * lks1[6] * lks2[2] +
			      prior22 * prior12 * lks1[7] * lks2[5] +
			      prior22 * prior22 * lks1[8] * lks2[8] +
			      1e-30);
	  }
	}
      }
    }
  }
  return likelihood;
}
*/

void mpuMixupLikelihood::GetPriors(double * priors, double freq, int i)
{
  if (mpu[i].sex == SEX_MALE)
    GetMalePriors(priors, freq);
  else
    GetFemalePriors(priors, freq);
}

void mpuMixupLikelihood::GetMalePriors(double * priors, double freq)
{
  double fpq = freq * (1.0-freq) * inbreedingCoeff;
  
  switch (chromosomeType)
    {
    case CT_AUTOSOME:
      priors[0] = freq * freq + fpq;
      priors[1] = 2 * (1. - freq) * freq * (1-fpq);
      priors[2] = (1. - freq) * (1. - freq) + fpq;
      break;
    case CT_CHRY: case CT_CHRX: case CT_MITO:
      priors[0] = freq;        /* would be zero for females */
      priors[1] = 0.0;
      priors[2] = 1. - freq;   /* would be zero for females */
      break;
    default:
      error("Cannot recognize chromosome type %d",chromosomeType);
    }
}

void mpuMixupLikelihood::GetFemalePriors(double * priors, double freq)
{
  double fpq = freq * (1.0-freq) * inbreedingCoeff;
  
  switch (chromosomeType)
    {
    case CT_AUTOSOME : case CT_CHRX:
      priors[0] = freq * freq + fpq;
      priors[1] = 2 * (1. - freq) * freq * (1. - fpq);
      priors[2] = (1. - freq) * (1. - freq) + fpq;
      break;
    case CT_CHRY: 
      priors[0] = 0.0;            /* would be freq for males */
      priors[1] = 0.0;
      priors[2] = 0.0;            /* would be 1. - freq for males */
      break;
    case CT_MITO :
      priors[0] = freq;
      priors[1] = 0;
      priors[2] = 1. - freq;
      break;
    default:
      error("Cannot recognize chromosome type %d",chromosomeType);
    }
}

// Optimize frequency 
void mpuMixupLikelihood::OptimizeFrequencyEM(bool hwe) {
  // fit AF with/without HWE
  double p = .05 + rand()/(RAND_MAX+1.)*.9; // start from random AF [.05-.95]
  double fracMale = .5; // assume 50/50 mixture when sex of contaminating genotype is unknown
  double q = 1.-p;
  int rounds, maxiter = 20;
  double llk, sum1, sum2, prior11, prior12, prior22, prior1, prior2, l, p1, p2, scale, delta, nref;


  if ( hwe ) { // when hwe is set, calculate the null likelihood
    nullLLK = 0;
    for (int i = 0; i < n; i++) {  // calculate likelihoods
      nullLLK += mpu[i].GLs[0]; // logLLK for no SNP
    }
    nullLLK *= log(10);
  }

  // start with HWE assumption from random guess
  prior11 = q * q;
  prior12 = p * q * 2.0;
  prior22 = p * p;
  prior1 = q;
  prior2 = p; 

  if ( chromosomeType == CT_MITO ) {
    prior11 = prior1;
    prior12 = 0.0;
    prior22 = prior2;
  }
  else if ( inbreedingCoeff > 0 ) {
    double fpq = p * q * inbreedingCoeff;
    prior11 += fpq;
    prior12 -= (2*fpq);
    prior22 += fpq;
  }
   
  for(rounds = 0; ; ++rounds) {
    llk = sum1 = sum2 = 0;
    for (int i = 0; i < n; i++) {  // calculate likelihoods
      if ( pairs[i] < 0 ) {
	l = mpu[i].computeGP(prior11, prior12, prior22, prior1, prior2, chromosomeType, fracMale);
      }
      else {
	l = mpu[i].computePairGP(mpu[pairs[i]], prior11, prior12, prior22, prior1, prior2, chromosomeType);
      }
      if ( pairs[i] < i ) llk += log(l + MINLK);
      sum1 += ( mpu[i].GPs[1] / l );
      sum2 += ( mpu[i].GPs[2] / l );
    }

    // update rule
    if ( hwe ) {
      p = (sum1 + 2. *sum2) / (2.0*n);
      if ( ( fabs(p + prior1 - 1.) < 1e-6 ) || ( rounds == maxiter ) ) { // converged
	hweAF = p;
	hweLLK = llk;

	/*
	Acount = 0.5;
	ABcount = 1.0;
	for (int i = 0; i < n; i++) {
	  scale = -10.0 * (mpu[i].GLs[8] + mpu[i].GLs[0] - 2 * mpu[i].GLs[4]) + 6 * mpu[i].nbase;
	  delta = -10.0 * (mpu[i].GLs[8] - mpu[i].GLs[0]);
	  if (scale < 4) scale = 4;
	  if (scale < fabs(delta)) scale = fabs(delta);
	  nref = 0.5 * mpu[i].nbase * (1.0 + delta / (scale + MINLK));
	  Acount += (mpu[i].GPs[1] * nref);
	  ABcount += (mpu[i].GPs[1] * mpu[i].nbase);
	  }*/
	return;
      }
      else {
	q = 1.-p;
	prior11 = q * q;
	prior12 = p * q * 2.0;
	prior22 = p * p;
	prior1 = q;
	prior2 = p; 
      }
    }
    else {
      p1 = sum1/n;
      p2 = sum2/n;
      if ( ( fabs( p1 - prior12 ) + fabs( p2 - prior22 ) < 1e-6 ) || ( rounds == maxiter ) ) { // converged
	hwdAF1 = p1;
	hwdAF2 = p2;
	hwdLLK = llk;

	Acount = 0.5;
	ABcount = 1.0;
	for (int i = 0; i < n; i++) {
	  scale = -10.0 * (mpu[i].GLs[8] + mpu[i].GLs[0] - 2 * mpu[i].GLs[4]) + 6 * mpu[i].nbase;
	  delta = -10.0 * (mpu[i].GLs[8] - mpu[i].GLs[0]);
	  if (scale < 4) scale = 4;
	  if (scale < fabs(delta)) scale = fabs(delta);
	  nref = 0.5 * mpu[i].nbase * (1.0 + delta / (scale + MINLK));
	  Acount += (mpu[i].GPs[1] * nref);
	  ABcount += (mpu[i].GPs[1] * mpu[i].nbase);
	}
	return;
      }
      else {
	prior11 = 1.-p1-p2;
	prior12 = p1;
	prior22 = p2;
	prior1 = p2/(1.0-p1+MINLK);
	prior2 = 1.0-prior1;
      }
    }
  }
}

/*
double mpuMixupLikelihood::OptimizeFrequency()
{
  for(int i=0; i < n; ++i) { 
    mpu[i].computeIBDGL(); 
    //if ( pairs[i] < 0 ) {
    //  mpu[i].computeMixGL(0); 
    //}
    //else {
    //  mpu[i].computeIBDGL(); 
    //}
  }

  a = 0.00001; fa = f(a);
  b = 0.4; fb = f(b);
  c = 0.99999; fc = f(c);

  if ( fc == fb ) { 
    min = 0.99999999;
    fmin = fb; 
    return min; 
  }
  
  Brent(0.00001);

  if ( hasContam ) {
    //double org = min;
    //notice("%s:%d\t%lf\t%d", mpu[0].chrom.c_str(), mpu[0].pos, min, 0);
    double prev;
    int iter = 0;
    do {
      prev = min;
      for(int i=0; i < n; ++i) { 
	if ( pairs[i] < 0 ) {
	  mpu[i].computeMixGL(1.-min); 
	}
	else {
	  mpu[i].computeIBDGL(); 
	}
      } 

      a = 0.00001; fa = f(a);
      b = min; fb = f(min);
      c = 0.99999; fc = f(c);

      Brent(0.00001);
      ++iter;
    } while ( ( iter < 20 ) && ( fabs(min - prev) > 0.0001 ) );
    //notice("%s:%d\t%lf\t%lf\t%d", mpu[0].chrom.c_str(), mpu[0].pos, org, min, iter);
  }

  return min;
}
*/

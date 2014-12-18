#ifndef __MPU_FILE_H__
#define __MPU_FILE_H__

#include "pFile.h"
#include "fVcf.h"
#include "PhredHelper.h"
#include "BaseAsciiMap.h"
#include "Constant.h"

#include <vector>
#include <cstdlib>
#include <cmath>
extern phredConverter phredConv;
#define ENDPOS 1000000000
#define MINGL -30.0
#define MINLK 1e-30

class mpuFile {
 public:
  pFile tf;
  int nMarkers;
  std::string chrom;
  int pos;
  char ref;
  char alt;
  int nbase;
  int maxDP;
  bool mixMode;
  char* pBases;
  char* pBQs;
  char* pMQs;
  char* pCYs;
  std::vector<uint8_t> bases;
  std::vector<uint8_t> bQs;
  std::vector<uint8_t> mQs;
  std::vector<uint8_t> cycles;
  std::vector<uint8_t> strands;
  double GLs[9];
  double LKs[9];
  double GPs[3];
  //double pairGPs[9];
  char* line;
  double fmix;
  bool nobgzf;
  char sex;
  int maxCY;

 mpuFile(int maxDepth = 255, bool _nobgzf = false) : fmix(0), nobgzf(_nobgzf), sex(SEX_FEMALE), maxCY(70) {
    setMaxDepth(maxDepth);
  }

  void setMaxDepth(int maxDepth) {
    maxDP = maxDepth;
    bases.resize(maxDP);
    bQs.resize(maxDP);
    mQs.resize(maxDP);
    cycles.resize(maxDP);
    strands.resize(maxDP);
  }

  bool load(const char* filename, const char* region = NULL, bool printHeader = false) {
    tf.load(filename, region, printHeader, nobgzf);
    return true;
  }

  static int GenotypeIndex(int base1, int base2)
  {
    return base1 < base2 ? (base1) *(9 - base1) / 2 + (base2 - base1) :
      (base2) *(9 - base2) / 2 + (base1 - base2);
  }

  static int GenotypeIndexChar(char base1, char base2) {
    return GenotypeIndex(BaseAsciiMap::base2int[(int)base1],BaseAsciiMap::base2int[(int)base2]);
  }
  
  /*
  void computeMixGL(double af = 0) {
    // use bases, bQs to compute GLs
    int i,j,k,l;
    for(i=0; i < 10; ++i) {
      GLs[i] = 0;
      LKs[i] = 1;
    }

    if ( nbase > 0 ) {
      if ( ( af == 0 ) || ( fmix == 0 ) ) { // becomes independent of af contamination
	//fprintf(stderr,"%d",nbase);
	//for(i=0; i < nbase; ++i) {
	//  fprintf(stderr," %d:%d",bases[i],bQs[i]);
	//}

	for(i=0; i < nbase; ++i) {
	  for(j=0, l=0; j < 4; ++j) {
	    for(k=j; k < 4; ++k, ++l) {
	      if ( j == bases[i] ) {
		if ( k == bases[i] ) { // both alleles match
		  GLs[l] += phredConv.phred2LogMat[bQs[i]];		  
		}
		else {  // only one allele match
		  GLs[l] += phredConv.phred2HalfLogMat3[bQs[i]];
		}
	      }
	      else if ( k == bases[i] ) { // only one allele match
		GLs[l] += phredConv.phred2HalfLogMat3[bQs[i]];
	      }
	      else {
		GLs[l] += (-0.1*bQs[i]-phredConv.log3);
	      }
	    }
	  }
	}
	double maxGL = GLs[0];
	for(i=1; i < 10; ++i) { 
	  if ( maxGL < GLs[i] ) maxGL = GLs[i];
	  //fprintf(stderr," %.2lf",GLs[i]);
	}

	for(i=0; i < 10; ++i) { 
	  GLs[i] -= maxGL;
	  if ( GLs[i] < MINGL ) GLs[i] = MINGL;
	  LKs[i] = pow(10,GLs[i]);
	}
	//fprintf(stderr,"\n");
      }
      else { // Matthew's part to fill in - currently only fills in relevant likelihoods, not all 10
	int nref = BaseAsciiMap::base2int[(int)ref];
	int nalt = BaseAsciiMap::base2int[(int)alt];
	int geno, cgeno;
	int igenos[3];
	igenos[0] = GenotypeIndex(nref,nref);
	igenos[1] = GenotypeIndex(nref,nalt);
	igenos[2] = GenotypeIndex(nalt,nalt);

	for(i=0; i < 10; ++i) {
	  GLs[i] = MINGL; // ignore non-ref non-alt
	  LKs[i] = MINLK;
	}

	for(geno=0; geno<3; ++geno) {
	  //geno=0 (ref-ref); geno=1 (ref-alt het), geno=2 (alt-alt)
	  //double clike = 0;
	  double llike[3];
	  double cgfreqs[3];
	  for(cgeno=0; cgeno<3; ++cgeno) {
	    double lprob = 0;
	    //double cprob = 1;
	    double cgfreq = (cgeno==0) ? (1-af)*(1-af) : ((cgeno==1) ? 2*af*(1-af) : (af)*(af));
	    for(i=0; i<nbase; ++i) {
	      //cprob *= ((1.-fmix) * probBaseGivenGeno(geno,nref,nalt,bQs[i], bases[i]) + (fmix) * probBaseGivenGeno(cgeno,nref,nalt,bQs[i], bases[i]));
	      lprob += log((1.-fmix) * probBaseGivenGeno(geno,nref,nalt,bQs[i], bases[i]) + (fmix) * probBaseGivenGeno(cgeno,nref,nalt,bQs[i], bases[i]));
	    }
	    //printf("%d lprob[%d %d] = %lg\n",pos,geno,cgeno,lprob);
	    //clike += cprob * cgfreq;
	    llike[cgeno] = lprob;
	    cgfreqs[cgeno] = cgfreq;
	  }
	  // cprob[0] * cgfreq[0] + cprob[1] * cgfreq[1] + cprob[2] * cgfreq[2]
	  double maxLL = llike[0];
	  if ( maxLL < llike[1] ) maxLL = llike[1];
	  if ( maxLL < llike[2] ) maxLL = llike[2];
	  GLs[igenos[geno]] = (maxLL + log( exp(llike[0]-maxLL)*cgfreqs[0] + exp(llike[1]-maxLL)*cgfreqs[1] + exp(llike[2]-maxLL)*cgfreqs[2]))/log(10.);
	  //LKs[igenos[geno]] = pow(10,GLs[igenos[geno]]);
	  //printf("%d\tllike[%d (%d)] = %lg\n",pos,geno,igenos[geno],GLs[igenos[geno]]);
	  //LKs[igenos[geno]] = clike;
	  //GLs[igenos[geno]] = log10(clike);
	}

	double maxGL = GLs[igenos[0]];
	if ( maxGL < GLs[igenos[1]] ) maxGL = GLs[igenos[1]];
	if ( maxGL < GLs[igenos[2]] ) maxGL = GLs[igenos[2]];
	
	//notice("maxGL=%lf",maxGL);
	GLs[igenos[0]] -= maxGL;
	GLs[igenos[1]] -= maxGL;
	GLs[igenos[2]] -= maxGL;
	if ( GLs[igenos[0]] < MINGL ) GLs[igenos[0]] = MINGL;
	if ( GLs[igenos[1]] < MINGL ) GLs[igenos[1]] = MINGL;
	if ( GLs[igenos[2]] < MINGL ) GLs[igenos[2]] = MINGL;
	LKs[igenos[0]] = pow(10,GLs[igenos[0]]);
	LKs[igenos[1]] = pow(10,GLs[igenos[1]]);
	LKs[igenos[2]] = pow(10,GLs[igenos[2]]);
	//printf("** %d\tllike[%d (%d)] = %lg\n",pos,geno,igenos[1],GLs[igenos[1]]);
      }
    }
  }
  */

  double computeGP(double prior11, double prior12, double prior22, double prior1, double prior2, int chromosomeType, double fracMale = 0.5) {
    if ( fmix == 0 ) {
      GPs[0] = prior11 * LKs[0];
      GPs[1] = prior12 * LKs[4];
      GPs[2] = prior22 * LKs[8];
    }
    else {
      switch(chromosomeType) {
      case CT_MITO: case CT_AUTOSOME:
	GPs[0] = prior11 * ( prior11 * LKs[0] + prior12 * LKs[1] + prior22 * LKs[2] );
	GPs[1] = prior12 * ( prior11 * LKs[3] + prior12 * LKs[4] + prior22 * LKs[5] );
	GPs[2] = prior22 * ( prior11 * LKs[6] + prior12 * LKs[7] + prior22 * LKs[8] );
	break;
      case CT_CHRY:
	if ( sex == SEX_MALE ) { // Assume 50/50 mixture  REF/. p REF/REF p2 REF/ALT ALT/. ALT/ALT ALT/ALT
	  GPs[0] = prior1 * ( ( fracMale * prior1 + (1.-fracMale) ) * LKs[0] + ( fracMale * prior2 ) * LKs[2] );
	  GPs[1] = 0;
	  GPs[2] = prior2 * ( ( fracMale * prior1 ) * LKs[6] + ( fracMale * prior2 + (1.-fracMale) ) * LKs[8] );
	}
	else {
	  GPs[0] = prior11;
	  GPs[1] = prior12;
	  GPs[2] = prior22;
	}
	break;
      case CT_CHRX:
	if ( sex == SEX_MALE ) {
	  GPs[0] = prior1 * ( ( fracMale * prior1 + (1.-fracMale) * prior11 ) * LKs[0] +    // REF/REF
			      + ( fracMale * prior1 + (1.-fracMale) * prior12 ) * LKs[1] +  // REF/HET
			      + ( fracMale * prior1 + (1.-fracMale) * prior22 ) * LKs[2] ); // REF/ALT
	  GPs[1] = 0;
	  GPs[2] = prior2 * ( ( fracMale * prior1 + (1.-fracMale) * prior11 ) * LKs[6] +    // ALT/REF
			      + ( fracMale * prior1 + (1.-fracMale) * prior12 ) * LKs[7] +  // ALT/HET
			      + ( fracMale * prior1 + (1.-fracMale) * prior22 ) * LKs[8] ); // ALT/ALT	  
	}
	else {
	  GPs[0] = prior11 * ( ( fracMale * prior1 + (1.-fracMale) * prior11 ) * LKs[0] +   // REF/REF
			       ( fracMale * prior1 + (1.-fracMale) * prior12 ) * LKs[1] +   // REF/HET
			       ( fracMale * prior1 + (1.-fracMale) * prior22 ) * LKs[2] );  // REF/ALT
	  GPs[1] = prior12 * ( ( fracMale * prior1 + (1.-fracMale) * prior11 ) * LKs[3] +   // HET/REF
			       ( fracMale * prior1 + (1.-fracMale) * prior12 ) * LKs[4] +   // HET/HET
			       ( fracMale * prior1 + (1.-fracMale) * prior22 ) * LKs[5] );  // HET/ALT
	  GPs[2] = prior22 * ( ( fracMale * prior1 + (1.-fracMale) * prior11 ) * LKs[6] +   // ALT/REF
			       ( fracMale * prior1 + (1.-fracMale) * prior12 ) * LKs[7] +   // ALT/HET
			       ( fracMale * prior1 + (1.-fracMale) * prior22 ) * LKs[8] ); // ALT/ALT
	}
      }
    }
    return GPs[0]+GPs[1]+GPs[2];
  }

  double computePairGP(mpuFile& pair, double prior11, double prior12, double prior22, double prior1, double prior2, int chromosomeType) {
    if ( fmix == 0 ) {
      GPs[0] = prior11 * LKs[0];
      GPs[1] = prior12 * LKs[4];
      GPs[2] = prior22 * LKs[8];
    }
    else {
      double* lks2 = pair.LKs;
      switch(chromosomeType) {
      case CT_MITO: case CT_AUTOSOME:
	GPs[0] = prior11 * ( prior11 * LKs[0] * lks2[0] + prior12 * LKs[1] * lks2[3] + prior22 * LKs[2] * lks2[6] );
	GPs[1] = prior12 * ( prior11 * LKs[3] * lks2[1] + prior12 * LKs[4] * lks2[4] + prior22 * LKs[5] * lks2[7] );
	GPs[2] = prior22 * ( prior11 * LKs[6] * lks2[2] + prior12 * LKs[7] * lks2[5] + prior22 * LKs[8] * lks2[8] );
	break;
      case CT_CHRY:
	if ( sex == SEX_MALE ) {
	  if ( pair.sex == SEX_MALE ) {
	    GPs[0] = prior1 * ( prior1 * LKs[0] * lks2[0] +  prior2 * LKs[2] * lks2[6] );  // REF/ALT
	    GPs[1] = 0;
	    GPs[2] = prior2 * ( prior1 * LKs[6] * lks2[2] +  prior2 * LKs[8] * lks2[8] );  // REF/ALT
	  }
	  else { 
	    GPs[0] = prior1 * LKs[0] * lks2[0]; // REF/.
	    GPs[1] = 0;
	    GPs[2] = prior2 * LKs[8] * lks2[8]; // ALT/.
	  }
	}
	else {
	  GPs[0] = prior1 * LKs[0] * lks2[0];   // ./REF
	  GPs[1] = 0;                           
	  GPs[2] = prior2 * LKs[8] * lks2[8];   // ./ALT
	}
	break;
      case CT_CHRX:
	if ( sex == SEX_MALE ) {
	  if ( pair.sex == SEX_MALE ) {
	    GPs[0] = prior1 * ( prior1 * LKs[0] * lks2[0] +  prior2 * LKs[2] * lks2[6] );  // REF/ALT
	    GPs[1] = 0;
	    GPs[2] = prior2 * ( prior1 * LKs[6] * lks2[2] +  prior2 * LKs[8] * lks2[8] );  // REF/ALT
	  }
	  else {
	    GPs[0] = prior1 * ( prior11 * LKs[0] * lks2[0] + prior12 * LKs[1] * lks2[3] + prior22 * LKs[2] * lks2[6] );
	    GPs[1] = 0;
	    GPs[2] = prior2 * ( prior11 * LKs[6] * lks2[2] + prior12 * LKs[7] * lks2[5] + prior22 * LKs[8] * lks2[8] );
	  }
	}
	else if ( pair.sex == SEX_MALE ) {
	  GPs[0] = prior11 * ( prior1 * LKs[0] * lks2[0] +  prior2 * LKs[2] * lks2[6] );  
	  GPs[1] = prior12 * ( prior1 * LKs[3] * lks2[1] +  prior2 * LKs[5] * lks2[7] );  
	  GPs[2] = prior22 * ( prior1 * LKs[6] * lks2[2] +  prior2 * LKs[8] * lks2[8] );  
	}
	else {
	  GPs[0] = prior11 * ( prior11 * LKs[0] * lks2[0] + prior12 * LKs[1] * lks2[3] + prior22 * LKs[2] * lks2[6] );
	  GPs[1] = prior12 * ( prior11 * LKs[3] * lks2[1] + prior12 * LKs[4] * lks2[4] + prior22 * LKs[5] * lks2[7] );
	  GPs[2] = prior22 * ( prior11 * LKs[6] * lks2[2] + prior12 * LKs[7] * lks2[5] + prior22 * LKs[8] * lks2[8] );
	}
      }
    }
    return GPs[0]+GPs[1]+GPs[2];
  }

  void computeIBDGL() {
    // use bases, bQs to compute GLs
    int i; //,j,k,l;
    for(i=0; i < 9; ++i) {
      GLs[i] = 0;
      LKs[i] = 1;
      //GPs[i] = 0;
    }

    if ( nbase > 0 ) {
      if ( fmix == 0 ) { // becomes independent of af contamination
	int nref = BaseAsciiMap::base2int[(int)ref];
	int nalt = BaseAsciiMap::base2int[(int)alt];

	for(i=0; i < nbase; ++i) {
	  if ( bases[i] == nref ) {
	    GLs[0] += phredConv.phred2LogMat[bQs[i]];
	    GLs[3] += phredConv.phred2HalfLogMat3[bQs[i]];
	    GLs[6] += (-0.1*bQs[i]-phredConv.log3);
	  }
	  else if ( bases[i] == nalt ) {
	    GLs[0] += (-0.1*bQs[i]-phredConv.log3);
	    GLs[3] += phredConv.phred2HalfLogMat3[bQs[i]];
	    GLs[6] += phredConv.phred2LogMat[bQs[i]];
	  }
	  else {
	    GLs[0] += (-0.1*bQs[i]-phredConv.log3);
	    GLs[3] += (-0.1*bQs[i]-phredConv.log3);
	    GLs[6] += (-0.1*bQs[i]-phredConv.log3);
	  }
	}
	GLs[2] = GLs[1] = GLs[0];
	GLs[5] = GLs[4] = GLs[3];
	GLs[8] = GLs[7] = GLs[6];

	double maxGL = GLs[0];
	if ( maxGL < GLs[3] ) maxGL = GLs[3];
	if ( maxGL < GLs[6] ) maxGL = GLs[6];

	for(i=0; i < 9; ++i) { 
	  GLs[i] -= maxGL;
	  if ( GLs[i] < MINGL ) GLs[i] = MINGL;
	  LKs[i] = pow(10,GLs[i]);
	}
      }
      else {
	int nref = BaseAsciiMap::base2int[(int)ref];
	int nalt = BaseAsciiMap::base2int[(int)alt];
	int geno, cgeno;

	for(i=0; i < 9; ++i) {
	  GLs[i] = MINGL;
	  LKs[i] = MINLK;
	}

	for(geno=0; geno<3; ++geno) {
	  for(cgeno=0; cgeno<3; ++cgeno) {
	    double lprob = 0;
	    for(i=0; i<nbase; ++i) {
	      lprob += log((1.-fmix) * probBaseGivenGeno(geno,nref,nalt,bQs[i], bases[i]) + (fmix) * probBaseGivenGeno(cgeno,nref,nalt,bQs[i], bases[i]));
	    }
	    GLs[geno*3+cgeno] = lprob/log(10);
	  }
	}

	double maxGL = GLs[0];
	for(i=1; i < 9; ++i) { 
	  if ( maxGL < GLs[i] ) maxGL = GLs[i];
	}
	for(i=0; i < 9; ++i) { 
	  GLs[i] -= maxGL;
	  if ( GLs[i] < MINGL ) GLs[i] = MINGL;
	  LKs[i] = pow(10,GLs[i]);
	}
      }
    }
    //notice("%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf",GLs[0],GLs[1],GLs[2],GLs[3],GLs[4],GLs[5],GLs[6],GLs[7],GLs[8]);
  }

  double probBaseGivenGeno(int geno, uint8_t nref, uint8_t nalt, uint8_t qual, uint8_t base) {
    //geno=0 (ref-ref); geno=1 (ref-alt het), geno=2 (alt-alt)
    if ((base != nref) && (base != nalt)) {
      //return 2/3.0 * phredConv.phred2Err[qual];
      return phredConv.phred2Err[qual]/3.0; // hyun
    }
    if (geno==1) {
      //return 1/2.0 + phredConv.phred2Err[qual]/3.0;
      return 1/2.0 - phredConv.phred2Err[qual]/3.0; // hyun
    }
    if ((geno==0 && base==nref) || (geno==2 && base==nalt)) {
        return phredConv.phred2Mat[qual];
    } else {
        return phredConv.phred2Err[qual]/3.0;
    }
  }

  bool advanceTo(const char* chr, int bp, char refB, char altB ) {
    //notice("foo %s %d %c %c",chr,bp,refB,altB);
    if ( ( chrom.compare(chr) != 0 ) && ( tf.getType() == 2 ) ) {
      char buf[256];
      sprintf(buf,"%s:%d",chr,bp);
      tf.updateRegion(buf);
      chrom = chr;
      pos = 0;
    }

    pMQs = pCYs = NULL;

    char *pch, *nch;
    int i,j,l;
    while( ( ( bp > pos ) && (line = (char*)tf.getLine()) != NULL ) ) {
      pch = line;
      for(i=0, j=0; pch != NULL; ++i) {
	nch = strchr(pch,'\t');
	switch(i) {
	case 0:
	  chrom.assign( pch, nch - pch );
	  break;
	case 1:
	  pos = atoi(pch);
	  break;
	case 2:
	  ref = pch[0];
	  break;
	case 3:
	  nbase = atoi(pch);
	  break;
	case 4:
	  pBases = pch;
	  break;
	case 5:
	  pBQs = pch;
	  break;
	case 6:
	  pMQs = pch;
	  break;
	case 7:
	  pCYs = pch;
	  break;
	default:
	  error("Unrecognized index %d, pch = %s",i,pch);
	}
	if ( nch == NULL ) { pch = NULL; }
	else { pch = nch+1; }
      }
      //notice("%d\t%s\t%d\t%c\t%d",bp,chrom.c_str(),pos,ref,nbase);
    }

    if ( line == NULL ) {
      pos = ENDPOS;
      nbase = 0;
      computeIBDGL();
      return false;
    }
    else if ( bp < pos ) {
      nbase = 0;
      computeIBDGL();
      return false;
    }
    else if ( bp == pos ) {
      //if ( refB != ref ) error("ref %c do not match with %c",refB);
      // parse 
      if ( nbase > maxDP ) nbase = maxDP-1;
      uint8_t nref = BaseAsciiMap::base2int[(int)ref];
      for(i=0, j=0; i < nbase; ++j) {
	switch( pBases[j] ) {
	case '$': // do nothing
	  break;
	case '^': // ignore the next character
	  ++j;
	  break;
	case '-': case '+':
	  l = atoi(&pBases[j+1]);
	  j += (l+(l<10 ? 1 : (l < 100 ? 2 : 3)));
	  break;
	case '.': case ',':
	  bases[i] = nref;
	  bQs[i] = pBQs[i] - 33;
	  if ( pMQs != NULL ) mQs[i] = pMQs[i] - 33;
	  if ( pCYs != NULL ) { cycles[i] = (uint8_t)strtol(pCYs, &pCYs, 10); ++pCYs; if ( maxCY > cycles[i] ) maxCY = cycles[i]; }
	  strands[i] = ( pBases[j] == '.' );
	  ++i;
	  break;
	case 'a': case 'A':
	  bases[i] = 0;
	  bQs[i] = pBQs[i] - 33;
	  if ( pMQs != NULL ) mQs[i] = pMQs[i] - 33;
	  if ( pCYs != NULL ) { cycles[i] = (uint8_t)strtol(pCYs, &pCYs, 10); ++pCYs; if ( maxCY > cycles[i] ) maxCY = cycles[i]; }
	  strands[i] = ( pBases[j] == 'A' );
	  ++i;
	  break;
	case 'c': case 'C':
	  bases[i] = 1;
	  bQs[i] = pBQs[i] - 33;
	  if ( pMQs != NULL ) mQs[i] = pMQs[i] - 33;
	  if ( pCYs != NULL ) { cycles[i] = (uint8_t)strtol(pCYs, &pCYs, 10); ++pCYs; if ( maxCY > cycles[i] ) maxCY = cycles[i]; }
	  strands[i] = ( pBases[j] == 'C' );
	  ++i;
	  break;
	case 'g': case 'G':
	  bases[i] = 2;
	  bQs[i] = pBQs[i] - 33;
	  if ( pMQs != NULL ) mQs[i] = pMQs[i] - 33;
	  if ( pCYs != NULL ) { cycles[i] = (uint8_t)strtol(pCYs, &pCYs, 10); ++pCYs; if ( maxCY > cycles[i] ) maxCY = cycles[i]; }
	  strands[i] = ( pBases[j] == 'G' );
	  ++i;
	  break;
	case 't': case 'T':
	  bases[i] = 3;
	  bQs[i] = pBQs[i] - 33;
	  if ( pMQs != NULL ) mQs[i] = pMQs[i] - 33;
	  if ( pCYs != NULL ) { cycles[i] = (uint8_t)strtol(pCYs, &pCYs, 10); ++pCYs; if ( maxCY > cycles[i] ) maxCY = cycles[i]; }
	  strands[i] = ( pBases[j] == 'T' );
	  ++i;
	  break;
	case '*': case 'N': case 'n':
	  bases[i] = 4;
	  bQs[i] = pBQs[i] - 33;
	  if ( pMQs != NULL ) mQs[i] = pMQs[i] - 33;
	  if ( pCYs != NULL ) { cycles[i] = (uint8_t)strtol(pCYs, &pCYs, 10); ++pCYs; if ( maxCY > cycles[i] ) maxCY = cycles[i]; }
	  strands[i] = ( pBases[j] == 'N' );
	  ++i;
	  break;
	default:
	  error("Cannot recognize %c in the pileups at i=%d, j= %d, nbase=%d, pBases=%s\n",pBases[j],i,j,nbase,pBases);
	}
      }
      alt = altB;
      computeIBDGL();
      return true;
    }
    else {
      error("ERROR: The base positions are not ordered?");
      return false;
    }
  }
};

#endif

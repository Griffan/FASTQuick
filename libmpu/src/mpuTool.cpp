#include "Parameters.h"
#include "Error.h"
#include "fVcf.h"
#include "wFile.h"
#include "BaseAsciiMap.h"
#include "mpuFile.h"
#include "mpuMixupLikelihoods.h"
#include "Constant.h"
#include "mpuVerify.h"
#include "MathGold.h"
#include "MathGenMin.h"

#include <iostream>
#include <string>
#include <map>
#include <ctime>
#include <climits>
#include <vector>
#include <getopt.h>
#include <utility>

void makeGeno(int g1, int g2, std::string& s) {
  char buf[255];
  sprintf(buf,"%d/%d",g1,g2);
  s = buf;
}

void int2str(int x, std::string& s) {
  char buf[255];
  sprintf(buf,"%d",x);
  s = buf;
}

std::pair<int,int> GetBestGenotype(mpuFile& mpu) {
  int best = 0;
  double sum = mpu.GPs[0];
  for(int i=1; i < 3; ++i) {
    sum += mpu.GPs[i];
    if ( mpu.GPs[i] > mpu.GPs[best] ) 
      best = i;
  }
  double error = (sum - mpu.GPs[best])/(sum+1e-30);
  int qual = (error < 3.16e-26) ? 255 : int(-log(error)/log(10)*10+0.5);
  return std::pair<int,int>(best,qual);
}

void GetGenotypeLikelihoods(mpuFile& mpu, double* priors, char* buf) {
  if ( mpu.fmix == 0 ) {
    sprintf(buf,":%.2lf,%.2lf,%.2lf",mpu.GLs[0],mpu.GLs[4],mpu.GLs[8]);
  }
  else {
    double l10 = log(10.);
    if ( priors[1] == 0 ) {
      if ( priors[2] == 0 ) { // monomorphic
	sprintf(buf,":%.2lf,%.2lf,%.2lf",mpu.GLs[0],mpu.GLs[4],mpu.GLs[8]);
      }
      else {
	double l0 = mpu.GPs[0] / priors[0];
	double l2 = mpu.GPs[2] / priors[2];
	if ( l0 > l2 ) {
	  l2 = log(l2/l0)/l10;
	  if ( l2 < MINGL ) l2 = MINGL;
	  sprintf(buf,":%.2lf,%.2lf,%.2lf",0.0,MINGL,l2);
	}
	else {
	  l0 = log(l0/l2)/l10;
	  if ( l0 < MINGL ) l0 = MINGL;
	  sprintf(buf,":%.2lf,%.2lf,%.2lf",l0,MINGL,0.0);
	}
      }
    }
    else {
      double l0 = mpu.GPs[0] / priors[0];
      double l1 = mpu.GPs[1] / priors[1];
      double l2 = mpu.GPs[2] / priors[2];
      if ( l0 > l1 ) {
	if ( l0 > l2 ) {
	  l1 = log(l1/l0)/l10;
	  if ( l1 < MINGL ) l1 = MINGL;
	  l2 = log(l2/l0)/l10;
	  if ( l2 < MINGL ) l2 = MINGL;
	  sprintf(buf,":%.2lf,%.2lf,%.2lf",0.0,l1,l2);
	}
	else {
	  l1 = log(l1/l2)/l10;
	  if ( l1 < MINGL ) l1 = MINGL;
	  l0 = log(l0/l2)/l10;
	  if ( l0 < MINGL ) l0 = MINGL;
	  sprintf(buf,":%.2lf,%.2lf,%.2lf",l0,l1,0.0);
	}
      }
      else if ( l1 > l2 ) {
	l2 = log(l2/l1)/l10;
	if ( l2 < MINGL ) l2 = MINGL;
	l0 = log(l0/l1)/l10;
	if ( l0 < MINGL ) l0 = MINGL;
	sprintf(buf,":%.2lf,%.2lf,%.2lf",l0,0.0,l2);
      }
      else {
	l1 = log(l1/l2)/l10;
	if ( l1 < MINGL ) l1 = MINGL;
	l0 = log(l0/l2)/l10;
	if ( l0 < MINGL ) l0 = MINGL;
	sprintf(buf,":%.2lf,%.2lf,%.2lf",l0,l1,0.0);
      }
    }
  }
}

void ReportDate(wFile& output) {
  time_t systemTime;
  time(&systemTime);
  
  tm * digestedTime;
  digestedTime = gmtime(&systemTime);
  
  output.printf("##filedate=%04d%02d%02d\n", digestedTime->tm_year + 1900,
		digestedTime->tm_mon + 1,
		digestedTime->tm_mday);
}

double ReportMixupGenotypes(mpuMixupLikelihood& lk, mpuFile * mpu, int n, const char* chrom, int position, int refAllele, int al1, int al2, std::string& info, std::string& genotypes, std::string& filter, bool reportPairGL) {
  info.clear();
  genotypes.clear();

  double priors[2][3];
  char buf[65536];

  int al3, al4;
  for(al3=0; al3 < 4; ++al3) { if ( ( al1 != al3 ) && ( al2 != al3 ) ) break; }
  for(al4=al3+1; al4 < 4; ++al4) { if ( ( al1 != al4 ) && ( al2 != al4 ) ) break; }

  lk.GetFemalePriors(priors[0], 1-lk.hweAF);
  lk.GetMalePriors(priors[1], 1-lk.hweAF);

  int label1 = al1 == refAllele ? 0 : 1;
  int label2 = al2 == refAllele ? 0 : al1 == al2 ? label1 : label1 + 1;
  
  //int genoRR = 0, genoR1 = 0, genoR2 = 0;
  
  if (label2 == 2) {
    error("Sorry, currently multi-allelic variants are not supported");    
    //genoRR = mpuLikelihood::GenotypeIndex(refAllele, refAllele);
    //genoR1 = mpuLikelihood::GenotypeIndex(refAllele, al1);
    //genoR2 = mpuLikelihood::GenotypeIndex(refAllele, al2);
   }

  std::string label11[2], label12[2], label22[2];

  if (lk.chromosomeType == CT_CHRY)
    label11[0] = label12[0] = label22[0] = ".";
  else if (lk.chromosomeType == CT_MITO) {
    int2str(label1, label11[0]);
    label12[0] = ".";
    int2str(label2, label22[0]);
  }
  else { /* CT_AUTO, CT_CHRX */
    makeGeno(label1, label1, label11[0]);
    makeGeno(label1, label2, label12[0]);
    makeGeno(label2, label2, label22[0]);
  }
  
  if (lk.chromosomeType != CT_AUTOSOME) {
    int2str(label1, label11[1]);
    label12[1] = ".";
    int2str(label2, label22[1]);
  }
  else {
    makeGeno(label1, label1, label11[1]);
    makeGeno(label1, label2, label12[1]);
    makeGeno(label2, label2, label22[1]);
  }

  // available INFO field for filtering
  // NS : Number of 

  int ns         = 0; // # of samples with coverage
  int dp         = 0; // total depth 
  int gc[3]      = {0, 0, 0}; // allele count for each genotypes
  int gd[3]      = {0, 0, 0}; // depth for each genotypes
  int gm[3]      = {0, 0, 0}; // MQ for each genotypes
  int sa[4] = {1,1,1,1};   // strand by allele contingency table
  int ta[4] = {1,1,1,1};   // within-5bp of tail by allele contingency table
  int sh[4] = {1,1,1,1};   // strand by allele contingency table only at het sites
  int th[4] = {1,1,1,1};   // within-5bp of tail by allele contingency only at het sites
  int q1 = 60, q2 = 1800, c1 = 100, c2 = 5000, a1 = 1, qa = 30, ca = 50, nr = 0, na = 0, ne = 0, m1 = 0, m2 = 0, ma1 = 0, ma2 = 0, lm = 0;
  // What statistics do we need?
  // Per individual
  // # of REF/ALT/OTH bases (for allele balance)
  double sgpb[9] = {0,0,0,0,0,0,0,0,0};
  int nMaxAlt = 0;
  double fMaxAlt = 0;
  int i, j, b, bq, s, c, m; // base, BQ, strand, cycle, MQ, offset  

  // iterate over n individuals
  for (i = 0; i < n; i++) {
    // find out the best guss genotypes
    int sex = mpu[i].sex == SEX_MALE ? 1 : 0;
    mpuFile& mp = mpu[i]; // assume that OptimizeFrequency() is already called
  
    // Report on the best genotype for the current SNP model
    std::pair<int,int> bestQual = GetBestGenotype(mp);
    //int quality = GetBestRatio(mp.LKs, priors[sex]);

    int quality = bestQual.second;
    int bestS = bestQual.first;
    std::string & labelS = bestS == 0 ? label11[sex] : bestS == 1 ? label12[sex] : label22[sex];
    bool nocall = (labelS[0] == '.');
    int  depth  = mp.nbase;
    dp += depth;
      
    sprintf(buf,"\t%s:%d:%d",labelS.c_str(), depth, nocall ? 0 : quality);
    genotypes += buf;
      
    if (labelS[0] != '.') {
      if ( depth > 0 ) {
	ns++;
	gd[bestS] += depth;
      }
      ++gc[bestS];

      if (label2 < 2) {
	GetGenotypeLikelihoods(mp, priors[sex], buf);
	genotypes += buf;
	if ( reportPairGL ) {
	  sprintf(buf,":%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf",mp.GLs[0],mp.GLs[1],mp.GLs[2],mp.GLs[3],mp.GLs[4],mp.GLs[5],mp.GLs[6],mp.GLs[7],mp.GLs[8]);
	  genotypes += buf;
	}
      }
      else {
	error("Currently multi-allelic variants are not supported");
      }
    }

    double gpb[9] = {0,0,0,0,0,0,0,0,0};
    for(j=0; j < (int)mp.nbase; ++j) {
      b = mp.bases[j];
      bq = mp.bQs[j];
      m = mp.mQs[j];
      c = mp.cycles[j];
      s = mp.strands[j];
      ++m1;    // # of MQs
      m2 += m; // sumsqMQ
      gm[bestS] += m;

      // b=0:ref, b=1:alt, b=3: others
      if ( b == al1 ) b = 0;
      else if ( b == al2 ) b = 1;
      else b = 3;

      if ( m <= 13 ) { ++lm; };   // low MQ
      if ( b <= 1 ) {
	++sa[s*2+b];              
	++ta[((c < 5) || (mp.maxCY-c < 5))*2+b];
	if ( bestS == 1 ) {
	  ++sh[s*2+b];              
	  ++th[((c < 5) || (mp.maxCY-c < 5))*2+b];
	}
	q1 += bq;
	q2 += (bq*bq);
	c1 += c;
	c2 += (c*c);
	a1 += b;
	qa += bq*b;
	ca += c*b;
	if ( b == 1 ) {
	  ++ma1;
	  ma2 += m;
	}
	++gpb[3*bestS+b];
	if ( bq >= 13 ) {
	  if ( b == 0 ) ++nr;
	  else ++na;
	}
      }
      else {
	++gpb[3*bestS+2];
	if ( bq >= 13 ) {
	  ++ne;
	}
      }
    }

    //if ( nMaxRef < gpb[0]+gpb[3]+gpb[6] ) nMaxRef = gpb[0]+gpb[3]+gpb[6];
    if ( nMaxAlt < gpb[1]+gpb[4]+gpb[7] ) {
      nMaxAlt = gpb[1]+gpb[4]+gpb[7];
    }
    double fAlt = (double)(gpb[1]+gpb[4]+gpb[7])/(double)(gpb[0]+gpb[1]+gpb[3]+gpb[4]+gpb[6]+gpb[7]+0.5);
    if ( fAlt > fMaxAlt ) fMaxAlt = fAlt;
    for(j=0; j < 9; ++j) { sgpb[j] += gpb[j]; }
  }
  
  //if (label1 == 0 && label2 == 0)
  //return;

  double qual = (lk.hweLLK - lk.nullLLK)/log(10);
  if ( qual < 0 ) qual = 0;
  
  int sasum = sa[0]+sa[1]+sa[2]+sa[3];
  int shsum = sh[0]+sh[1]+sh[2]+sh[3];
  double STR = (double)(sa[0]*sa[3]-sa[1]*sa[2])/sqrt((double)(sa[0]+sa[1])*(sa[2]+sa[3])*(sa[0]+sa[2])*(sa[1]+sa[3]));
  double STZ = STR * sqrt((double)sasum);
  double SHR = (double)(sh[0]*sh[3]-sh[1]*sh[2])/sqrt((double)(sh[0]+sh[1]+1e-6)*(sh[2]+sh[3]+1e-6)*(sh[0]+sh[2]+1e-6)*(sh[1]+sh[3]+1e-6));
  double SHZ = SHR * sqrt((double)shsum);
  //double BQR = (double)(qa - q1*a1/(double)sasum)/sqrt((double)(q2-q1*q1/(double)sasum)*(a1-a1*a1/(double)sasum));
  //double BQZ = BQR * sqrt((double)sasum);
  //double CBR = (double)(ca - c1*a1/(double)sasum)/sqrt((double)(c2-c1*c1/(double)sasum)*(a1-a1*a1/(double)sasum));
  //double CBZ = CBR * sqrt((double)sasum);
  double TBR = (double)(ta[1]*ta[2]-ta[0]*ta[3])/sqrt((double)(ta[0]+ta[1])*(ta[2]+ta[3])*(ta[0]+ta[2])*(ta[1]+ta[3]));
  double TBZ = TBR * sqrt((double)sasum);
  double THR = (double)(th[1]*th[2]-th[0]*th[3])/sqrt((double)(th[0]+th[1]+1e-6)*(th[2]+th[3]+1e-6)*(th[0]+th[2]+1e-6)*(th[1]+th[3]+1e-6));
  double THZ = THR * sqrt((double)shsum);
  double FOB = (double)(ne+.1)/(double)(na+ne+.1);
  int MMQ = (int)(m2/(m1+1e-6));
  int AMQ = (int)((ma2)/(ma1+1e-6));
  double LMQ = lm/(m1+1e-6);
  double ABE = (sgpb[3*1+0]+0.5)/(sgpb[3*1+0]+sgpb[3*1+1]+1.);
  double NRO = (sgpb[5]/(gc[1]+.01) + sgpb[8]/(gc[2]+.01) + 1)/(sgpb[2]/(gc[0]+.01) + 1);
  double ABL = lk.Acount/lk.ABcount;
  double FIC = 1.-lk.hwdAF1/(2*lk.hweAF*(1-lk.hweAF));
  double SLRT = (FIC > 0 ? 2 : -2 )*(lk.hwdLLK > lk.hweLLK ? lk.hwdLLK-lk.hweLLK : 0);
  double QD = qual / (gd[1]*0.5 + gd[2] + 0.5);
  double LQR = (double)(m1-nr-na-ne)/(double)m1;
  //sprintf(buf,"DP=%d;NS=%d", dp, ns); info += buf;
  //sprintf(buf,";AN=%d", ac[0] + ac[1] + ac[2] + ac[3]); info += buf;

  sprintf(buf,"DP=%d;NS=%d;AN=%d;AC=%d;AF=%.6lf;HWDAF=%.5lf,%.5lf;FIC=%.3lf;SLRT=%.3lf;ABL=%.4lf;STR=%.3lf;STZ=%.3lf;TBR=%.3lf;TBZ=%.3lf;FOB=%.3lf;MMQ=%d;AMQ=%d;LMQ=%.3lf;ABE=%.3lf;NRO=%.3lf;NMA=%d;FMA=%.3lf;QD=%.3lf;SHR=%.3lf;SHZ=%.3lf;THR=%.3lf;THZ=%.3lf;LQR=%.3lf",
	  dp, ns, 2*(gc[0]+gc[1]+gc[2]), gc[1]+gc[2]+gc[2], 
	  lk.hweAF, lk.hwdAF1, lk.hwdAF2, FIC, SLRT, ABL,
	  STR, STZ, TBR, TBZ,
	  FOB, MMQ, AMQ, LMQ, ABE, NRO, nMaxAlt, fMaxAlt,
	  QD, SHR, SHZ, THR, THZ, LQR//gc[1], gc[2],
	  //gd[0]/(gc[0]+1e-3),(gd[1]+gd[2])/(gc[1]+gc[2]+1e-3),
	  //int(gm[0]/(gd[0]+1e-3)),int((gm[1]+gm[2])/(gd[1]+gd[2]+1e-3))
	  );//,gpb[0]/ns[0],gpb[1]/ns[0],gpb[2]/ns[0],gpb[3]/ns[1],gpb[4]/ns[1],gpb[5]/ns[1],gpb[6]/ns[2],gpb[7]/ns[2],gpb[8]/ns[2]);
  info = buf;
  if ( gc[1]+gc[2] == 0 ) qual = 0;

  buf[0] = '\0';
  sprintf(buf,"%s%s%s%s%s%s%s%s%s%s%s%s%s",
	  qual < 1 ? ";q1" : "",
	  qual < 5 ? ";q5" : "",
	  ns + ns < n ? ";ns50" : "",
	  FIC < -0.1 ? ";fic-10" : "",
	  SLRT < -5 ? ";slrt-5" : "",
	  ABE > 0.67 ? ";AB67" : "",
	  (fabs(STR/0.1+STZ/20.) > 1) ? ";STR10STZ20" : "",
	  (THR/0.2 + THR/5. > 1) ? ";THR20THZ5" : "",
	  ((MMQ < 20) || (AMQ < 20)) ? ";mq20"  : "",
	  LMQ > 0.20 ? ";LMQ20" : "",
	  LQR > 0.20 ? ";LQR20" : "",
	  fMaxAlt < 0.3 ? ";fma30" : "",
	  nMaxAlt < 2 ? ";nma1" : "");
  if ( buf[0] == '\0' ) filter = "PASS";
  else filter = buf+1;
  return qual;
}

/*
void ReportSNP(wFile& baseCalls, mpuLikelihood & lk, int n, const char* chrom, int position, int refBase, int allele1, int allele2, double posterior ) {
  if (allele2 == refBase) {
    int swap = allele1;
    allele1 = allele2;
    allele2 = swap;
  }
  
  char alleles[] = { 'A', 'C', 'G', 'T', 0 };
  
  mpuFile* mpu = lk.mpu;
    
  // #CHROM\tPOS\tID
  baseCalls.printf("%s\t%d\t.\t", chrom, position+1);
  
  // REF
  int nalleles = 1;
  baseCalls.printf("%c\t", alleles[refBase]);
  
  // ALT
  if (allele1 != refBase)
    baseCalls.printf("%c", alleles[allele1]), nalleles++;

  if (allele2 != refBase && allele2 != allele1)
    baseCalls.printf("%s%c", nalleles > 1 ? "," : "", alleles[allele2]), nalleles++;

  if (nalleles == 1)
    baseCalls.printf(".");
  baseCalls.printf("\t");
  baseCalls.printf("100\t");
  baseCalls.printf("PASS\t");

  // Find best frequency
  lk.SetAlleles(allele1, allele2);
  lk.OptimizeFrequency();

  std::string info, genotypes;
  ReportGenotypes(lk, mpu, n, chrom, position, refBase, allele1, allele2, info, genotypes);
  
  baseCalls.printf("%s\t",info.c_str());
  baseCalls.printf("GT:GD:GQ:GL");
  baseCalls.printf("%s\n",genotypes.c_str());
}
*/

void ReportMixupSNP(wFile& baseCalls, mpuMixupLikelihood & lk, int n, const char* chrom, int position, int refBase, int allele1, int allele2, double posterior, bool printMono = false, bool reportPairGL = false) {
  if (allele2 == refBase) {
    int swap = allele1;
    allele1 = allele2;
    allele2 = swap;
  }
  
  char alleles[] = { 'A', 'C', 'G', 'T', 0 };
  
  mpuFile* mpu = lk.mpu;

  // Find best frequency
  lk.SetAlleles(allele1, allele2);
  lk.OptimizeFrequencyEM(false); // HWD
  lk.OptimizeFrequencyEM(true);  // HWE

  std::string info, genotypes, filter;
  double qual = ReportMixupGenotypes(lk, mpu, n, chrom, position, refBase, allele1, allele2, info, genotypes, filter, reportPairGL); 

  if ( ( qual == 0 ) && ( !printMono ) ) return;
    
  // #CHROM\tPOS\tID
  baseCalls.printf("%s\t%d\t.\t", chrom, position+1);
  
  // REF
  int nalleles = 1;
  baseCalls.printf("%c\t", alleles[refBase]);
  
  // ALT
  if (allele1 != refBase)
    baseCalls.printf("%c", alleles[allele1]), nalleles++;

  if (allele2 != refBase && allele2 != allele1)
    baseCalls.printf("%s%c", nalleles > 1 ? "," : "", alleles[allele2]), nalleles++;

  if (nalleles == 1)
    baseCalls.printf(".");


  baseCalls.printf("\t");
  baseCalls.printf("%.2lf\t",qual);
  baseCalls.printf("%s\t",filter.c_str());
  
  baseCalls.printf("%s\t",info.c_str());
  //baseCalls.printf("GT:GM:GD:GQ:GL:LM");
  baseCalls.printf("GT:DP:GQ:GL");
  if ( reportPairGL ) { baseCalls.printf(":LM"); }
  baseCalls.printf("%s\n",genotypes.c_str());
}

/*
int runGenotype(int argc, char ** argv) {
   printf("mpuGenotype -- Extract VCF based on pilups files\n");
   printf("(c) 2013 Hyun Min Kang, Matthew Flickinger, Goncalo Abecasis \n\n");

   std::string pedfile;
   std::string datfile;
   std::string colmpu("MPU");
   std::string colcon;
   std::string invcf;
   std::string callfile;
   std::string region;
   std::string rule;
   double minMAF = 0;
   double minCallRate = 0;
   bool ignoreFilter = false;
   bool verbose = false;
   bool nobgzf = false;
   double defaultContam = 0;
   double minContam = 0;
   int unit = 10000L;

   ParameterList pl;

   std::string xLabel("X"), yLabel("Y"), mitoLabel("MT");
   int    xStart = 2699520, xStop = 154931044;

   BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Pedigree File")
         LONG_STRINGPARAMETER("ped", &pedfile)
         LONG_STRINGPARAMETER("dat", &datfile)
         LONG_STRINGPARAMETER("col-mpu",&colmpu)
         LONG_STRINGPARAMETER("col-con",&colcon)
         LONG_DOUBLEPARAMETER("default-contam",&defaultContam)
         LONG_DOUBLEPARAMETER("thres-contam",&minContam)
      LONG_PARAMETER_GROUP("Input VCF")
         LONG_STRINGPARAMETER("invcf", &invcf)
         LONG_STRINGPARAMETER("region",&region)
         LONG_STRINGPARAMETER("rule",&rule)
         LONG_INTPARAMETER("unit",&unit)
         LONG_DOUBLEPARAMETER("minMAF",&minMAF)
         LONG_DOUBLEPARAMETER("minCallRate",&minCallRate)
         LONG_PARAMETER("ignoreFilter",&ignoreFilter)
         LONG_PARAMETER("notabix",&nobgzf)
      LONG_PARAMETER_GROUP("Chromosome Labels")
         LONG_STRINGPARAMETER("xChr", &xLabel)
         LONG_STRINGPARAMETER("yChr", &yLabel)
         LONG_STRINGPARAMETER("mito", &mitoLabel)
         LONG_INTPARAMETER("xStart", &xStart)
         LONG_INTPARAMETER("xStop", &xStop)
      LONG_PARAMETER_GROUP("Output")
         LONG_STRINGPARAMETER("out", &callfile)
         LONG_PARAMETER("verbose", &verbose)
   END_LONG_PARAMETERS();

   pl.Add(new LongParameters("Available Options", longParameters));
   pl.Read(argc, argv);
   pl.Status();

   // sanity check of input arguments
   if ( invcf.empty() || pedfile.empty() || callfile.empty()  ) {
     error("--invcf, --ped, --out are required parameters");
   }

   time_t t;
   time(&t);

   printf("Analysis started on %s\n", ctime(&t));
   fflush(stdout);

   fPed pedf(pedfile.c_str(), datfile.c_str());
   int n = pedf.ninds;

   mpuFile *mpu = new mpuFile[n];
   //int firstMpu = n;
   for (int i = n - 1; i >= 0; i--) {
     mpu[i].nobgzf = nobgzf;
     if (!mpu[i].load(pedf.getPheno(colmpu.c_str(),i).c_str()))
	 error("Failed to open MPU file [%s]", mpu[i].load(pedf.getPheno(colmpu.c_str(),i).c_str()));
     //else
     //firstMpu = i;
   }

   bool hasContam = false;
   if ( colcon.empty() ) {
     for(int i=0; i < n; ++i) { mpu[i].fmix = defaultContam; }
     if ( defaultContam > 0 ) hasContam = true;
   }
   else {
     for(int i=0; i < n; ++i) { 
       mpu[i].fmix = atof(pedf.getPheno(colcon.c_str(),i).c_str());
       if ( mpu[i].fmix < minContam ) {
	 mpu[i].fmix = defaultContam;
       }
       else {
	 hasContam = true;
       }
     }
     if ( defaultContam > 0 ) hasContam = true;
   }

   printf("Generating VCFs for files ...\n");
   for (int i = 0; i < n; i++)
     printf("%s\t%lf\n", pedf.getPheno(colmpu.c_str(),i).c_str(), mpu[i].fmix);
   printf("\n");

   wFile baseCalls(callfile.c_str());
   baseCalls.printf("##fileformat=VCFv4.1\n");
   ReportDate(baseCalls);
   baseCalls.printf("##source=mpuGenotype\n");
   baseCalls.printf("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth at Site\">\n");
   baseCalls.printf("##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root Mean Squared Mapping Quality\">\n");
   baseCalls.printf("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Coverage\">\n");
   baseCalls.printf("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of Alleles in Samples with Coverage\">\n");
   baseCalls.printf("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Alternate Allele Counts in Samples with Coverage\">\n");
   baseCalls.printf("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Alternate Allele Frequencies\">\n");
   baseCalls.printf("##INFO=<ID=AB,Number=1,Type=Float,Description=\"Allele Balance in Heterozygotes\">\n");
   baseCalls.printf("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
   baseCalls.printf("##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Genotype Likelihoods for Genotypes 0/0,0/1,1/1\">\n");
   baseCalls.printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
   for (int i = 0; i < n; i++) {
     baseCalls.printf("\t%s", pedf.getPheno("IND_ID",i).c_str());
   }
   baseCalls.printf("\n");

   std::string buffer;

   std::string curChrom;
   int curPos;
   int chromosomeType = CT_AUTOSOME;
   mpuLikelihood lkGeno(n, mpu, hasContam);

   fVcf tvcf;
   tvcf.load(invcf.c_str(), region.c_str(), "GT", rule.c_str(), !ignoreFilter, NULL);
   for (int i = 0; i < n; i++)
     lkGeno.sexes[i] = (atoi(pedf.getPheno("SEX",i).c_str()) == SEX_MALE ? SEX_MALE : SEX_FEMALE);

   int M, i, m;
   double af, maf;
   for(M=0; tvcf.readMarkers(unit); ) {
     M += tvcf.nMarkers;
     fprintf(stderr,"Processing %d markers..\n", M);
     for(i=0, m=0; i < tvcf.nMarkers; ++i) { // for each marker
       af = tvcf.alleleFreq(i);
       maf = af > 0.5 ? 1-af : af;
       if ( ( maf >= minMAF ) && 
	    ( tvcf.callRate(i) >= minCallRate ) )
	 { 
	   uint8_t refBase = BaseAsciiMap::base2int[(int)tvcf.refs[i][0]];
	   uint8_t altBase = BaseAsciiMap::base2int[(int)tvcf.alts[i][0]];
	   uint8_t allele1 = refBase;
	   uint8_t allele2 = altBase;
	   curChrom = tvcf.chroms[i];
	   curPos = tvcf.pos1s[i]-1;

	   if (curChrom == xLabel) chromosomeType = CT_CHRX;
	   if (curChrom == yLabel) chromosomeType = CT_CHRY;
	   if (curChrom == mitoLabel) chromosomeType = CT_MITO;

	   int     totalDepth = 0, nSamplesCovered = 0;
	   //double  rmsMapQuality = 0.0;
	   
	   for (int j = 0; j < n; ++j) {
	     int depth = 0;
	     if ( mpu[j].advanceTo(tvcf.chroms[i].c_str(), tvcf.pos1s[i], tvcf.refs[i][0], tvcf.alts[i][0]) ) {
	       depth = mpu[j].nbase;
	       totalDepth += depth;
	       ++nSamplesCovered;
	     }
	   }
	   
	   //lkGeno.position = curPos;
	   lkGeno.chromosomeType = chromosomeType != CT_CHRX ? chromosomeType :
	     curPos >= xStart && curPos <= xStop ? CT_CHRX : CT_AUTOSOME;
	   
	   //fprintf(stderr,"%s",curChrom.c_str());
	   int qual = 100;
	   ReportSNP(baseCalls, lkGeno, n, curChrom.c_str(), curPos, refBase, allele1, allele2, phredConv.phred2Mat[qual]);	 
	   ++m;
	 }
     }
   }
   tvcf.close();
   fprintf(stderr,"Finished processing total of %d markers\n", M);

   baseCalls.close();

   delete [] mpu;

   time(&t);
   printf("\nAnalysis completed on %s\n", ctime(&t));
   fflush(stdout);

   return 0;
}

int runMixGenotype(int argc, char ** argv) {
   printf("mpuMixup -- Extract VCF based on pilups files\n");
   printf("(c) 2013 Hyun Min Kang, Matthew Flickinger, Goncalo Abecasis \n\n");

   std::string pedfile;
   std::string datfile;
   std::string colmpu("MPU");
   std::string colcon;
   std::string colmixid("MIXID");
   std::string invcf;
   std::string callfile;
   std::string region;
   std::string rule;
   double minMAF = 0;
   double minCallRate = 0;
   double defaultContam = 0;
   double minContam = 0;
   bool ignoreFilter = false;
   bool verbose = false;
   bool nobgzf = false;
   int unit = 10000L;

   ParameterList pl;

   //std::string xLabel("X"), yLabel("Y"), mitoLabel("MT");
   //int    xStart = 2699520, xStop = 154931044;

   BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Pedigree File")
         LONG_STRINGPARAMETER("ped", &pedfile)
         LONG_STRINGPARAMETER("dat", &datfile)
         LONG_STRINGPARAMETER("col-mpu",&colmpu)
         LONG_STRINGPARAMETER("col-con",&colcon)
         LONG_STRINGPARAMETER("col-mixid",&colmixid)
         LONG_DOUBLEPARAMETER("default-contam",&defaultContam)
         LONG_DOUBLEPARAMETER("thres-contam",&minContam)
      LONG_PARAMETER_GROUP("Input VCF")
         LONG_STRINGPARAMETER("invcf", &invcf)
         LONG_STRINGPARAMETER("region",&region)
         LONG_STRINGPARAMETER("rule",&rule)
         LONG_INTPARAMETER("unit",&unit)
         LONG_DOUBLEPARAMETER("minMAF",&minMAF)
         LONG_DOUBLEPARAMETER("minCallRate",&minCallRate)
         LONG_PARAMETER("ignoreFilter",&ignoreFilter)
         LONG_PARAMETER("notabix",&nobgzf)
      LONG_PARAMETER_GROUP("Chromosome Labels")
         LONG_STRINGPARAMETER("xChr", &fVcf::xLabel)
         LONG_STRINGPARAMETER("yChr", &fVcf::yLabel)
         LONG_STRINGPARAMETER("mito", &fVcf::mitoLabel)
         LONG_INTPARAMETER("xStart", &fVcf::xStart)
         LONG_INTPARAMETER("xStop", &fVcf::xStop)
      LONG_PARAMETER_GROUP("Output")
         LONG_STRINGPARAMETER("out", &callfile)
         LONG_PARAMETER("verbose", &verbose)
   END_LONG_PARAMETERS();

   pl.Add(new LongParameters("Available Options", longParameters));
   pl.Read(argc, argv);
   pl.Status();

   // sanity check of input arguments
   if ( invcf.empty() || pedfile.empty() || callfile.empty()  ) {
     error("--invcf, --ped, --out are required parameters");
   }

   time_t t;
   time(&t);

   printf("Analysis started on %s\n", ctime(&t));
   fflush(stdout);

   fPed pedf(pedfile.c_str(), datfile.c_str());
   int n = pedf.ninds;

   mpuFile *mpu = new mpuFile[n];
   //int firstMpu = n;
   int* pairs = new int[n];

   for (int i = n - 1; i >= 0; i--) {
     mpu[i].nobgzf = nobgzf;
     if (!mpu[i].load(pedf.getPheno(colmpu.c_str(),i).c_str()))
	 error("Failed to open MPU file [%s]", mpu[i].load(pedf.getPheno(colmpu.c_str(),i).c_str()));
     //else
     //firstMpu = i;

     const std::string& mixId = pedf.getPheno(colmixid.c_str(),i);
     if ( pedf.iinds.find(mixId) != pedf.iinds.end() ) {
       pairs[i] = pedf.iinds[mixId];
     }
     else {
       pairs[i] = -1;
     }
   }

   bool hasContam = false;
   if ( colcon.empty() ) {
     for(int i=0; i < n; ++i) { mpu[i].fmix = defaultContam; }
     if ( defaultContam > 0 ) hasContam = true;
   }
   else {
     for(int i=0; i < n; ++i) { 
       mpu[i].fmix = atof(pedf.getPheno(colcon.c_str(),i).c_str());
       if ( mpu[i].fmix < minContam ) {
	 mpu[i].fmix = defaultContam;
       }
       else {
	 hasContam = true;
       }
     }
     if ( defaultContam > 0 ) hasContam = true;
   }

   printf("Generating VCFs for files ...\n");
   for (int i = 0; i < n; i++)
     printf("%s\t%lf\t%s\n", pedf.getPheno(colmpu.c_str(),i).c_str(), mpu[i].fmix,pairs[i] < 0 ? "." : pedf.getPheno("IND_ID",pairs[i]).c_str());
   printf("\n");

   wFile baseCalls(callfile.c_str());
   baseCalls.printf("##fileformat=VCFv4.1\n");
   ReportDate(baseCalls);
   baseCalls.printf("##source=mpuGenotype\n");
   baseCalls.printf("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth at Site\">\n");
   baseCalls.printf("##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root Mean Squared Mapping Quality\">\n");
   baseCalls.printf("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Coverage\">\n");
   baseCalls.printf("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of Alleles in Samples with Coverage\">\n");
   baseCalls.printf("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Alternate Allele Counts in Samples with Coverage\">\n");
   baseCalls.printf("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Alternate Allele Frequencies\">\n");
   baseCalls.printf("##INFO=<ID=AB,Number=1,Type=Float,Description=\"Allele Balance in Heterozygotes\">\n");
   baseCalls.printf("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
   baseCalls.printf("##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Genotype Likelihoods for Genotypes 0/0,0/1,1/1\">\n");
   baseCalls.printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
   for (int i = 0; i < n; i++) {
     baseCalls.printf("\t%s", pedf.getPheno("IND_ID",i).c_str());
   }
   baseCalls.printf("\n");

   std::string buffer;

   std::string curChrom;
   int curPos;
   int chromosomeType = CT_AUTOSOME;
   mpuMixupLikelihood lkGeno(n, mpu, pairs, hasContam);

   fVcf tvcf;
   tvcf.load(invcf.c_str(), region.c_str(), "GT", rule.c_str(), !ignoreFilter, NULL);
   for (int i = 0; i < n; i++)
     lkGeno.sexes[i] = (atoi(pedf.getPheno("SEX",i).c_str()) == SEX_MALE ? SEX_MALE : SEX_FEMALE);

   int M, i;
   double af, maf;
   for(M=0; tvcf.readMarkers(unit); ) {
     M += tvcf.nMarkers;
     fprintf(stderr,"Processing %d markers..\n", M);
     for(i=0; i < tvcf.nMarkers; ++i) { // for each marker
       af = tvcf.alleleFreq(i);
       maf = af > 0.5 ? 1-af : af;
       if ( ( maf >= minMAF ) && 
	    ( tvcf.callRate(i) >= minCallRate ) )
	 { 
	   uint8_t refBase = BaseAsciiMap::base2int[(int)tvcf.refs[i][0]];
	   uint8_t altBase = BaseAsciiMap::base2int[(int)tvcf.alts[i][0]];
	   uint8_t allele1 = refBase;
	   uint8_t allele2 = altBase;
	   curChrom = tvcf.chroms[i];
	   curPos = tvcf.pos1s[i]-1;

	   lkGeno.chromosomeType = chromosomeType = fVcf::getChromosomeType(curChrom.c_str(), curPos);
	   //if (curChrom == xLabel) chromosomeType = CT_CHRX;
	   //if (curChrom == yLabel) chromosomeType = CT_CHRY;
	   //if (curChrom == mitoLabel) chromosomeType = CT_MITO;

	   int     totalDepth = 0, nSamplesCovered = 0;
	   //double  rmsMapQuality = 0.0;
	   
	   for (int j = 0; j < n; ++j) {
	     int depth = 0;
	     if ( mpu[j].advanceTo(tvcf.chroms[i].c_str(), tvcf.pos1s[i], tvcf.refs[i][0], tvcf.alts[i][0]) ) {
	       depth = mpu[j].nbase;
	       totalDepth += depth;
	       ++nSamplesCovered;
	     }
	   }
	   
	   //lkGeno.chromosomeType = chromosomeType != CT_CHRX ? chromosomeType :
	   //  curPos >= xStart && curPos <= xStop ? CT_CHRX : CT_AUTOSOME;
	   
	   int qual = 100;
	   ReportMixupSNP(baseCalls, lkGeno, n, curChrom.c_str(), curPos, refBase, allele1, allele2, phredConv.phred2Mat[qual]);	 
	 }
     }
   }
   tvcf.close();
   fprintf(stderr,"Finished processing total of %d markers\n", M);

   baseCalls.close();

   delete [] mpu;
   delete [] pairs;

   time(&t);
   printf("\nAnalysis completed on %s\n", ctime(&t));
   fflush(stdout);

   return 0;
}
*/

int runVerify(int argc, char** argv) {
  printf("mpuVerify 1.0.0 -- verify identity and purity of sequence reads\n"
        "(c) 2013 Hyun Min Kang, Goo Jun, and Goncalo Abecasis\n\n");

  mpuVerifyArgs args;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Files")
    LONG_STRINGPARAMETER("vcf",&args.sVcfFile)
    LONG_STRINGPARAMETER("mpu",&args.sMpuFile)
    LONG_STRINGPARAMETER("subset",&args.sSubsetInds)
    LONG_STRINGPARAMETER("smID",&args.sSMID)
    LONG_STRINGPARAMETER("scanPair",&args.sScanPair)

    LONG_PARAMETER_GROUP("VCF analysis options")
    LONG_DOUBLEPARAMETER("genoError",&args.genoError)
    LONG_DOUBLEPARAMETER("minAF",&args.minAF)
    LONG_DOUBLEPARAMETER("minCallRate",&args.minCallRate)

    LONG_PARAMETER_GROUP("Individuals to compare with chip data")
    EXCLUSIVE_PARAMETER("site",&args.bSiteOnly)
    EXCLUSIVE_PARAMETER("self",&args.bSelfOnly)
    EXCLUSIVE_PARAMETER("best",&args.bFindBest)

    LONG_PARAMETER_GROUP("Chip-free optimization options")
    EXCLUSIVE_PARAMETER("free-none",&args.bFreeNone)
    EXCLUSIVE_PARAMETER("free-mix",&args.bFreeMixOnly)
    EXCLUSIVE_PARAMETER("free-refBias",&args.bFreeRefBiasOnly)
    EXCLUSIVE_PARAMETER("free-full",&args.bFreeFull)

    LONG_PARAMETER_GROUP("With-chip optimization options")
    EXCLUSIVE_PARAMETER("chip-none",&args.bChipNone)
    EXCLUSIVE_PARAMETER("chip-mix",&args.bChipMixOnly)
    EXCLUSIVE_PARAMETER("chip-refBias",&args.bChipRefBiasOnly)
    EXCLUSIVE_PARAMETER("chip-full",&args.bChipFull)

    LONG_PARAMETER_GROUP("MPU analysis options")
    LONG_PARAMETER("precise",&args.bPrecise)
    LONG_INTPARAMETER("minMapQ",&args.minMapQ)
    LONG_INTPARAMETER("maxDepth",&args.maxDepth)
    LONG_INTPARAMETER("minQ",&args.minQ)
    LONG_INTPARAMETER("maxQ",&args.maxQ)
    LONG_DOUBLEPARAMETER("grid",&args.grid)

    LONG_PARAMETER_GROUP("Modeling Reference Bias")
    LONG_DOUBLEPARAMETER("refRef",&args.pRefRef)
    LONG_DOUBLEPARAMETER("refHet",&args.pRefHet)
    LONG_DOUBLEPARAMETER("refAlt",&args.pRefAlt)

    LONG_PARAMETER_GROUP("Output options")
    LONG_STRINGPARAMETER("out",&args.sOutFile)
    LONG_PARAMETER("verbose",&args.bVerbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options",longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // check the validity of input files
  if ( args.sVcfFile.empty() ) {
    error("--vcf [vcf file] required");
  }

  if ( args.sMpuFile.empty() ) {
    error("--mpu [mpu file] is required");
  }

  if ( args.sOutFile.empty() ) {
    error("--out [output prefix] is required");
  }

  if ( args.sSMID.empty() ) {
    error("--smID [sample ID] is required");
  }

  if ( ! ( args.bSiteOnly || args.bSelfOnly || args.bFindBest ) ) {
    warning("--self option was autotomatically turned on by default. Specify --best option if you wanted to check across all possible samples in the VCF");
    args.bSelfOnly = true;
  }

  if ( !( args.bFreeNone || args.bFreeMixOnly || args.bFreeRefBiasOnly || args.bFreeFull ) ) {
    warning("--free-mix option is automatically turned on by default");
    args.bFreeMixOnly = true;
  }

  if ( !( args.bChipNone || args.bChipMixOnly || args.bChipRefBiasOnly || args.bChipFull ) ) {
    warning("--chip-mix option is automatically turned on by default");
    args.bChipMixOnly = true;
  }

  if ( ( args.maxDepth > 20 ) && ( !args.bPrecise ) ) {
    warning("--precise option is not turned on at --maxDepth %d : may be prone to precision errors",args.maxDepth);
  }

  if ( ( args.bChipRefBiasOnly ) && ( !args.bSelfOnly ) ) {
    error("--self must be set for --chip-refBias to work. Skipping..");
  }

  if ( ( !args.sScanPair.empty() ) && ( ( !args.bFreeNone ) || ( !args.bFindBest ) || ( !args.bChipMixOnly ) ) ) {
    error("--scanPair option is only compatible with --best, --free-none, and --chip-mix");
  }

  // check timestamp
  time_t t;
  time(&t);
  notice("Analysis started on %s",ctime(&t));

  // load arguments
  mpuVerify vbid(&args);

  // load input VCF and BAM files
  notice("Opening Input Files");
  vbid.loadFiles(args.sMpuFile.c_str(), args.sVcfFile.c_str());

  int selfIndex = -1;
  if ( !args.sScanPair.empty() ) { // if scanPair is on
    std::string indID1(args.sScanPair.c_str());
    for(int i=0; i < (int)vbid.pGenotypes->indids.size(); ++i) {
      if ( vbid.pGenotypes->indids[i] == indID1 ) {
	selfIndex = i;
      }
    }
    if ( selfIndex < 0 ) {
      error("Cannot find individual %s from the VCF file",indID1.c_str());
    }
  }

  // Check which genotype-free method is used
  if ( args.bFreeNone ) {  // if no genotype-free mode is tested. skip it
    // do nothing for genotype-free estimation
    notice("Skipping chip-free estimation of sample mixture");
  }
  else if ( args.bFreeMixOnly ) { // only mixture is estimated.
    // genotype-free method
    notice("Performing chip-free estimation of sample mixture at fixed reference bias parameters (%lf, %lf, %lf)",args.pRefRef,args.pRefHet,args.pRefAlt);

    mpuVerify::mixLLK mix(&vbid);
    mix.OptimizeLLK();
    notice("Optimal per-sample fMix = %lf, LLK0 = %lf, LLK1 = %lf\n",mix.fMix,mix.llk0,mix.llk1);
    vbid.mixOut.llk0 = mix.llk0;
    vbid.mixOut.llk1 = mix.llk1;
    vbid.mixOut.fMix = mix.fMix;
  }
  else if ( args.bFreeRefBiasOnly ) {
    notice("Performing chip-free estimation of reference-bias without sample mixture");
    mpuVerify::refBiasMixLLKFunc myFunc(&vbid);
    AmoebaMinimizer myMinimizer;
    Vector startingPoint(2);
    startingPoint[0] = 0;      // pRefHet = 0.5
    startingPoint[1] = -4.595; // pRefAlt = 0.01
    myMinimizer.func = &myFunc;
    myMinimizer.Reset(2);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(1e-6);
    double pRefHet = mpuVerify::invLogit(myMinimizer.point[0]);
    double pRefAlt = mpuVerify::invLogit(myMinimizer.point[1]);
    notice("Reference Bias Estimated as ( Pr[refBase|HET] = %lf, Pr[refBase|ALT] = %lf) with LLK = %lf",pRefHet,pRefAlt,myMinimizer.fmin);
    vbid.mixOut.llk0 = myFunc.llk0;
    vbid.mixOut.llk1 = myFunc.llk1;
    vbid.mixOut.refHet = myFunc.pRefHet;
    vbid.mixOut.refAlt = myFunc.pRefAlt;
  }
  else if ( args.bFreeFull ) {
    notice("Performing chip-free estimation of reference-bias and sample mixture together");
    mpuVerify::fullMixLLKFunc myFunc(&vbid);
    AmoebaMinimizer myMinimizer;
    Vector startingPoint(3);
    startingPoint[0] = -3.91;  // start with fMix = 0.01
    startingPoint[1] = 0;      // pRefHet = 0.5
    startingPoint[2] = -4.595; // pRefAlt = 0.01
    myMinimizer.func = &myFunc;
    myMinimizer.Reset(3);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(1e-6);
    double fMix = mpuVerify::invLogit(myMinimizer.point[0]);
    if ( fMix > 0.5 ) 
      fMix = 1.-fMix;
    double pRefHet = mpuVerify::invLogit(myMinimizer.point[1]);
    double pRefAlt = mpuVerify::invLogit(myMinimizer.point[2]);
    notice("Optimal per-sample fMix = %lf\n",fMix);
    notice("Reference Bias Estimated as ( Pr[refBase|HET] = %lf, Pr[refBase|ALT] = %lf) with LLK = %lf",pRefHet,pRefAlt,myMinimizer.fmin);

    vbid.mixOut.llk0 = myFunc.llk0;
    vbid.mixOut.llk1 = myFunc.llk1;
    vbid.mixOut.fMix = myFunc.fMix;
    vbid.mixOut.refHet = myFunc.pRefHet;
    vbid.mixOut.refAlt = myFunc.pRefAlt;
  }
  notice("calculating depth distribution");  
  vbid.calculateDepthDistribution(args.maxDepth, vbid.mixOut);

  notice("finished calculating depth distribution");  

  int bestInd = -1;
  int selfInd = -1;

  if ( args.bChipNone ) {
    // do nothing
    notice("Skipping with-chip estimation of sample mixture");
  }
  else if ( args.bChipMixOnly ) {
    notice("Performing with-chip estimation of sample mixture at fixed reference bias parameter (%lf, %lf, %lf)",args.pRefRef,args.pRefHet,args.pRefAlt);

    if ( args.sScanPair.empty() ) {
      double maxIBD = -1;
      mpuVerify::ibdLLK ibd(&vbid);
      for(int i=0; i < (int)vbid.pGenotypes->indids.size(); ++i) {
	double fIBD = ibd.OptimizeLLK(i);
	notice("Comparing with individual %s.. Optimal fIBD = %lf, LLK0 = %lf, LLK1 = %lf",vbid.pGenotypes->indids[i].c_str(),fIBD, ibd.llk0, ibd.llk1);
	if ( maxIBD < fIBD ) {
	  bestInd = i;
	  vbid.bestOut.llk0 = ibd.llk0;
	  vbid.bestOut.llk1 = ibd.llk1;
	  vbid.bestOut.fMix = 1-ibd.fIBD;
	  maxIBD = ibd.fIBD;
	}
	  
	if ( vbid.pPile->sSM == vbid.pGenotypes->indids[i] ) {
	  selfInd = i;
	  vbid.selfOut.llk0 = ibd.llk0;
	  vbid.selfOut.llk1 = ibd.llk1;
	  vbid.selfOut.fMix = 1-ibd.fIBD;
	}
      }
	
      if ( bestInd >= 0 ) {
	notice("Best Matching Individual is %s with IBD = %lf",vbid.pGenotypes->indids[bestInd].c_str(),maxIBD);
	vbid.calculateDepthByGenotype(bestInd,vbid.bestOut);
      }
      
      if ( selfInd >= 0 ) {
	notice("Self Individual is %s with IBD = %lf",vbid.pGenotypes->indids[selfInd].c_str(),vbid.selfOut.fMix);
	vbid.calculateDepthByGenotype(selfInd,vbid.selfOut);
      }
    }
    else {
      double maxLLK = 1e9;
      mpuVerify::pairLLK ibdPair(&vbid);
      for(int i=0; i < (int)vbid.pGenotypes->indids.size(); ++i) {
	double fIBD = ibdPair.OptimizeLLK(selfIndex, i);
	notice("Comparing with individual %s.. Optimal fIBD = %lf, LLK0 = %lf, LLK1 = %lf",vbid.pGenotypes->indids[i].c_str(),fIBD, ibdPair.llk0, ibdPair.llk1);
	double fLLK = ibdPair.llk1;
	if ( fLLK < maxLLK ) {
	  bestInd = i;
	  vbid.bestOut.llk0 = ibdPair.llk0;
	  vbid.bestOut.llk1 = ibdPair.llk1;
	  vbid.bestOut.fMix = 1-ibdPair.fIBD;
	  maxLLK = fLLK;
	}
	
	if (vbid.pPile->sSM == vbid.pGenotypes->indids[i] ) {
	  selfInd = i;
	  vbid.selfOut.llk0 = ibdPair.llk0;
	  vbid.selfOut.llk1 = ibdPair.llk1;
	  vbid.selfOut.fMix = 1-ibdPair.fIBD;
	}
      }
      
      if ( bestInd >= 0 ) {
	notice("Best Matching Individual is %s with LLK = %lf",vbid.pGenotypes->indids[bestInd].c_str(),maxLLK);
	vbid.calculateDepthByGenotype(bestInd,vbid.bestOut);
      }
      
      if ( selfInd >= 0 ) {
	notice("Self Individual is %s with IBD = %lf",vbid.pGenotypes->indids[selfInd].c_str(),vbid.selfOut.fMix);
	vbid.calculateDepthByGenotype(selfInd,vbid.selfOut);
      }
    }
  }
  else if ( args.bChipRefBiasOnly ) {
    notice("Performing with-chip estimation of reference-bias without sample mixture");
    if ( args.bSelfOnly ) {
      mpuVerify::refBiasIbdLLKFunc myFunc(&vbid);
      AmoebaMinimizer myMinimizer;
      Vector startingPoint(2);
      startingPoint[0] = 0;      // pRefHet = 0.5
      startingPoint[1] = -4.595; // pRefAlt = 0.01
      myMinimizer.func = &myFunc;
      myMinimizer.Reset(2);
      myMinimizer.point = startingPoint;
      myMinimizer.Minimize(1e-6);
      double pRefHet = mpuVerify::invLogit(myMinimizer.point[0]);
      double pRefAlt = mpuVerify::invLogit(myMinimizer.point[1]);
      notice("Reference Bias Estimated as ( Pr[refBase|HET] = %lf, Pr[refBase|ALT] = %lf) with LLK = %lf",pRefHet,pRefAlt,myMinimizer.fmin);
      //vbid.setRefBiasParams(1.0, pRefHet, pRefAlt);
      
      vbid.selfOut.llk0 = myFunc.llk0;
      vbid.selfOut.llk1 = myFunc.llk1;
      vbid.selfOut.refHet = myFunc.pRefHet;
      vbid.selfOut.refAlt = myFunc.pRefAlt;
      vbid.calculateDepthByGenotype(0,vbid.selfOut);
    }
    else {
      warning("--self must be set for --chip-refBias to work. Skipping..");
    }
  }
  else if ( args.bChipFull ) {
    notice("Performing with-chip estimation of reference-bias and sample mixture together");
    double maxIBD = -1;
    
    for(int i=0; i < (int)vbid.pGenotypes->indids.size(); ++i) {
      mpuVerify::fullIbdLLKFunc myFunc(&vbid,i);
      AmoebaMinimizer myMinimizer;
      Vector startingPoint(3);
      startingPoint[0] = 3.91;  // start with fIBD = 0.99
      startingPoint[1] = 0;      // pRefHet = 0.5
      startingPoint[2] = -4.595; // pRefAlt = 0.01
      myMinimizer.func = &myFunc;
      
      myFunc.indIdx = i;
      myMinimizer.Reset(3);
      myMinimizer.point = startingPoint;
      myMinimizer.Minimize(1e-6);
      double fIBD = mpuVerify::invLogit(myMinimizer.point[0]);
      double pRefHet = mpuVerify::invLogit(myMinimizer.point[1]);
      double pRefAlt = mpuVerify::invLogit(myMinimizer.point[2]);
      
      notice("Comparing with individual %s.. Optimal fIBD = %lf, LLK0 = %lf, LLK1 = %l",vbid.pGenotypes->indids[i].c_str(), fIBD, myFunc.llk0, myFunc.llk1);
      notice("Reference Bias Estimated as ( Pr[refBase|HET] = %lf, Pr[refBase|ALT] = %lf ) with LLK = %lf",pRefHet,pRefAlt,myMinimizer.fmin);
      if ( maxIBD < fIBD ) {
	bestInd = i;
	maxIBD = fIBD;
	vbid.bestOut.llk0 = myFunc.llk0;
	vbid.bestOut.llk1 = myFunc.llk1;
	vbid.bestOut.fMix = 1.-myFunc.fIBD;
	vbid.bestOut.refHet = myFunc.pRefHet;
	vbid.bestOut.refAlt = myFunc.pRefAlt;
      }

      if (vbid.pPile->sSM == vbid.pGenotypes->indids[i]) {
	selfInd = i;
	vbid.selfOut.llk0 = myFunc.llk0;
	vbid.selfOut.llk1 = myFunc.llk1;
	vbid.selfOut.fMix = 1.-myFunc.fIBD;
	vbid.selfOut.refHet = myFunc.pRefHet;
	vbid.selfOut.refAlt = myFunc.pRefAlt;
	vbid.calculateDepthByGenotype(i, vbid.selfOut);
      }
    }
    if ( bestInd >= 0 ) {
      notice("Best Matching Individual is %s with IBD = %lf",vbid.pGenotypes->indids[bestInd].c_str(),maxIBD);
      vbid.calculateDepthByGenotype(bestInd, vbid.bestOut);
    }
    
    if ( selfInd >= 0 ) {
      notice("Self Individual is %s with IBD = %lf",vbid.pGenotypes->indids[selfInd].c_str(),vbid.selfOut.fMix);
      vbid.calculateDepthByGenotype(selfInd,vbid.selfOut);
    }
  }

  // PRINT OUTPUT FILE - ".selfSM"
  // [SEQ_ID]  : SAMPLE ID in the sequence file
  // [CHIP_ID] : SAMPLE ID in the chip file (NA if not available)
  // [#SNPS] : Number of markers evaluated
  // [#READS]   : Number of reads evaluated
  // [AVG_DP]   : Mean depth
  // [FREEMIX]  : Chip-free estimated alpha (% MIX in 0-1 scale), NA if unavailable
  // [FREELK1]  : Chip-free log-likelihood at estimated alpha
  // [FREELK0]  : Chip-free log-likelihood at 0% contamination
  // [CHIPIBD]  : With-chip estimated alpha (% MIX in 0-1 scale)
  // [CHIPLK1]  : With-chip log-likelihood at estimated alpha
  // [CHIPLK0]  : With-chip log-likelihood at 0% contamination
  // [DPREF]    : Depth at reference site in the chip
  // [RDPHET]   : Relative depth at HET site in the chip
  // [RDPALT]   : Relative depth at HOMALT site in the chip
  // [FREE_RF]  : Pr(Ref|Ref) site estimated without chip data
  // [FREE_RH]  : Pr(Ref|Het) site estimated without chip data
  // [FREE_RA]  : Pr(Ref|Alt) site estimated without chip data
  // [CHIP_RF]  : Pr(Ref|Ref) site estimated with chip data
  // [CHIP_RH]  : Pr(Ref|Het) site estimated with chip data
  // [CHIP_RA]  : Pr(Ref|Alt) site estimated with chip data
  // [DPREF]    : Depth at reference alleles
  // [RDPHET]   : Relative depth at heterozygous alleles
  // [RDPALT]   : Relative depth at hom-alt alleles

  std::string soutfile = std::string(args.sOutFile.c_str());
  std::string selfSMFN = soutfile + ".selfSM";
  std::string bestSMFN = soutfile + ".bestSM";
  std::string dpSMFN = soutfile + ".depthSM";

  fprintf(stderr,"foo %s\n",selfSMFN.c_str());
  fprintf(stderr,"right before calling ifopen(%s)\n",selfSMFN.c_str());

  wFile selfSMF(selfSMFN.c_str());
  wFile bestSMF;
  if ( args.bFindBest ) bestSMF.open(bestSMFN.c_str());

  wFile dpSMF(dpSMFN.c_str());

  dpSMF.printf("#RG\tDEPTH\t#SNPs\t%%SNPs\t%%CUMUL\n");
  int nCumMarkers = 0;
  for(int i=args.maxDepth; i >= 0; --i) {
    nCumMarkers += vbid.mixOut.depths[i];
    dpSMF.printf("ALL\t%d\t%d\t%.5lf\t%.5lf\n",i, vbid.mixOut.depths[i],(double) vbid.mixOut.depths[i]/(double)vbid.nMarkers,(double)nCumMarkers/(double)vbid.nMarkers);
  }
  dpSMF.close();

  const char* headers[] = {"#SEQ_ID","RG","CHIP_ID","#SNPS","#READS","AVG_DP","FREEMIX","FREELK1","FREELK0","FREE_RH","FREE_RA","CHIPMIX","CHIPLK1","CHIPLK0","CHIP_RH","CHIP_RA","DPREF","RDPHET","RDPALT"};
  int nheaders = sizeof(headers)/sizeof(headers[0]);

  for(int i=0; i < nheaders; ++i) { selfSMF.printf("%s%s",i>0 ? "\t" : "",headers[i]); }
  selfSMF.printf("\n");
  selfSMF.printf("%s\tALL",vbid.pPile->sSM.c_str());
  selfSMF.printf("\t%s",selfInd >= 0 ? vbid.pGenotypes->indids[selfInd].c_str() : "NA");
  selfSMF.printf("\t%d\t%d\t%.2lf",vbid.nMarkers,vbid.mixOut.numReads[0],(double)vbid.mixOut.numReads[0]/(double)vbid.nMarkers);
  if ( args.bFreeNone ) { selfSMF.printf("\tNA\tNA\tNA\tNA\tNA"); }
  else if ( args.bFreeMixOnly ) { selfSMF.printf("\t%.5lf\t%.2lf\t%.2lf\tNA\tNA",vbid.mixOut.fMix,vbid.mixOut.llk1,vbid.mixOut.llk0); }
  else if ( args.bFreeRefBiasOnly ) { selfSMF.printf("\tNA\t%.2lf\t%.2lf\t%.5lf\t%.5lf",vbid.mixOut.llk1,vbid.mixOut.llk0,vbid.mixOut.refHet,vbid.mixOut.refAlt); }
  else if ( args.bFreeFull ) { selfSMF.printf("\t%.5lf\t%.2lf\t%.2lf\t%.5lf\t%.5lf",vbid.mixOut.fMix,vbid.mixOut.llk1,vbid.mixOut.llk0,vbid.mixOut.refHet,vbid.mixOut.refAlt); }
  else { error("Invalid option in handling bFree"); }

  if ( args.bChipNone || bestInd < 0 ) { selfSMF.printf("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"); }
  else if ( args.bChipMixOnly ) { selfSMF.printf("\t%.5lf\t%.2lf\t%.2lf\tNA\tNA\t%.3lf\t%.4lf\t%.4lf",vbid.selfOut.fMix,vbid.selfOut.llk1,vbid.selfOut.llk0,(double)vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[1], (double)vbid.selfOut.numReads[2]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[2], (double)vbid.selfOut.numReads[3]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[3]); }
  else if ( args.bChipMixOnly ) { selfSMF.printf("\tNA\t%.2lf\t%.2lf\t%.5lf\t%.5lf\t%.3lf\t%.4lf\t%.4lf",vbid.selfOut.llk1, vbid.selfOut.llk0, vbid.selfOut.refHet, vbid.selfOut.refAlt, (double)vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[1], (double)vbid.selfOut.numReads[2]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[2], (double)vbid.selfOut.numReads[3]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[3]); }
  else if ( args.bChipFull ) { selfSMF.printf("\t%.5lf\t%.2lf\t%.2lf\t%.5lf\t%.5lf\t%.3lf\t%.4lf\t%.4lf", vbid.selfOut.fMix, vbid.selfOut.llk1, vbid.selfOut.llk0, vbid.selfOut.refHet, vbid.selfOut.refAlt, (double)vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[1], (double)vbid.selfOut.numReads[2]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[2], (double)vbid.selfOut.numReads[3]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[3]); }
  else { error("Invalid option in handling bChip"); }
  selfSMF.printf("\n");
  selfSMF.close();

  if ( args.bFindBest ) {
    for(int i=0; i < nheaders; ++i) { bestSMF.printf("%s%s",i>0 ? "\t" : "",headers[i]); }
    bestSMF.printf("\n");
    bestSMF.printf("%s\tALL",vbid.pPile->sSM.c_str());
    bestSMF.printf("\t%s",bestInd >= 0 ? vbid.pGenotypes->indids[bestInd].c_str() : "NA");
    bestSMF.printf("\t%d\t%d\t%.2lf",vbid.nMarkers,vbid.mixOut.numReads[0],(double)vbid.mixOut.numReads[0]/(double)vbid.nMarkers);
    if ( args.bFreeNone ) { bestSMF.printf("\tNA\tNA\tNA\tNA\tNA"); }
    else if ( args.bFreeMixOnly ) { bestSMF.printf("\t%.5lf\t%.2lf\t%.2lf\tNA\tNA",vbid.mixOut.fMix,vbid.mixOut.llk1,vbid.mixOut.llk0); }
    else if ( args.bFreeRefBiasOnly ) { bestSMF.printf("\tNA\t%.2lf\t%.2lf\t%.5lf\t%.5lf",vbid.mixOut.llk1,vbid.mixOut.llk0,vbid.mixOut.refHet,vbid.mixOut.refAlt); }
    else if ( args.bFreeFull ) { bestSMF.printf("\t%.5lf\t%.2lf\t%.2lf\t%.5lf\t%.5lf",vbid.mixOut.fMix,vbid.mixOut.llk1,vbid.mixOut.llk0,vbid.mixOut.refHet,vbid.mixOut.refAlt); }
    else { error("Invalid option in handling bFree"); }
    
    if ( args.bChipNone || bestInd < 0 ) { bestSMF.printf("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"); }
    else if ( args.bChipMixOnly ) { bestSMF.printf("\t%.5lf\t%.2lf\t%.2lf\tNA\tNA\t%.3lf\t%.4lf\t%.4lf",vbid.bestOut.fMix,vbid.bestOut.llk1,vbid.bestOut.llk0,(double)vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[1], (double)vbid.bestOut.numReads[2]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[2], (double)vbid.bestOut.numReads[3]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[3]); }
    else if ( args.bChipMixOnly ) { bestSMF.printf("\tNA\t%.2lf\t%.2lf\t%.5lf\t%.5lf\t%.3lf\t%.4lf\t%.4lf",vbid.bestOut.llk1, vbid.bestOut.llk0, vbid.bestOut.refHet, vbid.bestOut.refAlt, (double)vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[1], (double)vbid.bestOut.numReads[2]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[2], (double)vbid.bestOut.numReads[3]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[3]); }
    else if ( args.bChipFull ) { bestSMF.printf("\t%.5lf\t%.2lf\t%.2lf\t%.5lf\t%.5lf\t%.3lf\t%.4lf\t%.4lf", vbid.bestOut.fMix, vbid.bestOut.llk1, vbid.bestOut.llk0, vbid.bestOut.refHet, vbid.bestOut.refAlt, (double)vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[1], (double)vbid.bestOut.numReads[2]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[2], (double)vbid.bestOut.numReads[3]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[3]); }
    else { error("Invalid option in handling bChip"); }
    bestSMF.printf("\n");
    bestSMF.close();
  }

  time(&t);
  notice("Analysis finished on %s",ctime(&t));

  return 0;
}

/*
int runMerge(int argc, char ** argv) {
   printf("mpuTool merge -- merge multiple pileups\n");
   printf("(c) 2013 Hyun Min Kang, Matthew Flickinger, Goncalo Abecasis \n\n");

   std::string pedfile;
   std::string datfile;
   std::string colmpu("MPU");
   std::string outfile;
   std::string region;
   bool verbose = false;
   bool nobgzf = false;
   bool single = false;
   int maxDP = DEFAULT_MAX_DP;

   ParameterList pl;

   BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Pedigree File")
         LONG_STRINGPARAMETER("ped", &pedfile)
         LONG_STRINGPARAMETER("dat", &datfile)
         LONG_STRINGPARAMETER("col-mpu",&colmpu)
      LONG_PARAMETER_GROUP("Input Options")
         LONG_STRINGPARAMETER("region",&region)
         LONG_PARAMETER("notabix",&nobgzf)
         LONG_INTPARAMETER("max-dp",&maxDP)
      LONG_PARAMETER_GROUP("Output")
         LONG_PARAMETER("single",&single)
         LONG_STRINGPARAMETER("out", &outfile)
         LONG_PARAMETER("verbose", &verbose)
   END_LONG_PARAMETERS();

   pl.Add(new LongParameters("Available Options", longParameters));
   int argstart = pl.ReadWithTrailer(argc, argv) + 1;
   pl.Status();

   // sanity check of input arguments
   if ( outfile.empty()  ) {
     error("--out are required parameters");
   }

   time_t t;
   time(&t);

   printf("Analysis started on %s\n", ctime(&t));
   fflush(stdout);

   int n = 0;
   int ns = 0;
   mpuSet *mpu;
   int i, j, k;
   if ( argc > argstart ) {
     if ( !pedfile.empty() ) error("--ped option cannot be used with trailing arguments of MPU files");
     n = argc-argstart;
     ns = n;
     argv += argstart;

     mpu = new mpuSet[n];
     for(i=0; i < n; ++i) {
       mpu[i].nobgzf = nobgzf;
       mpu[i].maxDP = maxDP;
       if ( !mpu[i].load(argv[i],region.empty() ? NULL : region.c_str()) )
	 error("Failed to open MPU file [%s]", argv[i]);
       mpu[i].inds[0] = argv[i];
     }
   }
   else {
     if ( pedfile.empty() ) error("--ped option is required without trailing arguments of MPU files");
     fPed pedf(pedfile.c_str(), datfile.c_str());
     n = pedf.ninds;

     mpu = new mpuSet[n];
     //int firstMpu = n;
     for (i = n - 1; i >= 0; i--) {
       mpu[i].nobgzf = nobgzf;
       mpu[i].maxDP = maxDP;
       if (!mpu[i].load(pedf.getPheno(colmpu.c_str(),i).c_str(),region.empty() ? NULL : region.c_str()))
	 error("Failed to open MPU file [%s]", pedf.getPheno(colmpu.c_str(),i).c_str());
       //else
       //firstMpu = i;
       
       if ( !mpu[i].hasHeader ) {
	 mpu[i].inds[0] = pedf.getPheno("IND_ID",i).c_str();
       }
       ns += mpu[i].nInds;
     }
   }

   printf("Generating Merged Output files across %d individuals...\n", ns);

   // write header 
   mpuSet ompu(ns, maxDP);   
   wFile wf(outfile.c_str());
   for(i=0, k=0; i < n; ++i) {
     for(j=0; j < mpu[i].nInds; ++j, ++k) {
       ompu.inds[k] = mpu[i].inds[j];
     }
   }
   if ( !single ) ompu.writeHeader(wf);

   std::vector<bool> adv(n,false);
   while (true) {
     // find the minimum position
     int imin = 0;
     int c;
     std::fill(adv.begin(), adv.end(), false);
     adv[0] = true;

     ompu.hasMQ = mpu[0].hasMQ;
     ompu.hasCY = mpu[0].hasCY;

     for(i=1; i < n; ++i) {
       c = mpu[imin].compare(mpu[i]);
       if ( c > 0 ) {
	 std::fill(adv.begin(), adv.begin()+i, false);
	 imin = i;
	 adv[i] = true;
	 ompu.hasMQ = mpu[i].hasMQ;
	 ompu.hasCY = mpu[i].hasCY;
       }
       else if ( c == 0 ) {
	 adv[i] = true;
	 ompu.hasMQ |= mpu[i].hasMQ;
	 ompu.hasCY |= mpu[i].hasCY;
       }
     }
     if ( mpu[imin].isEOF() ) break;

     // assign values in the marker
     ompu.chrom = mpu[imin].chrom;
     ompu.pos = mpu[imin].pos;
     ompu.id = mpu[imin].id;
     ompu.ref = mpu[imin].ref;
     ompu.alt = mpu[imin].alt;
     ompu.qual = mpu[imin].qual;
     ompu.filter = mpu[imin].filter;
     ompu.info = mpu[imin].info;
     for(i=0, j=0; i < n; ++i) {
       if ( adv[i] ) {
	 std::copy(mpu[i].nBases.begin(), mpu[i].nBases.end(), ompu.nBases.begin()+j);
	 std::copy( mpu[i].bases.begin(), mpu[i].bases.end(), ompu.bases.begin()+(j*(maxDP+1)));
	 std::copy( mpu[i].bQs.begin(), mpu[i].bQs.end(), ompu.bQs.begin()+(j*(maxDP+1)));
	 std::copy( mpu[i].mQs.begin(), mpu[i].mQs.end(), ompu.mQs.begin()+(j*(maxDP+1)));
	 std::copy( mpu[i].cycles.begin(), mpu[i].cycles.end(), ompu.cycles.begin()+(j*(maxDP+1)));
	 std::copy( mpu[i].strands.begin(), mpu[i].strands.end(), ompu.strands.begin()+(j*(maxDP+1)));
	 //std::copy( ompu.nBases.begin()+j, ompu.nBases.begin()+(j+mpu[i].nInds), mpu[i].nBases);
	 //std::copy( ompu.bases.begin()+j*(maxDP+1), ompu.bases.begin()+(j+mpu[i])*(maxDP+1).nInds, mpu[i].bases);
	 //std::copy( ompu.bQs.begin()+j*(maxDP+1), ompu.bQs.begin()+(j+mpu[i])*(maxDP+1).nInds, mpu[i].bQs);
	 //std::copy( ompu.mQs.begin()+j*(maxDP+1), ompu.mQs.begin()+(j+mpu[i])*(maxDP+1).nInds, mpu[i].mQs);
	 //std::copy( ompu.cycles.begin()+j*(maxDP+1), ompu.cycles.begin()+(j+mpu[i])*(maxDP+1).nInds, mpu[i].cycles);
	 //std::copy( ompu.strands.begin()+j*(maxDP+1), ompu.strands.begin()+(j+mpu[i])*(maxDP+1).nInds, mpu[i].strands);
       }
       else {
	 for(k=0; k < mpu[i].nInds; ++k) {
	   ompu.nBases[j+k] = 0;
	 }
       }
       j += mpu[i].nInds;
     }

     // print the marker
     ompu.writeMarker(wf,single);

     for(i=0; i < n; ++i) {
       if ( adv[i] ) mpu[i].next();
     }
   }
   delete [] mpu;

   time(&t);
   printf("\nAnalysis completed on %s\n", ctime(&t));
   fflush(stdout);

   return 0;
}
*/
/*
int runSetGenotype(int argc, char ** argv) {
   printf("mpuSetGenotype -- Create VCF based on pileups files\n");
   printf("(c) 2013 Hyun Min Kang, Matthew Flickinger, Goncalo Abecasis \n\n");

   std::string pedfile;
   std::string datfile;
   std::string msetfile;
   std::string colmix;
   std::string colpair;
   std::string callfile;
   std::string region;
   std::string invcf;
   std::string rule;
   double minMAF = 0;
   double minCallRate = 0;
   bool ignoreFilter = false;
   int unit = 10000L;
   bool verbose = false;

   ParameterList pl;

   std::string xLabel("X"), yLabel("Y"), mitoLabel("MT");
   int    xStart = 2699520, xStop = 154931044;

   BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Pileup File")
         LONG_STRINGPARAMETER("mset",&msetfile)
         LONG_STRINGPARAMETER("region",&region)
      LONG_PARAMETER_GROUP("Input VCF")
         LONG_STRINGPARAMETER("invcf", &invcf)
         LONG_STRINGPARAMETER("rule",&rule)
         LONG_INTPARAMETER("unit",&unit)
         LONG_DOUBLEPARAMETER("minMAF",&minMAF)
         LONG_DOUBLEPARAMETER("minCallRate",&minCallRate)
         LONG_PARAMETER("ignoreFilter",&ignoreFilter)
      LONG_PARAMETER_GROUP("Pedigree File")
         LONG_STRINGPARAMETER("ped", &pedfile)
         LONG_STRINGPARAMETER("dat", &datfile)
         LONG_STRINGPARAMETER("col-mix",&colmix)
         LONG_STRINGPARAMETER("col-pair",&colpair)
      LONG_PARAMETER_GROUP("Chromosome Labels")
         LONG_STRINGPARAMETER("xChr", &xLabel)
         LONG_STRINGPARAMETER("yChr", &yLabel)
         LONG_STRINGPARAMETER("mito", &mitoLabel)
         LONG_INTPARAMETER("xStart", &xStart)
         LONG_INTPARAMETER("xStop", &xStop)
      LONG_PARAMETER_GROUP("Output")
         LONG_STRINGPARAMETER("out", &callfile)
         LONG_PARAMETER("verbose", &verbose)
   END_LONG_PARAMETERS();

   pl.Add(new LongParameters("Available Options", longParameters));
   pl.Read(argc, argv);
   pl.Status();

   // sanity check of input arguments
   if ( callfile.empty() || msetfile.empty() ) {
     error("--out, --mset are required parameters");
   }

   time_t t;
   time(&t);

   printf("Analysis started on %s\n", ctime(&t));
   fflush(stdout);

   mpuSet mset(msetfile.c_str(),region.empty() ? NULL : region.c_str(),
	       region.empty() ? false : true);
   mpuSetLikelihood mlk(&mset, pedfile.empty() ? NULL : pedfile.c_str(),
			datfile.empty() ? NULL : datfile.c_str(),
			colmix.empty() ? NULL : colmix.c_str(),
			colpair.empty() ? NULL : colpair.c_str());

   int n = mset.nInds;

   mlk.hasMix  = !colmix.empty();
   mlk.hasPair = !colpair.empty();

   wFile baseCalls(callfile.c_str());
   baseCalls.printf("##fileformat=VCFv4.1\n");
   ReportDate(baseCalls);
   baseCalls.printf("##source=mpuTool\n");
   baseCalls.printf("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth at Site\">\n");
   baseCalls.printf("##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root Mean Squared Mapping Quality\">\n");
   baseCalls.printf("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Coverage\">\n");
   baseCalls.printf("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of Alleles in Samples with Coverage\">\n");
   baseCalls.printf("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Alternate Allele Counts in Samples with Coverage\">\n");
   baseCalls.printf("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Alternate Allele Frequencies\">\n");
   baseCalls.printf("##INFO=<ID=AB,Number=1,Type=Float,Description=\"Allele Balance in Heterozygotes\">\n");
   baseCalls.printf("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
   baseCalls.printf("##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Genotype Likelihoods for Genotypes 0/0,0/1,1/1\">\n");
   baseCalls.printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
   for (int i = 0; i < n; i++) {
     baseCalls.printf("\t%s", mset.inds[i].c_str());
   }
   baseCalls.printf("\n");

   if ( !invcf.empty() ) {
     fVcf tvcf;
     tvcf.load(invcf.c_str(), region.c_str(), "GT", rule.c_str(), !ignoreFilter, NULL);

     int M, i;
     double af, maf;
     for(M=0; tvcf.readMarkers(unit); ) {
       M += tvcf.nMarkers;
       fprintf(stderr,"Processing %d markers..\n", M);
       for(i=0, m=0; i < tvcf.nMarkers; ++i) { // for each marker
	 af = tvcf.alleleFreq(i);
	 maf = af > 0.5 ? 1-af : af;
	 if ( ( maf >= minMAF ) && 
	      ( tvcf.callRate(i) >= minCallRate ) )
	   { 
	     uint8_t refBase = BaseAsciiMap::base2int[(int)tvcf.refs[i][0]];
	     uint8_t altBase = BaseAsciiMap::base2int[(int)tvcf.alts[i][0]];
	     uint8_t allele1 = refBase;
	     uint8_t allele2 = altBase;

	     //notice("Processing %s:%d_%c/%c",tvcf.chroms[i].c_str(), tvcf.pos1s[i], tvcf.refs[i][0], tvcf.alts[i][0]);
	     if ( mset.advanceTo(tvcf.chroms[i].c_str(), tvcf.pos1s[i]) )
	       mlk.ReportSNP(baseCalls,allele1,allele2);
	     else 
	       notice("Skipping %s:%d_%c/%d",tvcf.chroms[i].c_str(), tvcf.pos1s[i], tvcf.refs[i][0], tvcf.alts[i][0]);
	   }
       }
     }
     tvcf.close();
   }
   else {
     for(int i=0; mset.next(); ++i) {
       if ( (i+1) % 1000 == 0 ) 
	 notice("Processing %d markers at %s:%d",i+1,mset.chrom.c_str(),mset.pos);
       mlk.ReportSNP(baseCalls);
     }
   }
   baseCalls.close();

   time(&t);
   printf("\nAnalysis completed on %s\n", ctime(&t));
   fflush(stdout);

   return 0;
}
*/

int runDeconvoluteGenotype(int argc, char ** argv) {
   printf("mpuMixup -- Extract VCF based on pilups files\n");
   printf("(c) 2013 Hyun Min Kang, Matthew Flickinger, Goncalo Abecasis \n\n");

   std::string pedfile; // PED file containing pileups
   std::string datfile; // DAT file if PED header is not self-contained
   std::string colmpu("MPU");  // COLUMN NAME containing pileup file
   std::string colcon;         // COLUMN NAME containing contamination estimate
   std::string colmixid;       // COLUMN NAME containing contaminating sample name (if known)
   std::string invcf;          // INPUT VCF containing the site list to genotype
   std::string callfile;       // OUTPUT VCF 
   std::string region;         // STRING with [chr]:[beg]-[end] format specifying region to call
   std::string rule;           // FILTER in the INFO field in the input VCF
   double minMAF = 0;          // MAF threshold from the input VCF
   double minCallRate = 0;     // CALL RATE threshold from the input VCF
   bool ignoreFilter = false;  // INCLUDE non-PASS variant if set

   double defaultContam = 0;   // DEFAULT CONTAMINATION ESTIMATE for every sample
   double offsetContam = 0;    // OFFSET OF CONTAMINATION ESTIMATE to be added to every sample
   double scaleContam = 1.;    // SCALE PARAMETER OF THE CONTAMINATION ESTIMATE to be multiplied 
   double thresContam = 0;     // THRESHOLD OF CONTAMINATION ESTIMATE TO BE IGNORED
   double inbreedingCoeff = 0; // INBREEDING COEFFICIENT for modeling structured data

   bool verbose = false;
   bool nobgzf = false;
   bool reportPairGL = false;
   bool printMono = false;
   int unit = 10000L;

   ParameterList pl;

   std::string xLabel("XT"), yLabel("YT"), mitoLabel("MT");
   int    xStart = 2699520, xStop = 154931044;

   BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Pedigree File")
         LONG_STRINGPARAMETER("ped", &pedfile)
         LONG_STRINGPARAMETER("dat", &datfile)
         LONG_STRINGPARAMETER("col-mpu",&colmpu)
         LONG_STRINGPARAMETER("col-con",&colcon)
         LONG_STRINGPARAMETER("col-mixid",&colmixid)
         LONG_DOUBLEPARAMETER("default-contam",&defaultContam)
         LONG_DOUBLEPARAMETER("thres-contam",&thresContam)
         LONG_DOUBLEPARAMETER("offset-contam",&offsetContam)
         LONG_DOUBLEPARAMETER("scale-contam",&scaleContam)
         LONG_DOUBLEPARAMETER("inbreeding-coeff",&inbreedingCoeff)
      LONG_PARAMETER_GROUP("Input VCF")
         LONG_STRINGPARAMETER("invcf", &invcf)
         LONG_STRINGPARAMETER("region",&region)
         LONG_STRINGPARAMETER("rule",&rule)
         LONG_INTPARAMETER("unit",&unit)
         LONG_DOUBLEPARAMETER("minMAF",&minMAF)
         LONG_DOUBLEPARAMETER("minCallRate",&minCallRate)
         LONG_PARAMETER("ignoreFilter",&ignoreFilter)
         LONG_PARAMETER("report-pair-GL",&reportPairGL)
         LONG_PARAMETER("print-mono",&printMono)
         LONG_PARAMETER("notabix",&nobgzf)
      LONG_PARAMETER_GROUP("Chromosome Labels")
         LONG_STRINGPARAMETER("xChr", &xLabel)
         LONG_STRINGPARAMETER("yChr", &yLabel)
         LONG_STRINGPARAMETER("mito", &mitoLabel)
         LONG_INTPARAMETER("xStart", &xStart)
         LONG_INTPARAMETER("xStop", &xStop)
      LONG_PARAMETER_GROUP("Output")
         LONG_STRINGPARAMETER("out", &callfile)
         LONG_PARAMETER("verbose", &verbose)
   END_LONG_PARAMETERS();

   pl.Add(new LongParameters("Available Options", longParameters));
   pl.Read(argc, argv);
   pl.Status();

   // sanity check of input arguments
   if ( invcf.empty() || pedfile.empty() || callfile.empty()  ) {
     error("--invcf, --ped, --out are required parameters");
   }

   time_t t;
   time(&t);

   printf("Analysis started on %s\n", ctime(&t));
   fflush(stdout);

   fPed pedf(pedfile.c_str(), datfile.c_str());
   int n = pedf.ninds;

   mpuFile *mpu = new mpuFile[n];
   int* pairs = new int[n];

   for (int i = n - 1; i >= 0; i--) {
     mpu[i].nobgzf = nobgzf;
     if (!mpu[i].load(pedf.getPheno(colmpu.c_str(),i).c_str()))
	 error("Failed to open MPU file [%s]", mpu[i].load(pedf.getPheno(colmpu.c_str(),i).c_str()));

     if ( colmixid.empty() ) {
       pairs[i] = -1;
     }
     else {
       const std::string& mixId = pedf.getPheno(colmixid.c_str(),i);
       if ( pedf.iinds.find(mixId) != pedf.iinds.end() ) {
	 pairs[i] = pedf.iinds[mixId];
       }
       else {
	 pairs[i] = -1;
       }
     }
   }

   bool hasContam = false;
   if ( colcon.empty() ) {
     for(int i=0; i < n; ++i) { 
       mpu[i].fmix = defaultContam * scaleContam + offsetContam ; 
       if ( mpu[i].fmix > 1 ) mpu[i].fmix = 1;
       else if ( mpu[i].fmix < 0 ) mpu[i].fmix = 0;
     }
   }
   else {
     for(int i=0; i < n; ++i) { 
       mpu[i].fmix = atof(pedf.getPheno(colcon.c_str(),i).c_str()) * scaleContam + offsetContam;
       if ( mpu[i].fmix < thresContam ) {
	 mpu[i].fmix = defaultContam * scaleContam + offsetContam;
       }
       else {
	 hasContam = true;
       }
       if ( mpu[i].fmix > 1 ) mpu[i].fmix = 1;
       else if ( mpu[i].fmix < 0 ) mpu[i].fmix = 0;
     }
   }
   if ( defaultContam * scaleContam + offsetContam > 0 ) hasContam = true;

   printf("Generating VCFs for files ...\n");
   for (int i = 0; i < n; i++)
     printf("%s\t%lf\t%s\n", pedf.getPheno(colmpu.c_str(),i).c_str(), mpu[i].fmix,pairs[i] < 0 ? "." : pedf.getPheno("IND_ID",pairs[i]).c_str());
   printf("\n");

   wFile baseCalls(callfile.c_str());
   baseCalls.printf("##fileformat=VCFv4.1\n");
   ReportDate(baseCalls);
   baseCalls.printf("##source=mpuDeconvoluteGenotype\n");
   baseCalls.printf("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth at Site\">\n");
   baseCalls.printf("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Coverage, Less is worse\">\n");
   baseCalls.printf("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of Alleles in Samples with Coverage, Less is worse\">\n");
   baseCalls.printf("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Alternate Allele Counts in Samples with Coverage\">\n");
   baseCalls.printf("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Alternate Allele Frequencies\">\n");
   baseCalls.printf("##INFO=<ID=HWDAF,Number=.,Type=Float,Description=\"Genotype Frequencies Under Hardy-Weinberg Disequilbrium\">\n");
   baseCalls.printf("##INFO=<ID=FIC,Number=1,Type=Float,Description=\"Inbreeding Coefficient Statistics. Negative is worse\">\n");
   baseCalls.printf("##INFO=<ID=SLRT,Number=1,Type=Float,Description=\"Signed likelihood-ratio test statistics testing Hardy-Weinberg Equilbrium. Negative is worse\">\n");
   baseCalls.printf("##INFO=<ID=ABL,Number=1,Type=Float,Description=\"Allele Balance in Heterozygotes Computed from Genotype Likelihood. High is worse\">\n");
   baseCalls.printf("##INFO=<ID=STR,Number=1,Type=Float,Description=\"Strand-bias Pearson's correlation. Extreme is worse\">\n");
   baseCalls.printf("##INFO=<ID=STZ,Number=1,Type=Float,Description=\"Strand-bias z-score. Extreme is worse\">\n");
   baseCalls.printf("##INFO=<ID=TBR,Number=1,Type=Float,Description=\"Tail-bias Pearson's correlation. Positive is worse\">\n");
   baseCalls.printf("##INFO=<ID=TBZ,Number=1,Type=Float,Description=\"Strand-bias z-score. Positive is worse\">\n");
   baseCalls.printf("##INFO=<ID=FOB,Number=1,Type=Float,Description=\"Fraction of non-ref, non-alt bases. High is worse\">\n");
   baseCalls.printf("##INFO=<ID=MMQ,Number=1,Type=Integer,Description=\"Mean mapping quality. Low is worse\">\n");
   baseCalls.printf("##INFO=<ID=AMQ,Number=1,Type=Integer,Description=\"Mean mapping quality for alternative alleles only. Low is worse\">\n");
   baseCalls.printf("##INFO=<ID=LMQ,Number=1,Type=Integer,Description=\"Fraction of low mapping quality reads. High is worse\">\n");
   baseCalls.printf("##INFO=<ID=ABL,Number=1,Type=Float,Description=\"Allele balance in heterozygotes computed from called genotypes. High iw worse\">\n");
   baseCalls.printf("##INFO=<ID=NRO,Number=1,Type=Float,Description=\"Enrichment of non-ref, non-alt bases in non-reference genotypes. High is worse\">\n");
   baseCalls.printf("##INFO=<ID=NMA,Number=1,Type=Integer,Description=\"Maximum non-ref allele count per sample. Low is worse\">\n");
   baseCalls.printf("##INFO=<ID=FMA,Number=1,Type=Integer,Description=\"Maximum fraction of non-ref alleles per sample. Low is worse\">\n");
   baseCalls.printf("##INFO=<ID=QD,Number=1,Type=Float,Description=\"Quality per depth. Low is worse\">\n");
   baseCalls.printf("##INFO=<ID=SHR,Number=1,Type=Float,Description=\"Strand bias correlation coefficient in heterozygotes only. Extreme is worse\">\n");
   baseCalls.printf("##INFO=<ID=SHZ,Number=1,Type=Float,Description=\"Strand bias z-score in heterozygotes only. Extreme is worse\">\n");
   baseCalls.printf("##INFO=<ID=THR,Number=1,Type=Float,Description=\"Tail bias correlation coefficient in heterozygotes only. High is worse\">\n");
   baseCalls.printf("##INFO=<ID=THZ,Number=1,Type=Float,Description=\"Tail bias z-score in heterozygotes only. High is worse\">\n");
   baseCalls.printf("##INFO=<ID=LQR,Number=1,Type=Float,Description=\"Ratio of low quality bases. High is worse\">\n");

   baseCalls.printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
   baseCalls.printf("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
   baseCalls.printf("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
   baseCalls.printf("##FORMAT=<ID=GL,Number=.,Type=Float,Description=\"Genotype Likelihoods for Genotypes 0/0,0/1,1/1\">\n");
   baseCalls.printf("##FORMAT=<ID=LM,Number=.,Type=Float,Description=\"Genotype Likelihoods for Pair of Genotypes\">\n");
   baseCalls.printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
   for (int i = 0; i < n; i++) {
     baseCalls.printf("\t%s", pedf.getPheno("IND_ID",i).c_str());
   }
   baseCalls.printf("\n");

   std::string buffer;

   std::string curChrom;
   int curPos;
   int chromosomeType = CT_AUTOSOME;
   mpuMixupLikelihood lkGeno(n, mpu, pairs, hasContam, inbreedingCoeff);

   fVcf tvcf;
   tvcf.load(invcf.c_str(), region.c_str(), "GT", rule.c_str(), !ignoreFilter, NULL);
   //for (int i = 0; i < n; i++)
   //lkGeno.sexes[i] = (atoi(pedf.getPheno("SEX",i).c_str()) == SEX_MALE ? SEX_MALE : SEX_FEMALE);

   int M, i, j;
   double af, maf;
   for(M=0; tvcf.readMarkers(unit); ) {
     M += tvcf.nMarkers;
     fprintf(stderr,"Processing %d markers..\n", M);
     for(i=0; i < tvcf.nMarkers; ++i) { // for each marker
       af = tvcf.alleleFreq(i);
       maf = af > 0.5 ? 1-af : af;
       if ( ( maf >= minMAF ) && 
	    ( tvcf.callRate(i) >= minCallRate ) &&
	    ( tvcf.refs[i].size() == 1 ) &&
	    ( tvcf.alts[i].size() == 1 ) 
	    )
	 { 
	   uint8_t refBase = BaseAsciiMap::base2int[(int)tvcf.refs[i][0]];
	   uint8_t altBase = BaseAsciiMap::base2int[(int)tvcf.alts[i][0]];
	   uint8_t allele1 = refBase;
	   uint8_t allele2 = altBase;
	   curChrom = tvcf.chroms[i];
	   curPos = tvcf.pos1s[i]-1;

	   if (curChrom == xLabel) chromosomeType = CT_CHRX;
	   if (curChrom == yLabel) chromosomeType = CT_CHRY;
	   if (curChrom == mitoLabel) chromosomeType = CT_MITO;

	   //int     totalDepth = 0, nSamplesCovered = 0;
	   //double  rmsMapQuality = 0.0;

	   //notice("foo");
	   
	   for (j = 0; j < n; ++j) {
	     //int depth = 0;
	     //notice("j=%d",j);
	     mpu[j].advanceTo(tvcf.chroms[i].c_str(), tvcf.pos1s[i], tvcf.refs[i][0], tvcf.alts[i][0]);
	     //if ( mpu[j].advanceTo(tvcf.chroms[i].c_str(), tvcf.pos1s[i], tvcf.refs[i][0], tvcf.alts[i][0]) ) {
	        //depth = mpu[j].nbase;
	       //totalDepth += depth;
	       //++nSamplesCovered;
	     //}
	     //if ( mpu[j].ref == 'N' ) break;
	   }

	   //if ( j < n ) continue;

	   //notice("bar");
	   
	   lkGeno.chromosomeType = chromosomeType != CT_CHRX ? chromosomeType :
	     curPos >= xStart && curPos <= xStop ? CT_CHRX : CT_AUTOSOME;
	   
	   int qual = 100;
	   ReportMixupSNP(baseCalls, lkGeno, n, curChrom.c_str(), curPos, refBase, allele1, allele2, phredConv.phred2Mat[qual], printMono, reportPairGL);	 
	 }
     }
   }
   tvcf.close();
   fprintf(stderr,"Finished processing total of %d markers\n", M);

   baseCalls.close();

   delete [] mpu;
   delete [] pairs;

   time(&t);
   printf("\nAnalysis completed on %s\n", ctime(&t));
   fflush(stdout);

   return 0;
}

/*
// new version of verifyBamID to detect contamination by population
int runDetectContamination(int argc, char** argv) {
  printf("mpuDetectContamination 1.0.0 -- verify identity and purity of sequence reads\n"
        "(c) 2013 Hyun Min Kang, Goo Jun, and Goncalo Abecasis\n\n");

  mpuVerifyArgs args;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Files")
    LONG_STRINGPARAMETER("vcf",&args.sVcfFile)
    LONG_STRINGPARAMETER("mpu",&args.sMpuFile)
    LONG_STRINGPARAMETER("subset",&args.sSubsetInds)
    LONG_STRINGPARAMETER("smID",&args.sSMID)
    LONG_STRINGPARAMETER("scanPair",&args.sScanPair)

    LONG_PARAMETER_GROUP("VCF analysis options")
    LONG_DOUBLEPARAMETER("genoError",&args.genoError)
    LONG_DOUBLEPARAMETER("minAF",&args.minAF)
    LONG_DOUBLEPARAMETER("minCallRate",&args.minCallRate)

    LONG_PARAMETER_GROUP("Individuals to compare with chip data")
    EXCLUSIVE_PARAMETER("site",&args.bSiteOnly)
    EXCLUSIVE_PARAMETER("self",&args.bSelfOnly)
    EXCLUSIVE_PARAMETER("best",&args.bFindBest)

    LONG_PARAMETER_GROUP("Chip-free optimization options")
    EXCLUSIVE_PARAMETER("free-none",&args.bFreeNone)
    EXCLUSIVE_PARAMETER("free-mix",&args.bFreeMixOnly)
    EXCLUSIVE_PARAMETER("free-refBias",&args.bFreeRefBiasOnly)
    EXCLUSIVE_PARAMETER("free-full",&args.bFreeFull)

    LONG_PARAMETER_GROUP("With-chip optimization options")
    EXCLUSIVE_PARAMETER("chip-none",&args.bChipNone)
    EXCLUSIVE_PARAMETER("chip-mix",&args.bChipMixOnly)
    EXCLUSIVE_PARAMETER("chip-refBias",&args.bChipRefBiasOnly)
    EXCLUSIVE_PARAMETER("chip-full",&args.bChipFull)

    LONG_PARAMETER_GROUP("MPU analysis options")
    LONG_PARAMETER("precise",&args.bPrecise)
    LONG_PARAMETER("ignore-sex",&args.bIgnoreSex)
    LONG_INTPARAMETER("minMapQ",&args.minMapQ)
    LONG_INTPARAMETER("maxDepth",&args.maxDepth)
    LONG_INTPARAMETER("minQ",&args.minQ)
    LONG_INTPARAMETER("maxQ",&args.maxQ)
    LONG_DOUBLEPARAMETER("grid",&args.grid)

    LONG_PARAMETER_GROUP("Modeling Reference Bias")
    LONG_DOUBLEPARAMETER("refRef",&args.pRefRef)
    LONG_DOUBLEPARAMETER("refHet",&args.pRefHet)
    LONG_DOUBLEPARAMETER("refAlt",&args.pRefAlt)

    LONG_PARAMETER_GROUP("Output options")
    LONG_STRINGPARAMETER("out",&args.sOutFile)
    LONG_PARAMETER("verbose",&args.bVerbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options",longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // check the validity of input files
  if ( args.sVcfFile.empty() ) {
    error("--vcf [vcf file] required");
  }

  if ( args.sMpuFile.empty() ) {
    error("--mpu [mpu file] is required");
  }

  if ( args.sOutFile.empty() ) {
    error("--out [output prefix] is required");
  }

  if ( args.sSMID.empty() ) {
    error("--smID [sample ID] is required");
  }

  if ( ! ( args.bSiteOnly || args.bSelfOnly || args.bFindBest ) ) {
    warning("--self option was autotomatically turned on by default. Specify --best option if you wanted to check across all possible samples in the VCF");
    args.bSelfOnly = true;
  }

  if ( !( args.bFreeNone || args.bFreeMixOnly || args.bFreeRefBiasOnly || args.bFreeFull ) ) {
    warning("--free-mix option is automatically turned on by default");
    args.bFreeMixOnly = true;
  }

  if ( !( args.bChipNone || args.bChipMixOnly || args.bChipRefBiasOnly || args.bChipFull ) ) {
    warning("--chip-mix option is automatically turned on by default");
    args.bChipMixOnly = true;
  }

  if ( ( args.maxDepth > 20 ) && ( !args.bPrecise ) ) {
    warning("--precise option is not turned on at --maxDepth %d : may be prone to precision errors",args.maxDepth);
  }

  if ( ( args.bChipRefBiasOnly ) && ( !args.bSelfOnly ) ) {
    error("--self must be set for --chip-refBias to work. Skipping..");
  }

  if ( ( !args.sScanPair.empty() ) && ( ( !args.bFreeNone ) || ( !args.bFindBest ) || ( !args.bChipMixOnly ) ) ) {
    error("--scanPair option is only compatible with --best, --free-none, and --chip-mix");
  }

  // check timestamp
  time_t t;
  time(&t);
  notice("Analysis started on %s",ctime(&t));

  // load arguments
  mpuVerify vbid(&args);

  // load input VCF and BAM files
  notice("Opening Input Files");
  vbid.loadFiles(args.sMpuFile.c_str(), args.sVcfFile.c_str());

  int selfIndex = -1;
  if ( !args.sScanPair.empty() ) { // if scanPair is on
    std::string indID1(args.sScanPair.c_str());
    for(int i=0; i < (int)vbid.pGenotypes->indids.size(); ++i) {
      if ( vbid.pGenotypes->indids[i] == indID1 ) {
	selfIndex = i;
      }
    }
    if ( selfIndex < 0 ) {
      error("Cannot find individual %s from the VCF file",indID1.c_str());
    }
  }

  // Check which genotype-free method is used
  if ( args.bFreeNone ) {  // if no genotype-free mode is tested. skip it
    // do nothing for genotype-free estimation
    notice("Skipping chip-free estimation of sample mixture");
  }
  else if ( args.bFreeMixOnly ) { // only mixture is estimated.
    // genotype-free method
    notice("Performing chip-free estimation of sample mixture at fixed reference bias parameters (%lf, %lf, %lf)",args.pRefRef,args.pRefHet,args.pRefAlt);

    mpuVerify::mixLLK mix(&vbid);
    mix.OptimizeLLK();
    notice("Optimal per-sample fMix = %lf, LLK0 = %lf, LLK1 = %lf\n",mix.fMix,mix.llk0,mix.llk1);
    vbid.mixOut.llk0 = mix.llk0;
    vbid.mixOut.llk1 = mix.llk1;
    vbid.mixOut.fMix = mix.fMix;
  }
  else if ( args.bFreeRefBiasOnly ) {
    notice("Performing chip-free estimation of reference-bias without sample mixture");
    mpuVerify::refBiasMixLLKFunc myFunc(&vbid);
    AmoebaMinimizer myMinimizer;
    Vector startingPoint(2);
    startingPoint[0] = 0;      // pRefHet = 0.5
    startingPoint[1] = -4.595; // pRefAlt = 0.01
    myMinimizer.func = &myFunc;
    myMinimizer.Reset(2);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(1e-6);
    double pRefHet = mpuVerify::invLogit(myMinimizer.point[0]);
    double pRefAlt = mpuVerify::invLogit(myMinimizer.point[1]);
    notice("Reference Bias Estimated as ( Pr[refBase|HET] = %lf, Pr[refBase|ALT] = %lf) with LLK = %lf",pRefHet,pRefAlt,myMinimizer.fmin);
    vbid.mixOut.llk0 = myFunc.llk0;
    vbid.mixOut.llk1 = myFunc.llk1;
    vbid.mixOut.refHet = myFunc.pRefHet;
    vbid.mixOut.refAlt = myFunc.pRefAlt;
  }
  else if ( args.bFreeFull ) {
    notice("Performing chip-free estimation of reference-bias and sample mixture together");
    mpuVerify::fullMixLLKFunc myFunc(&vbid);
    AmoebaMinimizer myMinimizer;
    Vector startingPoint(3);
    startingPoint[0] = -3.91;  // start with fMix = 0.01
    startingPoint[1] = 0;      // pRefHet = 0.5
    startingPoint[2] = -4.595; // pRefAlt = 0.01
    myMinimizer.func = &myFunc;
    myMinimizer.Reset(3);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(1e-6);
    double fMix = mpuVerify::invLogit(myMinimizer.point[0]);
    if ( fMix > 0.5 ) 
      fMix = 1.-fMix;
    double pRefHet = mpuVerify::invLogit(myMinimizer.point[1]);
    double pRefAlt = mpuVerify::invLogit(myMinimizer.point[2]);
    notice("Optimal per-sample fMix = %lf\n",fMix);
    notice("Reference Bias Estimated as ( Pr[refBase|HET] = %lf, Pr[refBase|ALT] = %lf) with LLK = %lf",pRefHet,pRefAlt,myMinimizer.fmin);

    vbid.mixOut.llk0 = myFunc.llk0;
    vbid.mixOut.llk1 = myFunc.llk1;
    vbid.mixOut.fMix = myFunc.fMix;
    vbid.mixOut.refHet = myFunc.pRefHet;
    vbid.mixOut.refAlt = myFunc.pRefAlt;
  }
  notice("calculating depth distribution");  
  vbid.calculateDepthDistribution(args.maxDepth, vbid.mixOut);

  notice("finished calculating depth distribution");  

  int bestInd = -1;
  int selfInd = -1;

  if ( args.bChipNone ) {
    // do nothing
    notice("Skipping with-chip estimation of sample mixture");
  }
  else if ( args.bChipMixOnly ) {
    notice("Performing with-chip estimation of sample mixture at fixed reference bias parameter (%lf, %lf, %lf)",args.pRefRef,args.pRefHet,args.pRefAlt);

    if ( args.sScanPair.empty() ) {
      double maxIBD = -1;
      mpuVerify::ibdLLK ibd(&vbid);
      for(int i=0; i < (int)vbid.pGenotypes->indids.size(); ++i) {
	double fIBD = ibd.OptimizeLLK(i);
	notice("Comparing with individual %s.. Optimal fIBD = %lf, LLK0 = %lf, LLK1 = %lf",vbid.pGenotypes->indids[i].c_str(),fIBD, ibd.llk0, ibd.llk1);
	if ( maxIBD < fIBD ) {
	  bestInd = i;
	  vbid.bestOut.llk0 = ibd.llk0;
	  vbid.bestOut.llk1 = ibd.llk1;
	  vbid.bestOut.fMix = 1-ibd.fIBD;
	  maxIBD = ibd.fIBD;
	}
	  
	if ( vbid.pPile->sSM == vbid.pGenotypes->indids[i] ) {
	  selfInd = i;
	  vbid.selfOut.llk0 = ibd.llk0;
	  vbid.selfOut.llk1 = ibd.llk1;
	  vbid.selfOut.fMix = 1-ibd.fIBD;
	}
      }
	
      if ( bestInd >= 0 ) {
	notice("Best Matching Individual is %s with IBD = %lf",vbid.pGenotypes->indids[bestInd].c_str(),maxIBD);
	vbid.calculateDepthByGenotype(bestInd,vbid.bestOut);
      }
      
      if ( selfInd >= 0 ) {
	notice("Self Individual is %s with IBD = %lf",vbid.pGenotypes->indids[selfInd].c_str(),vbid.selfOut.fMix);
	vbid.calculateDepthByGenotype(selfInd,vbid.selfOut);
      }
    }
    else {
      double maxLLK = 1e9;
      mpuVerify::pairLLK ibdPair(&vbid);
      for(int i=0; i < (int)vbid.pGenotypes->indids.size(); ++i) {
	double fIBD = ibdPair.OptimizeLLK(selfIndex, i);
	notice("Comparing with individual %s.. Optimal fIBD = %lf, LLK0 = %lf, LLK1 = %lf",vbid.pGenotypes->indids[i].c_str(),fIBD, ibdPair.llk0, ibdPair.llk1);
	double fLLK = ibdPair.llk1;
	if ( fLLK < maxLLK ) {
	  bestInd = i;
	  vbid.bestOut.llk0 = ibdPair.llk0;
	  vbid.bestOut.llk1 = ibdPair.llk1;
	  vbid.bestOut.fMix = 1-ibdPair.fIBD;
	  maxLLK = fLLK;
	}
	
	if (vbid.pPile->sSM == vbid.pGenotypes->indids[i] ) {
	  selfInd = i;
	  vbid.selfOut.llk0 = ibdPair.llk0;
	  vbid.selfOut.llk1 = ibdPair.llk1;
	  vbid.selfOut.fMix = 1-ibdPair.fIBD;
	}
      }
      
      if ( bestInd >= 0 ) {
	notice("Best Matching Individual is %s with LLK = %lf",vbid.pGenotypes->indids[bestInd].c_str(),maxLLK);
	vbid.calculateDepthByGenotype(bestInd,vbid.bestOut);
      }
      
      if ( selfInd >= 0 ) {
	notice("Self Individual is %s with IBD = %lf",vbid.pGenotypes->indids[selfInd].c_str(),vbid.selfOut.fMix);
	vbid.calculateDepthByGenotype(selfInd,vbid.selfOut);
      }
    }
  }
  else if ( args.bChipRefBiasOnly ) {
    notice("Performing with-chip estimation of reference-bias without sample mixture");
    if ( args.bSelfOnly ) {
      mpuVerify::refBiasIbdLLKFunc myFunc(&vbid);
      AmoebaMinimizer myMinimizer;
      Vector startingPoint(2);
      startingPoint[0] = 0;      // pRefHet = 0.5
      startingPoint[1] = -4.595; // pRefAlt = 0.01
      myMinimizer.func = &myFunc;
      myMinimizer.Reset(2);
      myMinimizer.point = startingPoint;
      myMinimizer.Minimize(1e-6);
      double pRefHet = mpuVerify::invLogit(myMinimizer.point[0]);
      double pRefAlt = mpuVerify::invLogit(myMinimizer.point[1]);
      notice("Reference Bias Estimated as ( Pr[refBase|HET] = %lf, Pr[refBase|ALT] = %lf) with LLK = %lf",pRefHet,pRefAlt,myMinimizer.fmin);
      //vbid.setRefBiasParams(1.0, pRefHet, pRefAlt);
      
      vbid.selfOut.llk0 = myFunc.llk0;
      vbid.selfOut.llk1 = myFunc.llk1;
      vbid.selfOut.refHet = myFunc.pRefHet;
      vbid.selfOut.refAlt = myFunc.pRefAlt;
      vbid.calculateDepthByGenotype(0,vbid.selfOut);
    }
    else {
      warning("--self must be set for --chip-refBias to work. Skipping..");
    }
  }
  else if ( args.bChipFull ) {
    notice("Performing with-chip estimation of reference-bias and sample mixture together");
    double maxIBD = -1;
    
    for(int i=0; i < (int)vbid.pGenotypes->indids.size(); ++i) {
      mpuVerify::fullIbdLLKFunc myFunc(&vbid,i);
      AmoebaMinimizer myMinimizer;
      Vector startingPoint(3);
      startingPoint[0] = 3.91;  // start with fIBD = 0.99
      startingPoint[1] = 0;      // pRefHet = 0.5
      startingPoint[2] = -4.595; // pRefAlt = 0.01
      myMinimizer.func = &myFunc;
      
      myFunc.indIdx = i;
      myMinimizer.Reset(3);
      myMinimizer.point = startingPoint;
      myMinimizer.Minimize(1e-6);
      double fIBD = mpuVerify::invLogit(myMinimizer.point[0]);
      double pRefHet = mpuVerify::invLogit(myMinimizer.point[1]);
      double pRefAlt = mpuVerify::invLogit(myMinimizer.point[2]);
      
      notice("Comparing with individual %s.. Optimal fIBD = %lf, LLK0 = %lf, LLK1 = %l",vbid.pGenotypes->indids[i].c_str(), fIBD, myFunc.llk0, myFunc.llk1);
      notice("Reference Bias Estimated as ( Pr[refBase|HET] = %lf, Pr[refBase|ALT] = %lf ) with LLK = %lf",pRefHet,pRefAlt,myMinimizer.fmin);
      if ( maxIBD < fIBD ) {
	bestInd = i;
	maxIBD = fIBD;
	vbid.bestOut.llk0 = myFunc.llk0;
	vbid.bestOut.llk1 = myFunc.llk1;
	vbid.bestOut.fMix = 1.-myFunc.fIBD;
	vbid.bestOut.refHet = myFunc.pRefHet;
	vbid.bestOut.refAlt = myFunc.pRefAlt;
      }

      if (vbid.pPile->sSM == vbid.pGenotypes->indids[i]) {
	selfInd = i;
	vbid.selfOut.llk0 = myFunc.llk0;
	vbid.selfOut.llk1 = myFunc.llk1;
	vbid.selfOut.fMix = 1.-myFunc.fIBD;
	vbid.selfOut.refHet = myFunc.pRefHet;
	vbid.selfOut.refAlt = myFunc.pRefAlt;
	vbid.calculateDepthByGenotype(i, vbid.selfOut);
      }
    }
    if ( bestInd >= 0 ) {
      notice("Best Matching Individual is %s with IBD = %lf",vbid.pGenotypes->indids[bestInd].c_str(),maxIBD);
      vbid.calculateDepthByGenotype(bestInd, vbid.bestOut);
    }
    
    if ( selfInd >= 0 ) {
      notice("Self Individual is %s with IBD = %lf",vbid.pGenotypes->indids[selfInd].c_str(),vbid.selfOut.fMix);
      vbid.calculateDepthByGenotype(selfInd,vbid.selfOut);
    }
  }

  // PRINT OUTPUT FILE - ".selfSM"
  // [SEQ_ID]  : SAMPLE ID in the sequence file
  // [CHIP_ID] : SAMPLE ID in the chip file (NA if not available)
  // [#SNPS] : Number of markers evaluated
  // [#READS]   : Number of reads evaluated
  // [AVG_DP]   : Mean depth
  // [FREEMIX]  : Chip-free estimated alpha (% MIX in 0-1 scale), NA if unavailable
  // [FREELK1]  : Chip-free log-likelihood at estimated alpha
  // [FREELK0]  : Chip-free log-likelihood at 0% contamination
  // [CHIPIBD]  : With-chip estimated alpha (% MIX in 0-1 scale)
  // [CHIPLK1]  : With-chip log-likelihood at estimated alpha
  // [CHIPLK0]  : With-chip log-likelihood at 0% contamination
  // [DPREF]    : Depth at reference site in the chip
  // [RDPHET]   : Relative depth at HET site in the chip
  // [RDPALT]   : Relative depth at HOMALT site in the chip
  // [FREE_RF]  : Pr(Ref|Ref) site estimated without chip data
  // [FREE_RH]  : Pr(Ref|Het) site estimated without chip data
  // [FREE_RA]  : Pr(Ref|Alt) site estimated without chip data
  // [CHIP_RF]  : Pr(Ref|Ref) site estimated with chip data
  // [CHIP_RH]  : Pr(Ref|Het) site estimated with chip data
  // [CHIP_RA]  : Pr(Ref|Alt) site estimated with chip data
  // [DPREF]    : Depth at reference alleles
  // [RDPHET]   : Relative depth at heterozygous alleles
  // [RDPALT]   : Relative depth at hom-alt alleles

  std::string soutfile = std::string(args.sOutFile.c_str());
  std::string selfSMFN = soutfile + ".selfSM";
  std::string bestSMFN = soutfile + ".bestSM";
  std::string dpSMFN = soutfile + ".depthSM";

  fprintf(stderr,"foo %s\n",selfSMFN.c_str());
  fprintf(stderr,"right before calling ifopen(%s)\n",selfSMFN.c_str());

  wFile selfSMF(selfSMFN.c_str());
  wFile bestSMF;
  if ( args.bFindBest ) bestSMF.open(bestSMFN.c_str());

  wFile dpSMF(dpSMFN.c_str());

  dpSMF.printf("#RG\tDEPTH\t#SNPs\t%%SNPs\t%%CUMUL\n");
  int nCumMarkers = 0;
  for(int i=args.maxDepth; i >= 0; --i) {
    nCumMarkers += vbid.mixOut.depths[i];
    dpSMF.printf("ALL\t%d\t%d\t%.5lf\t%.5lf\n",i, vbid.mixOut.depths[i],(double) vbid.mixOut.depths[i]/(double)vbid.nMarkers,(double)nCumMarkers/(double)vbid.nMarkers);
  }
  dpSMF.close();

  const char* headers[] = {"#SEQ_ID","RG","CHIP_ID","#SNPS","#READS","AVG_DP","FREEMIX","FREELK1","FREELK0","FREE_RH","FREE_RA","CHIPMIX","CHIPLK1","CHIPLK0","CHIP_RH","CHIP_RA","DPREF","RDPHET","RDPALT"};
  int nheaders = sizeof(headers)/sizeof(headers[0]);

  for(int i=0; i < nheaders; ++i) { selfSMF.printf("%s%s",i>0 ? "\t" : "",headers[i]); }
  selfSMF.printf("\n");
  selfSMF.printf("%s\tALL",vbid.pPile->sSM.c_str());
  selfSMF.printf("\t%s",selfInd >= 0 ? vbid.pGenotypes->indids[selfInd].c_str() : "NA");
  selfSMF.printf("\t%d\t%d\t%.2lf",vbid.nMarkers,vbid.mixOut.numReads[0],(double)vbid.mixOut.numReads[0]/(double)vbid.nMarkers);
  if ( args.bFreeNone ) { selfSMF.printf("\tNA\tNA\tNA\tNA\tNA"); }
  else if ( args.bFreeMixOnly ) { selfSMF.printf("\t%.5lf\t%.2lf\t%.2lf\tNA\tNA",vbid.mixOut.fMix,vbid.mixOut.llk1,vbid.mixOut.llk0); }
  else if ( args.bFreeRefBiasOnly ) { selfSMF.printf("\tNA\t%.2lf\t%.2lf\t%.5lf\t%.5lf",vbid.mixOut.llk1,vbid.mixOut.llk0,vbid.mixOut.refHet,vbid.mixOut.refAlt); }
  else if ( args.bFreeFull ) { selfSMF.printf("\t%.5lf\t%.2lf\t%.2lf\t%.5lf\t%.5lf",vbid.mixOut.fMix,vbid.mixOut.llk1,vbid.mixOut.llk0,vbid.mixOut.refHet,vbid.mixOut.refAlt); }
  else { error("Invalid option in handling bFree"); }

  if ( args.bChipNone || bestInd < 0 ) { selfSMF.printf("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"); }
  else if ( args.bChipMixOnly ) { selfSMF.printf("\t%.5lf\t%.2lf\t%.2lf\tNA\tNA\t%.3lf\t%.4lf\t%.4lf",vbid.selfOut.fMix,vbid.selfOut.llk1,vbid.selfOut.llk0,(double)vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[1], (double)vbid.selfOut.numReads[2]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[2], (double)vbid.selfOut.numReads[3]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[3]); }
  else if ( args.bChipMixOnly ) { selfSMF.printf("\tNA\t%.2lf\t%.2lf\t%.5lf\t%.5lf\t%.3lf\t%.4lf\t%.4lf",vbid.selfOut.llk1, vbid.selfOut.llk0, vbid.selfOut.refHet, vbid.selfOut.refAlt, (double)vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[1], (double)vbid.selfOut.numReads[2]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[2], (double)vbid.selfOut.numReads[3]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[3]); }
  else if ( args.bChipFull ) { selfSMF.printf("\t%.5lf\t%.2lf\t%.2lf\t%.5lf\t%.5lf\t%.3lf\t%.4lf\t%.4lf", vbid.selfOut.fMix, vbid.selfOut.llk1, vbid.selfOut.llk0, vbid.selfOut.refHet, vbid.selfOut.refAlt, (double)vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[1], (double)vbid.selfOut.numReads[2]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[2], (double)vbid.selfOut.numReads[3]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[3]); }
  else { error("Invalid option in handling bChip"); }
  selfSMF.printf("\n");
  selfSMF.close();

  if ( args.bFindBest ) {
    for(int i=0; i < nheaders; ++i) { bestSMF.printf("%s%s",i>0 ? "\t" : "",headers[i]); }
    bestSMF.printf("\n");
    bestSMF.printf("%s\tALL",vbid.pPile->sSM.c_str());
    bestSMF.printf("\t%s",bestInd >= 0 ? vbid.pGenotypes->indids[bestInd].c_str() : "NA");
    bestSMF.printf("\t%d\t%d\t%.2lf",vbid.nMarkers,vbid.mixOut.numReads[0],(double)vbid.mixOut.numReads[0]/(double)vbid.nMarkers);
    if ( args.bFreeNone ) { bestSMF.printf("\tNA\tNA\tNA\tNA\tNA"); }
    else if ( args.bFreeMixOnly ) { bestSMF.printf("\t%.5lf\t%.2lf\t%.2lf\tNA\tNA",vbid.mixOut.fMix,vbid.mixOut.llk1,vbid.mixOut.llk0); }
    else if ( args.bFreeRefBiasOnly ) { bestSMF.printf("\tNA\t%.2lf\t%.2lf\t%.5lf\t%.5lf",vbid.mixOut.llk1,vbid.mixOut.llk0,vbid.mixOut.refHet,vbid.mixOut.refAlt); }
    else if ( args.bFreeFull ) { bestSMF.printf("\t%.5lf\t%.2lf\t%.2lf\t%.5lf\t%.5lf",vbid.mixOut.fMix,vbid.mixOut.llk1,vbid.mixOut.llk0,vbid.mixOut.refHet,vbid.mixOut.refAlt); }
    else { error("Invalid option in handling bFree"); }
    
    if ( args.bChipNone || bestInd < 0 ) { bestSMF.printf("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"); }
    else if ( args.bChipMixOnly ) { bestSMF.printf("\t%.5lf\t%.2lf\t%.2lf\tNA\tNA\t%.3lf\t%.4lf\t%.4lf",vbid.bestOut.fMix,vbid.bestOut.llk1,vbid.bestOut.llk0,(double)vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[1], (double)vbid.bestOut.numReads[2]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[2], (double)vbid.bestOut.numReads[3]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[3]); }
    else if ( args.bChipMixOnly ) { bestSMF.printf("\tNA\t%.2lf\t%.2lf\t%.5lf\t%.5lf\t%.3lf\t%.4lf\t%.4lf",vbid.bestOut.llk1, vbid.bestOut.llk0, vbid.bestOut.refHet, vbid.bestOut.refAlt, (double)vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[1], (double)vbid.bestOut.numReads[2]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[2], (double)vbid.bestOut.numReads[3]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[3]); }
    else if ( args.bChipFull ) { bestSMF.printf("\t%.5lf\t%.2lf\t%.2lf\t%.5lf\t%.5lf\t%.3lf\t%.4lf\t%.4lf", vbid.bestOut.fMix, vbid.bestOut.llk1, vbid.bestOut.llk0, vbid.bestOut.refHet, vbid.bestOut.refAlt, (double)vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[1], (double)vbid.bestOut.numReads[2]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[2], (double)vbid.bestOut.numReads[3]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[3]); }
    else { error("Invalid option in handling bChip"); }
    bestSMF.printf("\n");
    bestSMF.close();
  }

  time(&t);
  notice("Analysis finished on %s",ctime(&t));

  return 0;
}
*/

/*
// new version of verifyBamID to detect contamination by population
int runLikelihoodPCA(int argc, char** argv) {
  printf("mpuInferPopulationAndSex 1.0.0 -- \n"
        "(c) 2013 Hyun Min Kang, Goo Jun, and Goncalo Abecasis\n\n");

  mpuVerifyArgs args;
  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Files")
    LONG_STRINGPARAMETER("vcf",&args.sVcfFile)
    LONG_STRINGPARAMETER("mpu",&args.sMpuFile)
    LONG_STRINGPARAMETER("subset",&args.sSubsetInds)
    LONG_STRINGPARAMETER("smID",&args.sSMID)
    LONG_STRINGPARAMETER("scanPair",&args.sScanPair)

    LONG_PARAMETER_GROUP("VCF analysis options")
    LONG_DOUBLEPARAMETER("genoError",&args.genoError)
    LONG_DOUBLEPARAMETER("minAF",&args.minAF)
    LONG_DOUBLEPARAMETER("minCallRate",&args.minCallRate)

    LONG_PARAMETER_GROUP("Individuals to compare with chip data")
    EXCLUSIVE_PARAMETER("site",&args.bSiteOnly)
    EXCLUSIVE_PARAMETER("self",&args.bSelfOnly)
    EXCLUSIVE_PARAMETER("best",&args.bFindBest)

    LONG_PARAMETER_GROUP("Chip-free optimization options")
    EXCLUSIVE_PARAMETER("free-none",&args.bFreeNone)
    EXCLUSIVE_PARAMETER("free-mix",&args.bFreeMixOnly)
    EXCLUSIVE_PARAMETER("free-refBias",&args.bFreeRefBiasOnly)
    EXCLUSIVE_PARAMETER("free-full",&args.bFreeFull)

    LONG_PARAMETER_GROUP("With-chip optimization options")
    EXCLUSIVE_PARAMETER("chip-none",&args.bChipNone)
    EXCLUSIVE_PARAMETER("chip-mix",&args.bChipMixOnly)
    EXCLUSIVE_PARAMETER("chip-refBias",&args.bChipRefBiasOnly)
    EXCLUSIVE_PARAMETER("chip-full",&args.bChipFull)

    LONG_PARAMETER_GROUP("MPU analysis options")
    LONG_PARAMETER("precise",&args.bPrecise)
    LONG_PARAMETER("ignore-sex",&args.bIgnoreSex)
    LONG_INTPARAMETER("minMapQ",&args.minMapQ)
    LONG_INTPARAMETER("maxDepth",&args.maxDepth)
    LONG_INTPARAMETER("minQ",&args.minQ)
    LONG_INTPARAMETER("maxQ",&args.maxQ)
    LONG_DOUBLEPARAMETER("grid",&args.grid)

    LONG_PARAMETER_GROUP("Modeling Reference Bias")
    LONG_DOUBLEPARAMETER("refRef",&args.pRefRef)
    LONG_DOUBLEPARAMETER("refHet",&args.pRefHet)
    LONG_DOUBLEPARAMETER("refAlt",&args.pRefAlt)

    LONG_PARAMETER_GROUP("Output options")
    LONG_STRINGPARAMETER("out",&args.sOutFile)
    LONG_PARAMETER("verbose",&args.bVerbose)
  END_LONG_PARAMETERS();

  pl.Add(new LongParameters("Available Options",longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // check the validity of input files
  if ( args.sVcfFile.empty() ) {
    error("--vcf [vcf file] required");
  }

  if ( args.sMpuFile.empty() ) {
    error("--mpu [mpu file] is required");
  }

  if ( args.sOutFile.empty() ) {
    error("--out [output prefix] is required");
  }

  if ( args.sSMID.empty() ) {
    error("--smID [sample ID] is required");
  }

  if ( ! ( args.bSiteOnly || args.bSelfOnly || args.bFindBest ) ) {
    warning("--self option was autotomatically turned on by default. Specify --best option if you wanted to check across all possible samples in the VCF");
    args.bSelfOnly = true;
  }

  if ( !( args.bFreeNone || args.bFreeMixOnly || args.bFreeRefBiasOnly || args.bFreeFull ) ) {
    warning("--free-mix option is automatically turned on by default");
    args.bFreeMixOnly = true;
  }

  if ( !( args.bChipNone || args.bChipMixOnly || args.bChipRefBiasOnly || args.bChipFull ) ) {
    warning("--chip-mix option is automatically turned on by default");
    args.bChipMixOnly = true;
  }

  if ( ( args.maxDepth > 20 ) && ( !args.bPrecise ) ) {
    warning("--precise option is not turned on at --maxDepth %d : may be prone to precision errors",args.maxDepth);
  }

  if ( ( args.bChipRefBiasOnly ) && ( !args.bSelfOnly ) ) {
    error("--self must be set for --chip-refBias to work. Skipping..");
  }

  if ( ( !args.sScanPair.empty() ) && ( ( !args.bFreeNone ) || ( !args.bFindBest ) || ( !args.bChipMixOnly ) ) ) {
    error("--scanPair option is only compatible with --best, --free-none, and --chip-mix");
  }

  // check timestamp
  time_t t;
  time(&t);
  notice("Analysis started on %s",ctime(&t));

  // load arguments
  mpuVerify vbid(&args);

  // load input VCF and BAM files
  notice("Opening Input Files");
  vbid.loadFiles(args.sMpuFile.c_str(), args.sVcfFile.c_str());

  int selfIndex = -1;
  if ( !args.sScanPair.empty() ) { // if scanPair is on
    std::string indID1(args.sScanPair.c_str());
    for(int i=0; i < (int)vbid.pGenotypes->indids.size(); ++i) {
      if ( vbid.pGenotypes->indids[i] == indID1 ) {
	selfIndex = i;
      }
    }
    if ( selfIndex < 0 ) {
      error("Cannot find individual %s from the VCF file",indID1.c_str());
    }
  }

  // Check which genotype-free method is used
  if ( args.bFreeNone ) {  // if no genotype-free mode is tested. skip it
    // do nothing for genotype-free estimation
    notice("Skipping chip-free estimation of sample mixture");
  }
  else if ( args.bFreeMixOnly ) { // only mixture is estimated.
    // genotype-free method
    notice("Performing chip-free estimation of sample mixture at fixed reference bias parameters (%lf, %lf, %lf)",args.pRefRef,args.pRefHet,args.pRefAlt);

    mpuVerify::mixLLK mix(&vbid);
    mix.OptimizeLLK();
    notice("Optimal per-sample fMix = %lf, LLK0 = %lf, LLK1 = %lf\n",mix.fMix,mix.llk0,mix.llk1);
    vbid.mixOut.llk0 = mix.llk0;
    vbid.mixOut.llk1 = mix.llk1;
    vbid.mixOut.fMix = mix.fMix;
  }
  else if ( args.bFreeRefBiasOnly ) {
    notice("Performing chip-free estimation of reference-bias without sample mixture");
    mpuVerify::refBiasMixLLKFunc myFunc(&vbid);
    AmoebaMinimizer myMinimizer;
    Vector startingPoint(2);
    startingPoint[0] = 0;      // pRefHet = 0.5
    startingPoint[1] = -4.595; // pRefAlt = 0.01
    myMinimizer.func = &myFunc;
    myMinimizer.Reset(2);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(1e-6);
    double pRefHet = mpuVerify::invLogit(myMinimizer.point[0]);
    double pRefAlt = mpuVerify::invLogit(myMinimizer.point[1]);
    notice("Reference Bias Estimated as ( Pr[refBase|HET] = %lf, Pr[refBase|ALT] = %lf) with LLK = %lf",pRefHet,pRefAlt,myMinimizer.fmin);
    vbid.mixOut.llk0 = myFunc.llk0;
    vbid.mixOut.llk1 = myFunc.llk1;
    vbid.mixOut.refHet = myFunc.pRefHet;
    vbid.mixOut.refAlt = myFunc.pRefAlt;
  }
  else if ( args.bFreeFull ) {
    notice("Performing chip-free estimation of reference-bias and sample mixture together");
    mpuVerify::fullMixLLKFunc myFunc(&vbid);
    AmoebaMinimizer myMinimizer;
    Vector startingPoint(3);
    startingPoint[0] = -3.91;  // start with fMix = 0.01
    startingPoint[1] = 0;      // pRefHet = 0.5
    startingPoint[2] = -4.595; // pRefAlt = 0.01
    myMinimizer.func = &myFunc;
    myMinimizer.Reset(3);
    myMinimizer.point = startingPoint;
    myMinimizer.Minimize(1e-6);
    double fMix = mpuVerify::invLogit(myMinimizer.point[0]);
    if ( fMix > 0.5 ) 
      fMix = 1.-fMix;
    double pRefHet = mpuVerify::invLogit(myMinimizer.point[1]);
    double pRefAlt = mpuVerify::invLogit(myMinimizer.point[2]);
    notice("Optimal per-sample fMix = %lf\n",fMix);
    notice("Reference Bias Estimated as ( Pr[refBase|HET] = %lf, Pr[refBase|ALT] = %lf) with LLK = %lf",pRefHet,pRefAlt,myMinimizer.fmin);

    vbid.mixOut.llk0 = myFunc.llk0;
    vbid.mixOut.llk1 = myFunc.llk1;
    vbid.mixOut.fMix = myFunc.fMix;
    vbid.mixOut.refHet = myFunc.pRefHet;
    vbid.mixOut.refAlt = myFunc.pRefAlt;
  }
  notice("calculating depth distribution");  
  vbid.calculateDepthDistribution(args.maxDepth, vbid.mixOut);

  notice("finished calculating depth distribution");  

  int bestInd = -1;
  int selfInd = -1;

  if ( args.bChipNone ) {
    // do nothing
    notice("Skipping with-chip estimation of sample mixture");
  }
  else if ( args.bChipMixOnly ) {
    notice("Performing with-chip estimation of sample mixture at fixed reference bias parameter (%lf, %lf, %lf)",args.pRefRef,args.pRefHet,args.pRefAlt);

    if ( args.sScanPair.empty() ) {
      double maxIBD = -1;
      mpuVerify::ibdLLK ibd(&vbid);
      for(int i=0; i < (int)vbid.pGenotypes->indids.size(); ++i) {
	double fIBD = ibd.OptimizeLLK(i);
	notice("Comparing with individual %s.. Optimal fIBD = %lf, LLK0 = %lf, LLK1 = %lf",vbid.pGenotypes->indids[i].c_str(),fIBD, ibd.llk0, ibd.llk1);
	if ( maxIBD < fIBD ) {
	  bestInd = i;
	  vbid.bestOut.llk0 = ibd.llk0;
	  vbid.bestOut.llk1 = ibd.llk1;
	  vbid.bestOut.fMix = 1-ibd.fIBD;
	  maxIBD = ibd.fIBD;
	}
	  
	if ( vbid.pPile->sSM == vbid.pGenotypes->indids[i] ) {
	  selfInd = i;
	  vbid.selfOut.llk0 = ibd.llk0;
	  vbid.selfOut.llk1 = ibd.llk1;
	  vbid.selfOut.fMix = 1-ibd.fIBD;
	}
      }
	
      if ( bestInd >= 0 ) {
	notice("Best Matching Individual is %s with IBD = %lf",vbid.pGenotypes->indids[bestInd].c_str(),maxIBD);
	vbid.calculateDepthByGenotype(bestInd,vbid.bestOut);
      }
      
      if ( selfInd >= 0 ) {
	notice("Self Individual is %s with IBD = %lf",vbid.pGenotypes->indids[selfInd].c_str(),vbid.selfOut.fMix);
	vbid.calculateDepthByGenotype(selfInd,vbid.selfOut);
      }
    }
    else {
      double maxLLK = 1e9;
      mpuVerify::pairLLK ibdPair(&vbid);
      for(int i=0; i < (int)vbid.pGenotypes->indids.size(); ++i) {
	double fIBD = ibdPair.OptimizeLLK(selfIndex, i);
	notice("Comparing with individual %s.. Optimal fIBD = %lf, LLK0 = %lf, LLK1 = %lf",vbid.pGenotypes->indids[i].c_str(),fIBD, ibdPair.llk0, ibdPair.llk1);
	double fLLK = ibdPair.llk1;
	if ( fLLK < maxLLK ) {
	  bestInd = i;
	  vbid.bestOut.llk0 = ibdPair.llk0;
	  vbid.bestOut.llk1 = ibdPair.llk1;
	  vbid.bestOut.fMix = 1-ibdPair.fIBD;
	  maxLLK = fLLK;
	}
	
	if (vbid.pPile->sSM == vbid.pGenotypes->indids[i] ) {
	  selfInd = i;
	  vbid.selfOut.llk0 = ibdPair.llk0;
	  vbid.selfOut.llk1 = ibdPair.llk1;
	  vbid.selfOut.fMix = 1-ibdPair.fIBD;
	}
      }
      
      if ( bestInd >= 0 ) {
	notice("Best Matching Individual is %s with LLK = %lf",vbid.pGenotypes->indids[bestInd].c_str(),maxLLK);
	vbid.calculateDepthByGenotype(bestInd,vbid.bestOut);
      }
      
      if ( selfInd >= 0 ) {
	notice("Self Individual is %s with IBD = %lf",vbid.pGenotypes->indids[selfInd].c_str(),vbid.selfOut.fMix);
	vbid.calculateDepthByGenotype(selfInd,vbid.selfOut);
      }
    }
  }
  else if ( args.bChipRefBiasOnly ) {
    notice("Performing with-chip estimation of reference-bias without sample mixture");
    if ( args.bSelfOnly ) {
      mpuVerify::refBiasIbdLLKFunc myFunc(&vbid);
      AmoebaMinimizer myMinimizer;
      Vector startingPoint(2);
      startingPoint[0] = 0;      // pRefHet = 0.5
      startingPoint[1] = -4.595; // pRefAlt = 0.01
      myMinimizer.func = &myFunc;
      myMinimizer.Reset(2);
      myMinimizer.point = startingPoint;
      myMinimizer.Minimize(1e-6);
      double pRefHet = mpuVerify::invLogit(myMinimizer.point[0]);
      double pRefAlt = mpuVerify::invLogit(myMinimizer.point[1]);
      notice("Reference Bias Estimated as ( Pr[refBase|HET] = %lf, Pr[refBase|ALT] = %lf) with LLK = %lf",pRefHet,pRefAlt,myMinimizer.fmin);
      //vbid.setRefBiasParams(1.0, pRefHet, pRefAlt);
      
      vbid.selfOut.llk0 = myFunc.llk0;
      vbid.selfOut.llk1 = myFunc.llk1;
      vbid.selfOut.refHet = myFunc.pRefHet;
      vbid.selfOut.refAlt = myFunc.pRefAlt;
      vbid.calculateDepthByGenotype(0,vbid.selfOut);
    }
    else {
      warning("--self must be set for --chip-refBias to work. Skipping..");
    }
  }
  else if ( args.bChipFull ) {
    notice("Performing with-chip estimation of reference-bias and sample mixture together");
    double maxIBD = -1;
    
    for(int i=0; i < (int)vbid.pGenotypes->indids.size(); ++i) {
      mpuVerify::fullIbdLLKFunc myFunc(&vbid,i);
      AmoebaMinimizer myMinimizer;
      Vector startingPoint(3);
      startingPoint[0] = 3.91;  // start with fIBD = 0.99
      startingPoint[1] = 0;      // pRefHet = 0.5
      startingPoint[2] = -4.595; // pRefAlt = 0.01
      myMinimizer.func = &myFunc;
      
      myFunc.indIdx = i;
      myMinimizer.Reset(3);
      myMinimizer.point = startingPoint;
      myMinimizer.Minimize(1e-6);
      double fIBD = mpuVerify::invLogit(myMinimizer.point[0]);
      double pRefHet = mpuVerify::invLogit(myMinimizer.point[1]);
      double pRefAlt = mpuVerify::invLogit(myMinimizer.point[2]);
      
      notice("Comparing with individual %s.. Optimal fIBD = %lf, LLK0 = %lf, LLK1 = %l",vbid.pGenotypes->indids[i].c_str(), fIBD, myFunc.llk0, myFunc.llk1);
      notice("Reference Bias Estimated as ( Pr[refBase|HET] = %lf, Pr[refBase|ALT] = %lf ) with LLK = %lf",pRefHet,pRefAlt,myMinimizer.fmin);
      if ( maxIBD < fIBD ) {
	bestInd = i;
	maxIBD = fIBD;
	vbid.bestOut.llk0 = myFunc.llk0;
	vbid.bestOut.llk1 = myFunc.llk1;
	vbid.bestOut.fMix = 1.-myFunc.fIBD;
	vbid.bestOut.refHet = myFunc.pRefHet;
	vbid.bestOut.refAlt = myFunc.pRefAlt;
      }

      if (vbid.pPile->sSM == vbid.pGenotypes->indids[i]) {
	selfInd = i;
	vbid.selfOut.llk0 = myFunc.llk0;
	vbid.selfOut.llk1 = myFunc.llk1;
	vbid.selfOut.fMix = 1.-myFunc.fIBD;
	vbid.selfOut.refHet = myFunc.pRefHet;
	vbid.selfOut.refAlt = myFunc.pRefAlt;
	vbid.calculateDepthByGenotype(i, vbid.selfOut);
      }
    }
    if ( bestInd >= 0 ) {
      notice("Best Matching Individual is %s with IBD = %lf",vbid.pGenotypes->indids[bestInd].c_str(),maxIBD);
      vbid.calculateDepthByGenotype(bestInd, vbid.bestOut);
    }
    
    if ( selfInd >= 0 ) {
      notice("Self Individual is %s with IBD = %lf",vbid.pGenotypes->indids[selfInd].c_str(),vbid.selfOut.fMix);
      vbid.calculateDepthByGenotype(selfInd,vbid.selfOut);
    }
  }

  // PRINT OUTPUT FILE - ".selfSM"
  // [SEQ_ID]  : SAMPLE ID in the sequence file
  // [CHIP_ID] : SAMPLE ID in the chip file (NA if not available)
  // [#SNPS] : Number of markers evaluated
  // [#READS]   : Number of reads evaluated
  // [AVG_DP]   : Mean depth
  // [FREEMIX]  : Chip-free estimated alpha (% MIX in 0-1 scale), NA if unavailable
  // [FREELK1]  : Chip-free log-likelihood at estimated alpha
  // [FREELK0]  : Chip-free log-likelihood at 0% contamination
  // [CHIPIBD]  : With-chip estimated alpha (% MIX in 0-1 scale)
  // [CHIPLK1]  : With-chip log-likelihood at estimated alpha
  // [CHIPLK0]  : With-chip log-likelihood at 0% contamination
  // [DPREF]    : Depth at reference site in the chip
  // [RDPHET]   : Relative depth at HET site in the chip
  // [RDPALT]   : Relative depth at HOMALT site in the chip
  // [FREE_RF]  : Pr(Ref|Ref) site estimated without chip data
  // [FREE_RH]  : Pr(Ref|Het) site estimated without chip data
  // [FREE_RA]  : Pr(Ref|Alt) site estimated without chip data
  // [CHIP_RF]  : Pr(Ref|Ref) site estimated with chip data
  // [CHIP_RH]  : Pr(Ref|Het) site estimated with chip data
  // [CHIP_RA]  : Pr(Ref|Alt) site estimated with chip data
  // [DPREF]    : Depth at reference alleles
  // [RDPHET]   : Relative depth at heterozygous alleles
  // [RDPALT]   : Relative depth at hom-alt alleles

  std::string soutfile = std::string(args.sOutFile.c_str());
  std::string selfSMFN = soutfile + ".selfSM";
  std::string bestSMFN = soutfile + ".bestSM";
  std::string dpSMFN = soutfile + ".depthSM";

  fprintf(stderr,"foo %s\n",selfSMFN.c_str());
  fprintf(stderr,"right before calling ifopen(%s)\n",selfSMFN.c_str());

  wFile selfSMF(selfSMFN.c_str());
  wFile bestSMF;
  if ( args.bFindBest ) bestSMF.open(bestSMFN.c_str());

  wFile dpSMF(dpSMFN.c_str());

  dpSMF.printf("#RG\tDEPTH\t#SNPs\t%%SNPs\t%%CUMUL\n");
  int nCumMarkers = 0;
  for(int i=args.maxDepth; i >= 0; --i) {
    nCumMarkers += vbid.mixOut.depths[i];
    dpSMF.printf("ALL\t%d\t%d\t%.5lf\t%.5lf\n",i, vbid.mixOut.depths[i],(double) vbid.mixOut.depths[i]/(double)vbid.nMarkers,(double)nCumMarkers/(double)vbid.nMarkers);
  }
  dpSMF.close();

  const char* headers[] = {"#SEQ_ID","RG","CHIP_ID","#SNPS","#READS","AVG_DP","FREEMIX","FREELK1","FREELK0","FREE_RH","FREE_RA","CHIPMIX","CHIPLK1","CHIPLK0","CHIP_RH","CHIP_RA","DPREF","RDPHET","RDPALT"};
  int nheaders = sizeof(headers)/sizeof(headers[0]);

  for(int i=0; i < nheaders; ++i) { selfSMF.printf("%s%s",i>0 ? "\t" : "",headers[i]); }
  selfSMF.printf("\n");
  selfSMF.printf("%s\tALL",vbid.pPile->sSM.c_str());
  selfSMF.printf("\t%s",selfInd >= 0 ? vbid.pGenotypes->indids[selfInd].c_str() : "NA");
  selfSMF.printf("\t%d\t%d\t%.2lf",vbid.nMarkers,vbid.mixOut.numReads[0],(double)vbid.mixOut.numReads[0]/(double)vbid.nMarkers);
  if ( args.bFreeNone ) { selfSMF.printf("\tNA\tNA\tNA\tNA\tNA"); }
  else if ( args.bFreeMixOnly ) { selfSMF.printf("\t%.5lf\t%.2lf\t%.2lf\tNA\tNA",vbid.mixOut.fMix,vbid.mixOut.llk1,vbid.mixOut.llk0); }
  else if ( args.bFreeRefBiasOnly ) { selfSMF.printf("\tNA\t%.2lf\t%.2lf\t%.5lf\t%.5lf",vbid.mixOut.llk1,vbid.mixOut.llk0,vbid.mixOut.refHet,vbid.mixOut.refAlt); }
  else if ( args.bFreeFull ) { selfSMF.printf("\t%.5lf\t%.2lf\t%.2lf\t%.5lf\t%.5lf",vbid.mixOut.fMix,vbid.mixOut.llk1,vbid.mixOut.llk0,vbid.mixOut.refHet,vbid.mixOut.refAlt); }
  else { error("Invalid option in handling bFree"); }

  if ( args.bChipNone || bestInd < 0 ) { selfSMF.printf("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"); }
  else if ( args.bChipMixOnly ) { selfSMF.printf("\t%.5lf\t%.2lf\t%.2lf\tNA\tNA\t%.3lf\t%.4lf\t%.4lf",vbid.selfOut.fMix,vbid.selfOut.llk1,vbid.selfOut.llk0,(double)vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[1], (double)vbid.selfOut.numReads[2]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[2], (double)vbid.selfOut.numReads[3]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[3]); }
  else if ( args.bChipMixOnly ) { selfSMF.printf("\tNA\t%.2lf\t%.2lf\t%.5lf\t%.5lf\t%.3lf\t%.4lf\t%.4lf",vbid.selfOut.llk1, vbid.selfOut.llk0, vbid.selfOut.refHet, vbid.selfOut.refAlt, (double)vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[1], (double)vbid.selfOut.numReads[2]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[2], (double)vbid.selfOut.numReads[3]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[3]); }
  else if ( args.bChipFull ) { selfSMF.printf("\t%.5lf\t%.2lf\t%.2lf\t%.5lf\t%.5lf\t%.3lf\t%.4lf\t%.4lf", vbid.selfOut.fMix, vbid.selfOut.llk1, vbid.selfOut.llk0, vbid.selfOut.refHet, vbid.selfOut.refAlt, (double)vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[1], (double)vbid.selfOut.numReads[2]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[2], (double)vbid.selfOut.numReads[3]*vbid.selfOut.numGenos[1]/vbid.selfOut.numReads[1]/vbid.selfOut.numGenos[3]); }
  else { error("Invalid option in handling bChip"); }
  selfSMF.printf("\n");
  selfSMF.close();

  if ( args.bFindBest ) {
    for(int i=0; i < nheaders; ++i) { bestSMF.printf("%s%s",i>0 ? "\t" : "",headers[i]); }
    bestSMF.printf("\n");
    bestSMF.printf("%s\tALL",vbid.pPile->sSM.c_str());
    bestSMF.printf("\t%s",bestInd >= 0 ? vbid.pGenotypes->indids[bestInd].c_str() : "NA");
    bestSMF.printf("\t%d\t%d\t%.2lf",vbid.nMarkers,vbid.mixOut.numReads[0],(double)vbid.mixOut.numReads[0]/(double)vbid.nMarkers);
    if ( args.bFreeNone ) { bestSMF.printf("\tNA\tNA\tNA\tNA\tNA"); }
    else if ( args.bFreeMixOnly ) { bestSMF.printf("\t%.5lf\t%.2lf\t%.2lf\tNA\tNA",vbid.mixOut.fMix,vbid.mixOut.llk1,vbid.mixOut.llk0); }
    else if ( args.bFreeRefBiasOnly ) { bestSMF.printf("\tNA\t%.2lf\t%.2lf\t%.5lf\t%.5lf",vbid.mixOut.llk1,vbid.mixOut.llk0,vbid.mixOut.refHet,vbid.mixOut.refAlt); }
    else if ( args.bFreeFull ) { bestSMF.printf("\t%.5lf\t%.2lf\t%.2lf\t%.5lf\t%.5lf",vbid.mixOut.fMix,vbid.mixOut.llk1,vbid.mixOut.llk0,vbid.mixOut.refHet,vbid.mixOut.refAlt); }
    else { error("Invalid option in handling bFree"); }
    
    if ( args.bChipNone || bestInd < 0 ) { bestSMF.printf("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"); }
    else if ( args.bChipMixOnly ) { bestSMF.printf("\t%.5lf\t%.2lf\t%.2lf\tNA\tNA\t%.3lf\t%.4lf\t%.4lf",vbid.bestOut.fMix,vbid.bestOut.llk1,vbid.bestOut.llk0,(double)vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[1], (double)vbid.bestOut.numReads[2]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[2], (double)vbid.bestOut.numReads[3]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[3]); }
    else if ( args.bChipMixOnly ) { bestSMF.printf("\tNA\t%.2lf\t%.2lf\t%.5lf\t%.5lf\t%.3lf\t%.4lf\t%.4lf",vbid.bestOut.llk1, vbid.bestOut.llk0, vbid.bestOut.refHet, vbid.bestOut.refAlt, (double)vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[1], (double)vbid.bestOut.numReads[2]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[2], (double)vbid.bestOut.numReads[3]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[3]); }
    else if ( args.bChipFull ) { bestSMF.printf("\t%.5lf\t%.2lf\t%.2lf\t%.5lf\t%.5lf\t%.3lf\t%.4lf\t%.4lf", vbid.bestOut.fMix, vbid.bestOut.llk1, vbid.bestOut.llk0, vbid.bestOut.refHet, vbid.bestOut.refAlt, (double)vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[1], (double)vbid.bestOut.numReads[2]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[2], (double)vbid.bestOut.numReads[3]*vbid.bestOut.numGenos[1]/vbid.bestOut.numReads[1]/vbid.bestOut.numGenos[3]); }
    else { error("Invalid option in handling bChip"); }
    bestSMF.printf("\n");
    bestSMF.close();
  }

  time(&t);
  notice("Analysis finished on %s",ctime(&t));

  return 0;
}
*/
/*
int main(int argc, char** argv) {
  if ( argc < 2 ) {
    printf("mpuTool -- Fast analytic tools for analyzing and manipulating pileups\n");
    printf("Copyright (c) 2013 Hyun Min Kang, Goo Jun, Matthew Flicinger, and Goncalo Abecasis\n");
    printf("Usage : %s [command] [options]\n",argv[0]);
    printf("\tType one of the following commands below to get detailed usage\n");
    printf("\t%s extract      [options] : extract GLs from a single VCF\n",argv[0]);
    printf("\t%s genotype     [options] : genotype across multiple pileups in known sites\n",argv[0]);
    printf("\t%s mixup        [options] : genotype across multiple mixed up pileups in known sites\n",argv[0]);
    printf("\t%s verify       [options] : detect contamination based on genotype data and population allele frequencies\n",argv[0]);
    printf("\t%s set-genotype [options] : genotype\n",argv[0]);
  }
  else {
    std::string cmd(argv[1]);
    //if ( cmd == "extract" ) {
      //return runExtract(argc-1,argv+1);
    //}
    if ( cmd == "genotype" ) {
      //return runGenotype(argc-1,argv+1);
      return runDeconvoluteGenotype(argc-1,argv+1);
    }
    //else if ( cmd == "mixup" ) {
      //return runMixGenotype(argc-1,argv+1);
    //}
    else if ( cmd == "verify" ) {
      return runVerify(argc-1,argv+1);
    }
    //else if ( cmd == "merge" ) {
    //return runMerge(argc-1,argv+1);
    //}
    //else if ( cmd == "set-genotype" ) {
      //return runSetGenotype(argc-1,argv+1);
    //}
    else {
      error("Unrecognized command %s\n",argv[0]);
    }
  }
  return 0;
}*/


#ifndef __TABIXED_VCF_H
#define __TABIXED_VCF_H

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cstdio>
#include <ctime>
//#include <cmath>
#include <math.h>
#include <set>
#include "pFile.h"
#include "wFile.h"
#include "Error.h"
#include "boolParser.h"
#include "PhredHelper.h"
#include "Constant.h"

class fVcf {
public:
  // member variables
  int nInds;           // number of individuals subselected
  int nMarkers;        // number of markers
  std::string key;     // key string
  boolParser parser;   // per-marker pattern matching parser
  bool passOnly;       // filtering option
  bool infoAF;         // extract site information only
  bool hardGenotype;   // hard/soft genotypes
  pFile tf;            // file handle
  std::vector<std::string> inds;     // individual IDs
  std::vector<std::string> markers;  // marker IDs
  std::vector<std::string> chroms;   // chromosome names
  std::vector<int> pos1s;            // marker positions
  std::vector<std::string> refs;     // reference alleles
  std::vector<std::string> alts;     // non-reference alleles
  std::vector<float> genos;          // GT : genotypes
  std::vector<uint8_t> phases;       // haplotype phases unphased/noinfo/0|1/1|0
  std::vector<double> GPs;           // GP : genotype probabilities
  std::vector<uint8_t> PLs;          // GL or PL : genotype likelihoods
  std::vector<int> depths;           // GD or DP : genotype depth
  std::vector<int> numAlleles;       // AN 
  std::vector<float> sumAlleles;     // AC
  std::vector<float> sumsqAlleles;   // squared sum of AC
  std::vector<int> icols;            // individual indices to subset
  std::vector<float> AFs;            // allele frequency estimates

  std::set<std::string> markerSet;
  std::string groupID, groupChrom;
  int groupBeg, groupEnd;
  std::vector<std::string> headers;

  bool glFlag;
  bool phredFlag;
  bool gpFlag;
  bool ignoreEmpty;
  std::vector<double> p2e;
  int nwarns;

  static const float NAN_FLT;
  static const double NAN_DBL;
  static std::string xLabel;
  static std::string yLabel;
  static std::string mitoLabel;
  static int xStart;
  static int xStop;

  static int getChromosomeType(const char* s, int pos) {
    if ( xLabel == s ) { 
      return ( ( pos >= xStart ) && ( pos <= xStop ) ) ? CT_CHRX : CT_AUTOSOME;
    }
    else if ( yLabel == s ) return CT_CHRY; 
    else if ( mitoLabel == s ) return CT_MITO;
    else return CT_AUTOSOME;
  }

  static int parseMarkerID(const char* markerID, std::string& chrom, std::string& beg, std::string& end, std::string& name, std::string& anno) {
    const char* pp = markerID;
    const char* pn = markerID;
    size_t found;
    int step = 0; // 0 : CHROM[:], 1 : BEG [-_], 2 : END [_], 3 : NAME [:]
    char buf[255];
    anno.clear();
    while( *pn != '\0' ) {
      switch(step) {
      case 0: // CHROM[:]
	if ( *pn == ':' ) { // END PARSING
	  chrom.assign( pp, pn - pp ); // copy chrom
	  step = 1;
	  pp = pn+1;
	}
	// otherwise, just skip else {}
	break;
      case 1: // BEG
	if ( *pn == '-' ) { // INTERVAL is given
	  beg.assign( pp, pn - pp );
	  step = 2;
	  pp = pn+1;
	}
	else if ( *pn == '_' ) { // BEG==END
	  beg.assign( pp, pn - pp );
	  end = beg;
	  step = 3;
	  pp = pn+1;
	}
	break;
      case 2: // END
	if ( *pn == '_' ) { // BEG==END
	  end.assign( pp, pn - pp );
	  step = 3;
	  pp = pn+1;
	}
	break;
      case 3: // NAME
	if ( *pn == '_' ) {
	  name.assign( pp, pn - pp );
	  step = 4;
	  pp = pn+1;
	}
	break;
      }
      ++pn;
    }

    if ( pn > pp ) {
      switch(step) {
      case 0:
	chrom.assign( pp, pn - pp );
	break;
      case 1:
	beg.assign( pp, pn - pp );
	end = beg;
	step = 2;
	break;
      case 2:
	end.assign( pp, pn - pp );
	break;
      case 3:
	name.assign( pp, pn - pp );
	found = name.find('/');
	if ( ( found != std::string::npos ) && ( found != 1 ) ) {
	  sprintf(buf, "%d", atoi(beg.c_str()) + (int)found - 1);
	  end = buf;
	}
	break;
      case 4:
	anno.assign( pp, pn - pp );
	break;
      }
    }
    return step;
  }

  void init(const char* vcf, const char* region = NULL, const char* _key = "GT", const char* rule = NULL, bool pass = true) {
    nInds = 0;
    nMarkers = 0;
    key = _key;
    parser = rule;
    passOnly = pass;
    hardGenotype = ( key == "GT" ? true : false);
    tf.load(vcf, region, true);

    if ( key == "GL" ) {
      glFlag = true;
      phredFlag = false;
    }
    else if ( key == "PL" ) {
      glFlag = phredFlag = true;
    }
    else {
      glFlag = phredFlag = false;
    }
    gpFlag = false;
  }

  void load(const char* vcf, const char* region, const char* _key, const char* rule, bool pass, const char* indf) {
    // read indf files if needed
    std::set<std::string> idset;
    char* line = NULL;
    if ( indf != NULL ) {
      pFile tind(indf);
      while( (line = (char*)tind.getLine()) != NULL ) {
	if ( line[0] == '#' ) continue;
	char* p = line;
	while( ( *p != ' ' ) && ( *p & 0xf0 ) ) { ++p; }
	std::string id(line,p-line);
	if ( idset.find(id) != idset.end() ) {
	  error("ERROR: Duplicate individual ID %s in %s",id.c_str(),indf);
	}
	idset.insert(id);
      }
    }
    load(vcf, region, _key, rule, pass, idset);
  }

  void load(const char* vcf, const char* region, const char* _key, const char* rule, bool pass, std::set<std::string>& idset) {
    init(vcf, region, _key, rule, pass);
    
    // read VCF header
    char* line;
    //int nc = 0;
    while ( (line = (char*)tf.getLine()) != NULL ) {
      //nc = tf.getLength();
      if ( line[0] == '#' ) {
	if ( strncmp(line,"#CHROM",6) == 0 ) {
	  nInds = parseInds(line, idset);
	  break;
	}
	else if ( line[1] == '#' ) {
	  // store meta lines;
	  headers.push_back(line);
	}
      }
      else {
	error("Non-header line %s is observed before #CHROM in %s",line,vcf);
      }
    }
   
    if ( ( idset.size() > 0 ) && ( idset.size() != icols.size() ) ) {
      warning("Identified %d individuals from index file, and only %d overlaps with VCF",(int)idset.size(),(int)icols.size());
    }
  }

  void load(const char* vcf, const char* region, const char* _key, const char* rule, bool pass, std::vector<int>& subcols) {
    init(vcf, region, _key, rule, pass);
    
    // read VCF header
    char* line;
    //int nc;
    while ( (line = (char*)tf.getLine()) != NULL ) {
      //nc = tf.getLength();
      if ( line[0] == '#' ) {
	if ( strncmp(line,"#CHROM",6) == 0 ) {
	  nInds = parseInds(line, subcols);
	  break;
	}
	else if ( line[1] == '#' ) {
	  // parse meta line???
	}
      }
      else {
	error("Non-header line %s is observed before #CHROM in %s",line,vcf);
      }
    }
    
    if ( ( subcols.size() > 0 ) && ( (int)subcols.size() != nInds ) ) {
      warning("Identified %d individuals from index file, and nInds = %d",(int)subcols.size(),nInds);
    }
  }

 fVcf() : nInds(0), nMarkers(0), passOnly(false), infoAF(false), ignoreEmpty(true), nwarns(0) {}

  void clear() {
    markers.clear();
    chroms.clear();
    pos1s.clear();
    refs.clear();
    alts.clear();
    genos.clear();
    phases.clear();
    GPs.clear();
    PLs.clear();
    depths.clear();
    numAlleles.clear();
    sumAlleles.clear();
    sumsqAlleles.clear();
    //icols.clear();
    AFs.clear();
    markerSet.clear();
    nMarkers = 0;
  }

  float callRate(int m) {
    return nInds > 0 ? ((float)numAlleles[m] / (float)nInds / 2.) : 1.;
  }

  float alleleVar(int m, bool allN = false) {
    if ( allN ) {
      return(sumsqAlleles[m]/nInds - sumAlleles[m]*sumAlleles[m]/nInds/nInds);
    }
    else if ( numAlleles[m] > 0 ) {
      return(sumsqAlleles[m]/numAlleles[m] - sumAlleles[m]*sumAlleles[m]/numAlleles[m]/numAlleles[m]);
    }
    else {
      return 1e-6;
    }
  }

  float alleleSD(int m, bool allN = false) {
    return sqrtf(alleleVar(m,allN));
  }

  float alleleFreq(int m) {
    if ( ( (int)AFs.size() > m ) && ( AFs[m] >= 0 ) ) {
      return AFs[m];
    }
    else if ( glFlag ) {
      if ( (int)AFs.size() < m + 1 ) {
	AFs.resize(m+1,-1);
      }
      if ( AFs[m] < 0 ) {
	AFs[m] = (float)(emAF(m,1e-6).first);
      }
      return AFs[m];
    }
    else {
      return ( numAlleles[m] > 0 ) ? (sumAlleles[m] / (float)numAlleles[m]) : 0;
    }
  }

  int MAC(int m) {
    return (int)floor(((sumAlleles[m] < numAlleles[m]-sumAlleles[m]) ? sumAlleles[m] : numAlleles[m]-sumAlleles[m])+.5);
  }

  // works only for discrete genotypes
  // numAlleles[m] = r + h + a
  // sumAlleles[m] = h + 2a
  // sumsqAlleles[m] = h + 4a
  // a = ((h+4a)-(h+2a))/2
  // r = (r+h+a) - (h+2a) + ((h+4a)-(h+2a))/2
  //   = (r+h+a) - 1.5(h+2a) + 0.5(h+4a)
  int HOMMINC(int m) {
    if ( hardGenotype ) {
      return ((int)floor( ((sumAlleles[m] < numAlleles[m]-sumAlleles[m]) ? (sumsqAlleles[m]-sumAlleles[m])/2. : ( numAlleles[m] - 1.5 * sumAlleles[m] + 0.5 * sumsqAlleles[m] )) + .5));
    }
    else {
      float g;
      int cnts[3] = {0,0,0};
      for(int j=0; j < nInds; ++j) {
	g = genos[m*nInds + j];
	if ( !isnan(g) ) { // do not count missing at any place
	  if ( g < 0.5 ) ++cnts[0];
	  else if ( g >= 1.5 ) ++cnts[2];
	  else ++cnts[1];
	}
      }
      return cnts[2];
    }
  }

  // works only 
  // numAlleles[m] = 2*(r + h + a)
  // sumAlleles[m] = h + 2a
  // sumsqAlleles[m] = h + 4a
  // a = ((h+4a)-(h+2a))/2
  // r = (r+h+a) - (h+2a) + ((h+4a)-(h+2a))/2
  //   = (r+h+a) - 1.5(h+2a) + 0.5(h+4a)
  // h = (r+h+a) - a - r
  void GENOCNT(int m, int* cnts) {
    if ( hardGenotype ) {
      cnts[0] = (int)(( numAlleles[m]/2 - 1.5 * sumAlleles[m] + 0.5 * sumsqAlleles[m] ) + 0.5);
      cnts[2] = (int)((sumsqAlleles[m]-sumAlleles[m])/2.+0.5);
      cnts[1] = numAlleles[m]/2 - cnts[0] - cnts[2];
    }
    else {
      float g;
      cnts[0] = cnts[1] = cnts[2] = 0;
      for(int j=0; j < nInds; ++j) {
	g = genos[m*nInds + j];
	if ( !isnan(g) ) { // do not count missing at any place
	  if ( g < 0.5 ) ++cnts[0];
	  else if ( g >= 1.5 ) ++cnts[2];
	  else ++cnts[1];
	}
      }
    }
  }

  void CASECTRLCNT(int m, int* cnts, std::vector<bool>& isCases) {
    float g;
    cnts[0] = cnts[1] = cnts[2] = cnts[3] = cnts[4] = cnts[5] = 0;
    for(int j=0; j < nInds; ++j) {
      g = genos[m*nInds + j];
      if ( !isnan(g) ) { // do not count missing at any place
	if ( g < 0.5 ) ++cnts[0+isCases[j]*3];
	else if ( g >= 1.5 ) ++cnts[2+isCases[j]*3];
	else ++cnts[1+isCases[j]*3];
      }
    }
  }

  float MAF(int m) {
    float f = alleleFreq(m);
    return (f > .5 ? 1.-f : f);
  }

  float RSQ(int m) {
    int n = numAlleles[m]/2;
    if ( n > 1 ) {
      float s = sumAlleles[m];
      float varObs = (sumsqAlleles[m]-s*s/n)/(n-1);     // var(x)
      float varExp = 2. * s * ( n + n - s ) / ( n * n ); // 2pq
      return (varObs > varExp) ? 1.0 : varObs/varExp;
    }
    else {
      return 0;
    }
  }

  void print(FILE* fp) {
    float v;
    //fprintf(stderr,"%d %d %s\n",nInds,nMarkers,inds[nInds-1].c_str());
    if ( inds.size() > 0 ) {
      fprintf(fp, "#MARKER");
      for(int i=0; i < nInds; ++i) {
	if ( inds.empty() ) {
	  fprintf(fp, "\tID%d", i+1);
	}
	else {
	  fprintf(fp, "\t%s", inds[i].c_str());
	}
      }
      fprintf(fp,"\n");
    }
    for(int i=0; i < nMarkers; ++i) {
      fprintf(fp, "%s", markers[i].c_str());
      for(int j=0; j < nInds; ++j) {
	v = genos[i*nInds + j];
	if ( isnan(v) ) {
	  fprintf(fp,"\tNA");
	}
	else {
	  fprintf(fp, "\t%.4f", genos[i*nInds + j]);
	}
      }
      fprintf(fp,"\n");
    }
  }

  int parseMarkers(char* line, int startIdx = 9) {
    //notice("fVcf::parseMarkers() called");
    //, std::vector<std::string>& m, std::vector<float>& v, int startIdx = 9, const char* key = "GT", const char** infoSubstrs = NULL, int nSubstr = 0) {
    const char* p;
    const char* n;
    const char* pch = line;
    const char* nch = NULL;
    std::string chrom;
    std::string pos;
    std::string ref;
    std::string alt;
    std::string markerId;
    std::string markerKey;
    std::string s;
    int keyIdx = 0;
    int i, j, k;
    int lKey = (int)key.size();
    int AN = 0;
    float AC = 0, sqAC = 0;
    float f;

    for(i=0, j=0; pch != NULL; ++i) {
      nch = strchr(pch, '\t');
      if ( i < startIdx ) {
	if ( i < 7 ) {
	  if ( nch == NULL ) s.assign( pch );
	  else s.assign( pch, nch - pch );

	  switch(i) {
	  case 0:                // copy chromosomes
	    chrom = s; break;
	  case 1:                // copy pos
	    pos = s; break;
	  case 2:                // marker_id
	    markerId = s; break;
	  case 3:                // reference allele
	    ref = s; break;
	  case 4:                // non-reference alleles
	    alt = s;  // if it is not in the predefined markerSet, skip the marker
	    if ( ! markerSet.empty() ) {
	      markerKey = chrom+":"+pos+"_"+ref+"/"+alt;
	      if ( markerSet.find(markerKey) == markerSet.end() ) {
		return -1;
	      }
	    }
	    break;
	  case 6:
	    if ( ( passOnly ) && ( s != "PASS" ) ) return -1;
	  }
	}
	else if ( i == 7 ) { // parse INFO field
	  if ( ( nch == NULL ) || ( infoAF ) ) {
	    // If the VCF is site only, parse AC & AN information
	    // from the INFO field, if exists
	    if ( nch != NULL ) {
	      s.assign( pch, nch - pch );
	      pch = s.c_str();
	    }

	    const char* ic;
	    AC = AN = -1;
	    double AF = -1;
	    if ( (ic = strstr(pch, "AC=")) != NULL ) {
	      AC = (float)atoi(ic+3);
	    }
	    if ( (ic = strstr(pch, "AN=")) != NULL ) {
	      AN = atoi(ic+3);
	    }
	    if ( (ic = strstr(pch, "AF=")) != NULL ) {
	      AF = atof(ic+3);
	    }

	    if ( ( AN < 0 ) || ( AC < 0 ) ) {
	      if ( AF >= 0 ) {
		if ( AFs.size() == markers.size() ) {
		  AFs.push_back(AF);
		}
		else {
		  error("AFs.size() != markers.size() after %s",markers.back().c_str());
		}
	      }
	      else {
		if ( nwarns < 1 ) {
		  warning("AN/AC or AF is not observed in INFO field. Assuming AC=1, AN=100");
		  ++nwarns;
		}

		AC = 1;
		AN = 100;
		AFs.push_back(0.01);
		sqAC = floor(2*AN*0.01*(1+0.01)+0.5);
	      }
	    }
	    else {
	      double af = (double)AC/(double)AN;
	      if ( AFs.size() == markers.size() ) {
		AFs.push_back(AF);
	      }
	      else {
		error("AFs.size() != markers.size() after %s",markers.back().c_str());
	      }
	      sqAC = floor(2*AN*af*(1.+af)+.5); // assumes HWE sqAC = p^2 * 4 + 2p(1-p) = 2p(1+p)
	    }

	    if ( !parser.parse( pch ) ) return -1; 
	  }
	  else { if ( !parser.parse( pch, nch - pch ) ) return -1; }
	}
	else if ( i == 8 ) { // parse FORMAT field
	  if ( pch[0] == key[0] && pch[1] == key[1] ) { // comparing two characters are not exactly right
	    keyIdx = 0;
	  }
	  else if ( nch == NULL ) {
	    error("VCF has FORMAT field but does not have any genotype");
	  }
	  else {
	    k = 0;
	    keyIdx = 0;
	    p = pch;
	    while( p < nch ) {
	      if ( *p == ':' ) { 
		if ( k >= lKey ) {
		  break;
		}
		else {
		  ++keyIdx;
		  k = 0;
		}
	      }
	      else {
		if ( ( k == 2 ) || ( key[k] == *p ) ) {
		  //if ( key[k] == *p ) {
		  ++k;
		}
		else {
		  k = 0;
		}
	      }
	      ++p;
	    }
	    if ( ( p == nch ) && ( k != lKey ) ) {
	      if ( nwarns < 100 ) {
		warning("Cannot find %s in the FORMAT field at marker %s:%s_%s/%s .. Skipping",key.c_str(),chrom.c_str(),pos.c_str(),ref.c_str(),alt.c_str());
		++nwarns;
	      }
	      return -1;
	    }
	  }
	}
      }
      else {
	if ( ( icols.empty() && ignoreEmpty ) || ( ( j < (int)icols.size() ) && ( icols[j] == i - startIdx ) ) ) {
	  p = pch;
	  n = NULL;
	  
	  // reach to the key index
	  if ( keyIdx > 0 ) {
	    for(int i=0; (i < keyIdx) && (p != NULL); ++i) {
	      n = strchr(p, ':');
	      p = (n == NULL) ? NULL : n+1;
	    }
	  }
	  
	  // if the field contains '/' or '|', add two values
	  // if the field is '.', return NA
	  // otherwise, convert the field into float and return
	  if ( ( p == NULL ) || ( p[0] == '.' ) ) { // missing
	    if ( glFlag ) {
	      PLs.push_back(0);
	      PLs.push_back(0);
	      PLs.push_back(0);
	    }
	    else {
	      genos.push_back(NAN_FLT);
	      phases.push_back(0); // unphased
	    }
	  }
	  else if ( ( p[1] == '/' ) || ( p[1] == '|' ) ) {
	    f = (float)((p[0] - '0') + ( p[2] - '0' )); // ignore phase info
	    genos.push_back(f);
	    phases.push_back((p[1] == '|') ? ( (p[0] == p[2]) ? 1 : (p[0] < p[2] ? 2 : 3) ) : 0);
	    AN += 2;
	    AC += f;
	    sqAC += f*f;
	  }
	  else { // Genotypes are not in 0|0 format
	    if ( glFlag ) { // search for three commas
	      const char* c1 = strchr(p,',');
	      if ( c1 == NULL ) error("Cannot parse %s field (currently suppports only biallelic autosomal variants",key.c_str());
	      const char* c2 = strchr(c1+1,',');
	      if ( c2 == NULL ) error("Cannot parse %s field (currently suppports only biallelic autosomal variants",key.c_str());
	      if ( phredFlag ) {
		int pl = atoi(p);
		if ( pl > 255 ) pl = 255;
		PLs.push_back(pl);
		pl = atoi(c1+1);
		if ( pl > 255 ) pl = 255;
		PLs.push_back(pl);
		pl = atoi(c2+1);
		if ( pl > 255 ) pl = 255;
		PLs.push_back(pl);
	      }
	      else {
		float gl = strtof(p,NULL);
		if ( gl < -25.5 ) PLs.push_back(255);
		else PLs.push_back((int)(-10*gl+.5));

		gl = strtof(c1+1,NULL);
		if ( gl < -25.5 ) PLs.push_back(255);
		else PLs.push_back((int)(-10*gl+.5));

		gl = strtof(c2+1,NULL);
		if ( gl < -25.5 ) PLs.push_back(255);
		else PLs.push_back((int)(-10*gl+.5));
	      }
	    }
	    else {  // take genos as dosages
	      f = strtof(p,NULL);
	      genos.push_back(f);
	      phases.push_back(0);
	      AN += 2;
	      AC += f;
	      sqAC += f*f;
	    }
	  }
	  ++j;
	}
	else {
	  //std::cout << "Skipping " << i << "\t" << j << std::endl;
	}
      }
      pch = ( nch == NULL ) ? NULL : nch + 1;
    }

    if ( (nInds == 0) && ( icols.empty() && ignoreEmpty ) && ((int)genos.size() == j) ) nInds = j;
    
    if ( nInds != j ) {
      fprintf(stderr,"i=%d, j=%d, nInds=%d, icols.size()=%d, genos.back()=%d\n",i,j,nInds,(int)icols.size(),(int)genos.back());
      abort();
    }

    // if GL or PL flag is set, automatically put dosage as quantitative genotypes
    if ( glFlag ) {
      //notice("fVcf::parseMarkers() - glFlag is on");

      // set allele frequency
      int m = markers.size();  //notice("fVcf::parseMarkers() - m = %d",m);
      float af = (float)(emAF(m,1e-6).first); //notice("fVcf::parseMarkers() - af = %f",af);
      if ( (int)AFs.size() > m ) { AFs[m] = af; }
      else if ( (int)AFs.size() == m ) { AFs.push_back(af); }
      else { error("fVcf::parseMarkers() -- AFs.size() < m"); }

      double p0,p1,p2;
      int kos;
      for(int k=0; k < j; ++k) {
	// calculate genotype dosages under HWE
	// Pr(G|D) = Pr(D|G)Pr(G)
	kos = 3*(m*nInds+k);
	p0 = phredConv.phred2Err[PLs[kos+0]] * (1.-AFs[m]) * (1.-AFs[m]);
	p1 = phredConv.phred2Err[PLs[kos+1]] * 2. * AFs[m] * (1.-AFs[m]);
	p2 = phredConv.phred2Err[PLs[kos+2]] * AFs[m] * AFs[m];
	if ( p0+p1+p2 > 0 ) {
	  genos.push_back((p1+2*p2)/(p0+p1+p2));
	  phases.push_back(0);
	}
	else {
	  genos.push_back(AFs[m]*2.);
	  phases.push_back(0);
	}
	if ( gpFlag ) { GPs.push_back(p0); GPs.push_back(p1); GPs.push_back(p2); }

	AN += 2;
	AC += genos.back();
	sqAC += genos.back()*genos.back();
      }
    }

    //notice("AC=%f, AN=%d, sqAC=%f",AC,AN,sqAC);

    if ( markerId == "." ) {
      markers.push_back(chrom+":"+pos+"_"+ref+"/"+alt);
    }
    else {
      markers.push_back(chrom+":"+pos+"_"+ref+"/"+alt+"_"+markerId);
    }

    chroms.push_back(chrom);
    pos1s.push_back(atoi(pos.c_str()));
    refs.push_back(ref);
    alts.push_back(alt);
    numAlleles.push_back(AN);
    sumAlleles.push_back(AC);
    sumsqAlleles.push_back(sqAC);

    //fprintf(stderr,"%s\n",markers.back().c_str());

    return j;
  }

  // Parse GT : GD/DP : GL/PL
  int fullParseMarkers(char* line, int startIdx = 9) {
    char* p;
    char* n;
    char* pch = line;
    char* nch = NULL;
    std::string chrom;
    std::string pos;
    std::string ref;
    std::string alt;
    std::string markerId;
    std::string markerKey;
    std::string s;

    int GTidx = -1, PLidx = -1, GLidx = -1, DPidx = -1, GDidx = -1;

    int i, j;
    int AN = 0;
    float AC = 0, sqAC = 0;
    float gt, gl;
    int dp, pl;
    int pls[3];
    uint8_t phase = 0;
    char *c1, *c2;

    for(i=0, j=0; pch != NULL; ++i) {
      nch = strchr(pch, '\t');
      if ( i < startIdx ) {
	if ( i < 7 ) {
	  if ( nch == NULL ) s.assign( pch );
	  else s.assign( pch, nch - pch );
	  switch(i) {
	  case 0:                // copy chromosomes
	    chrom = s; break;
	  case 1:                // copy pos
	    pos = s; break;
	  case 2:                // marker_id
	    markerId = s; break;
	  case 3:                // reference allele
	    ref = s; break;
	  case 4:                // non-reference alleles
	    alt = s;  // if it is not in the predefined markerSet, skip the marker
	    if ( ! markerSet.empty() ) {
	      markerKey = chrom+":"+pos+"_"+ref+"/"+alt;
	      if ( markerSet.find(markerKey) == markerSet.end() ) {
		return -1;
	      }
	    }
	    break;
	  case 6:
	    if ( ( passOnly ) && ( s != "PASS" ) ) return -1;
	  }
	}
	else if ( i == 7 ) { // parse INFO field
	  if ( nch == NULL ) { if ( !parser.parse( pch ) ) return -1; }
	  else { if ( !parser.parse( pch, nch - pch ) ) return -1; }
	}
	else if ( i == 8 ) { // parse FORMAT field
	  if ( nch == NULL ) {
	    error("VCF has FORMAT field but does not have any genotype");
	  }
	  else {
	    int ncolons = 0;
	    p = pch-1;
	    while(p < nch) {
	      if ( ( p < pch ) || ( *p == ':' ) ) {
		if ( ( p+3 >= nch ) || ( p[3] == ':' ) ) {
		  switch(p[1]) {
		  case 'G':
		    if ( p[2] == 'T' ) GTidx = ncolons;
		    else if ( p[2] == 'D' ) GDidx = ncolons;
		    else if ( p[2] == 'L' ) GLidx = ncolons;
		    break;
		  case 'P':
		    if ( p[2] == 'L' ) PLidx = ncolons;
		    break;
		  case 'D':
		    if ( p[2] == 'P' ) DPidx = ncolons;
		    break;
		  }
		}
		++ncolons;
	      }
	      ++p;
	    }
	  }
	}
      }
      else {
	if ( ( icols.empty() && ignoreEmpty ) || ( ( j < (int)icols.size() ) && ( icols[j] == i - startIdx ) ) ) {
	  p = pch;
	  n = NULL;

	  gt = NAN_FLT;
	  phase = 0;
	  dp = 0;
	  pls[0] = pls[1] = pls[2] = 0;

	  for(int i2=0; ((nch == NULL ) || (p < nch)) && (p != NULL); ++i2) {
	    if ( GTidx == i2 ) {
	      if ( ( p[1] == '/' ) || ( p[1] == '|' ) ) {
		gt = (float)((p[0] - '0') + ( p[2] - '0' )); // ignore phase info
		phase = ((p[1] == '|') ? ( (p[0] == p[2]) ? 1 : ((p[0] < p[2]) ? 2 : 3) ) : 0);
	      }
	    }
	    else if ( ( GDidx == i2 ) || ( DPidx == i2 ) ) {
	      dp = atoi(p);
	    }
	    else if ( GLidx == i2 ) {
	      c1 = strchr(p,',');
	      if ( c1 == NULL ) error("Cannot parse %s field (currently suppports only biallelic autosomal variants",key.c_str());
	      c2 = strchr(c1+1,',');
	      if ( c2 == NULL ) error("Cannot parse %s field (currently suppports only biallelic autosomal variants",key.c_str());
	      
	      gl = strtof(p,NULL);
	      if ( gl < -25.5 ) pls[0] = 255;
	      else pls[0] = (int)(-10*gl+.5);
	      
	      gl = strtof(c1+1,NULL);
	      if ( gl < -25.5 ) pls[1] = 255;
	      else pls[1] = (int)(-10*gl+.5);
	      
	      gl = strtof(c2+1,NULL);
	      if ( gl < -25.5 ) pls[2] = 255;
	      else pls[2] = (int)(-10*gl+.5);
	    }
	    else if ( PLidx == i2 ) {
	      c1 = strchr(p,',');
	      if ( c1 == NULL ) error("Cannot parse %s field (currently suppports only biallelic autosomal variants",key.c_str());
	      c2 = strchr(c1+1,',');
	      if ( c2 == NULL ) error("Cannot parse %s field (currently suppports only biallelic autosomal variants",key.c_str());
	      
	      pl = atoi(p);
	      if ( pl > 255 ) pl = 255;
	      pls[0] = pl;

	      pl = atoi(c1+1);
	      if ( pl > 255 ) pl = 255;
	      pls[1] = pl;

	      pl = atoi(c2+1);
	      if ( pl > 255 ) pl = 255;
	      pls[2] = pl;
	    }
	    n = strchr(p, ':');
	    p = (n == NULL) ? NULL : n+1;
	  }

	  genos.push_back(gt);
	  phases.push_back(phase);
	  PLs.push_back(pls[0]); 	  
	  PLs.push_back(pls[1]); 
	  PLs.push_back(pls[2]);
	  depths.push_back(dp);
	  if ( !isnan(gt) ) {
	    AN += 2;
	    AC += gt;
	    sqAC += (gt*gt);
	  }
	  ++j;
	}
	else {
	  //std::cout << "Skipping " << i << "\t" << j << std::endl;
	}
      }
      pch = ( nch == NULL ) ? NULL : nch + 1;
    }

    //error("depths.size() = %d, PLs.size() = %d / %f %d %d %d %d %d / %d",(int)depths.size(),(int)PLs.size(), genos.back(), depths.back(), (int)PLs.back(), (int)pls[0], (int)pls[1], (int)pls[2], j);

    if ( (nInds == 0) && ( icols.empty() && ignoreEmpty ) && ((int)genos.size() == j) ) nInds = j;

    if ( nInds != j ) {
      fprintf(stderr,"i=%d, j=%d, nInds=%d, icols.size()=%d, genos.back()=%d\n",i,j,nInds,(int)icols.size(),(int)genos.back());
      abort();
    }

    if ( markerId == "." ) {
      markers.push_back(chrom+":"+pos+"_"+ref+"/"+alt);
    }
    else {
      markers.push_back(chrom+":"+pos+"_"+ref+"/"+alt+"_"+markerId);
    }
    chroms.push_back(chrom);
    pos1s.push_back(atoi(pos.c_str()));
    refs.push_back(ref);
    alts.push_back(alt);
    numAlleles.push_back(AN);
    sumAlleles.push_back(AC);
    sumsqAlleles.push_back(sqAC);

    //fprintf(stderr,"%s\n",markers.back().c_str());

    return j;
  }

  // Parse GT : GD/DP : GL/PL
  int writeSubsetMarker(wFile& wf, const char* line, int startIdx = 9) {
    // this assumes that GT field was calculated already using parseMarker()
    int i,j,k;
    const char *p;
    const char* pch = line;
    const char* nch = NULL;
    for(i=0, j=0; pch != NULL; ++i) {
      nch = strchr(pch, '\t');
      if ( nch == NULL ) { nch = pch + strlen(pch); }
      if ( i < startIdx ) {
	if ( i > 0 ) wf.putc('\t');
	if ( ( i < 7 ) || (i == 8 ) ) {   // copy until FILTER field
	  for(p=pch; p < nch; ++p) { wf.putc(*p); }
	}
	else if ( i == 7 ) { // special handling for INFO field
	  // remove AC, AN, AF columns if exists
	  std::vector<std::string> tokens;
	  std::string info( pch, nch - pch );
	  pFile::tokenizeLine(info.c_str(),";",tokens);
	  wf.printf("AC=%d;AN=%d;NS=%d;AF=%.5f",(int)sumAlleles.back(),(int)numAlleles.back(),(int)numAlleles.back()/2,( numAlleles.back() > 0 ) ? (sumAlleles.back() / (float)numAlleles.back()) : 0);
	  for(k=0; k < (int)tokens.size(); ++k) {
	    if ( ( tokens[k].compare(0,3,"AC=") == 0 ) ||
		 ( tokens[k].compare(0,3,"AN=") == 0 ) ||
		 ( tokens[k].compare(0,3,"NS=") == 0 ) ||
		 ( tokens[k].compare(0,3,"AF=") == 0 ) ) {
	      // skip
	    }
	    else { 
	      wf.printf(";%s",tokens[k].c_str());
	    }
	  }
	}
      }
      else {
	if ( ( icols.empty() && ignoreEmpty ) || ( ( j < (int)icols.size() ) && ( icols[j] == i - startIdx ) ) ) {
	  wf.putc('\t');
	  for(p=pch; p < nch; ++p) { wf.putc(*p); }
	  ++j;
	}
      }
      pch = ( *nch == '\0' ) ? NULL : nch+1;
    }
    wf.putc('\n');
    //error("foo");
    return j;
  }

  int parseInds(char* line, std::set<std::string>& idset, int startIdx = 9) {
    char* pch = line;
    char* nch = NULL;
    int i, j;

    icols.clear();
    for(i=0, j=0; pch != NULL; ++i) {
      nch = strchr(pch, '\t');
      if ( i >= startIdx ) {
	std::string id = (nch == NULL) ? std::string(pch) : std::string(pch,nch-pch);
	if ( ( idset.empty() && ignoreEmpty ) || (idset.find(id) != idset.end()) ) {
	  icols.push_back(i - startIdx);
	  inds.push_back(id);
	  ++j;
	}
      }
      pch = ( nch == NULL ) ? NULL : nch + 1;
    }
    //notice("icols.size() = %d",(int)icols.size());
    return j;
  }

  int parseInds(char* line, std::vector<std::string>& subids, int startIdx = 9) {
    char* pch = line;
    char* nch = NULL;
    int i, j;

    icols.clear();
    for(i=0, j=0; pch != NULL; ++i) {
      nch = strchr(pch, '\t');
      if ( i >= startIdx ) {
	std::string id = (nch == NULL) ? std::string(pch) : std::string(pch,nch-pch);
	if ( subids[(int)icols.size()] == id ) {
	  icols.push_back(i - startIdx);
	  inds.push_back(id);
	  ++j;
	}
      }
      pch = ( nch == NULL ) ? NULL : nch + 1;
    }
    if ( j != (int)subids.size() ) {
      error("ERROR in fVcf::parseInds() - not all subids matches");
    }
    return j;
  }

  int parseInds(char* line, std::vector<int>& subcols, int startIdx = 9) {
    //notice("fVcf::parseInds(char*, std::vector<int>&) called");

    char* pch = line;
    char* nch = NULL;
    int i, j;

    icols.clear();
    for(i=0, j=0; pch != NULL; ++i) {
      nch = strchr(pch, '\t');
      if ( i >= startIdx ) {
	if ( ( subcols.empty() && ignoreEmpty ) || subcols[j] == i - startIdx ) {
	  icols.push_back(i - startIdx);
	  std::string id = (nch == NULL) ? std::string(pch) : std::string(pch,nch-pch);
	  inds.push_back(id);
	  ++j;
	}
	else { // Skip the individual
	  //std::cerr << "Skipping " << i << "\t" << j << std::endl;
	  //abort();
	}
      }
      pch = ( nch == NULL ) ? NULL : nch + 1;
    }
    return j;
  }

  // read markers until reach the end of file
  int readMarkerGroup(const char** pMarkers, int nMarkers, int startIndex = 0, bool sepchr = false) {
    std::vector<std::string> markerIDs;
    for(int i=0; i < nMarkers; ++i) {
      markerIDs.push_back(pMarkers[i]);
    }
    return readMarkerGroup(markerIDs,startIndex,sepchr);
  }

  int readMarkerGroup(std::vector<std::string>& markerIDs, 
		      int startIndex = 0, bool sepchr = false) 
  {
    //int nc;
    char* line = NULL;

    std::string curChrom;
    int beg = 1000000000, end = 0;

    std::vector<std::string> tokens;
    for(int i=startIndex; i < (int)markerIDs.size(); ++i) {
      pFile::tokenizeLine(markerIDs[i].c_str(),":_/",tokens);
      if ( tokens.size() != 4 ) {
	error("Cannot parse markerID %s",markerIDs[i].c_str());
      }
      std::string region(tokens[0]);
      region += ":";
      region += tokens[1];
      region += "-";
      region += tokens[1];

      if ( i == startIndex ) {
	curChrom = tokens[0];
	beg = end = atoi(tokens[1].c_str());
      }
      else {
	if ( curChrom != tokens[0] ) { 
	  curChrom = "multichrs";
	  beg = end = 0;
	  //error("Currently group across different chromosomes are not supported"); 
	}
	else {
	  int bp = atoi(tokens[1].c_str());
	  if ( beg > bp ) { beg = bp; }
	  if ( end < bp ) { end = bp; }
	}
      }

      tf.updateRegion(region.c_str(),sepchr);

      markerSet.clear();
      markerSet.insert(markerIDs[i]);
      while ( (line = (char*)tf.getLine()) != NULL ) {
	//nc = tf.getLength();
	int cols2 = parseMarkers(line); 
	if ( cols2 >= 0 ) {
	  if ( nInds != cols2 ) {
	    error("fVcf::readMarkers() : Column size does not match : %d vs %d at marker %s", nInds, cols2, markers[markers.size()-1].c_str() );
	    abort();
	  }
	  ++nMarkers;
	}
      }
    }

    if ( startIndex > 0 ) {
      char tmp[1024];
      sprintf(tmp,"%s:%d-%d_%s",curChrom.c_str(),beg,end,markerIDs[0].c_str());
      //markerIDs[0] = tmp;
      groupID = tmp;
      groupChrom = curChrom;
      groupBeg = beg;
      groupEnd = end;
    }
    return nMarkers;
  }

  // read markers until reach the end of file
  int readMarkers(int m = 0, bool del = true) {
    //int nc;
    char* line = NULL;
    
    if ( del ) clear();
    
    int mStart = nMarkers;

    //fprintf(stderr,"fVcf::readMarkers(%d) called",m);

    while ( (line = (char*)tf.getLine()) != NULL ) {
      //notice("bar");
      //nc = tf.getLength();
      //fprintf(stderr,"nMarkers=%d, line[0] = %c, icols.size() = %d, nc=%d\n",nMarkers,line[0],(int)icols.size(), nc);
      if ( line[0] == '#' ) {
	if ( strncmp(line,"#CHROM",6) == 0 ) {
	  nInds = parseInds(line, icols);
	}
      }
      else {
	int cols2 = parseMarkers(line); 
	//notice("cols2 = %d",cols2);
	if ( cols2 >= 0 ) {
	  if ( nInds != cols2 ) {
	    error("fVcf::readMarkers() : Column size does not match : %d vs %d at marker %s", nInds, cols2, markers[markers.size()-1].c_str() );
	    abort();
	  }
	  ++nMarkers;
	  if ( ( m > 0 ) && ( m <= nMarkers - mStart ) ) break;
	}
      }
    }
    //notice("Returning %d",nMarkers-mStart);
    return ( nMarkers-mStart );
  };

  // read markers until reach the end of file
  int fullReadMarkers(int m = 0, bool del = true) {
    //int nc;
    char* line = NULL;
    
    if ( del ) clear();
    
    int mStart = nMarkers;

    //notice("fVcf::readMarkers(%d) called",m);

    while ( (line = (char*)tf.getLine()) != NULL ) {
      //notice("bar");
      //nc = tf.getLength();
      //notice("nMarkers=%d, line[0] = %c, icols.size() = %d, nc=%d",nMarkers,line[0],(int)icols.size(), nc);
      if ( line[0] == '#' ) {
	if ( strncmp(line,"#CHROM",6) == 0 ) {
	  nInds = parseInds(line, icols);
	}
      }
      else {
	int cols2 = fullParseMarkers(line); 
	//notice("cols2 = %d",cols2);
	if ( cols2 >= 0 ) {
	  if ( nInds != cols2 ) {
	    error("fVcf::readMarkers() : Column size does not match : %d vs %d at marker %s", nInds, cols2, markers[markers.size()-1].c_str() );
	    abort();
	  }
	  ++nMarkers;
	  if ( ( m > 0 ) && ( m <= nMarkers - mStart ) ) break;
	}
      }
    }
    //notice("Returning %d",nMarkers-mStart);
    return ( nMarkers-mStart );
  };

  void close() {
    tf.close();
  }

  /*
  static FILE* openVCF(const char* filename, const char* region = NULL, const char* tabix = "/usr/cluster/bin/tabix" ) {
    // check if the file exists
    if ( exists(filename) == 0 ) {
      std::cerr << "Cannot open file " << filename << " for reading" << std::endl;
      return NULL;
    }
    
    std::string fn(filename);
    FILE* fp;
    if ( fn.substr(fn.size()-3) == ".gz" ) {
      std::string cmd = tabix;
      if ( region == NULL ) {
	cmd = "zcat ";
      }
      else if ( exists(tabix) == 0 ) {
	std::cerr << "Cannot find tabix binary " << tabix << ". Failed opening with region specified";
	return NULL;
      }
      else {
	cmd += " -h ";
      }
      cmd += fn;
      if ( region != NULL ) {
	cmd += " ";
	cmd += region;
      }
      //std::cout << cmd << std::endl;
      fp = popen(cmd.c_str(),"r");
      return fp;
    }
    else {
      fp = fopen(filename,"r");
      return fp;
    }
  }
  */

  /**
     Computes HWE allele frequencies using EM algorithm
     Input:
     1)Genotype Likelihoods (qscores)
     
     Output: 
     1) estimated HWE allele frequencies
     2) sample size
  */
  std::pair<double,double> emAF(int m, double eps) {
    // initial step : start from randomly assigned allele frequencies
    std::vector<int> indices;
    for(int i=0; i < nInds; ++i) {
      indices.push_back(i);
    }
    return emAF(m, indices,eps);
  }

  std::pair<double,double> emAF(int m, std::vector<int>& indices, double eps) {
    // initialization : pick an arbitrary AF from 
    //notice("fVcf::emAF(%d, %d, %lf), p2e.size() = %d",m,(int)indices.size(),eps,(int)p2e.size());

    std::vector<int>::iterator it;
    int i, c, n, r;
    double p = .5 + rand()/(RAND_MAX+1.)*.3; // random AF
    double q;
    double f0,f1,f2, fsum, sum;
    double llk = 0;
    std::vector<double> post;
    n = (int)indices.size();
    post.resize(n*3);

    for(r = 0; r < 100; ++r) {
      sum = 0;
      c = 0;
      q = 1.-p;
      for(it = indices.begin(); it != indices.end(); ++it) {
	i = 3*m*nInds + (*it)*3;
	f0 = q * q * phredConv.phred2Err[PLs[i]];
	f1 = 2. * p * q * phredConv.phred2Err[PLs[i+1]];
	f2 = p * p * phredConv.phred2Err[PLs[i+2]];
	fsum = f0+f1+f2;
	post[c++] = f0/fsum;
	sum += (post[c++] = f1/fsum);
	sum += (2 * (post[c++] = f2/fsum));
      }
      p = sum / (2*n);
      if ( fabs(p + q - 1.) < eps ) break;
    }

    //notice("fVcf::emAF() - p = %lf",p);

    // Pr(Data|AF) = \sum_g Pr(Data|g)Pr(g|AF)
    q = 1.-p;
    for(it = indices.begin(); it != indices.end(); ++it) {
      i = 3*m*nInds + (*it)*3;
      f0 = q * q * phredConv.phred2Err[PLs[i]];
      f1 = 2 * p * q * phredConv.phred2Err[PLs[i+1]];
      f2 = p * p * phredConv.phred2Err[PLs[i+2]];
      fsum = f0+f1+f2;
      llk += log(fsum);
    }
    //notice("fVcf::emAF() finished - returning (%lf, %lf)",p,llk);
    return std::pair<double,double>(p,llk);
  }

  double LRT(int m, std::vector<int>& g1, std::vector<int>& g2, double& af1, double& af2) {
    std::vector<int> g;
    std::pair<double,double> p;
    for(int i=0; i < (int)g1.size(); ++i) { g.push_back(g1[i]); }
    for(int i=0; i < (int)g2.size(); ++i) { g.push_back(g2[i]); }
    double llk0 = emAF(m, g, 1e-6).second;
    p = emAF(m, g1, 1e-6);
    af1 = p.first;
    double llk1 = p.second;
    p = emAF(m, g2, 1e-6);
    af2 = p.first;
    llk1 += p.second;

    //printf("%d\t%lf\t%lf\t%lf\n",m,llk0,llk1,llk1-llk0);
    
    return( 2*(llk1 - llk0) );
  }
};

#endif // __TABIXED_VCF_H

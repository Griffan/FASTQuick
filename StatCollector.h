/*
 * StatCollector.h
 *
 *  Created on: 2014Äê7ÔÂ20ÈÕ
 *      Author: Administrator
 */
#include "Utility.h"
#include "cmath"
#include "fstream"
#include "./libStatGen-master/vcf/VcfFileReader.h"
#include "./libStatGen-master/vcf/VcfHeader.h"
#include "./libStatGen-master/vcf/VcfRecord.h"
#include "libbwa/bwtaln.h"
#ifndef STATCOLLECTOR_H_
#define STATCOLLECTOR_H_

using namespace std;
typedef string Seq;
typedef string Qual;
typedef unordered_map<unsigned int, unsigned int  >  SeqPairIndex;
typedef unordered_map<string, SeqPairIndex> string_map;
//typedef vector<VcfRecord> VcfRecVec;

class StatCollector
{
private:
	string_map PositionTable;
	unsigned int index;
	vector<string> SeqVec,QualVec;
	vector<vector<unsigned int> > CycleVec;

	vector<VcfRecord> VcfRecVec;
	string_map VcfTable;
	vector<vector<unsigned char> > QualDist;//40*40
	vector<vector<unsigned char> > CycleDist;//100*40

	unordered_map<string, unordered_map<unsigned int, bool> > VariantProxyTable;
	uint64_t total_base;
	uint64_t total_region_size;
public:
	StatCollector();
	StatCollector(const string & OutFile);
	int  addAlignment(const string & PosName, const string & seq, const string & qual, const int & n_cigar, const bwa_cigar_t * cigar, const unsigned int & pos, const gap_opt_t* opt);
	int restoreVcfSites(const string & VcfPath,const gap_opt_t* opt);
	int getDepthDist(const string & outputPath,const vector<uint64_t> & DepthDist);
	int getEmpRepDist(const string & outputPath,const vector<uint64_t> & EmpRepDist,const vector<uint64_t> &misEmpRepDist);
	int getEmpCycleDist(const string & outputPath,const vector<uint64_t> & EmpCycleDist, const vector<uint64_t> misEmpCycleDist);
	int outputPileup(const string & statPrefix);
	int processCore(const string & statPrefix);
	int getGenoLikelihood(const string & statPrefix);
	virtual ~StatCollector();
};

#endif /* STATCOLLECTOR_H_ */

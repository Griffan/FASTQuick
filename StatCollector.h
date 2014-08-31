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
#include "libbwa/bntseq.h"
#ifndef STATCOLLECTOR_H_
#define STATCOLLECTOR_H_

using namespace std;
typedef string Seq;
typedef string Qual;
typedef unordered_map< int, unsigned int  >  SeqPairIndex;
typedef unordered_map<string, SeqPairIndex> string_map;
typedef unordered_map<string, unordered_map< int, bool> > bool_table;
//typedef vector<VcfRecord> VcfRecVec;

class StatCollector
{
private:
	string_map PositionTable;
	unsigned int index;
	vector<string> SeqVec,QualVec;
	vector<vector<unsigned int> > CycleVec;
	vector<vector<unsigned char> > MaqVec;
	vector<vector<bool> >StrandVec;

	vector<VcfRecord*> VcfRecVec;
	string_map VcfTable;
	vector<vector<unsigned char> > QualDist;//40*40
	vector<vector<unsigned char> > CycleDist;//100*40

	vector<uint64_t> DepthDist;
	vector<uint64_t> EmpRepDist;
	vector<uint64_t> misEmpRepDist;
	vector<uint64_t> EmpCycleDist;
	vector<uint64_t> misEmpCycleDist;
	vector<uint64_t> GCDist;
	string_map GC;

	bool_table VariantProxyTable;
	uint64_t total_base;
	uint64_t total_region_size;
public:
	StatCollector();
	StatCollector(const string & OutFile);
	//int  addAlignment(const string & PosName, const string & seq, const string & qual, const int & n_cigar, const bwa_cigar_t * cigar, const unsigned int & pos, const gap_opt_t* opt);
	int addAlignment(const bntseq_t *bns, bwa_seq_t *p,const gap_opt_t* opt);
	int restoreVcfSites(const string & VcfPath,const gap_opt_t* opt);
	int getDepthDist(const string & outputPath,const vector<uint64_t> & DepthDist);
	int getGCDist(const string & outputPath,const vector<uint64_t> & DepthDist, const vector<uint64_t> & GCDist);
	int getEmpRepDist(const string & outputPath,const vector<uint64_t> & EmpRepDist,const vector<uint64_t> &misEmpRepDist);
	int getEmpCycleDist(const string & outputPath,const vector<uint64_t> & EmpCycleDist, const vector<uint64_t> misEmpCycleDist);
	int outputPileup(const string & statPrefix);
	int processCore(const string & statPrefix);
	int getGenoLikelihood(const string & statPrefix);

	inline bool ConstructFakeSeqQual(const string &seq, const string &qual,const int & n_cigar, const bwa_cigar_t* cigar,string & newSeq, string &newQual)
	{

		int last=0;
		for( int k=0;k< n_cigar;++k)
		{
				int cl= __cigar_len(cigar[k]);
				int cop= "MIDS"[__cigar_op(cigar[k])];
				switch (cop) {
				 case 'M':
					 newSeq+=seq.substr(last,cl);
					 newQual+=qual.substr(last,cl);
					 break;
				 case 'D':
					 newSeq+=string('N',cl);
					 newQual+=string('!',cl);
					 break;
				 case 'I':
					 break;
				 default:
					 break;
				}

	/*	string MD(md);
		last=0;


		int total_len=0;
		for(int i=0;i!=MD.size();++i)
			if(isdigit(MD[i])) continue;
			else if(MD[i]=='^')
			{
				i++;
				while(!isdigit(MD[i]))
				{
					i++;
					total_len++;
				}
				last=i;

			}else
			{
				impeccable=false;
				int len=atoi(MD.substr(last,i-last).c_str())+1;
				total_len+=len;
				count+=newQual[total_len-1]-33;
				last=i+1;
			}

	}

		if(impeccable||count >50) return false;
		else return true;
		*/
		}
		return true;
	}

	virtual ~StatCollector();
};

#endif /* STATCOLLECTOR_H_ */

/*
 * StatCollector.h
 *
 *  Created on: 2014Äê7ÔÂ20ÈÕ
 *      Author: Administrator
 */
#include "Utility.h"
#include "cmath"
#include "fstream"
#include "VcfFileReader.h"
#include "VcfHeader.h"
#include "VcfRecord.h"
#include "libbwa/bwtaln.h"
#include "libbwa/bntseq.h"
#ifndef STATCOLLECTOR_H_
#define STATCOLLECTOR_H_


typedef std::string Seq;
typedef std::string Qual;
typedef std::unordered_map< int, unsigned int  >  SeqPairIndex;
typedef std::unordered_map<std::string, SeqPairIndex> string_map;
typedef std::unordered_map<std::string, std::unordered_map< int, bool> > bool_table;
//typedef vector<VcfRecord> VcfRecVec;

class StatCollector
{
private:
	string_map PositionTable;
	unsigned int index;
	//nsigned int vcf_index;
	vector<uint32_t> DepthVec;
	/****for SNP site******/
	vector<std::string> SeqVec,QualVec;
	vector<vector<unsigned int> > CycleVec;
	vector<vector<unsigned char> > MaqVec;
	vector<vector<bool> >StrandVec;
	/************************/
	vector<VcfRecord*> VcfRecVec;
	string_map VcfTable;
	//string_map VcfObTable;//actually covered by reads
	vector<vector<unsigned char> > QualDist;//40*40
	vector<vector<unsigned char> > CycleDist;//100*40

	vector<int> DepthDist;
	vector<int> EmpRepDist;
	vector<int> misEmpRepDist;
	vector<int> EmpCycleDist;
	vector<int> misEmpCycleDist;
	vector<int> GCDist;
	vector<int> InsertSizeDist;
	vector<int> MaxInsertSizeDist;
	string_map GC;

	bool_table VariantProxyTable;


	uint64_t total_base;
	uint64_t total_region_size;

	std::unordered_map<std::string,bool> duplicateTable;



public:
	StatCollector();
	StatCollector(const std::string & OutFile);
	int addAlignment(const bntseq_t *bns, bwa_seq_t *p,const gap_opt_t* opt);
	int addPairAlignment(const bntseq_t *bns, bwa_seq_t *p, bwa_seq_t *q,const gap_opt_t* opt , std::ofstream & fout,int &);
	int IsDuplicated(const bntseq_t *bns, const bwa_seq_t *p, const bwa_seq_t *q,const gap_opt_t* opt, int type, std::ofstream & fout);
	int restoreVcfSites(const std::string & VcfPath,const gap_opt_t* opt);
	int getDepthDist(const std::string & outputPath);
	int getGCDist(const std::string & outputPath,const vector<int> & PosNum);
	int getEmpRepDist(const std::string & outputPath);
	int getEmpCycleDist(const std::string & outputPath);
	int getInsertSizeDist(const std::string & outputPath);
	int outputPileup(const std::string & statPrefix);

	int processCore(const std::string & statPrefix);
	int getGenoLikelihood(const std::string & statPrefix);
	inline int isPartialAlign(const bwa_seq_t * q)
	{
		for (int k = 0; k < q->n_cigar; ++k)
			{

				//int cl = __cigar_len(q->cigar[k]);
				int cop = "MIDS"[__cigar_op(q->cigar[k])];
				switch (cop)
				{
				case 'S':
					return 1;
				case 'H':
					return 1;
				default:
					break;
				}
			}
		return 0;
	}
	inline bool ConstructFakeSeqQual(const std::string &seq, const std::string &qual,const int & n_cigar, const bwa_cigar_t* cigar,std::string & newSeq, std::string &newQual)
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
					 newSeq+=std::string('N',cl);
					 newQual+=std::string('!',cl);
					 break;
				 case 'I':
					 break;
				 default:
					 break;
				}

	/*	std::string MD(md);
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

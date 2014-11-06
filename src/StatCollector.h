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
#include "../libbwa/bwtaln.h"
#include "../libbwa/bntseq.h"
#include "SamRecord.h"
#include "BamInterface.h"
#include "SamFile.h"
#ifndef STATCOLLECTOR_H_
#define STATCOLLECTOR_H_


//typedef std::string Seq;
//typedef std::string Qual;

typedef std::unordered_map<std::string, std::unordered_map< int, unsigned int  > > unsort_map;
typedef std::map<std::string, std::map< int, unsigned int  > > sort_map;
typedef std::unordered_map<std::string, std::unordered_map< int, bool> > bool_table;
//typedef vector<VcfRecord> VcfRecVec;
class FileStatCollector
{
public:
	long long  NumRead;
	long long  NumBase;
	std::string FileName1,FileName2;
	FileStatCollector():NumRead(0),NumBase(0),FileName1(0),FileName2(0)
	{}
	FileStatCollector(const char* filename):NumRead(0),NumBase(0),FileName1(filename),FileName2(0)
	{
	}
	FileStatCollector(const char* filename1, const char* filename2):NumRead(0),NumBase(0),FileName1(filename1),FileName2(filename2)
	{
	}
};
class StatCollector
{
private:
	unsort_map PositionTable;
	unsigned int index;
	//nsigned int vcf_index;
	std::vector<uint32_t> DepthVec;
	std::vector<uint32_t> Q20DepthVec;
	std::vector<uint32_t> Q30DepthVec;
	/****for SNP site******/
	std::vector<std::string> SeqVec,QualVec;
	std::vector<std::vector<unsigned int> > CycleVec;
	std::vector<std::vector<unsigned char> > MaqVec;
	std::vector<std::vector<bool> >StrandVec;
	/************************/
	std::vector<VcfRecord*> VcfRecVec;
	sort_map VcfTable;
	unsort_map dbSNPTable;
	//string_map VcfObTable;//actually covered by reads
	std::vector<std::vector<unsigned char> > QualDist;//40*40
	std::vector<std::vector<unsigned char> > CycleDist;//100*40


	std::vector<int> DepthDist;
	std::vector<int> EmpRepDist;
	std::vector<int> misEmpRepDist;
	std::vector<int> EmpCycleDist;
	std::vector<int> misEmpCycleDist;
	std::vector<int> GCDist;
	std::vector<int> InsertSizeDist;
	std::vector<int> MaxInsertSizeDist;
	unsort_map GC;

	bool_table VariantProxyTable;


	uint64_t total_base;
	uint64_t total_region_size;
	uint64_t ref_genome_size;

	std::unordered_map<std::string,bool> duplicateTable;

	std::unordered_map<std::string,ContigStatus> contigStatusTable;

	std::vector<FileStatCollector> FSCVec;

private:
	int addSingleAlignment(const bntseq_t *bns, bwa_seq_t *p,const gap_opt_t* opt);
	int addSingleAlignment( SamRecord& p,const gap_opt_t* opt);
public:
	StatCollector();
	StatCollector(const std::string & OutFile);

	int addAlignment(const bntseq_t *bns, bwa_seq_t *p, bwa_seq_t *q,const gap_opt_t* opt , std::ofstream & fout,int &);
	int IsDuplicated(const bntseq_t *bns, const bwa_seq_t *p, const bwa_seq_t *q,const gap_opt_t* opt, int type, std::ofstream & fout);
//overload functions for direct bam reading
	int addAlignment(SamFileHeader & SFH, SamRecord * p, SamRecord* q, const gap_opt_t* opt, std::ofstream & fout, int &);
	int IsDuplicated(  SamFileHeader& SFH, SamRecord& p, SamRecord& q, const gap_opt_t* opt, int type, std::ofstream & fout);

	int ReadAlignmentFromBam( const gap_opt_t* opt, /*SamFileHeader& SFH, SamFile& BamIO, */const char * BamFile, std::ofstream & fout,int & total_add);

	int restoreVcfSites(const std::string & VcfPath,const gap_opt_t* opt);
	int getDepthDist(const std::string & outputPath,const gap_opt_t* opt);
	int getGCDist(const std::string & outputPath,const std::vector<int> & PosNum);
	int getEmpRepDist(const std::string & outputPath);
	int getEmpCycleDist(const std::string & outputPath);
	int getInsertSizeDist(const std::string & outputPath);
	int getSexChromInfo(const std::string & outputPath);
	int outputPileup(const std::string & statPrefix);

	int processCore(const std::string & statPrefix, const gap_opt_t*opt);
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
	inline int isPartialAlign( SamRecord& q)
	{
		if(std::string(q.getCigar()).find('S')!= std::string::npos)
			return 1;
		else
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

	int addFSC(FileStatCollector a);
	int getGenomeSize(std::string RefPath);
	int SummaryOutput(const std::string & outputPath,const gap_opt_t* opt);
	virtual ~StatCollector();
};

#endif /* STATCOLLECTOR_H_ */

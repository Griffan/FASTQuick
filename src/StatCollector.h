/*The MIT License (MIT)
Copyright (c) 2017 Fan Zhang, Hyun Min Kang
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */
/* Contact: Fan Zhang <fanzhang@umich.edu> */


#include "Utility.h"
#include "cmath"
#include "fstream"
#include "../misc/vcf/VcfFileReader.h"
#include "../misc/vcf/VcfHeader.h"
#include "../misc/vcf/VcfRecord.h"
#include "../misc/bam/SamRecord.h"
#include "../misc/bam/BamInterface.h"
#include "../misc/bam/SamFile.h"
#include "../misc/bam/SamValidation.h"
#include "../libbwa/bwtaln.h"
#include "../libbwa/bntseq.h"

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
	FileStatCollector(const char* filename):NumRead(0),NumBase(0),FileName1(filename),FileName2("")
	{
	}
	FileStatCollector(const char* filename1, const char* filename2):NumRead(0),NumBase(0),FileName1(filename1),FileName2(filename2)
	{
	}
};
class StatCollector
{
private:

	/*Expected Info*/
	uint64_t total_base;
	uint64_t total_region_size;
	uint64_t ref_genome_size;

	/*Actual Statistics*/
	uint64_t NumPCRDup;
	uint64_t NumBaseMapped;
	uint64_t NumPositionCovered;//position with depth larger than 0
	//uint64_t NumPositionCovered1;
	uint64_t NumPositionCovered2;//larger than 1
	uint64_t NumPositionCovered5;//larger than 4
	uint64_t NumPositionCovered10;// larger than 9
	unsort_map PositionTable;//chrom pos absolute index of the site
	unsigned int index;
	//nsigned int vcf_index;
	std::vector<uint32_t> DepthVec;
	std::vector<uint32_t> Q20DepthVec;
	std::vector<uint32_t> Q30DepthVec;
	/****for SNP site******/
	std::vector<std::string> SeqVec,QualVec;
	std::vector<std::vector<int> > CycleVec;
	std::vector<std::vector<unsigned char> > MaqVec;
	std::vector<std::vector<bool> >StrandVec;
	/************************/
	std::vector<VcfRecord*> VcfRecVec;
	sort_map VcfTable;
	unsort_map dbSNPTable;
	//string_map VcfObTable;//actually covered by reads
	std::vector<size_t> QualDist;//40*40
	std::vector<size_t> CycleDist;//100*40


	std::vector<size_t> DepthDist;
	std::vector<size_t> EmpRepDist;
	std::vector<size_t> misEmpRepDist;
	std::vector<size_t> EmpCycleDist;
	std::vector<size_t> misEmpCycleDist;
	std::vector<size_t> GCDist;
	std::vector<size_t> InsertSizeDist;
	std::vector<size_t> MaxInsertSizeDist;
	unsort_map GC;

	bool_table VariantProxyTable;

	std::unordered_map<std::string,bool> duplicateTable;

	std::unordered_map<std::string,ContigStatus> contigStatusTable;

	std::vector<FileStatCollector> FSCVec;

private:
	bool addSingleAlignment(const bntseq_t *bns, bwa_seq_t *p, const gap_opt_t *opt);

	bool addSingleAlignment(SamRecord &p, const gap_opt_t *opt);

	void StatVecDistUpdate(const std::string &qual, int left_to_right_coord,
						   unsigned int tmp_index, const std::string &RefSeq, const std::string &seq,
						   int tmpCycle);

	void AddBaseInfoToNewCoord(const std::string &Chrom, int i,
							   const std::string &qual, int left_to_right_coord,
							   const std::string &RefSeq, const std::string &seq, int tmpCycle);

	void UpdateInfoVecAtMarker(int tmpCycleVcfTable,
							   int tmpCycle, int tmp_left_to_right_coord,
							   int left_to_right_coord, int realCoord, int cl,
							   const char *sign, bool strand, const std::string &Chrom,
							   unsigned int tmp_index, const std::string &seq, const std::string &qual,
							   u_char mapQ);

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
	int releaseVcfSites();
	int getDepthDist(const std::string & outputPath,const gap_opt_t* opt);
	int getGCDist(const std::string & outputPath,const std::vector<int> & PosNum);
	int getEmpRepDist(const std::string & outputPath);
	int getEmpCycleDist(const std::string & outputPath);
	int getInsertSizeDist(const std::string &outputPath);
	int getSexChromInfo(const std::string & outputPath);
	int outputPileup(const std::string & statPrefix, const gap_opt_t* opt);

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
		}
		return true;
	}
	int updateRefseqByMD(std::string&RefSeq, std::string&MD)
	{
		int last(0), total_len(0)/*pos on seq*/;
		for (uint32_t i = 0; i != MD.size(); ++i)//remember we are making reference sequence
			if (isdigit(MD[i]))
				continue;
			else if (MD[i] == '^')
			{
				int len = atoi(MD.substr(last, i - last).c_str());//len from last to current deletion, '^' not included
				total_len += len;
				int start_on_read = total_len;//1 based
				i++;
				std::string tmp;
				while (!isdigit(MD[i])) // we don't need to take care of Deletion
				{
					i++;
					total_len++;
					tmp += MD[i];
				}
				std::string left = RefSeq.substr(0, start_on_read);
				std::string right = RefSeq.substr(start_on_read, RefSeq.length() - start_on_read + 1);
				RefSeq = left + tmp + right;
				total_len--;
				last = i;
			}
			else
			{
				int len = atoi(MD.substr(last, i - last).c_str()) + 1;//len from last to current mismatch, include mismatch
				total_len += len;
				//if (strcmp(p.getReadName(),"WTCHG_8105:7:66:5315:89850#0")==0){ fprintf(stderr, "we see strange i:%dth out of %d char of MD tag %s, %d out of %d on read:%s", i,MD.length(),MD.c_str(),total_len-1,RefSeq.length(), p.getReadName()); exit(EXIT_FAILURE); }
				RefSeq[total_len - 1] = MD[i];
				last = i + 1;
			}
		return 0;
	}

	int addFSC(FileStatCollector a);
	int getGenomeSize(std::string RefPath);
	int SummaryOutput(const std::string & outputPath,const gap_opt_t* opt);

	double Q20AvgDepth();
	double Q30AvgDepth();
	size_t MIS500();
	size_t MIS300();
	double Q20BaseFraction();
	double Q30BaseFraction();
	virtual ~StatCollector();
};

#endif /* STATCOLLECTOR_H_ */

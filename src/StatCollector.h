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
	std::vector<uint32_t> DepthVec;//depth at each site
	std::vector<uint32_t> Q20DepthVec;
	std::vector<uint32_t> Q30DepthVec;
	/****for SNP site******/
	std::vector<std::string> SeqVec,QualVec;
	std::vector<std::vector<int> > CycleVec;
	std::vector<std::vector<unsigned char> > MaqVec;
	std::vector<std::vector<bool> >StrandVec;
	/************************/
	std::vector<VcfRecord*> VcfRecVec;
	sort_map VcfTable;// (chr,pos) -> index

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
//	std::vector<size_t> MaxInsertSizeDist;
	unsort_map GC;

	bool_table VariantProxyTable;

	std::unordered_map<std::string,bool> duplicateTable;

	std::unordered_map<std::string,ContigStatus> contigStatusTable;

	std::vector<FileStatCollector> FSCVec;

private:
	bool AddSingleAlignment(const bntseq_t *bns, bwa_seq_t *p, const gap_opt_t *opt);

	bool AddSingleAlignment(SamRecord &p, const gap_opt_t *opt);

	void StatVecDistUpdate(const std::string &qual, unsigned int tmpIndex, const std::string &refSeq,
                               const std::string &seq, int tmpCycle, int relativeCoordOnRead,
                               int relativeCoordOnRef);

	void AddBaseInfoToNewCoord(const std::string &chrom, int i, const std::string &qual,
                                   const std::string &refSeq, const std::string &seq, int tmpCycle,
                                   int relativeCoordOnRead, int relativeCoordOnRef);

	void UpdateInfoVecAtMarker(int tmpCycle, int absoluteSite, int cl, const char *sign, bool strand,
                               const std::string &chrom, const std::string &seq, const std::string &qual,
                               u_char mapQ, int relativeCoordOnRead);

public:
	StatCollector();
	StatCollector(const std::string & OutFile);

	int AddAlignment(const bntseq_t *bns, bwa_seq_t *p, bwa_seq_t *q, const gap_opt_t *opt, std::ofstream &fout,
                     long &total_add_failed);
	int IsDuplicated(const bntseq_t *bns, const bwa_seq_t *p, const bwa_seq_t *q,const gap_opt_t* opt, int type, std::ofstream & fout);
//overload functions for direct bam reading
	int AddAlignment(SamFileHeader &SFH, SamRecord *p, SamRecord *q, const gap_opt_t *opt, std::ofstream &fout,
					 long &total_add);
	int IsDuplicated(  SamFileHeader& SFH, SamRecord& p, SamRecord& q, const gap_opt_t* opt, int type, std::ofstream & fout);

	int ReadAlignmentFromBam( const gap_opt_t* opt, /*SamFileHeader& SFH, SamFile& BamIO, */const char * BamFile, std::ofstream & fout,int & total_add);

	int RestoreVcfSites(const std::string &RefPath, const gap_opt_t *opt);
	int ReleaseVcfSites();
	int GetDepthDist(const std::string &outputPath, const gap_opt_t *opt);
	int GetGCDist(const std::string &outputPath, const std::vector<int> &PosNum);
	int GetEmpRepDist(const std::string &outputPath);
	int GetEmpCycleDist(const std::string &outputPath);
	int GetInsertSizeDist(const std::string &outputPath);
	int GetSexChromInfo(const std::string &outputPath);
	int GetPileup(const std::string &statPrefix, const gap_opt_t *opt);

	int ProcessCore(const std::string &statPrefix, const gap_opt_t *opt);
	int GetGenoLikelihood(const std::string &statPrefix);
	inline int IsPartialAlign(const bwa_seq_t *q)
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
	inline int IsPartialAlign(SamRecord &q)
	{
		if(std::string(q.getCigar()).find('S')!= std::string::npos)
			return 1;
		else
			return 0;
	}

	std::string RecoverRefseqByMDandCigar(const std::string &readSeq, std::string MD, const std::string &cigarString);
    std::string RecoverRefseqByMDandCigar(const std::string &readSeq, std::string MD, const bwa_cigar_t * cigar, int n_cigar);

	int AddFSC(FileStatCollector a);
	int GetGenomeSize(std::string RefPath);
	int SummaryOutput(const std::string & outputPath,const gap_opt_t* opt);

	double Q20AvgDepth();
	double Q30AvgDepth();
	size_t MIS500();
	size_t MIS300();
	double Q20BaseFraction();
	double Q30BaseFraction();
	virtual ~StatCollector();

	int AddMatchBaseInfo(const gap_opt_t *opt, const std::string &seq, const std::string &qual,
                         const std::string &refSeq, const std::string &chr, int readRealStart,
                         int refRealStart, int refRealEnd, const char *sign, bool strand, u_char mapQ,
                         int matchLen, int tmpCycle, int relativeCoordOnRead, int relativeCoordOnRef);

	int UpdateInfoVecAtRegularSite(const gap_opt_t *opt, const std::string &seq, const std::string &qual,
                                   const std::string &refSeq, const std::string &chr, int readRealStart,
                                   int refRealStart, int refRealEnd, const char *sign, bool strand,
                                   int matchLen, int tmpCycle, int relativeCoordOnRead,
                                   int relativeCoordOnRef);
private:
//    std::vector<std::string> DebugSeqVec,DebugQualVec;
//    std::vector<std::vector<int> > DebugCycleVec;
//    std::vector<std::vector<unsigned char> > DebugMaqVec;
//    std::vector<std::vector<bool> >DebugStrandVec;
};

#endif /* STATCOLLECTOR_H_ */

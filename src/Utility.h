/* The MIT License

   Copyright (c) 2009 Genome Research Ltd (GRL), 2010 Broad Institute

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Fan Zhang <fanzhang@umich.edu> */
#ifndef UTILITY_
#define UTILITY_

#include  <string>
#include  <unordered_map>
#include <vector>
#include <fstream>
//#define SEQ_INIT_Len 128
//class FPC_seq_t{
//public:
//	char *name;
//	ubyte_t *seq, *rseq, *qual;
//	//char* original_seq;
//	uint32_t len : 20, strand : 1, type : 2, filtered : 1, extra_flag : 8;//change dummy bit into filtered bit
//	uint32_t n_mm : 8, n_gapo : 8, n_gape : 8, mapQ : 8;
//	int score;
//	int clip_len;
//	// alignments in SA coordinates
//	int n_aln;
//	bwt_aln1_t *aln;
//	// multiple hits
//	int n_multi;
//	bwt_multi1_t *multi;
//	// alignment information
//	bwtint_t sa, pos;
//	uint64_t c1 : 28, c2 : 28, seQ : 8; // number of top1 and top2 hits; single-end mapQ
//	int n_cigar;
//	bwa_cigar_t *cigar;
//	// for multi-threading only
//	int tid;
//	// barcode
//	char bc[16]; // null terminated; up to 15 bases
//	// NM and MD tags
//	uint32_t full_len : 20, nm : 12;
//	char *md;
//	//int  count;
//
//	FPC_seq_t()
//	{
//		name = new char[SEQ_INIT_Len * 2];
//		seq = new ubyte_t[SEQ_INIT_Len];
//		qual = new ubyte_t[SEQ_INIT_Len];
//		//rseq = new ubyte_t[SEQ_INIT_Len];
//		len = 0;
//		strand = 0;
//		type = 0;
//		filtered = 0;
//		extra_flag = 0;
//		n_mm = 0;
//		n_gape = 0;
//		n_gapo = 0;
//		mapQ = 0;
//		score = 0;
//		clip_len = 0;
//		n_aln = 0;
//		aln = 0;
//		n_multi = 0;
//		multi = 0;
//		sa = 0;
//		pos = 0;
//		c1 = 0;
//		c2 = 0;
//		seQ = 0;
//		n_cigar = 0;
//		cigar = 0;
//		tid = -1;
//		full_len = 0;
//		nm = 0;
//		md = 0;
//	}
//	~FPC_seq_t()
//	{
//		delete[] name;
//		delete[] seq;
//		delete[] qual;
//		//delete[] rseq;
//	}
//};
//static std::unordered_map<std::string, bool> VariantLongTable;
class _GCstruct
{
public:
	//char *Chrom;
	//char ChrNameLength;
	//unsigned int* pos;
	unsigned char* GC;
	unsigned int len;
	_GCstruct( unsigned int length)
	{
		len=length;
		//ChrNameLength=ChrName.size();
		//Chrom=new char[ChrNameLength+1];
	//	strcpy(Chrom,ChrName.c_str());
		//Chrom[ChrNameLength]='\0';
		//pos=new uint32_t [len];
		GC= new unsigned char [len];
	}
	_GCstruct( )
	{
		len=0;
		//ChrNameLength=ChrName.size();
		//Chrom=new char[ChrNameLength+1];
	//	strcpy(Chrom,ChrName.c_str());
		//Chrom[ChrNameLength]='\0';
		//pos=new uint32_t [len];
		GC=0;
	}
	~_GCstruct()
	{
		//delete Chrom;
		//delete pos;
		if(GC!=0)
		{
		delete [] GC;
		}
	}
	void write(std::ofstream& fout)
	{
		fout.write((char*)this,sizeof(_GCstruct));
		//fout.write((char*) &this->ChrNameLength,1);
		//fout.write((char*)&this->Chrom,this->ChrNameLength);
		//fout.write((char*) &this->len, sizeof(unsigned int));
		//fout.write((char*) &this->pos, this->len*sizeof(unsigned int));
		fout.write((char*) this->GC, this->len*sizeof(unsigned char));
	}
	void read(std::ifstream & fin)
	{
		fin.read((char*)this,sizeof(_GCstruct));
		//fin.read((char*) &this->ChrNameLength,1);
		//this->Chrom= new char[this->ChrNameLength+1];
		//fin.read((char*)&this->Chrom,this->ChrNameLength);
		//this->Chrom[this->ChrNameLength+1]='\0';
		//fin.read((char*) &this->len, sizeof(unsigned int));
		//this->pos = new unsigned int [this->len];
		this->GC = new unsigned char [this->len];
		//fin.read((char*) &this->pos, this->len*sizeof(unsigned int));
		fin.read((char*) this->GC, this->len*sizeof(unsigned char));
	}
};



class ContigStatus
{
public:
	 ContigStatus()
	{
		name="default";
		length=0;
		numOverlappedReads=0;
		numPairOverlappedReads=0;
		numFullyIncludedReads=0;
		numFullyIncludedPairedReads=0;
	}
	 ContigStatus(std::string chr, int len)
	{
		name=chr;
		length=len;
		numOverlappedReads=0;
		numPairOverlappedReads=0;
		numFullyIncludedReads=0;
		numFullyIncludedPairedReads=0;
	}
	 ~ContigStatus()
	 {

	 }
//	 SamRecord bwtSeq2SamRecord(const bntseq_t *bns, bwa_seq_t *p)
//	 {
//		 SamRecord tmp;
//
//	 }
//	 int addPair(SamRecord& p, SamRecord&q)
//	 {
//		 return 0;
//	 }
//	 int addSingle(SamRecord& p)
//	 {
//			if(std::string(p.getReferenceName()).find(name)!=std::string::npos||std::string(p.getReferenceName()).find("Y")!=std::string::npos)
//			{
//				addNumOverlappedReads();
//			}
//
//		 return 0;
//	 }
//	 int addPair(const bntseq_t *bns, bwa_seq_t *p, bwa_seq_t *q)
//	 {
//		 return 0;
//	 }
//	 int addSingle(const bntseq_t *bns, bwa_seq_t *p)
//	 {
//		 return 0;
//	 }
	 int addNumOverlappedReads()
	 {
		 numOverlappedReads++;
		 return 0;
	 }
	 int addNumPairOverlappedReads()
	 {
		 numPairOverlappedReads++;
		 return 0;
	 }
	 int addNumFullyIncludedReads()
	 {
		 numFullyIncludedReads++;
		 return 0;
	 }
	 int addNumFullyIncludedPairedReads()
	 {
		 numFullyIncludedPairedReads++;
		 return 0;
	 }
	 int getNumOverlappedReads()
	 {

		 return 		 numOverlappedReads;
	 }
	 int getNumPairOverlappedReads()
	 {

		 return 		 numPairOverlappedReads++;
	 }
	 int getNumFullyIncludedReads()
	 {

		 return 		 numFullyIncludedReads++;
	 }
	 int getNumFullyIncludedPairedReads()
	 {

		 return 		 numFullyIncludedPairedReads++;
	 }
	 std::string getName()
	 {
		 return name;
	 }
	 int getLength()
	 {
		 return length;
	 }
private:
		std::string name;
		int length;
		int numOverlappedReads;
		int numPairOverlappedReads;//same contig
		int numFullyIncludedReads;
		int numFullyIncludedPairedReads;//same contig
};

#endif /* UTILITY_ */

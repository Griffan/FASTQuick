/*
 * Utility
 *
 *  Created on: 2014Äê7ÔÂ8ÈÕ
 *      Author: Administrator
 *      For variables that need to be shared between indexing and mapping stages.
 */

#ifndef UTILITY_
#define UTILITY_

#include  <string>
#include  <unordered_map>
#include <vector>
#include <fstream>
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

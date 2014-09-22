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
		delete GC;
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

#endif /* UTILITY_ */

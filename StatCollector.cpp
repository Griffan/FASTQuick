/*
 * StatCollector.cpp
 *
 *  Created on: 2014Äê7ÔÂ20ÈÕ
 *      Author: Administrator
 */

#include "StatCollector.h"
#include "libbwa/bwase.h"
#include "libbwa/bamlite.h"


StatCollector::StatCollector()
{
	// TODO Auto-generated constructor stub

}
StatCollector::StatCollector(const string & OutFile)
{
	PositionTable.clear();
	VcfRecVec.clear();
}
int StatCollector::addAlignment(const string & PosName, const string & seq, const string & qual, const int & n_cigar, const bwa_cigar_t * cigar, const unsigned int & pos, const gap_opt_t* opt)
{
	size_t colonPos=PosName.find(":");
	size_t atPos = PosName.find("@");
	string Chrom=PosName.substr(0,colonPos);
	unsigned int refCoord= atoi(PosName.substr(colonPos+1,atPos-colonPos+1).c_str());// coordinate of variant site
	size_t slashPos = PosName.find("/");
	string ref = PosName.substr(atPos+1,slashPos-atPos-1);//ref allele
	unsigned int realCoord = refCoord - opt->flank_len +pos-1;//real coordinate of current reads on reference
	//unsigned int tmpCoord=realCoord;
	unsigned int tmpCycle=0;


	//cigar string
	for( int k=0;k< n_cigar;++k)
	{
	//int cop =cigar[k] & BAM_CIGAR_MASK; // operation
	//int cl = cigar[k] >> BAM_CIGAR_SHIFT; // length
	int cl= __cigar_len(cigar[k]);
	int cop= "MIDS"[__cigar_op(cigar[k])];
	switch (cop) {
	 case 'M':
	      //printf("M");
	      //printf("[%d-%d]", pos, pos + cl - 1);
		 for(int i=realCoord;i!=realCoord+cl-1;++i)
		 {
			 //tmpCycle=i-tmpCoord;
			 if( PositionTable[Chrom].find(i)!= PositionTable[Chrom].end())
			 {
				 index++;
				 SeqVec.push_back(string(""));
				 QualVec.push_back(string(""));
				 CycleVec.push_back(vector<unsigned int> (0));
				 SeqVec[index]+=seq[tmpCycle];
				 QualVec[index]+=qual[tmpCycle];
				 CycleVec[index].push_back(tmpCycle);
			 PositionTable[Chrom][i]=index;
			 }
			 else
			 {
				 unsigned int tmp_index=PositionTable[Chrom][i];
				 SeqVec[tmp_index]+=seq[tmpCycle];
				 QualVec[tmp_index]+=qual[tmpCycle];
				 CycleVec[tmp_index].push_back(tmpCycle);
			 }
			 tmpCycle++;
		 }
	      realCoord+=cl;
	      break;

//	 case BAM_CHARD_CLIP:
//	      printf("H");
//	      /* printf("[%d]", pos);  // No coverage */
//	      /* pos is not advanced by this operation */
//	      break;

	 case 'S':
//	      printf("S");
	      /* printf("[%d]", pos);  // No coverage */
	      /* pos is not advanced by this operation */
		 tmpCycle+=cl;
	      break;

	 case 'D':
//	      printf("D");
	      /* printf("[%d-%d]", pos, pos + cl - 1);  // Spans positions, No Coverage */
	     //tmpCoord+=cl;
		 realCoord+=cl;

	      break;

//	 case BAM_CPAD:
//	      printf("P");
	      /* printf("[%d-%d]", pos, pos + cl - 1); // Spans positions, No Coverage */
//	      pos+=cl;
//	      break;

	 case 'I':
//	      printf("I");
	      /* printf("[%d]", pos); // Special case - adds <cl> bp "throughput", but not genomic position "coverage" */
	      /* How you handle this is application dependent */
	      /* pos is not advanced by this operation */
	     //tmpCoord-=cl;
	     tmpCycle+=cl;
	      break;

//	 case BAM_CREF_SKIP:
//	      printf("S");
	      /* printf("[%d-%d]", pos, pos + cl - 1); /* Spans positions, No Coverage */
//	      pos+=cl;
//	      break;

	 default:
	      fprintf(stderr, "Unhandled cigar_op %d:%d\n", cop, cl);
	      printf("?");
	}
	}
	return 1;
}
int StatCollector::restoreVcfSites(const string & VcfPath,const gap_opt_t* opt)
{
	VcfHeader header;
	VcfFileReader reader;
	string SelectedSite=VcfPath+".SelectedSite.vcf";
	if(!reader.open(SelectedSite.c_str(),header))
	{
		reader.open(VcfPath.c_str(),header);
	}
	while(!reader.isEOF())
	{
		VcfRecord VcfLine;
		reader.readRecord(VcfLine);
		VcfRecVec.push_back(VcfLine);

		string chr(VcfLine.getChromStr());
		int pos=VcfLine.get1BasedPosition();
		VcfTable[chr][pos]=VcfRecVec.size()-1;

		int start=pos-opt->flank_len;
		int end=pos+opt->flank_len;
		for(int j=start;j!=end;++j)
		{
			VariantProxyTable[chr][j]=1;
		}
	}
	return 0;
}
int StatCollector::getDepthDist(const string & outputPath, const vector<uint64_t> & DepthDist)
{

	ofstream fout(outputPath+".DepthDist");
	for(int i=0;i!=DepthDist.size();++i)
	{
		fout<<i<<"\t"<<DepthDist[i]<<endl;
	}
	fout.close();
	return 0;
}
#define PHRED(x)	(-10)*log10(x)
int StatCollector::getEmpRepDist(const string & outputPath,const vector<uint64_t> & EmpRepDist,const vector<uint64_t> &misEmpRepDist)
{
	ofstream fout(outputPath+".EmpRepDist");
	for(int i=0;i!=EmpRepDist.size();++i)
	{
		fout<<i<<"\t"<<PHRED((double)(misEmpRepDist[i]+1)/(EmpRepDist[i]+2))<<endl;
	}
	fout.close();
	return 0;
}
int getEmpCycleDist(const string & outputPath,const vector<uint64_t> & EmpCycleDist, const vector<uint64_t> misEmpCycleDist)
{
	ofstream fout(outputPath+".EmpCycleDist");
	for(int i=0;i!=EmpCycleDist.size();++i)
	{
		fout<<i<<"\t"<<PHRED((double)(misEmpCycleDist[i]+1)/(EmpCycleDist[i]+2))<<endl;
	}
	fout.close();
	return 0;
}
int StatCollector::processCore(const string & statPrefix)
{
	vector<uint64_t> DepthDist(255);
	vector<uint64_t> EmpRepDist(255);
	vector<uint64_t> misEmpRepDist(255);
	vector<uint64_t> EmpCycleDist(255);
	vector<uint64_t> misEmpCycleDist(255);
	unsigned int VcfIndex=0;
	unsigned int PosIndex=0;
	for(string_map::iterator i=PositionTable.begin();i!=PositionTable.end();++i)//each chr
	{
		for(SeqPairIndex::iterator j=i->second.begin();j!=i->second.end();++j)//each site
		{
			VcfIndex=VcfTable[i->first][j->first];
			PosIndex=j->second;
			string RefStr(VcfRecVec[VcfIndex].getRefStr());
			if(RefStr.size()>1) continue;//only take care of snp for now

			for(int k=0;k!=SeqVec[PosIndex].size();++i)//every base in particular site
			{
			/************EmpRep**************************************************************/
			EmpRepDist[QualVec[PosIndex][k]]++;
			if(RefStr[0]!=SeqVec[PosIndex][k])
			{
				misEmpRepDist[QualVec[PosIndex][k]]++;
				misEmpCycleDist[CycleVec[PosIndex][k]]++;
			}
			/************EmpCycle**************************************************************/
			EmpCycleDist[CycleVec[PosIndex][k]]++;
			}
			/************DepthDist**************************************************************/
			if(VariantProxyTable[i->first][j->first])//if this is in a proxy region
			{
				if(SeqVec[j->second].size()>255)
					DepthDist[255]++;
				else
				DepthDist[SeqVec[j->second].size()]++;
			}
		}
	}


	getDepthDist(statPrefix,DepthDist);
	getEmpRepDist(statPrefix,EmpRepDist,misEmpRepDist);
	getEmpCycleDist(statPrefix,EmpCycleDist, misEmpCycleDist);
	return 0;
}
int StatCollector::outputPileup(const string & outputPath)
{
	ofstream fout(outputPath+".Pileup");
	for(string_map::iterator i=PositionTable.begin();i!=PositionTable.end();++i)//each chr
	{
		for(SeqPairIndex::iterator j=i->second.begin();j!=i->second.end();++j)//each site
		{
			fout<<i->first<<"\t"<<j->first<<"\t"<<SeqVec[j->second]<<"\t"<<QualVec[j->second]<<"\t";
					for(int k=0;k!=CycleVec[j->second].size();++k)
					fout<<CycleVec[j->second][k];
			fout <<endl;
		}
	}
	fout.close();
	return 0;
}

int findMaxAllele(size_t * a, size_t len)
{
	int max=0;
	int maxIndex=0;
	for(int i=0;i!=len;++i)
	{
		if(a[i]>max)
			{
			max=a[i];
			maxIndex=i;
			}
	}
	return maxIndex;
}
int countAllele(size_t*a, const string & seq)
{
	for(int i=0;i!=seq.size();++i)
	{
		if(seq[i]=='A')
			{
			a[0]++;continue;
			}
		if(seq[i]=='C')
			{
			a[1]++;continue;
			}
		if(seq[i]=='G')
			{
			a[2]++;continue;
			}
		if(seq[i]=='T')
			{
			a[3]++;continue;
			}
	}
	return 0;
}
#define REV_PHRED(x)	pow(10.0,x/(-10))
double calLikelihood(const string & seq, const string & qual,const char& maj, const char& min)
{
	double lik(0);
	if(maj==min)
	for(int i=0;i!=seq.size();++i)
	{
		if(seq[i]==maj)
		{
			lik+=log10(1-REV_PHRED(qual[i]));
		}
		else
		{
			lik+=log10(REV_PHRED(qual[i])/3);
		}
	}
	else
	{
		for(int i=0;i!=seq.size();++i)
		{
			if(seq[i]==maj|| seq[i]==min)
			{
				lik+=log10(1/2-REV_PHRED(qual[i])/3);
			}
			else
			{
				lik+=log10(REV_PHRED(qual[i])/3);
			}
		}
	}
	return lik;
}
int StatCollector::getGenoLikelihood(const string & outputPath)
{
	ofstream fout(outputPath+".likelihood");
	size_t numAllele[4]={0};
	char majAllele,minAllele;
	int maxIndex=0;
	for(string_map::iterator i=PositionTable.begin();i!=PositionTable.end();++i)//each chr
	{
		for(SeqPairIndex::iterator j=i->second.begin();j!=i->second.end();++j)//each site
		{
			countAllele(numAllele, SeqVec[j->second]);
			maxIndex=findMaxAllele(numAllele,4);
			majAllele="ACGT"[maxIndex];
			numAllele[maxIndex]=0;
			maxIndex=findMaxAllele(numAllele,4);
			minAllele="ACGT"[maxIndex];
			fout<<i->first<<"\t"<<j->first<<"\t"<<calLikelihood(SeqVec[j->second], QualVec[j->second],majAllele,minAllele)<<endl;

		}
	}
	fout.close();
	return 0;
}
StatCollector::~StatCollector()
{
	// TODO Auto-generated destructor stub
}


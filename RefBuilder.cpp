/*
 * RefBuilder.cpp
 *
 *  Created on: 2014Äê7ÔÂ9ÈÕ
 *      Author: Administrator
 */
#include "algorithm"
#include "RefBuilder.h"
#include "Utility.h"
#include "./samtools-faidx/faidx.h"
#include "./libStatGen-master/vcf/VcfFileReader.h"
#include "./libStatGen-master/vcf/VcfHeader.h"
#include "./libStatGen-master/vcf/VcfRecord.h"
using namespace std;

#define DEBUG 0

RefBuilder::RefBuilder()
{
	// TODO Auto-generated constructor stub

}

RefBuilder::RefBuilder(string VcfPath,string RefPath, string MaskPath, const gap_opt_t* opt)
{
	cerr<<"Initialization of RefBwt..."<<endl;
		//read in ref.fa and ref.fai
		//string RefFaiPath=RefPath+".fai";
		faidx_t * seq;
		seq=fai_load(RefPath.c_str());
		cerr<<"Loading Ref fai file done!\n";

		faidx_t * FastaMask;
		if(MaskPath!="Empty")
		{
		FastaMask=fai_load(MaskPath.c_str());
		cerr<<"Loading Mask fai file done!\n";
		}
		//read in vcf, hm3 sites
		VcfHeader header;
		VcfFileReader reader;

		string SelectedSite=VcfPath+".SelectedSite.vcf";
		InputFile FoutSelectedSite(SelectedSite.c_str(),"w");


		reader.open(VcfPath.c_str(),header);
		char region[128];
		unsigned int nseqs(0);
		unsigned int nmarker(0);
		while(!reader.isEOF())
		{
			if(MaskPath!="Empty"&&nmarker>100000)
			{
				break;
			}
			VcfRecord VcfLine;
			//vector<char> SeqVecIdx;
			reader.readRecord(VcfLine);
			if(VcfLine.getNumRefBases()  != 1)// filtering abnormal sites
				continue;

			string Chrom(VcfLine.getChromStr());
			int Position=VcfLine.get1BasedPosition();
			int dummy;
			sprintf(region, "%s:%d-%d", Chrom.c_str(),Position-opt->flank_len,Position+opt->flank_len);
			if(MaskPath!="Empty")
			{
			string MaskSeq(fai_fetch(FastaMask, region, &dummy));
			size_t n = std::count(MaskSeq.begin(), MaskSeq.end(), 'P');
			if(double(n)/MaskSeq.size()<0.9) continue;
			VcfLine.write(&FoutSelectedSite,1);
			}

			string FetchedSeq(fai_fetch(seq, region, &dummy));

			//vector<string> SeqVec;


			SeqVec.push_back(FetchedSeq);
			//SeqVecIdx.push_back(nseqs);


			//sprintf(region, "%s:%d@%s", Chrom.c_str(),Position,VcfLine.getAlleles(0));
			//string Coord(region);
			//string RefAllele(VcfLine.getAlleles(0));
			//RefTableIndex.insert ( make_pair(Coord,nseqs));
			//nseqs++;

			if(DEBUG)
			{
				string AltAllele(VcfLine.getAlleles(0));
				string RefSeq=FetchedSeq.substr(0,opt->flank_len)+AltAllele+FetchedSeq.substr(opt->flank_len+1,opt->flank_len);
				if(RefSeq != FetchedSeq)
				{
					cerr<<"Coordinate problem!!!!!!"<<endl<<"Number of Alt:"<<VcfLine.getNumAlts()<<endl;
					//cerr<<"Region: "<<Coord<<endl;
					cerr<<VcfLine.getAlleles(0)<<endl;
					cerr<<VcfLine.getAlleles(1)<<endl;
					cerr<<RefSeq<<endl;
					cerr<<FetchedSeq<<endl;
					exit(1);
				}
			}

			/*if(VcfLine.getNumAlts()>=1)
			for(unsigned int i=1;i<VcfLine.getNumAlts()+1;i++)
			{
				string AltAllele(VcfLine.getAlleles(i));
				SeqVec.push_back(FetchedSeq.substr(0,opt->flank_len)+AltAllele+FetchedSeq.substr(opt->flank_len+1,opt->flank_len));// position 500 was replaced by AltAllele

				//SeqVecIdx.push_back(nseqs);
				sprintf(region, "%s:%d@%s", Chrom.c_str(),Position,AltAllele.c_str());

				string Coord(region);
				RefTableIndex.insert ( make_pair(Coord,nseqs));
				nseqs++;
				//cerr<<"Insert alt  allele at:"<<Coord<<endl;
				//cout<<Coord<<"\t"<<SeqVec[nseqs-1]<<endl;
			}
			else
			{
				//cerr<<"Only ref allele at:"<<Coord<<endl;
			}*/
			SeqVec.push_back(FetchedSeq.substr(0,opt->flank_len)+string("N")+FetchedSeq.substr(opt->flank_len+1,opt->flank_len));// position 500 was replaced by AltAllele
			sprintf(region, "%s:%d@", Chrom.c_str(),Position);
			string AltAllele(VcfLine.getAlleles(0));
			sprintf(region, "%s%s",region,AltAllele.c_str());
			for(unsigned int i=1;i!=VcfLine.getNumAlts();i++)
			{
			sprintf(region, "%s/%s", region,string((VcfLine.getAlleles(i))).c_str());
			}
			//string Coord(region);
			RefTableIndex.insert ( make_pair(string(region),nseqs));
			nseqs++;
			nmarker++;


		}
		FoutSelectedSite.ifclose();
}

RefBuilder::~RefBuilder()
{
	// TODO Auto-generated destructor stub
	SeqVec.clear();
	RefTableIndex.clear();
}


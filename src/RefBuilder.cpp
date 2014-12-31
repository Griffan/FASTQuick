/*
 * RefBuilder.cpp
 *
 *  Created on: 2014Äê7ÔÂ9ÈÕ
 *      Author: Administrator
 */
#include "algorithm"
#include "RefBuilder.h"
#include "Utility.h"
#include "../misc/faidx.h"
#include "../libmpu/Error.h"
#include "VcfFileReader.h"
#include "VcfHeader.h"
#include "VcfRecord.h"
#include <fstream>
using namespace std;

#define DEBUG 0
//extern string Prefix;
extern void notice(const char*,...);
extern void warning(const char*,...);
extern void error(const char*,...);
RefBuilder::RefBuilder()
{
  // TODO Auto-generated constructor stub

}

RefBuilder::RefBuilder(const string& VcfPath, const string& RefPath, const string& DBsnpPath, const string& MaskPath, const gap_opt_t* opt)//, unordered_map<string,bool>& longRefTable)
{
  notice("Initialization of RefBwt...\n");
  //read in ref.fa and ref.fai
  //string RefFaiPath=RefPath+".fai";
  faidx_t * seq;
  seq=fai_load(RefPath.c_str());
  cerr<<endl;
  notice("Loading Ref fai file done!\n");

  faidx_t * FastaMask=0;
  if(MaskPath!="Empty")
    {
      FastaMask=fai_load(MaskPath.c_str());
      notice("Loading Mask fai file done!\n");
    }
  //read in vcf, hm3 sites
  VcfHeader header;
  VcfFileReader reader;

  string SelectedSite = RefPath + ".SelectedSite.vcf";
  InputFile FoutSelectedSite(SelectedSite.c_str(),"w");
  string GCpath = RefPath + ".gc";
  ofstream FGC(GCpath,ios_base::binary);
  string BedPath = RefPath + ".bed";
  ofstream BedFile(BedPath);
  //FGC.write((char*)opt->num_variant_short,sizeof(int));
  //_GCstruct * GCstruct=new _GCstruct [opt->num_variant_long+opt->num_variant_short];

  reader.open(VcfPath.c_str(),header);
  header.write(&FoutSelectedSite);
  char region[128];
  unsigned int nseqs(0);
  unsigned int nmarker(0);
  int last_pos=0;
  string last_chr;

  //int num_so_far=0;
  while(!reader.isEOF())
    {
      if(nmarker>=opt->num_variant_short)// for short region
        {
          //cerr<<"the nmarker:"<<nmarker<<endl;
          break;
        }
      VcfRecord VcfLine;
      //vector<char> SeqVecIdx;
      reader.readRecord(VcfLine);
      if(VcfLine.getNumRefBases()  != 1)// filtering abnormal sites
        continue;

      string Chrom(VcfLine.getChromStr());
      if(Chrom=="X"||Chrom=="x"||Chrom=="chrX"||Chrom=="chrx"||Chrom=="Y"||Chrom=="y"||Chrom=="chrY"||Chrom=="chry"||Chrom=="MT"||Chrom=="mt")
        continue;
      int Position=VcfLine.get1BasedPosition();
      if(Chrom==last_chr&&abs(Position-last_pos)<300)
        continue;
      int dummy;
      sprintf(region, "%s:%d-%d", Chrom.c_str(),Position-opt->flank_len,Position+opt->flank_len);
      if(MaskPath!="Empty")
        {
          //cerr<<"region:"<<region<<endl;
          string MaskSeq(fai_fetch(FastaMask, region, &dummy));
          size_t n = std::count(MaskSeq.begin(), MaskSeq.end(), 'P');
          if(double(n)/MaskSeq.size()<0.9)
            continue;

          if(!VcfLine.write(&FoutSelectedSite,1))
            {
              warning("Writing retained sites failed!\n");
              exit(1);
            }
        }

      string FetchedSeq(fai_fetch(seq, region, &dummy));

      {//Calculate GC content
        _GCstruct GCstruct(opt->flank_len*2+1);
        for(int i=Position-opt->flank_len,t=0;i!=Position+opt->flank_len+1;++i,++t)
          {
            sprintf(region, "%s:%d-%d",Chrom.c_str(), i-50,i+49);
            string Window(fai_fetch(seq, region, &dummy));
            int total=0;
            for(unsigned int j=0;j!=Window.size();++j)
              {
                if(Window[j]=='G'||Window[j]=='C'||Window[j]=='g'||Window[j]=='c')
                  total++;
              }

            //strcpy(GCstruct[num_so_far].chrom,Chrom.c_str());//potential overflow
            //GCstruct.pos[t]=i;
            GCstruct.GC[t]=total;
            //FGC.write((char*)&GCstruct,sizeof(GCstruct));
          }
        //cerr<<Chrom<<"\t"<<Position<<"\t"<<num_so_far<<endl;
        //num_so_far++;
        GCstruct.write(FGC);
      }


      if(DEBUG)
        {
          string AltAllele(VcfLine.getAlleles(0));
          string RefSeq=FetchedSeq.substr(0,opt->flank_len)+AltAllele+FetchedSeq.substr(opt->flank_len+1,opt->flank_len);
          if(RefSeq != FetchedSeq)
            {
              cerr<<"WARNING:Coordinate problem!!!!!!"<<endl<<"Number of Alt:"<<VcfLine.getNumAlts()<<endl;
              //cerr<<"Region: "<<Coord<<endl;
              cerr<<VcfLine.getAlleles(0)<<endl;
              cerr<<VcfLine.getAlleles(1)<<endl;
              cerr<<RefSeq<<endl;
              cerr<<FetchedSeq<<endl;
              exit(1);
            }
        }

      SeqVec.push_back(FetchedSeq.substr(0,opt->flank_len)+string("N")+FetchedSeq.substr(opt->flank_len+1,opt->flank_len));// position 500 was replaced by AltAllele
      sprintf(region, "%s:%d@%s/%s", Chrom.c_str(),Position, VcfLine.getRefStr(),VcfLine.getAltStr());
      RefTableIndex.insert( make_pair(string(region),nseqs));
      sprintf(region, "%s\t%d\t%d", Chrom.c_str(),Position-opt->flank_len, Position+opt->flank_len);
      BedFile<<region<<endl;
      nseqs++;
      nmarker++;
      last_pos=Position;
      last_chr=Chrom;
    }

  //cerr<<"******************************************************************"<<endl;
  //reader.open(VcfPath.c_str(),header);
  // below is for long ref
  nmarker=0;

  while(!reader.isEOF())
    {
      if(nmarker>=opt->num_variant_long)// for long region
        {
          break;
        }
      VcfRecord VcfLine;
      //vector<char> SeqVecIdx;
      reader.readRecord(VcfLine);
      if(VcfLine.getNumRefBases()  != 1)// filtering abnormal sites
        continue;

      string Chrom(VcfLine.getChromStr());
      if(Chrom=="X"||Chrom=="x"||Chrom=="chrX"||Chrom=="chrx"||Chrom=="Y"||Chrom=="y"||Chrom=="chrY"||Chrom=="chry"||Chrom=="MT"||Chrom=="mt")
        continue;
      int Position=VcfLine.get1BasedPosition();
      if(Chrom==last_chr&&abs(Position-last_pos)<300)
        continue;
      int dummy;
      sprintf(region, "%s:%d-%d", Chrom.c_str(),Position-opt->flank_long_len,Position+opt->flank_long_len);
      if(MaskPath!="Empty")
        {
          //	cerr<<"region:"<<region<<endl;
          string MaskSeq(fai_fetch(FastaMask, region, &dummy));
          size_t n = std::count(MaskSeq.begin(), MaskSeq.end(), 'P');
          if(double(n)/MaskSeq.size()<0.9)
            continue;

          if(!VcfLine.write(&FoutSelectedSite,1))
            {
              warning("Writing retained sites failed!\n");
              exit(1);
            }
        }

      string FetchedSeq(fai_fetch(seq, region, &dummy));

      {//Calculate GC content
        _GCstruct GCstruct(opt->flank_long_len*2+1);
        for(int i=Position-opt->flank_long_len,t=0;i!=Position+opt->flank_long_len+1;++i,++t)
          {
            sprintf(region, "%s:%d-%d",Chrom.c_str(), i-50,i+49);
            string Window(fai_fetch(seq, region, &dummy));
            int total=0;
            for(unsigned int j=0;j!=Window.size();++j)
              {
                if(Window[j]=='G'||Window[j]=='C'||Window[j]=='g'||Window[j]=='c')
                  total++;
              }
            //cerr<<Chrom<<"\t"<<i<<"\t"<<total<<endl;
            //strcpy(GCstruct[num_so_far].chrom,Chrom.c_str());//potential overflow
            //GCstruct.pos[t]=i;
            GCstruct.GC[t]=total;
          }
        GCstruct.write(FGC);
      }

      SeqVec.push_back(FetchedSeq.substr(0,opt->flank_long_len)+string("N")+FetchedSeq.substr(opt->flank_long_len+1,opt->flank_long_len));// position 500 was replaced by AltAllele
      sprintf(region, "%s:%d@%s/%s|L", Chrom.c_str(),Position, VcfLine.getRefStr(),VcfLine.getAltStr());
      RefTableIndex.insert( make_pair(string(region),nseqs));
      //sprintf(region, "%s:%d@%s/%s|L", Chrom.c_str(),Position, VcfLine.getRefStr(),VcfLine.getAltStr());
      //longRefTable.insert(make_pair(string(region),true));
      //VariantLongTable[string(region)]=1;
      sprintf(region, "%s\t%d\t%d", Chrom.c_str(),Position-opt->flank_long_len, Position+opt->flank_long_len);
      BedFile<<region<<endl;
      nseqs++;
      nmarker++;
      last_pos=Position;
      last_chr=Chrom;

    }

  int max_XorYmarker(0);
  if(opt->num_variant_short>=100000)
    max_XorYmarker=3000;
  else if(opt->num_variant_short>=10000)
    max_XorYmarker=300;
  else
    max_XorYmarker=100;
  // choosing chrX and chrY
  int Ynmarker=0,Xnmarker=0, chr_flag(-1);

  while(!reader.isEOF())
    {
      if(Ynmarker>=max_XorYmarker&&Xnmarker>=max_XorYmarker)// for long region
        {
          break;
        }
      VcfRecord VcfLine;
      //vector<char> SeqVecIdx;
      reader.readRecord(VcfLine);
      if(VcfLine.getNumRefBases()  != 1)// filtering abnormal sites
        continue;

      string Chrom(VcfLine.getChromStr());
      if(Chrom=="X"||Chrom=="x"||Chrom=="chrX"||Chrom=="chrx")
        {
          if(Xnmarker>=max_XorYmarker)
            continue;
          chr_flag=0;
        }
      else if(Chrom=="Y"||Chrom=="y"||Chrom=="chrY"||Chrom=="chry")
        {
          if(Ynmarker>=max_XorYmarker)
            continue;
          chr_flag=1;
        }
      else
        continue;


      int Position=VcfLine.get1BasedPosition();
      if(Chrom==last_chr&&abs(Position-last_pos)<300)
        continue;
      int dummy;
      sprintf(region, "%s:%d-%d", Chrom.c_str(),Position-250,Position+250);
      if(MaskPath!="Empty")
        {
          //	cerr<<"region:"<<region<<endl;
          string MaskSeq(fai_fetch(FastaMask, region, &dummy));
          size_t n = std::count(MaskSeq.begin(), MaskSeq.end(), 'P');
          if(double(n)/MaskSeq.size()<0.9)
            continue;

          if(!VcfLine.write(&FoutSelectedSite,1))
            {
              warning("Writing retained sites failed!\n");
              exit(1);
            }
        }

      string FetchedSeq(fai_fetch(seq, region, &dummy));

      {//Calculate GC content
        _GCstruct GCstruct(250*2+1);
        for( int i=Position-250,t=0;i!=Position+250+1;++i,++t)
          {
            sprintf(region, "%s:%d-%d",Chrom.c_str(), i-50,i+49);
            string Window(fai_fetch(seq, region, &dummy));
            int total=0;
            for(unsigned int j=0;j!=Window.size();++j)
              {
                if(Window[j]=='G'||Window[j]=='C'||Window[j]=='g'||Window[j]=='c')
                  total++;
              }
            //FGC<<Chrom<<"\t"<<i<<"\t"<<total<<endl;
            //strcpy(GCstruct[num_so_far].chrom,Chrom.c_str());//potential overflow
            //GCstruct.pos[t]=i;
            GCstruct.GC[t]=total;
          }
        GCstruct.write(FGC);
      }

      SeqVec.push_back(FetchedSeq.substr(0,250)+string("N")+FetchedSeq.substr(250+1,250));// position 500 was replaced by AltAllele
      sprintf(region, "%s:%d@%s/%s|L", Chrom.c_str(),Position, VcfLine.getRefStr(),VcfLine.getAltStr());
      RefTableIndex.insert( make_pair(string(region),nseqs));
      //sprintf(region, "%s:%d@%s/%s|L", Chrom.c_str(),Position, VcfLine.getRefStr(),VcfLine.getAltStr());
      //longRefTable.insert(make_pair(string(region),true));
      //VariantLongTable[string(region)]=1;
      sprintf(region, "%s\t%d\t%d", Chrom.c_str(),Position-250, Position+250);
      BedFile<<region<<endl;
      nseqs++;
      if(chr_flag==0)
        Xnmarker++;
      else
        Ynmarker++;

      last_pos=Position;
      last_chr=Chrom;

    }
  //cerr<<"the total gc sites:"<<num_so_far<<endl;
  //FGC.write((char*)GCstruct,(opt->num_variant_long*(4*opt->flank_len+1)+opt->num_variant_short*(2*opt->flank_len+1))*sizeof(_GCstruct));
  //int total=(opt->num_variant_long*(4*opt->flank_len+1)+opt->num_variant_short*(2*opt->flank_len+1))*sizeof(_GCstruct);
  //FGC.write((char*) &nmarker, sizeof( unsigned int));
  FGC.close();
  BedFile.close();
  char cmdline[2048];
  sprintf(cmdline, "tabix  -h -B %s %s  > %s.dpSNP.subset.vcf", DBsnpPath.c_str(), BedPath.c_str(), RefPath.c_str());
  int ret=system(cmdline);
  if(ret!=0)
    {
      warning("Building dbsnp subset.vcf failed!\n");
      exit(1);
    }
  FoutSelectedSite.ifclose();

  sprintf(cmdline, "(grep ^# %s.SelectedSite.vcf && grep -v  ^# %s.SelectedSite.vcf|sort -k1,1 -k2,2n) >%s.SelectedSite.vcf.tmp", RefPath.c_str(), RefPath.c_str(), RefPath.c_str());
  if(system(cmdline)!=0)
    {
      warning("Call command line:\n%s\nfailed!\n",cmdline);
	  exit(EXIT_FAILURE);
    }
  sprintf(cmdline, "mv %s.SelectedSite.vcf.tmp %s.SelectedSite.vcf", RefPath.c_str(), RefPath.c_str());
  if(system(cmdline)!=0)
    {
      warning("Call command line:\n%s\nfailed!\n",cmdline);
	  exit(EXIT_FAILURE);
    }
  sprintf(cmdline, "bgzip -f %s.SelectedSite.vcf", RefPath.c_str());
  if(system(cmdline)!=0)
    {
      warning("Call command line:\n%s\nfailed!\n",cmdline);
	  exit(EXIT_FAILURE);
    }
  sprintf(cmdline, "tabix -pvcf %s.SelectedSite.vcf.gz", RefPath.c_str());
  if(system(cmdline)!=0)
    {
      warning("Call command line:\n%s\nfailed!\n",cmdline);
	  exit(EXIT_FAILURE);
    }
  reader.close();

}

RefBuilder::~RefBuilder()
{
  // TODO Auto-generated destructor stub
  SeqVec.clear();
  RefTableIndex.clear();
}


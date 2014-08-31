/*
 * StatCollector.cpp
 *
 *  Created on: 2014Äê7ÔÂ20ÈÕ
 *      Author: Administrator
 */

#include "StatCollector.h"
#include "libbwa/bwase.h"
#include "libbwa/bamlite.h"
#include <fstream>
#include <iostream>
#include <sstream>

StatCollector::StatCollector()
{
	// TODO Auto-generated constructor stub
	cerr << "Using default initializer..." << endl;
	PositionTable.clear();
	VcfRecVec.clear();
	index = 0;
	DepthDist = vector<uint64_t>(256, 0);
	GCDist = vector<uint64_t>(256, 0);
	EmpRepDist = vector<uint64_t>(256, 0);
	misEmpRepDist = vector<uint64_t>(256, 0);
	EmpCycleDist = vector<uint64_t>(256, 0);
	misEmpCycleDist = vector<uint64_t>(256, 0);
}
StatCollector::StatCollector(const string & OutFile)
{
	PositionTable.clear();
	VcfRecVec.clear();
	index = 0;
	DepthDist = vector<uint64_t>(256, 0);
	GCDist = vector<uint64_t>(256, 0);
	EmpRepDist = vector<uint64_t>(256, 0);
	misEmpRepDist = vector<uint64_t>(256, 0);
	EmpCycleDist = vector<uint64_t>(256, 0);
	misEmpCycleDist = vector<uint64_t>(256, 0);
}
//int StatCollector::addAlignment(const string & PosName, const string & seq, const string & qual, const int & n_cigar, const bwa_cigar_t * cigar, const unsigned int & pos, const gap_opt_t* opt)
int StatCollector::addAlignment(const bntseq_t *bns, bwa_seq_t *p,
		const gap_opt_t* opt)
{
	int seqid(0), j(0);

	if (p->type == BWA_TYPE_NO_MATCH)
	{
		//j = 1;
		return false;
	}
	else
	{
		j = pos_end(p) - p->pos;
	}
	bns_coor_pac2real(bns, p->pos, j, &seqid);
	if (p->type != BWA_TYPE_NO_MATCH
			&& p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len)
	{
		return false; //this alignment bridges two adjacent reference sequences
	}
	if (/*p->type == BWA_TYPE_NO_MATCH || */(p->count) > 60)
	{
		return false;
	}

	string seq, qual, newSeq, newQual;

	if (p->strand == 0)
		for (j = 0; j != p->full_len; ++j)
		{
			seq += "ACGTN"[(int) (p)->seq[j]];
			qual += (char) p->qual[j];
		}
	else
		for (j = 0; j != p->full_len; ++j)
		{
			seq += "TGCAN"[(int) (p)->seq[p->full_len - 1 - j]];
			qual += (char) p->qual[p->full_len - 1 - j];
		}
	//if (p->strand) seq_reverse(p->len, p->qual, 0);
	//fprintf(stderr,"%s\t%s\t%d\n",p->name,p->md,p->count);

	//ConstructFakeSeqQual(seq,qual,p->n_cigar,p->cigar,newSeq,newQual);
	string RefSeq(seq), MD(p->md);
	int last(0), total_len(0);
	for (int i = 0; i != MD.size(); ++i)
		if (isdigit(MD[i]))
			continue;
		else if (MD[i] == '^')
		{
			i++;
			while (!isdigit(MD[i]))
			{
				i++;
				total_len++;
			}
			last = i;
		}
		else
		{
			int len = atoi(MD.substr(last, i - last).c_str()) + 1;
			total_len += len;
			RefSeq[total_len - 1] = MD[i];
			last = i + 1;
		}

	string PosName = string(bns->anns[seqid].name);
	int n_cigar = p->n_cigar;
	bwa_cigar_t* cigar = p->cigar;
	int pos = (int) (p->pos - bns->anns[seqid].offset + 1);

	size_t colonPos = PosName.find(":");
	size_t atPos = PosName.find("@");
	string Chrom = PosName.substr(0, colonPos);
	unsigned int refCoord = atoi(
			PosName.substr(colonPos + 1, atPos - colonPos + 1).c_str()); // coordinate of variant site
	size_t slashPos = PosName.find("/");
	string ref = PosName.substr(atPos + 1, slashPos - atPos - 1); //ref allele
	unsigned int realCoord = refCoord - opt->flank_len + pos - 1; //real coordinate of current reads on reference

	unsigned int tmpCycle = 0;

	unsigned int tmp_index = 0;
	//cigar string
	for (int k = 0; k < n_cigar; ++k)
	{

		int cl = __cigar_len(cigar[k]);
		int cop = "MIDS"[__cigar_op(cigar[k])];
		switch (cop)
		{
		case 'M':
			//printf("M");
			//printf("[%d-%d]", pos, pos + cl - 1);
			if (PositionTable.find(Chrom) != PositionTable.end()) //chrom exists
				for (int i = realCoord; i != realCoord + cl - 1 + 1;
						++i, tmpCycle++)
				{
					std::pair<unordered_map<int, bool>::iterator, bool> ret =
							VariantProxyTable[Chrom].insert(make_pair(i, 1));
					if (!ret.second) //if insert failed, this coord exists
					{
						tmp_index = PositionTable[Chrom][i];
						SeqVec[tmp_index] += seq[tmpCycle];
						QualVec[tmp_index] += qual[tmpCycle];
						CycleVec[tmp_index].push_back(tmpCycle);
						MaqVec[tmp_index].push_back(p->mapQ);
						StrandVec[tmp_index].push_back(p->strand);
						EmpRepDist[qual[tmpCycle]]++;
						if (RefSeq[tmpCycle] != seq[tmpCycle])
						{
							misEmpRepDist[qual[tmpCycle]]++;
							misEmpCycleDist[tmpCycle]++;
						}
						/************EmpCycle**************************************************************/
						EmpCycleDist[tmpCycle]++;
					}
					else //coord not exists
					{
						SeqVec.push_back(string(""));
						QualVec.push_back(string(""));
						CycleVec.push_back(vector<unsigned int>(0));
						MaqVec.push_back(vector<unsigned char>(0));
						StrandVec.push_back(vector<bool>(0));
						SeqVec[index] += seq[tmpCycle];
						QualVec[index] += qual[tmpCycle];
						CycleVec[index].push_back(tmpCycle);
						MaqVec[index].push_back(p->mapQ);
						StrandVec[index].push_back(p->strand);
						PositionTable[Chrom][i] = index;
						index++;
						EmpRepDist[qual[tmpCycle]]++;
						if (RefSeq[tmpCycle] != seq[tmpCycle])
						{
							misEmpRepDist[qual[tmpCycle]]++;
							misEmpCycleDist[tmpCycle]++;
						}
						/************EmpCycle**************************************************************/
						EmpCycleDist[tmpCycle]++;
					}
				}
			else // chrom not exists
			{
				for (int i = realCoord; i != realCoord + cl - 1 + 1;
						++i, tmpCycle++)
				{
					VariantProxyTable[Chrom][i] = 1;
					SeqVec.push_back(string(""));
					QualVec.push_back(string(""));
					CycleVec.push_back(vector<unsigned int>(0));
					MaqVec.push_back(vector<unsigned char>(0));
					StrandVec.push_back(vector<bool>(0));

					//tmp_index = PositionTable[Chrom][i];
					//cerr << p->name<<"\trealCoord:"<<i<<"\ttmp index :" << tmp_index <<"\ttmp cycle:"<<tmpCycle << endl;
					SeqVec[index] += seq[tmpCycle];
					QualVec[index] += qual[tmpCycle];
					CycleVec[index].push_back(tmpCycle);
					MaqVec[index].push_back(p->mapQ);
					StrandVec[index].push_back(p->strand);
					PositionTable[Chrom][i] = index;
					index++;

					EmpRepDist[qual[tmpCycle]]++;
					if (RefSeq[tmpCycle] != seq[tmpCycle])
					{
						misEmpRepDist[qual[tmpCycle]]++;
						misEmpCycleDist[tmpCycle]++;
					}
					/************EmpCycle**************************************************************/
					EmpCycleDist[tmpCycle]++;
				}
			}

		realCoord += cl;
		break;

//	 case BAM_CHARD_CLIP:
//	      printf("H");
//	      /* printf("[%d]", pos);  // No coverage
//	      /* pos is not advanced by this operation
//	      break;

		case 'S':
//	      printf("S");
//	      /* printf("[%d]", pos);  // No coverage
//	      /* pos is not advanced by this operation
		tmpCycle += cl;
		break;

		case 'D':
//	      printf("D");
		/* printf("[%d-%d]", pos, pos + cl - 1);  // Spans positions, No Coverage*/

		realCoord += cl;

		break;

//	 case BAM_CPAD:
//	      printf("P");
		/* printf("[%d-%d]", pos, pos + cl - 1); // Spans positions, No Coverage*/
//	      pos+=cl;
//	      break;
		case 'I':
//	      printf("I");
		/* printf("[%d]", pos); // Special case - adds <cl> bp "throughput", but not genomic position "coverage"*/
		/* How you handle this is application dependent*/
		/* pos is not advanced by this operation*/

		tmpCycle += cl;
		break;

//	 case BAM_CREF_SKIP:
//	      printf("S");
		//   printf("[%d-%d]", pos, pos + cl - 1); /* Spans positions, No Coverage*/
//	      pos+=cl;
//	      break;

		default:
		fprintf(stderr, "Unhandled cigar_op %d:%d\n", cop, cl);
		printf("?");
	}
}
return true;
}
int StatCollector::restoreVcfSites(const string & VcfPath, const gap_opt_t* opt)
{
VcfHeader header;
VcfFileReader reader;
string SelectedSite = VcfPath + ".SelectedSite.vcf";
if (!reader.open(SelectedSite.c_str(), header))
{
	reader.open(VcfPath.c_str(), header);
}
while (!reader.isEOF())
{
	VcfRecord* VcfLine = new VcfRecord;

	reader.readRecord(*VcfLine);
	VcfRecVec.push_back(VcfLine);

	string chr(VcfLine->getChromStr());
	int pos = VcfLine->get1BasedPosition();
	VcfTable[chr][pos] = VcfRecVec.size() - 1;
/*
	int start = pos - opt->flank_len;
	int end = pos + opt->flank_len;

	if (VariantProxyTable.find(chr) != VariantProxyTable.end())
		for (int j = start; j != end + 1; ++j)
		{
			std::pair<unordered_map<int, bool>::iterator, bool> ret =
					VariantProxyTable[chr].insert(make_pair(j, 1));
			if (!ret.second)
				continue;
			SeqVec.push_back(string(""));
			QualVec.push_back(string(""));
			CycleVec.push_back(vector<unsigned int>(0));
			MaqVec.push_back(vector<unsigned char>(0));
			StrandVec.push_back(vector<bool>(0));
			PositionTable[chr][j] = index;
			index++;
		}
	else
	{
		for (int j = start; j != end + 1; ++j)
		{
			VariantProxyTable[chr][j] = 1;
			SeqVec.push_back(string(""));
			QualVec.push_back(string(""));
			CycleVec.push_back(vector<unsigned int>(0));
			MaqVec.push_back(vector<unsigned char>(0));
			StrandVec.push_back(vector<bool>(0));
			PositionTable[chr][j] = index;
			index++;
		}
	}*/

}
string GCpath=VcfPath+".gc";
ifstream FGC(GCpath);
string line,Chrom,PosStr,GCStr;
while(getline(FGC,line))
{

	stringstream ss(line);
	ss>>Chrom>>PosStr>>GCStr;
	//cerr<<Chrom<<"\t"<<PosStr<<"\t"<<GCStr<<endl;
	GC[Chrom][atoi(PosStr.c_str())]=(unsigned int)floor(atof(GCStr.c_str())*100+0.5);
}
FGC.close();
return 0;
}
int StatCollector::getDepthDist(const string & outputPath,
	const vector<uint64_t> & DepthDist)
{

ofstream fout(outputPath + ".DepthDist");
for (int i = 0; i != DepthDist.size(); ++i)
{
	fout << i << "\t" << DepthDist[i] << endl;
}
fout.close();
return 0;
}

int StatCollector::getGCDist(const string & outputPath,const vector<uint64_t> & PosNum, const vector<uint64_t> & GCDist)
{
	ofstream fout(outputPath + ".GCDist");
	for (int i = 0; i != GCDist.size(); ++i)
	{
		fout << i << "\t";
		fout<< PosNum[i]==0?0:(GCDist[i]/double(PosNum[i]) );
		fout<< endl;
	}
	fout.close();
	return 0;
}

#define PHRED(x)	(-10)*log10(x)
int StatCollector::getEmpRepDist(const string & outputPath,
	const vector<uint64_t> & EmpRepDist, const vector<uint64_t> &misEmpRepDist)
{
ofstream fout(outputPath + ".EmpRepDist");
for (int i = 0; i != EmpRepDist.size(); ++i)
{
	fout << i << "\t" << (misEmpRepDist[i]) << "\t" << (EmpRepDist[i]) << "\t"
			<< PHRED((double )(misEmpRepDist[i] + 1) / (EmpRepDist[i] + 2))
			<< endl;
}
fout.close();
return 0;
}
int StatCollector::getEmpCycleDist(const string & outputPath,
	const vector<uint64_t> & EmpCycleDist,
	const vector<uint64_t> misEmpCycleDist)
{
ofstream fout(outputPath + ".EmpCycleDist");
for (int i = 0; i != EmpCycleDist.size(); ++i)
{
	fout << i << "\t"
			<< PHRED((double )(misEmpCycleDist[i] + 1) / (EmpCycleDist[i] + 2))
			<< endl;
}
fout.close();
return 0;
}
int StatCollector::processCore(const string & statPrefix)
{

unsigned int VcfIndex = 0;
unsigned int PosIndex = 0;
vector<uint64_t> PosNum(256,0);
for (string_map::iterator i = PositionTable.begin(); i != PositionTable.end();
		++i) //each chr
{
	for (SeqPairIndex::iterator j = i->second.begin(); j != i->second.end();
			++j) //each site
	{

		/************DepthDist**************************************************************/
		//cerr<<"We got here:"<<i->first<<"\t"<<j->first<<"\tindex:"<<j->second<<"\t"<<VariantProxyTable[i->first][j->first]<<"\t"<<SeqVec[j->second] <<"\tend"<<endl;
	//	if (VariantProxyTable[i->first][j->first] == 1) //if this is in a proxy region
		{
			if (SeqVec[j->second].size() > 255)
				DepthDist[255]++;
			else
				DepthDist[SeqVec[j->second].size()]++;
		}
		GCDist[GC[i->first][j->first]]+=SeqVec[j->second].size();
		PosNum[GC[i->first][j->first]]++;
		//cerr<<"We got here:"<<i->first<<"\t"<<j->first<<"\t"<<GC[i->first][j->first]<<endl;
	}
}

getDepthDist(statPrefix, DepthDist);
getGCDist(statPrefix, PosNum, GCDist);
getEmpRepDist(statPrefix, EmpRepDist, misEmpRepDist);
getEmpCycleDist(statPrefix, EmpCycleDist, misEmpCycleDist);
return 0;
}
int StatCollector::outputPileup(const string & outputPath)
{
ofstream fout(outputPath + ".Pileup");
for (string_map::iterator i = PositionTable.begin(); i != PositionTable.end();
		++i) //each chr
{
	for (SeqPairIndex::iterator j = i->second.begin(); j != i->second.end();
			++j) //each site
	{
		fout << i->first << "\t" << j->first << "\t" << SeqVec[j->second]
				<< "\t" << QualVec[j->second] << "\t";
		for (int k = 0; k != CycleVec[j->second].size(); ++k)
			fout << CycleVec[j->second][k];
		fout << endl;
	}
}
fout.close();
return 0;
}

int findMaxAllele(size_t * a, size_t len)
{
int max = 0;
int maxIndex = 0;
for (int i = 0; i != len; ++i)
{
	if (a[i] > max)
	{
		max = a[i];
		maxIndex = i;
	}
}
return maxIndex;
}
int countAllele(size_t*a, const string & seq)
{
for (int i = 0; i != seq.size(); ++i)
{
	if (seq[i] == 'A')
	{
		a[0]++;
		continue;
	}
	if (seq[i] == 'C')
	{
		a[1]++;
		continue;
	}
	if (seq[i] == 'G')
	{
		a[2]++;
		continue;
	}
	if (seq[i] == 'T')
	{
		a[3]++;
		continue;
	}
}
return 0;
}
#define REV_PHRED(x)	pow(10.0,x/(-10))
double calLikelihood(const string & seq, const string & qual, const char& maj,
	const char& min)
{
double lik(0);
if (maj == min)
	for (int i = 0; i != seq.size(); ++i)
	{
		if (seq[i] == maj)
		{
			lik += log10(1 - REV_PHRED(qual[i]));
		}
		else
		{
			lik += log10(REV_PHRED(qual[i]) / 3);
		}
	}
else
{
	for (int i = 0; i != seq.size(); ++i)
	{
		if (seq[i] == maj || seq[i] == min)
		{
			lik += log10(1 / 2 - REV_PHRED(qual[i]) / 3);
		}
		else
		{
			lik += log10(REV_PHRED(qual[i]) / 3);
		}
	}
}
return lik;
}
int StatCollector::getGenoLikelihood(const string & outputPath)
{
ofstream fout(outputPath + ".likelihood");
size_t numAllele[4] =
{ 0 };
char majAllele, minAllele;
int maxIndex = 0;
for (string_map::iterator i = PositionTable.begin(); i != PositionTable.end();
		++i) //each chr
{
	for (SeqPairIndex::iterator j = i->second.begin(); j != i->second.end();
			++j) //each site
	{
		countAllele(numAllele, SeqVec[j->second]);
		maxIndex = findMaxAllele(numAllele, 4);
		majAllele = "ACGT"[maxIndex];
		numAllele[maxIndex] = 0;
		maxIndex = findMaxAllele(numAllele, 4);
		minAllele = "ACGT"[maxIndex];
		fout << i->first << "\t" << j->first << "\t"
				<< calLikelihood(SeqVec[j->second], QualVec[j->second],
						majAllele, minAllele);
		fout << endl;
	}
}
fout.close();
return 0;
}
StatCollector::~StatCollector()
{
// TODO Auto-generated destructor stub
for (int i = 0; i != VcfRecVec.size(); ++i)
	delete VcfRecVec[i];
}


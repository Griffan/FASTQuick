/*
 * StatCollector.cpp
 *
 *  Created on: 2014Äê7ÔÂ20ÈÕ
 *      Author: Administrator
 */

#include "StatCollector.h"
#include "../libbwa/bwase.h"
#include "../libbwa/bamlite.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
using namespace std;
extern string Prefix;
extern void notice(const char*,...);
extern void warning(const char*,...);
extern void error(const char*,...);
StatCollector::StatCollector()
{
	// TODO Auto-generated constructor stub
	//cerr << "NOTE:Using default initializer..." << endl;
	PositionTable.clear();
	VcfRecVec.clear();
	duplicateTable.clear();
	VariantProxyTable.clear();
	GC.clear();
	index = 0;
	total_base = 0;
	total_region_size = 0;
	DepthDist = vector<int>(256, 0);
	GCDist = vector<int>(256, 0);
	EmpRepDist = vector<int>(256, 0);
	misEmpRepDist = vector<int>(256, 0);
	EmpCycleDist = vector<int>(256, 0);
	misEmpCycleDist = vector<int>(256, 0);
	InsertSizeDist = vector<int>(2048, 0);
	MaxInsertSizeDist = vector<int>(2048, 0);
}
StatCollector::StatCollector(const string & OutFile)
{
	PositionTable.clear();
	VcfRecVec.clear();
	duplicateTable.clear();
	VariantProxyTable.clear();
	GC.clear();
	index = 0;
	total_base = 0;
	total_region_size = 0;
	DepthDist = vector<int>(256, 0);
	GCDist = vector<int>(256, 0);
	EmpRepDist = vector<int>(256, 0);
	misEmpRepDist = vector<int>(256, 0);
	EmpCycleDist = vector<int>(256, 0);
	misEmpCycleDist = vector<int>(256, 0);
	InsertSizeDist = vector<int>(2048, 0);
	MaxInsertSizeDist = vector<int>(2048, 0);
}

int StatCollector::addSingleAlignment(const bntseq_t *bns, bwa_seq_t *p,
		const gap_opt_t* opt) // TODO: cycle of reverse read need to be fixed
{
	int seqid(0), j(0);

	if (p->type == BWA_TYPE_NO_MATCH)
	{
		//j = 1;
		return false;
	}
	else
	{
		j = pos_end(p) - p->pos; //length of read
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
	for (uint32_t i = 0; i != MD.size(); ++i)
		if (isdigit(MD[i]))
			continue;
		else if (MD[i] == '^')
		{
			i++;
			while (!isdigit(MD[i])) // we don't need to take care of Deletion
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

	//cerr << p->name << "\t" << RefSeq << "\t" << seq << "\t" << MD << endl;

	string PosName = string(bns->anns[seqid].name);

	int pos = (int) (p->pos - bns->anns[seqid].offset + 1);

	size_t colonPos = PosName.find(":");
	size_t atPos = PosName.find("@");
	string Chrom = PosName.substr(0, colonPos);
	unsigned int refCoord = atoi(
			PosName.substr(colonPos + 1, atPos - colonPos + 1).c_str()); // coordinate of variant site
	//size_t slashPos = PosName.find("/");
	//string ref = PosName.substr(atPos + 1, slashPos - atPos - 1); //ref allele
	unsigned int realCoord(0);
	unsigned int RefRealStart(0), RefRealEnd(0);
	if (PosName[PosName.size() - 1] == 'L')
	{
		realCoord = refCoord - opt->flank_long_len + pos - 1; //real coordinate of current reads on reference
		RefRealStart = refCoord - opt->flank_long_len;
		RefRealEnd = refCoord + opt->flank_long_len;
	}
	else
	{
		realCoord = refCoord - opt->flank_len + pos - 1; //real coordinate of current reads on reference
		RefRealStart = refCoord - opt->flank_len;
		RefRealEnd = refCoord + opt->flank_len;
	}
	realCoord += 100; //
	unsigned int tmpCycle(0), tmpCycleVcfTable(0), left_to_right_coord(0),
			tmp_left_to_right_coord(0);
	char sign[2] =
	{ 1, -1 };
	if (p->strand != 0)
	{
		tmpCycle = p->len - 1;
	}

	unsigned int tmp_index = 0;
	//cerr << p->name << "\t" << PosName<<"\t"<<pos<<"\t"<<Chrom<< "\t" << realCoord << "\t" << MD <<"\t";
	//cigar string
	if (p->cigar)
	{
		int n_cigar = p->n_cigar;
		bwa_cigar_t* cigar = p->cigar;
		for (int k = 0; k < n_cigar; ++k)
		{

			int cl = __cigar_len(cigar[k]);
			int cop = "MIDS"[__cigar_op(cigar[k])];
			switch (cop)
			{
			case 'M':
				//cerr<<"M";
				//cerr<< realCoord<<"-"<< realCoord + cl - 1;
				/*************************variant site****************************************/

				if (VcfTable.find(Chrom) != VcfTable.end())
				{
					tmpCycleVcfTable = tmpCycle;
					tmp_left_to_right_coord = left_to_right_coord;
					for (uint32_t i = realCoord; i != realCoord + cl - 1 + 1;
							++i, tmpCycleVcfTable += 1 * sign[p->strand], ++tmp_left_to_right_coord)
					{
						if (VcfTable[Chrom].find(i) != VcfTable[Chrom].end())// actual snp site
						{
							tmp_index = VcfTable[Chrom][i];
							SeqVec[tmp_index] += seq[tmp_left_to_right_coord];
							QualVec[tmp_index] += qual[tmp_left_to_right_coord];
							CycleVec[tmp_index].push_back(tmpCycleVcfTable);
							MaqVec[tmp_index].push_back(p->mapQ + 33);
							StrandVec[tmp_index].push_back(p->strand);
						}
					}
				}

				/*****************************************************************************/
				if (PositionTable.find(Chrom) != PositionTable.end()) //chrom exists
					for (uint32_t i = realCoord; i != realCoord + cl - 1 + 1;
							++i, tmpCycle += 1 * sign[p->strand], ++left_to_right_coord)
					{
						if (i < RefRealStart)
							continue;
						if (i > RefRealEnd)
							break;
						//std::pair<unordered_map<int, bool>::iterator, bool> ret =
						//		PositionTable[Chrom].insert(make_pair(i, 1));
						if (PositionTable[Chrom].find(i)
								!= PositionTable[Chrom].end())//!ret.second) //if insert failed, this coord exists
						{
							//cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
							tmp_index = PositionTable[Chrom][i];
							DepthVec[tmp_index]++;
							if (dbSNPTable[Chrom].find(i)
									== dbSNPTable[Chrom].end())	//not in dbsnp table
							{
								EmpRepDist[qual[left_to_right_coord]]++;
								if(qual[left_to_right_coord]-33>=20)
								{
									Q20DepthVec[tmp_index]++;
									if(qual[left_to_right_coord]-33>=30)
									{
										Q30DepthVec[tmp_index]++;
									}
								}
								if (RefSeq[left_to_right_coord]
										!= seq[left_to_right_coord])
								{
									misEmpRepDist[qual[left_to_right_coord]]++;
									misEmpCycleDist[tmpCycle]++;
								}
								/************EmpCycle**************************************************************/
								EmpCycleDist[tmpCycle]++;
							}
						}
						else //coord not exists
						{
							//cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
							DepthVec.push_back(0);
							Q30DepthVec.push_back(0);
							Q20DepthVec.push_back(0);

							DepthVec[index]++;
							PositionTable[Chrom][i] = index;
							index++;
							if (dbSNPTable[Chrom].find(i)
									== dbSNPTable[Chrom].end()) //not in dbsnp table
							{
								EmpRepDist[qual[left_to_right_coord]]++;
								if(qual[left_to_right_coord]-33>=20)
								{
									Q20DepthVec[tmp_index]++;
									if(qual[left_to_right_coord]-33>=30)
									{
										Q30DepthVec[tmp_index]++;
									}
								}
								if (RefSeq[left_to_right_coord]
										!= seq[left_to_right_coord])
								{
									misEmpRepDist[qual[left_to_right_coord]]++;
									misEmpCycleDist[tmpCycle]++;
								}
								/************EmpCycle**************************************************************/
								EmpCycleDist[tmpCycle]++;
							}
						}
					}
				else // chrom not exists
				{
					for (uint32_t i = realCoord; i != realCoord + cl - 1 + 1;
							++i, tmpCycle += 1 * sign[p->strand], ++left_to_right_coord)
					{
						if (i < RefRealStart)
							continue;
						if (i > RefRealEnd)
							break;
						//cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
						DepthVec.push_back(0);
						Q30DepthVec.push_back(0);
						Q20DepthVec.push_back(0);
						DepthVec[index]++;
						PositionTable[Chrom][i] = index;
						index++;
						if (dbSNPTable[Chrom].find(i)
								== dbSNPTable[Chrom].end()) //not in dbsnp table
						{
							EmpRepDist[qual[left_to_right_coord]]++;
							if(qual[left_to_right_coord]-33>=20)
							{
								Q20DepthVec[tmp_index]++;
								if(qual[left_to_right_coord]-33>=30)
								{
									Q30DepthVec[tmp_index]++;
								}
							}
							if (RefSeq[left_to_right_coord]
									!= seq[left_to_right_coord])
							{
								misEmpRepDist[qual[left_to_right_coord]]++;
								misEmpCycleDist[tmpCycle]++;
							}
							/************EmpCycle**************************************************************/
							EmpCycleDist[tmpCycle]++;
						}
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
				tmpCycle += cl * sign[p->strand];
				left_to_right_coord += cl;
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

				tmpCycle += cl * sign[p->strand];
				left_to_right_coord += cl;
				break;

//	 case BAM_CREF_SKIP:
//	      printf("S");
				//   printf("[%d-%d]", pos, pos + cl - 1); /* Spans positions, No Coverage*/
//	      pos+=cl;
//	      break;

			default:
				warning("Unhandled cigar_op %d:%d\n", cop, cl);
				printf("?");
			}
		}			//end for
	}
	else
	{
		/*************************variant site****************************************/

		if (VcfTable.find(Chrom) != VcfTable.end())
		{
			tmpCycleVcfTable = tmpCycle;
			tmp_left_to_right_coord = left_to_right_coord;
			for (uint32_t i = realCoord; i != realCoord + p->len - 1 + 1;
					++i, tmpCycleVcfTable += 1 * sign[p->strand], ++tmp_left_to_right_coord)
				if (VcfTable[Chrom].find(i) != VcfTable[Chrom].end())// actual snp site
				{
					tmp_index = VcfTable[Chrom][i];
					SeqVec[tmp_index] += seq[tmp_left_to_right_coord];
					QualVec[tmp_index] += qual[tmp_left_to_right_coord];
					CycleVec[tmp_index].push_back(tmpCycleVcfTable);
					MaqVec[tmp_index].push_back(p->mapQ + 33);
					StrandVec[tmp_index].push_back(p->strand);
				}
		}
		/************************************************************************************/
		if (PositionTable.find(Chrom) != PositionTable.end()) //chrom exists
			for (uint32_t i = realCoord; i != realCoord + p->len - 1 + 1;
					++i, tmpCycle += 1 * sign[p->strand], ++left_to_right_coord)
			{
				if (i < RefRealStart)
					continue;
				if (i > RefRealEnd)
					break;
				//std::pair<unordered_map<int, bool>::iterator, bool> ret =
				//		PositionTable[Chrom].insert(make_pair(i, 1));
				if (PositionTable[Chrom].find(i) != PositionTable[Chrom].end())	//!ret.second) //if insert failed, this coord exists
				{
					//cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
					tmp_index = PositionTable[Chrom][i];
					DepthVec[tmp_index]++;
					if (dbSNPTable[Chrom].find(i) == dbSNPTable[Chrom].end())//not in dbsnp table
					{
						EmpRepDist[qual[left_to_right_coord]]++;
						if(qual[left_to_right_coord]-33>=20)
						{
							Q20DepthVec[tmp_index]++;
							if(qual[left_to_right_coord]-33>=30)
							{
								Q30DepthVec[tmp_index]++;
							}
						}
						if (RefSeq[left_to_right_coord]
								!= seq[left_to_right_coord])
						{
							misEmpRepDist[qual[left_to_right_coord]]++;
							misEmpCycleDist[tmpCycle]++;
						}
						/************EmpCycle**************************************************************/
						EmpCycleDist[tmpCycle]++;
					}
				}
				else //coord not exists
				{
					//cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
					DepthVec.push_back(0);
					Q30DepthVec.push_back(0);
					Q20DepthVec.push_back(0);
					DepthVec[index]++;
					PositionTable[Chrom][i] = index;
					index++;
					if (dbSNPTable[Chrom].find(i) == dbSNPTable[Chrom].end()) //not in dbsnp table
					{
						EmpRepDist[qual[left_to_right_coord]]++;
						if(qual[left_to_right_coord]-33>=20)
						{
							Q20DepthVec[tmp_index]++;
							if(qual[left_to_right_coord]-33>=30)
							{
								Q30DepthVec[tmp_index]++;
							}
						}
						if (RefSeq[left_to_right_coord]
								!= seq[left_to_right_coord])
						{
							misEmpRepDist[qual[left_to_right_coord]]++;
							misEmpCycleDist[tmpCycle]++;
						}
						/************EmpCycle**************************************************************/
						EmpCycleDist[tmpCycle]++;
					}
				}
			}
		else // chrom not exists
		{
			for (uint32_t i = realCoord; i != realCoord + p->len - 1 + 1;
					++i, tmpCycle += 1 * sign[p->strand], ++left_to_right_coord)
			{
				if (i < RefRealStart)
					continue;
				if (i > RefRealEnd)
					break;
				//cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
				DepthVec.push_back(0);
				Q30DepthVec.push_back(0);
				Q20DepthVec.push_back(0);
				DepthVec[index]++;
				PositionTable[Chrom][i] = index;
				index++;
				if (dbSNPTable[Chrom].find(i) == dbSNPTable[Chrom].end()) //not in dbsnp table
				{
					EmpRepDist[qual[left_to_right_coord]]++;
					if(qual[left_to_right_coord]-33>=20)
					{
						Q20DepthVec[tmp_index]++;
						if(qual[left_to_right_coord]-33>=30)
						{
							Q30DepthVec[tmp_index]++;
						}
					}
					if (RefSeq[left_to_right_coord] != seq[left_to_right_coord])
					{
						misEmpRepDist[qual[left_to_right_coord]]++;
						misEmpCycleDist[tmpCycle]++;
					}
					/************EmpCycle**************************************************************/
					EmpCycleDist[tmpCycle]++;
				}
			}
		}
	}

	return true;
}

int StatCollector::addSingleAlignment(SamRecord& p, const gap_opt_t* opt) // TODO: cycle of reverse read need to be fixed
{
	//int seqid(0), j(0);

	if (p.getFlag() & SAM_FSU)
	{
		return false;
	}

//	if (/*p->type == BWA_TYPE_NO_MATCH || */(p->count) > 60)
//	{
//		return false;
//	}
//
	string seq(p.getSequence()), qual(p.getQuality()), newSeq, newQual;

//	if (p->strand == 0)
//		for (j = 0; j != p->full_len; ++j)
//		{
//			seq += "ACGTN"[(int) (p)->seq[j]];
//			qual += (char) p->qual[j];
//		}
//	else
//		for (j = 0; j != p->full_len; ++j)
//		{
//			seq += "TGCAN"[(int) (p)->seq[p->full_len - 1 - j]];
//			qual += (char) p->qual[p->full_len - 1 - j];
//		}
	//if (p->strand) seq_reverse(p->len, p->qual, 0);
	//fprintf(stderr,"%s\t%s\t%d\n",p->name,p->md,p->count);

	//ConstructFakeSeqQual(seq,qual,p->n_cigar,p->cigar,newSeq,newQual);
	string RefSeq(seq), MD(p.getStringTag("MD")->c_str());
	int last(0), total_len(0);
	for (uint32_t i = 0; i != MD.size(); ++i)
		if (isdigit(MD[i]))
			continue;
		else if (MD[i] == '^')
		{
			i++;
			while (!isdigit(MD[i])) // we don't need to take care of Deletion
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

	//cerr << p->name << "\t" << RefSeq << "\t" << seq << "\t" << MD << endl;

	string PosName = p.getReferenceName();

	int pos = p.get1BasedPosition();

	size_t colonPos = PosName.find(":");
	size_t atPos = PosName.find("@");
	string Chrom = PosName.substr(0, colonPos);
	unsigned int refCoord = atoi(
			PosName.substr(colonPos + 1, atPos - colonPos + 1).c_str()); // coordinate of variant site
	//size_t slashPos = PosName.find("/");
	//string ref = PosName.substr(atPos + 1, slashPos - atPos - 1); //ref allele
	unsigned int realCoord(0);
	unsigned int RefRealStart(0), RefRealEnd(0);
	if (PosName[PosName.size() - 1] == 'L')
	{
		realCoord = refCoord - opt->flank_long_len + pos - 1; //real coordinate of current reads on reference
		RefRealStart = refCoord - opt->flank_long_len;
		RefRealEnd = refCoord + opt->flank_long_len;
	}
	else
	{
		realCoord = refCoord - opt->flank_len + pos - 1; //real coordinate of current reads on reference
		RefRealStart = refCoord - opt->flank_len;
		RefRealEnd = refCoord + opt->flank_len;
	}
	realCoord += 100; //
	unsigned int tmpCycle(0), tmpCycleVcfTable(0), left_to_right_coord(0),
			tmp_left_to_right_coord(0);
	char sign[2] =
	{ 1, -1 };

	bool strand = (p.getFlag() & SAM_FSR);
	int mapQ = p.getMapQuality();
	if (strand != 0)
	{
		tmpCycle = p.getReadLength() - 1;
	}

	unsigned int tmp_index = 0;

	if (strcmp(p.getCigar(), "*") != 0)
	{

		string cigar = p.getCigar();
		int n_cigar = cigar.size();
		int digitNow(0), digitLast(0);
		//for (int k = 0; k < n_cigar; ++k)
		while (digitNow != n_cigar )
		{
			while (isdigit(cigar[digitNow]))
				digitNow++;
			int cl = atoi(
					cigar.substr(digitLast, digitNow - digitLast).c_str());//digitNow is "SMID"
			digitLast = digitNow + 1;
			char cop = cigar[digitNow];
			digitNow++;
			//std::cerr<<cigar<<":\t"<<cop<<"\t"<<cl<<std::endl;
			switch (cop)
			{
			case 'M':
				/*************************variant site****************************************/

				if (VcfTable.find(Chrom) != VcfTable.end())
				{
					tmpCycleVcfTable = tmpCycle;
					tmp_left_to_right_coord = left_to_right_coord;
					for (uint32_t i = realCoord; i != realCoord + cl - 1 + 1;
							++i, tmpCycleVcfTable += 1 * sign[strand], ++tmp_left_to_right_coord)
					{
						if (VcfTable[Chrom].find(i) != VcfTable[Chrom].end()) // actual snp site
						{
							tmp_index = VcfTable[Chrom][i];
							SeqVec[tmp_index] += seq[tmp_left_to_right_coord];
							QualVec[tmp_index] += qual[tmp_left_to_right_coord];
							CycleVec[tmp_index].push_back(tmpCycleVcfTable);
							MaqVec[tmp_index].push_back(mapQ + 33);
							StrandVec[tmp_index].push_back(strand);
						}
					}
				}
				/*****************************************************************************/
				if (PositionTable.find(Chrom) != PositionTable.end()) //chrom exists
					for (uint32_t i = realCoord; i != realCoord + cl - 1 + 1;
							++i, tmpCycle += 1 * sign[strand], ++left_to_right_coord)
					{
						if (i < RefRealStart)
							continue;
						if (i > RefRealEnd)
							break;
						if (PositionTable[Chrom].find(i)
								!= PositionTable[Chrom].end()) //!ret.second) //if insert failed, this coord exists
						{
							tmp_index = PositionTable[Chrom][i];
							DepthVec[tmp_index]++;
							if (dbSNPTable[Chrom].find(i)
									== dbSNPTable[Chrom].end())	//not in dbsnp table
							{
								EmpRepDist[qual[left_to_right_coord]]++;
								if(qual[left_to_right_coord]-33>=20)
								{
									Q20DepthVec[tmp_index]++;
									if(qual[left_to_right_coord]-33>=30)
									{
										Q30DepthVec[tmp_index]++;
									}
								}
								if (RefSeq[left_to_right_coord]
										!= seq[left_to_right_coord])
								{
									misEmpRepDist[qual[left_to_right_coord]]++;
									misEmpCycleDist[tmpCycle]++;
								}
								/************EmpCycle**************************************************************/
								EmpCycleDist[tmpCycle]++;
							}
						}
						else //coord not exists
						{
							DepthVec.push_back(0);
							Q30DepthVec.push_back(0);
							Q20DepthVec.push_back(0);
							DepthVec[index]++;
							PositionTable[Chrom][i] = index;
							index++;
							if (dbSNPTable[Chrom].find(i)
									== dbSNPTable[Chrom].end()) //not in dbsnp table
							{
								EmpRepDist[qual[left_to_right_coord]]++;
								if(qual[left_to_right_coord]-33>=20)
								{
									Q20DepthVec[tmp_index]++;
									if(qual[left_to_right_coord]-33>=30)
									{
										Q30DepthVec[tmp_index]++;
									}
								}
								if (RefSeq[left_to_right_coord]
										!= seq[left_to_right_coord])
								{
									misEmpRepDist[qual[left_to_right_coord]]++;
									misEmpCycleDist[tmpCycle]++;
								}
								/************EmpCycle**************************************************************/
								EmpCycleDist[tmpCycle]++;
							}
						}
					}
				else // chrom not exists
				{
					for (uint32_t i = realCoord; i != realCoord + cl - 1 + 1;
							++i, tmpCycle += 1 * sign[strand], ++left_to_right_coord)
					{
						if (i < RefRealStart)
							continue;
						if (i > RefRealEnd)
							break;
						DepthVec.push_back(0);
						Q30DepthVec.push_back(0);
						Q20DepthVec.push_back(0);
						DepthVec[index]++;
						PositionTable[Chrom][i] = index;
						index++;
						if (dbSNPTable[Chrom].find(i)
								== dbSNPTable[Chrom].end()) //not in dbsnp table
						{
							EmpRepDist[qual[left_to_right_coord]]++;
							if(qual[left_to_right_coord]-33>=20)
							{
								Q20DepthVec[tmp_index]++;
								if(qual[left_to_right_coord]-33>=30)
								{
									Q30DepthVec[tmp_index]++;
								}
							}
							if (RefSeq[left_to_right_coord]
									!= seq[left_to_right_coord])
							{
								misEmpRepDist[qual[left_to_right_coord]]++;
								misEmpCycleDist[tmpCycle]++;
							}
							/************EmpCycle**************************************************************/
							EmpCycleDist[tmpCycle]++;
						}
					}
				}
				realCoord += cl;
				break;
			case 'S':
				tmpCycle += cl * sign[strand];
				left_to_right_coord += cl;
				break;
			case 'D':
				realCoord += cl;
				break;
			case 'I':
				tmpCycle += cl * sign[strand];
				left_to_right_coord += cl;
				break;
			default:
				warning("Unhandled cigar_op %d:%d\n", cop, cl);
				printf("?");
			}
		}			//end for
	}
	else
	{
		/*************************variant site****************************************/

		if (VcfTable.find(Chrom) != VcfTable.end())
		{
			tmpCycleVcfTable = tmpCycle;
			tmp_left_to_right_coord = left_to_right_coord;
			for (uint32_t i = realCoord;
					i != realCoord + p.getReadLength() - 1 + 1;
					++i, tmpCycleVcfTable += 1 * sign[strand], ++tmp_left_to_right_coord)
				if (VcfTable[Chrom].find(i) != VcfTable[Chrom].end())// actual snp site
				{
					tmp_index = VcfTable[Chrom][i];
					SeqVec[tmp_index] += seq[tmp_left_to_right_coord];
					QualVec[tmp_index] += qual[tmp_left_to_right_coord];
					CycleVec[tmp_index].push_back(tmpCycleVcfTable);
					MaqVec[tmp_index].push_back(mapQ + 33);
					StrandVec[tmp_index].push_back(strand);
				}
		}
		/************************************************************************************/
		if (PositionTable.find(Chrom) != PositionTable.end()) //chrom exists
			for (uint32_t i = realCoord;
					i != realCoord + p.getReadLength() - 1 + 1;
					++i, tmpCycle += 1 * sign[strand], ++left_to_right_coord)
			{
				if (i < RefRealStart)
					continue;
				if (i > RefRealEnd)
					break;
				//std::pair<unordered_map<int, bool>::iterator, bool> ret =
				//		PositionTable[Chrom].insert(make_pair(i, 1));
				if (PositionTable[Chrom].find(i) != PositionTable[Chrom].end())	//!ret.second) //if insert failed, this coord exists
				{
					//cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
					tmp_index = PositionTable[Chrom][i];
					DepthVec[tmp_index]++;
					if (dbSNPTable[Chrom].find(i) == dbSNPTable[Chrom].end())//not in dbsnp table
					{
						EmpRepDist[qual[left_to_right_coord]]++;
						if(qual[left_to_right_coord]-33>=20)
						{
							Q20DepthVec[tmp_index]++;
							if(qual[left_to_right_coord]-33>=30)
							{
								Q30DepthVec[tmp_index]++;
							}
						}
						if (RefSeq[left_to_right_coord]
								!= seq[left_to_right_coord])
						{
							misEmpRepDist[qual[left_to_right_coord]]++;
							misEmpCycleDist[tmpCycle]++;
						}
						/************EmpCycle**************************************************************/
						EmpCycleDist[tmpCycle]++;
					}
				}
				else //coord not exists
				{
					//cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
					DepthVec.push_back(0);
					Q30DepthVec.push_back(0);
					Q20DepthVec.push_back(0);
					DepthVec[index]++;
					PositionTable[Chrom][i] = index;
					index++;
					if (dbSNPTable[Chrom].find(i) == dbSNPTable[Chrom].end()) //not in dbsnp table
					{
						EmpRepDist[qual[left_to_right_coord]]++;
						if(qual[left_to_right_coord]-33>=20)
						{
							Q20DepthVec[tmp_index]++;
							if(qual[left_to_right_coord]-33>=30)
							{
								Q30DepthVec[tmp_index]++;
							}
						}
						if (RefSeq[left_to_right_coord]
								!= seq[left_to_right_coord])
						{
							misEmpRepDist[qual[left_to_right_coord]]++;
							misEmpCycleDist[tmpCycle]++;
						}
						/************EmpCycle**************************************************************/
						EmpCycleDist[tmpCycle]++;
					}
				}
			}
		else // chrom not exists
		{
			for (uint32_t i = realCoord;
					i != realCoord + p.getReadLength() - 1 + 1;
					++i, tmpCycle += 1 * sign[strand], ++left_to_right_coord)
			{
				if (i < RefRealStart)
					continue;
				if (i > RefRealEnd)
					break;
				//cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
				DepthVec.push_back(0);
				Q30DepthVec.push_back(0);
				Q20DepthVec.push_back(0);
				DepthVec[index]++;
				PositionTable[Chrom][i] = index;
				index++;
				if (dbSNPTable[Chrom].find(i) == dbSNPTable[Chrom].end()) //not in dbsnp table
				{
					EmpRepDist[qual[left_to_right_coord]]++;
					if(qual[left_to_right_coord]-33>=20)
					{
						Q20DepthVec[tmp_index]++;
						if(qual[left_to_right_coord]-33>=30)
						{
							Q30DepthVec[tmp_index]++;
						}
					}
					if (RefSeq[left_to_right_coord] != seq[left_to_right_coord])
					{
						misEmpRepDist[qual[left_to_right_coord]]++;
						misEmpCycleDist[tmpCycle]++;
					}
					/************EmpCycle**************************************************************/
					EmpCycleDist[tmpCycle]++;
				}
			}
		}
	}

	return true;
}

int StatCollector::IsDuplicated(const bntseq_t *bns, const bwa_seq_t *p,
		const bwa_seq_t *q, const gap_opt_t* opt, int type, ofstream & fout)
{
	int MaxInsert(-1), MaxInsert2(-1);
	;
	int j(0), seqid_p(-1), seqid_q(-1);
	//cerr<<"Duplicate function entered"<<endl;
	if (type == 1) //q is aligned only
	{
		j = pos_end(q) - q->pos; //length of read
		bns_coor_pac2real(bns, q->pos, j, &seqid_q);
		MaxInsert = q->pos + j - bns->anns[seqid_q].offset;
		//cerr<<"Duplicate function exit single q"<<endl;
		fout << q->name << "\t" << MaxInsert << "\t" << 0 << "\t" << "*" << "\t"
				<< "*" << "\t" << bns->anns[seqid_q].name << "\t"
				<< q->pos - bns->anns[seqid_q].offset + 1 << "\tRevOnly"
				<< endl;
		return 0;
	}
	else if (type == 3) //p is aligned only
	{
		j = pos_end(p) - p->pos; //length of read
		bns_coor_pac2real(bns, p->pos, j, &seqid_p);
		MaxInsert = bns->anns[seqid_p].offset + bns->anns[seqid_p].len - p->pos;
		//cerr<<"Duplicate function exit single p"<<endl;
		fout << p->name << "\t" << MaxInsert << "\t" << 0 << "\t"
				<< bns->anns[seqid_p].name << "\t"
				<< p->pos - bns->anns[seqid_p].offset + 1 << "\t" << "*" << "\t"
				<< "*" << "\tFwdOnly" << endl;
		return 0;
	}
	else if (type == 2) // both aligned
	{
		if (p->pos < q->pos)
		{
			j = pos_end(q) - q->pos; //length of read
			bns_coor_pac2real(bns, q->pos, j, &seqid_q);
			MaxInsert = q->pos + j - bns->anns[seqid_q].offset;
			j = pos_end(p) - p->pos; //length of read
			bns_coor_pac2real(bns, p->pos, j, &seqid_p);
			MaxInsert2 = bns->anns[seqid_p].offset + bns->anns[seqid_p].len
					- p->pos;
		}
		else
		{
			j = pos_end(q) - q->pos; //length of read
			bns_coor_pac2real(bns, q->pos, j, &seqid_q);
			MaxInsert = bns->anns[seqid_q].offset + bns->anns[seqid_q].len
					- q->pos;
			j = pos_end(p) - p->pos; //length of read
			bns_coor_pac2real(bns, p->pos, j, &seqid_p);
			MaxInsert2 = p->pos + j - bns->anns[seqid_p].offset;
		}
	}
	else
	{
		warning("Alignment status fatal error!\n");
		exit(1);
	}

	if (MaxInsert2 > MaxInsert)
		MaxInsert = MaxInsert2;
	if (MaxInsert > 2047)
		MaxInsert = 2047;
	MaxInsertSizeDist[MaxInsert]++;
	if (seqid_p != seqid_q && seqid_p != -1 && seqid_q != -1)
	{
		//int ActualInsert(-1);
		//	ActualInsert=0;
		InsertSizeDist[0]++;
		//cerr<<"Duplicate function exit from diff chrom"<<endl;
		fout << p->name << "\t" << MaxInsert << "\t" << 0 << "\t"
				<< bns->anns[seqid_p].name << "\t"
				<< p->pos - bns->anns[seqid_p].offset + 1 << "\t"
				<< bns->anns[seqid_q].name << "\t"
				<< q->pos - bns->anns[seqid_q].offset + 1 << "\tDiffChrom"
				<< endl;
		return 0;
	}
	if (p->mapQ >= 20 && q->mapQ >= 20)
	{
		int ActualInsert(-1), start(0), end(0);

		if (p->pos < q->pos)
		{
			start = p->pos;
			end = q->pos + q->len;
			ActualInsert = end - start;
		}
		else
		{
			start = q->pos;
			end = p->pos + p->len;
			ActualInsert = end - start;
		}
		if (ActualInsert < 100000)
		{
			InsertSizeDist[ActualInsert]++;
			fout << p->name << "\t" << MaxInsert << "\t" << ActualInsert << "\t"
					<< bns->anns[seqid_p].name << "\t"
					<< p->pos - bns->anns[seqid_p].offset + 1 << "\t"
					<< bns->anns[seqid_q].name << "\t"
					<< q->pos - bns->anns[seqid_q].offset + 1 << "\tPropPair"
					<< endl;
			char start_end[256];
			sprintf(start_end, "%d:%d", start, end);
			pair<unordered_map<string, bool>::iterator, bool> iter =
					duplicateTable.insert(make_pair(string(start_end), true));
			if (!iter.second) //insert failed, duplicated
			{
				//	cerr<<"Duplicate function exit from duplicate"<<endl;
				return 1;
			}
		}
		else
		{
			//cerr<<"exit from Inf"<<endl;
			fout << p->name << "\t" << MaxInsert << "\t" << "-1" << "\t"
					<< bns->anns[seqid_p].name << "\t"
					<< p->pos - bns->anns[seqid_p].offset + 1 << "\t"
					<< bns->anns[seqid_q].name << "\t"
					<< q->pos - bns->anns[seqid_q].offset + 1 << "\tAbnormal"
					<< endl;
		}
	}
	else
	{
		//cerr<<"exit from LowQualf"<<endl;
		fout << p->name << "\t" << MaxInsert << "\t" << 0 << "\t"
				<< bns->anns[seqid_p].name << "\t"
				<< p->pos - bns->anns[seqid_p].offset + 1 << "\t"
				<< bns->anns[seqid_q].name << "\t"
				<< q->pos - bns->anns[seqid_q].offset + 1 << "\tLowQual"
				<< endl;
		return 2; //low quality
	}
	return 0;
}
//overload of function IsDuplicated for direct bam reading
int StatCollector::IsDuplicated(SamFileHeader& SFH,SamRecord& p, SamRecord& q,
		const gap_opt_t* opt, int type, ofstream & fout)
{
	int MaxInsert(-1), MaxInsert2(-1);

	if (type == 1) //q is aligned only
	{
		if(q.getFlag()&SAM_FR1)//if it's read one
			MaxInsert = atoi(SFH.getSQTagValue("LN",q.getReferenceName()))- q.get1BasedPosition() + 1;
		else
			MaxInsert = q.get1BasedPosition() + q.getReadLength();
		fout << q.getReadName() << "\t" << MaxInsert << "\t" << 0 << "\t" << "*"
				<< "\t" << "*" << "\t" << q.getReferenceName() << "\t"
				<< q.get1BasedPosition() << "\tRevOnly" << endl;
		return 0;
	}
	else if (type == 3) //p is aligned only
	{
		if(p.getFlag()&SAM_FR1)//if it's read one
		MaxInsert =  atoi(SFH.getSQTagValue("LN",p.getReferenceName()))
				- p.get1BasedPosition() + 1;
		else
			MaxInsert = p.get1BasedPosition() + p.getReadLength();
		fout << p.getReadName() << "\t" << MaxInsert << "\t" << 0 << "\t"
				<< p.getReferenceName() << "\t" << p.get1BasedPosition() << "\t"
				<< "*" << "\t" << "*" << "\tFwdOnly" << endl;
		return 0;
	}
	else if (type == 2) // both aligned
	{
		if (p.get1BasedPosition() < q.get1BasedPosition())
		{
			//j = pos_end(q) - q->pos; //length of read
			//bns_coor_pac2real(bns, q->pos, j, &seqid_q);
			MaxInsert = q.get1BasedPosition() + q.getReadLength();
			//j = pos_end(p) - p->pos; //length of read
			//bns_coor_pac2real(bns, p->pos, j, &seqid_p);
			MaxInsert2 = atoi(SFH.getSQTagValue("LN",p.getReferenceName()))
					- p.get1BasedPosition() + 1;
			;
		}
		else
		{
			//j = pos_end(q) - q->pos; //length of read
			//bns_coor_pac2real(bns, q->pos, j, &seqid_q);
			MaxInsert = atoi(SFH.getSQTagValue("LN",q.getReferenceName()))
					- q.get1BasedPosition() + 1;
			//j = pos_end(p) - p->pos; //length of read
			//bns_coor_pac2real(bns, p->pos, j, &seqid_p);
			MaxInsert2 = p.get1BasedPosition() + p.getReadLength();
		}
	}
	else
	{
		warning("Alignment status fatal error!\n");
		exit(1);
	}

	if (MaxInsert2 > MaxInsert)
		MaxInsert = MaxInsert2;
	if (MaxInsert > 2047)
		MaxInsert = 2047;
	MaxInsertSizeDist[MaxInsert]++;
	if (strcmp(p.getReferenceName() , q.getReferenceName())!=0)
	{

		InsertSizeDist[0]++;
		//cerr<<"Duplicate function exit from diff chrom"<<endl;
		fout << p.getReadName() << "\t" << MaxInsert << "\t" << 0 << "\t"
				<< p.getReferenceName() << "\t" << p.get1BasedPosition() << "\t"
				<< q.getReferenceName() << "\t" << q.get1BasedPosition()
				<< "\tDiffChrom" << endl;
		return 0;
	}
	if (p.getMapQuality() >= 20 && q.getMapQuality() >= 20)
	{
		int ActualInsert(-1), start(0), end(0);

		if (p.get1BasedPosition() < q.get1BasedPosition())
		{
			start = p.get1BasedPosition();
			end = q.get1BasedPosition() + q.getReadLength();
			ActualInsert = end - start;
		}
		else
		{
			start = q.get1BasedPosition();
			end = p.get1BasedPosition() + p.getReadLength();
			ActualInsert = end - start;
		}
		if (ActualInsert < 100000)
		{
			InsertSizeDist[ActualInsert]++;
			fout << p.getReadName() << "\t" << MaxInsert << "\t" << ActualInsert
					<< "\t" << p.getReferenceName() << "\t"
					<< p.get1BasedPosition() << "\t" << q.getReferenceName()
					<< "\t" << q.get1BasedPosition() << "\tPropPair" << endl;
			char start_end[256];
			sprintf(start_end, "%d:%d", start, end);
			pair<unordered_map<string, bool>::iterator, bool> iter =
					duplicateTable.insert(make_pair(string(start_end), true));
			if (!iter.second) //insert failed, duplicated
			{
				//	cerr<<"Duplicate function exit from duplicate"<<endl;
				return 1;
			}
		}
		else
		{
			//cerr<<"exit from Inf"<<endl;
			fout << p.getReadName() << "\t" << MaxInsert << "\t" << 0 << "\t"
					<< p.getReferenceName() << "\t" << p.get1BasedPosition()
					<< "\t" << q.getReferenceName() << "\t"
					<< q.get1BasedPosition() << "\tAbnormal" << endl;
		}
	}
	else
	{
		//cerr<<"exit from LowQualf"<<endl;
		fout << p.getReadName() << "\t" << MaxInsert << "\t" << 0 << "\t"
				<< p.getReferenceName() << "\t" << p.get1BasedPosition() << "\t"
				<< q.getReferenceName() << "\t" << q.get1BasedPosition()
				<< "\tLowQual" << endl;
		return 2; //low quality
	}
	return 0;
}
//return value: 0 failed, 1 add pair success, 2 add single success by using add pair interface
int StatCollector::addAlignment(const bntseq_t *bns, bwa_seq_t *p,
		bwa_seq_t *q, const gap_opt_t* opt, ofstream & fout, int &total_add) //TODO: resolve duplicate output option
{
	int seqid(0), seqid2(0), j(0), j2(0);
	if (p==0||p->type == BWA_TYPE_NO_MATCH)
	{
		if(p==0&&q!=0&&addSingleAlignment(bns,q, opt))//adding single via pair interface
		{
			bns_coor_pac2real(bns, q->pos, j2, &seqid2);
			string qname(bns->anns[seqid2].name);
			if (string(qname).find("chrY") != string::npos
								|| string(qname).find("Y") != string::npos
								|| string(qname).find("chrX") != string::npos
								|| string(qname).find("X") != string::npos)
			{
				 if (!isPartialAlign(q))
				 {
						contigStatusTable[qname].addNumOverlappedReads();
						contigStatusTable[qname].addNumFullyIncludedReads();
				 }
				 else
					 contigStatusTable[qname].addNumOverlappedReads();
			}
			total_add++;
			return 2;
		}
		j = 1;
		return 0;
	}
	else
	{
		j = pos_end(p) - p->pos; //length of read
	}
	bns_coor_pac2real(bns, p->pos, j, &seqid);
	string pname(bns->anns[seqid].name);
	if (q==0||q->type == BWA_TYPE_NO_MATCH)
	{
		if(q==0&&p!=0&&addSingleAlignment(bns,p, opt))//adding single via pair interface
		{
			if (string(pname).find("chrY") != string::npos
								|| string(pname).find("Y") != string::npos
								|| string(pname).find("chrX") != string::npos
								|| string(pname).find("X") != string::npos)
			{
				 if (!isPartialAlign(p))
				 {
						contigStatusTable[pname].addNumOverlappedReads();
						contigStatusTable[pname].addNumFullyIncludedReads();
				 }
				 else
					 contigStatusTable[pname].addNumOverlappedReads();
			}
			total_add++;
			return 2;
		}
		j2 = 1;
		return 0;
	}
	else
	{
		j2 = pos_end(q) - q->pos; //length of read
	}
	bns_coor_pac2real(bns, q->pos, j2, &seqid2);
	string qname(bns->anns[seqid2].name);
//until now both reads are aligned
if (isPartialAlign(p)) //p is partially aligned
	{
			if (string(qname).find("chrY") != string::npos
					|| string(qname).find("Y") != string::npos
					|| string(qname).find("chrX") != string::npos
					|| string(qname).find("X") != string::npos)
			{
				if (isPartialAlign(q))
				{
					contigStatusTable[qname].addNumOverlappedReads();
				}
				else //q is perfectly aligned
				{
					contigStatusTable[qname].addNumOverlappedReads();
					contigStatusTable[qname].addNumFullyIncludedReads();
				}
				if (pname == qname)
				{
					contigStatusTable[qname].addNumPairOverlappedReads();
				}
			}
			if (string(pname).find("chrY") != string::npos
					|| string(pname).find("Y") != string::npos
					|| string(pname).find("chrX") != string::npos
					|| string(pname).find("X") != string::npos)
			{
				contigStatusTable[pname].addNumOverlappedReads();
			}

			if (IsDuplicated(bns, p, q, opt, 2, fout) != 1 || opt->cal_dup)
			{
				if (addSingleAlignment(bns, p, opt) && addSingleAlignment(bns, q, opt))
					total_add += 2;
			}
			return 0;
	}
	else //p is perfectly aligned
	{
		 // p and q are both aligned
			if (string(qname).find("chrY") != string::npos
					|| string(qname).find("Y") != string::npos
					|| string(qname).find("chrX") != string::npos
					|| string(qname).find("X") != string::npos)
			{
				if (isPartialAlign(q)) //q is partially aligned
				{
					contigStatusTable[qname].addNumOverlappedReads();
					if (pname == qname)
					{
						contigStatusTable[qname].addNumPairOverlappedReads();
					}
				}
				else //q is perfectly aligned
				{
					contigStatusTable[qname].addNumOverlappedReads();
					contigStatusTable[qname].addNumFullyIncludedReads();
					if (pname == qname)
					{
						contigStatusTable[qname].addNumPairOverlappedReads();
						contigStatusTable[qname].addNumFullyIncludedPairedReads();
					}
				}

			}
			if (string(pname).find("chrY") != string::npos
					|| string(pname).find("Y") != string::npos
					|| string(pname).find("chrX") != string::npos
					|| string(pname).find("X") != string::npos)
			{
				contigStatusTable[pname].addNumOverlappedReads();
				contigStatusTable[pname].addNumFullyIncludedReads();
			}

			if (IsDuplicated(bns, p, q, opt, 2, fout) != 1 || opt->cal_dup)
			{
				if (addSingleAlignment(bns, p, opt) && addSingleAlignment(bns, q, opt))
					total_add += 2;
			}
			else
				return 0;
	}
	//cerr << "currently added reads " << total_add << endl;
	return 1;
}
//overload of addAlignment function for direct bam reading
//return value: 0 failed, 1 add pair success, 2 add single success by using add pair interface
int StatCollector::addAlignment(SamFileHeader & SFH, SamRecord * p, SamRecord* q,
		const gap_opt_t* opt, ofstream & fout, int &total_add)
{

//	if (p->type == BWA_TYPE_NO_MATCH)
	if (p==0||(p->getFlag() & SAM_FSU))
	{
		if (q==0||(q->getFlag() & SAM_FSU)) //both end are not mapped
			return 0;
		else if (isPartialAlign(*q))  //only one pair partially mapped
		{
			string qname(q->getReferenceName());
			if (string(qname).find("chrY") != string::npos
					|| string(qname).find("Y") != string::npos
					|| string(qname).find("chrX") != string::npos
					|| string(qname).find("X") != string::npos)
			{
				contigStatusTable[qname].addNumOverlappedReads();
			}
			if(p==0&&addSingleAlignment(*q, opt))//adding single via pair interface
			{
				total_add++;
				return 2;
			}
			return 0;
		}
		else // q is perfectly aligned, p is not
		{
			string qname(q->getReferenceName());
			if (string(qname).find("chrY") != string::npos
					|| string(qname).find("Y") != string::npos
					|| string(qname).find("chrX") != string::npos
					|| string(qname).find("X") != string::npos)
			{
				contigStatusTable[qname].addNumOverlappedReads();
				contigStatusTable[qname].addNumFullyIncludedReads();
			}

			if(p==0&&addSingleAlignment(*q, opt))//adding single via pair interface
			{
				total_add++;
				return 2;
			}
			else
				return 0;
		}
	}
	else if (isPartialAlign(*p)) //p is partially aligned
	{
		if (q==0||(q->getFlag() & SAM_FSU)) //q is not aligned
		{
			string pname(p->getReferenceName());
			if (string(pname).find("chrY") != string::npos
					|| string(pname).find("Y") != string::npos
					|| string(pname).find("chrX") != string::npos
					|| string(pname).find("X") != string::npos)
			{
				contigStatusTable[pname].addNumOverlappedReads();
			}
			if(q==0&&addSingleAlignment(*p, opt))//adding single via pair interface
			{
				total_add++;
				return 2;
			}
			return 0;
		}
		else // p is partially aligned q is aligned too
		{
			string pname(p->getReferenceName()),qname(q->getReferenceName());
			if (string(qname).find("chrY") != string::npos
					|| string(qname).find("Y") != string::npos
					|| string(qname).find("chrX") != string::npos
					|| string(qname).find("X") != string::npos)
			{
				if (isPartialAlign(*q))
				{
					contigStatusTable[qname].addNumOverlappedReads();
				}
				else //q is perfectly aligned
				{
					contigStatusTable[qname].addNumOverlappedReads();
					contigStatusTable[qname].addNumFullyIncludedReads();
				}
				if (pname==qname)
				{
					contigStatusTable[qname].addNumPairOverlappedReads();
				}
			}
			if (string(pname).find("chrY") != string::npos
					|| string(pname).find("Y") != string::npos
					|| string(pname).find("chrX") != string::npos
					|| string(pname).find("X") != string::npos)
			{
				contigStatusTable[pname].addNumOverlappedReads();
			}

			if (IsDuplicated(SFH,*p, *q, opt, 2, fout) != 1 || opt->cal_dup)
			{
				if (addSingleAlignment(*p, opt) && addSingleAlignment(*q, opt))
					total_add += 2;
			}
			return 0;
		}
	}
	else //p is perfectly aligned
	{
		if (q==0||(q->getFlag() & SAM_FSU)) // p is perfectly aligned, q is not
		{
			string pname(p->getReferenceName());
			if (string(pname).find("chrY") != string::npos
					|| string(pname).find("Y") != string::npos
					|| string(pname).find("chrX") != string::npos
					|| string(pname).find("X") != string::npos)
			{
				contigStatusTable[pname].addNumOverlappedReads();
				contigStatusTable[pname].addNumFullyIncludedReads();
			}
			if(q==0&&addSingleAlignment(*p, opt))//adding single via pair interface
			{
				total_add++;
				return 2;
			}
			else
				return 0;
		}
		else // p and q are both aligned
		{
			string pname(p->getReferenceName()),qname(q->getReferenceName());
			if (string(qname).find("chrY") != string::npos
					|| string(qname).find("Y") != string::npos
					|| string(qname).find("chrX") != string::npos
					|| string(qname).find("X") != string::npos)
			{
				if (isPartialAlign(*q)) //q is partially aligned
				{
					contigStatusTable[qname].addNumOverlappedReads();
					if (pname == qname)
					{
						contigStatusTable[qname].addNumPairOverlappedReads();
					}
				}
				else //q is perfectly aligned
				{
					contigStatusTable[qname].addNumOverlappedReads();
					contigStatusTable[qname].addNumFullyIncludedReads();
					if (pname==qname)
					{
						contigStatusTable[qname].addNumPairOverlappedReads();
						contigStatusTable[qname].addNumFullyIncludedPairedReads();
					}
				}

			}
			if (string(pname).find("chrY") != string::npos
					|| string(pname).find("Y") != string::npos
					|| string(pname).find("chrX") != string::npos
					|| string(pname).find("X") != string::npos)
			{
				contigStatusTable[pname].addNumOverlappedReads();
				contigStatusTable[pname].addNumFullyIncludedReads();
			}

			if (IsDuplicated(SFH, *p, *q, opt, 2, fout) != 1 || opt->cal_dup)
			{
				if (addSingleAlignment(*p, opt) && addSingleAlignment(*q, opt))
					total_add += 2;
			}
			else
				return 0;
		}
	}
	//cerr << "currently added reads " << total_add << endl;
	return 1;
}

int StatCollector::ReadAlignmentFromBam(const gap_opt_t* opt,const char * BamFile,std::ofstream & fout, int & total_add)
{
    SamFileHeader SFH;
    SamFile SFIO;
	if(!SFIO.OpenForRead(BamFile, &SFH))
	{
		cerr << SFIO.GetStatusMessage() << endl;
		warning("Reading Bam Header Failed!\n");
		exit(1);
	}
	FileStatCollector FSC(BamFile);
	unordered_map<string, SamRecord*> pairBuffer;
	while (1)
	{
		SamRecord*  SR= new SamRecord;
		//SFIO.ReadRecord( SFH, *SR);
		if (!SFIO.ReadRecord( SFH, *SR))
		{
			cerr << SFIO.GetStatusMessage() << endl;
			notice("End Bam File Reading...\n");
			delete SR;
			break;
		}
		FSC.NumRead++;
		FSC.NumBase+=SR->getReadLength();
		string readName;
		if(SR->getReadName()[SR->getReadNameLength()-2]=='\\')
		readName=string(SR->getReadName()).substr(0,SR->getReadNameLength() - 3);
		else
			readName=string(SR->getReadName());
		if (!(SR->getFlag() &SAM_FPP)) // read is not mapped in pair
		{
			addAlignment(SFH, SR, 0,opt, fout, total_add);
			delete SR;
			continue;
		}
		else if (pairBuffer.find(readName) == pairBuffer.end()) // mate not existed
		{
			pairBuffer.insert(make_pair(readName, SR));
		}
		else
		{
			SamRecord* SRtmp = pairBuffer[readName];
			addAlignment(SFH,SR, SRtmp, opt, fout, total_add);
			delete SR;
			delete SRtmp;
			pairBuffer.erase(readName);
		}
	}//end of while
	addFSC(FSC);
	for (unordered_map<string, SamRecord*>::iterator iter = pairBuffer.begin();
			iter != pairBuffer.end(); ++iter)
	{
		addAlignment(SFH,iter->second,0,opt, fout, total_add);
		delete iter->second;
	}
	return 0;
}

int StatCollector::restoreVcfSites(const string & VcfPath, const gap_opt_t* opt)
{

	VcfHeader header;
	VcfFileReader reader;
	string SelectedSite = Prefix + ".SelectedSite.vcf.gz";
	if (!reader.open(SelectedSite.c_str(), header))
	{
		warning("File open failed: %s\n",SelectedSite.c_str());
	}
	string GCpath = Prefix + ".gc";
	ifstream FGC(GCpath, ios_base::binary);
	//int num_so_far(0);
	while (!reader.isEOF())
	{
		VcfRecord* VcfLine = new VcfRecord;

		reader.readRecord(*VcfLine);
		VcfRecVec.push_back(VcfLine);

		string chr(VcfLine->getChromStr());
		int pos = VcfLine->get1BasedPosition();
		VcfTable[chr][pos] = VcfRecVec.size() - 1;
		_GCstruct GCstruct;
		GCstruct.read(FGC);
		int tmp_pos = pos - (GCstruct.len - 1) / 2;
		for (uint32_t i = 0; i != GCstruct.len; ++i)
		{
			GC[chr][tmp_pos + i] = GCstruct.GC[i];
		}
		//cerr<<chr<<"\t"<<pos<<"\t"<<num_so_far<<endl;
		//num_so_far++;

		SeqVec.push_back(string(""));
		QualVec.push_back(string(""));
		CycleVec.push_back(vector<unsigned int>(0));
		MaqVec.push_back(vector<unsigned char>(0));
		StrandVec.push_back(vector<bool>(0));
		//VcfTable[Chrom][i]=vcf_index;
		//vcf_index++;

	}
	string BedFile = Prefix + ".subset.vcf";
	if (!reader.open(BedFile.c_str(), header))
	{
		notice("Open %s failed!\n",BedFile.c_str());
		exit(1);
	}
	while (!reader.isEOF())
	{
		VcfRecord VcfLine;
		;
		reader.readRecord(VcfLine);
		string chr(VcfLine.getChromStr());
		int pos = VcfLine.get1BasedPosition();
		dbSNPTable[chr][pos] = 1;
	}

	//string line, Chrom, PosStr, GCStr;
	//_GCstruct * GCstruct = new _GCstruct [opt->num_variant_long*(4*opt->flank_len+1)+opt->num_variant_short*(2*opt->flank_len+1)];
	//FGC.read((char*)GCstruct,(opt->num_variant_long*(4*opt->flank_len+1)+opt->num_variant_short*(2*opt->flank_len+1))*sizeof(_GCstruct));
//	for(uint32_t i=0;i!=VcfRecVec.size();++i)
//	{
	//	GC[string(GCstruct[i].chrom)][GCstruct[i].pos] = GCstruct[i].GC;
	//cerr<<GCstruct[i].chrom<<"\t"<<GCstruct[i].pos<<"\t"<<(int)GCstruct[i].GC<<endl;
	//}
	//int total=0;
	//FGC.read((char*) &total, sizeof(int));
	//cerr<<"the total bytes of gc is :"<<total<<endl;
	/*
	 while (getline(FGC, line))
	 {

	 stringstream ss(line);
	 ss >> Chrom >> PosStr >> GCStr;
	 //cerr<<Chrom<<"\t"<<PosStr<<"\t"<<GCStr<<endl;
	 GC[Chrom][atoi(PosStr.c_str())] = atoi(GCStr.c_str());
	 }*/
	//delete GCstruct;
	FGC.close();
	return 0;
}
int StatCollector::getDepthDist(const string & outputPath, const gap_opt_t* opt)
{

	ofstream fout(outputPath + ".DepthDist");
	int sum = std::accumulate(DepthDist.begin(), DepthDist.end(), 0);
	int all = (opt->flank_len * 2 + 1) * opt->num_variant_short
			+ (opt->flank_long_len * 2 + 1) * opt->num_variant_long;
	fout << 0 << "\t" << all - sum << endl;
	for (uint32_t i = 1; i != DepthDist.size(); ++i)
	{
		fout << i << "\t" << DepthDist[i] << endl;
	}
	fout.close();
	return 0;
}

int StatCollector::getGCDist(const string & outputPath,
		const vector<int> & PosNum)
{
	ofstream fout(outputPath + ".GCDist");
	for (uint32_t i = 0; i != GCDist.size(); ++i)
	{
		fout << i << "\t" << GCDist[i] << "\t" << PosNum[i] << "\t";
		if (PosNum[i] == 0)
		{
			fout << 0;
		}
		else
		{
			fout << double(GCDist[i]) / PosNum[i];
		}
		fout << endl;
	}
	fout.close();
	return 0;
}

#define PHRED(x)	(-10)*log10(x)
int StatCollector::getEmpRepDist(const string & outputPath)
{
	ofstream fout(outputPath + ".EmpRepDist");
	for (uint32_t i = 33; i != EmpRepDist.size(); ++i)
	{
		fout << i << "\t" << (misEmpRepDist[i]) << "\t" << (EmpRepDist[i])
				<< "\t"
				<< PHRED((double )(misEmpRepDist[i] + 1) / (EmpRepDist[i] + 2))
				<< endl;
	}
	fout.close();
	return 0;
}
int StatCollector::getEmpCycleDist(const string & outputPath)
{
	ofstream fout(outputPath + ".EmpCycleDist");
	for (uint32_t i = 0; i != EmpCycleDist.size(); ++i)
	{
		fout << i << "\t" << misEmpCycleDist[i] << "\t" << EmpCycleDist[i]
				<< "\t"
				<< PHRED(
						(double )(misEmpCycleDist[i] + 1)
								/ (EmpCycleDist[i] + 2)) << endl;
	}
	fout.close();
	return 0;
}
int StatCollector::getInsertSizeDist(const string & outputPath)
{
	ofstream fout(outputPath + ".InsertSizeDist");
	for (uint32_t i = 0; i != InsertSizeDist.size(); ++i)
	{
		fout << i << "\t" << InsertSizeDist[i] << endl;
	}
	fout.close();
	return 0;
}
int StatCollector::getSexChromInfo(const string & outputPath)
{
	ofstream fout(outputPath + ".SexChromInfo");
	for (unordered_map<string,ContigStatus>::iterator iter = contigStatusTable.begin(); iter != contigStatusTable.end(); ++iter)
	{
		fout << iter->first<< "\t" <<  iter->second.getNumOverlappedReads()<<"\t"<<iter->second.getNumFullyIncludedReads()<<"\t"<<iter->second.getNumPairOverlappedReads()<<"\t"<<iter->second.getNumFullyIncludedPairedReads()<< endl;
	}
	fout.close();
	return 0;
}
int StatCollector::processCore(const string & statPrefix, const gap_opt_t* opt)
{
	vector<int> PosNum(256, 0);
	//int total(0);
	for (unsort_map::iterator i = PositionTable.begin();
			i != PositionTable.end(); ++i) //each chr
	{
		for (std::unordered_map<int, unsigned int>::iterator j =
				i->second.begin(); j != i->second.end(); ++j) //each site
		{
			/************DepthDist**************************************************************/
			{
				if (DepthVec[j->second] > 255)
					DepthDist[255]++;
				else
					DepthDist[DepthVec[j->second]]++;
			}
			GCDist[GC[i->first][j->first]] += DepthVec[j->second];
			PosNum[GC[i->first][j->first]]++;
		}
	}
	getDepthDist(statPrefix, opt);
	getGCDist(statPrefix, PosNum);
	getEmpRepDist(statPrefix);
	getEmpCycleDist(statPrefix);
	getInsertSizeDist(statPrefix);
	getSexChromInfo(statPrefix);
	outputPileup(statPrefix);
	return 0;
}
int StatCollector::outputPileup(const string & outputPath)
{
	ofstream fout(outputPath + ".Pileup");
	for (sort_map::iterator i = VcfTable.begin(); i != VcfTable.end(); ++i) //each chr
	{
		for (std::map<int, unsigned int>::iterator j = i->second.begin();
				j != i->second.end(); ++j) //each site
		{
			if (SeqVec[j->second].size() == 0)
				continue;
			fout << i->first << "\t" << j->first << "\t.\t"
					<< StrandVec[j->second].size() << "\t";
			for (uint32_t k = 0; k != StrandVec[j->second].size(); ++k)
			{
				if (StrandVec[j->second][k])
					fout << (char) toupper(SeqVec[j->second][k]);
				else
					fout << (char) tolower(SeqVec[j->second][k]);
			}
			fout << "\t" << QualVec[j->second] << "\t";
			for (uint32_t k = 0; k != MaqVec[j->second].size(); ++k)
			{
				fout << MaqVec[j->second][k];
			}
			fout << "\t";
			for (uint32_t k = 0; k != CycleVec[j->second].size(); ++k)
			{
				fout << CycleVec[j->second][k];
				if (k != CycleVec[j->second].size() - 1)
					fout << ",";
			}
			fout << endl;
		}
	}
	fout.close();
	char cmdline[2048];
	sprintf(cmdline,"bgzip -f %s.Pileup",outputPath.c_str());
	if(system(cmdline)!=0)
	{
		warning("Call command line:\n%s\nfailed!\n",cmdline);
	}
	sprintf(cmdline,"tabix  -s 1 -b 2 -e 2 %s.Pileup.gz",outputPath.c_str());
	if(system(cmdline)!=0)
	{
		warning("Call command line:\n%s\nfailed!\n",cmdline);
	}
	return 0;
}

int findMaxAllele(size_t * a, size_t len)
{
	uint32_t max = 0;
	uint32_t maxIndex = 0;
	for (uint32_t i = 0; i != len; ++i)
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
	for (uint32_t i = 0; i != seq.size(); ++i)
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
		for (uint32_t i = 0; i != seq.size(); ++i)
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
		for (uint32_t i = 0; i != seq.size(); ++i)
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
	for (unsort_map::iterator i = PositionTable.begin();
			i != PositionTable.end(); ++i) //each chr
	{
		for (std::unordered_map<int, unsigned int>::iterator j =
				i->second.begin(); j != i->second.end(); ++j) //each site
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
int StatCollector::addFSC(FileStatCollector a)
{
	FSCVec.push_back(a);
	return 0;
}
int StatCollector::SummaryOutput(const string & outputPath,const gap_opt_t* opt)
{
	ofstream fout(outputPath+".summary");
	 int max_XorYmarker(0);
	  if(opt->num_variant_short>=100000)
	    max_XorYmarker=3000;
	  else if(opt->num_variant_short>=10000)
	    max_XorYmarker=300;
	  else
	    max_XorYmarker=100;
	int total_genome_size= (opt->flank_len * 2 + 1) * opt->num_variant_short
			+ (opt->flank_long_len * 2 + 1) * opt->num_variant_long
			+ (501)*max_XorYmarker;
	int total_base(0);
	int total_reads(0);


	fout<<"|:FILE 1:|:FILE 2:|:# Reads:|:Average Length:|"<<endl;
	fout<<"|:------------------------------------------:|:------------------------------------------:|:---------------------:|:---------------------:|"<<endl;
	for(int i=0;i!=FSCVec.size();++i)
	{
		fout<<"|"<<FSCVec[i].FileName1<<"|"<<FSCVec[i].FileName2<<"|"<<FSCVec[i].NumRead<<"|"<<FSCVec[i].NumBase/(double) FSCVec[i].NumRead<<"|"<<endl;
		total_base+=FSCVec[i].NumBase;
		total_reads+=FSCVec[i].NumRead;
	}
	fout<<"|:All:"<<"|:-:|"<<total_reads<<"|"<<total_base/(double)total_reads<<"|"<<endl;
	fout<<endl;
	fout<<"Expected Read Depth = "<<(double)total_base/total_genome_size<<" ["<<total_base<<"/"<<total_genome_size<<"]"<<endl;
	fout<<"Estimated AvgDepth ="<<[&DepthDist](int total_genome_size)->double{int tmp(0);for(int i=0;i!=DepthDist.size();++i) tmp+=i*DepthDist[i]; return double(tmp)/total_genome_size; }<<endl;
	fout<<"Estimated percentage of accessible genome covered ="<<(1-(double)DepthDist[0]/total_genome_size)*100<<"/100"<<endl;
	fout<<"Estimated AvgDepth for Q20 bases ="<<[&Q20DepthVec](int total_genome_size)->double{int tmp(0);for(int i=0;i!=Q20DepthVec.size();++i) tmp+=Q20DepthVec[i]; return double(tmp)/total_genome_size; }<<endl;
	fout<<"Estimated AvgDepth for Q30 bases ="<<[&Q30DepthVec](int total_genome_size)->double{int tmp(0);for(int i=0;i!=Q30DepthVec.size();++i) tmp+=Q30DepthVec[i]; return double(tmp)/total_genome_size; }<<endl;
	fout<<"Median Insert Size(>=500bp) ="<<[&InsertSizeDist]()->double{int tmp(0),total(0);for(int i=500;i!=InsertSizeDist.size();++i) total+=InsertSizeDist[i]; for(int i=500;i!=InsertSizeDist.size();++i) {tmp+=InsertSizeDist[i];if(tmp>total/2) return i;} }<<endl;
	fout<<"Median Insert Size(>=300bp) ="<<[&InsertSizeDist]()->double{int tmp(0),total(0);for(int i=300;i!=InsertSizeDist.size();++i) total+=InsertSizeDist[i]; for(int i=300;i!=InsertSizeDist.size();++i) {tmp+=InsertSizeDist[i];if(tmp>total/2) return i;} }<<endl;
	return 0;
}
StatCollector::~StatCollector()
{
// TODO Auto-generated destructor stub
	for (uint32_t i = 0; i != VcfRecVec.size(); ++i)
		delete VcfRecVec[i];

}


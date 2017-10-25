
#include "StatCollector.h"
#include "InsertSizeEstimator.h"
#include "../libbwa/bwase.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;

//extern string Prefix;
extern void notice(const char *, ...);

extern void warning(const char *, ...);

extern void error(const char *, ...);

static std::string cigar_output(int n_cigar, const bwa_cigar_t *cigar, int len) {
    std::string tmp;
    int cl;
    char cop;
    for (int k = 0; k < n_cigar; ++k) {
        cl = __cigar_len(cigar[k]);
        cop = "MIDS"[__cigar_op(cigar[k])];
        tmp += to_string(cl);
        tmp.push_back(cop);
    }
    if (!cigar)
        tmp = to_string(len) + "M";
    return tmp;
}


int findMaxAllele(const size_t *a, size_t len) {
    uint32_t max = 0;
    uint32_t maxIndex = 0;
    for (uint32_t i = 0; i != len; ++i) {
        if (a[i] > max) {
            max = a[i];
            maxIndex = i;
        }
    }
    return maxIndex;
}

int countAllele(size_t *a, const string &seq) {
    for (uint32_t i = 0; i != seq.size(); ++i) {
        if (seq[i] == 'A') {
            a[0]++;
            continue;
        }
        if (seq[i] == 'C') {
            a[1]++;
            continue;
        }
        if (seq[i] == 'G') {
            a[2]++;
            continue;
        }
        if (seq[i] == 'T') {
            a[3]++;
            continue;
        }
    }
    return 0;
}

#define INSERT_SIZE_LIMIT 4096
#define REV_PHRED(x)    pow(10.0,(x/(-10.0)))
#define PHRED(x)    (-10)*log10(x)

StatCollector::StatCollector() {
    //cerr << "NOTE:Using default initializer..." << endl;
    PositionTable.clear();
    VcfRecVec.clear();
    duplicateTable.clear();
    VariantProxyTable.clear();
    GC.clear();
    index = 0;
    total_base = 0;
    total_region_size = 0;
    ref_genome_size = 0;
    NumPCRDup = 0;
    NumBaseMapped = 0;
    NumPositionCovered = 0;//position with depth larger than 0
    NumPositionCovered2 = 0;//larger than 1
    NumPositionCovered5 = 0;//larger than 4
    NumPositionCovered10 = 0;// larger than 9
    DepthDist = vector<size_t>(1024, 0);//up to depth 1024X
    CycleDist = vector<size_t>(512, 0);
    GCDist = vector<size_t>(256, 0);
    EmpRepDist = vector<size_t>(256, 0);
    misEmpRepDist = vector<size_t>(256, 0);
    EmpCycleDist = vector<size_t>(256, 0);
    misEmpCycleDist = vector<size_t>(256, 0);
    InsertSizeDist = vector<size_t>(INSERT_SIZE_LIMIT, 0);
//    MaxInsertSizeDist = vector<size_t>(INSERT_SIZE_LIMIT, 0);
}

StatCollector::StatCollector(const string &OutFile) {
    PositionTable.clear();
    VcfRecVec.clear();
    duplicateTable.clear();
    VariantProxyTable.clear();
    GC.clear();
    index = 0;
    total_base = 0;
    total_region_size = 0;
    ref_genome_size = 0;
    NumPCRDup = 0;
    NumBaseMapped = 0;
    NumPositionCovered = 0;//position with depth larger than 0
    NumPositionCovered2 = 0;//larger than 1
    NumPositionCovered5 = 0;//larger than 4
    NumPositionCovered10 = 0;// larger than 9
    DepthDist = vector<size_t>(1024, 0);//up to depth 1024X
    CycleDist = vector<size_t>(512, 0);
    GCDist = vector<size_t>(256, 0);
    EmpRepDist = vector<size_t>(256, 0);
    misEmpRepDist = vector<size_t>(256, 0);
    EmpCycleDist = vector<size_t>(256, 0);
    misEmpCycleDist = vector<size_t>(256, 0);
    InsertSizeDist = vector<size_t>(INSERT_SIZE_LIMIT, 0);
//    MaxInsertSizeDist = vector<size_t>(INSERT_SIZE_LIMIT, 0);
}

void StatCollector::StatVecDistUpdate(const string &qual,
                                      int left_to_right_coord, unsigned int tmp_index,
                                      const string &RefSeq, const string &seq, int tmpCycle) {
    EmpRepDist[qual[left_to_right_coord]]++;
    if (qual[left_to_right_coord] >= 20) {
        Q20DepthVec[tmp_index]++;
        if (qual[left_to_right_coord] >= 30) {
            Q30DepthVec[tmp_index]++;
        }
    }
    if (RefSeq[left_to_right_coord] != seq[left_to_right_coord]) {
        misEmpRepDist[qual[left_to_right_coord]]++;
        misEmpCycleDist[tmpCycle]++;
    }
    /************EmpCycle**************************************************************/
    EmpCycleDist[tmpCycle]++;
}

void StatCollector::AddBaseInfoToNewCoord(const string &Chrom, int i,
                                          const string &qual, int left_to_right_coord,
                                          const string &RefSeq, const string &seq, int tmpCycle) {
    //cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
    DepthVec.push_back(0);
    Q30DepthVec.push_back(0);
    Q20DepthVec.push_back(0);
    DepthVec[index]++;
    PositionTable[Chrom][i] = index;
    if (dbSNPTable[Chrom].find(i) == dbSNPTable[Chrom].end()) //not in dbsnp table
    {
        StatVecDistUpdate(qual, left_to_right_coord, index, RefSeq, seq,
                          tmpCycle);
    }
    index++;
}

void StatCollector::UpdateInfoVecAtMarker(int tmpCycleVcfTable,
                                          int tmpCycle, int tmp_left_to_right_coord,
                                          int left_to_right_coord, int realCoord, int cl,
                                          const char *sign, bool strand, const string &Chrom, unsigned int tmp_index,
                                          const string &seq, const string &qual, u_char mapQ) {
    tmpCycleVcfTable = tmpCycle;
    tmp_left_to_right_coord = left_to_right_coord;
    for (int i = realCoord; i != realCoord + cl - 1 + 1;
         ++i, tmpCycleVcfTable += 1 * sign[strand], ++tmp_left_to_right_coord) {
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

bool StatCollector::addSingleAlignment(const bntseq_t *bns, bwa_seq_t *p, const gap_opt_t *opt) //
{
    int seqid(0);
    int j(0);

    if (p->type == BWA_TYPE_NO_MATCH) {
        return false;
    }


    j = static_cast<int>(pos_end(p) - p->pos); //length of read
    bns_coor_pac2real(bns, p->pos, j, &seqid);
    if (p->type != BWA_TYPE_NO_MATCH
        && p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len) {
        //if (strcmp(p->name, "ERR018525.33349") == 0)
        //{
        //	fprintf(stderr, "ERR018525.33349 died here 2\n");
        //}
        return false; //this alignment bridges two adjacent reference sequences
    }

    string seq, qual, newSeq, newQual;

    if (p->strand == 0)
        for (j = 0; j != p->full_len; ++j) {
            seq += "ACGTN"[(int) (p)->seq[j]];
            qual += (char) (p->qual[j] - 33);//change from ASCII to Phred
        }
    else
        for (j = 0; j != p->full_len; ++j) {
            seq += "TGCAN"[(int) (p)->seq[p->full_len - 1 - j]];
            qual += (char) (p->qual[p->full_len - 1 - j] - 33);//0 based
        }
    //if (p->strand) seq_reverse(p->len, p->qual, 0);
    //fprintf(stderr,"%s\t%s\t%d\n",p->name,p->md,p->count);

    //ConstructFakeSeqQual(seq,qual,p->n_cigar,p->cigar,newSeq,newQual);
    string RefSeq(seq), MD(p->md);
    updateRefseqByMD(RefSeq, MD);

    string PosName = string(bns->anns[seqid].name);

    int pos = static_cast<int>(p->pos - bns->anns[seqid].offset + 1);

    size_t colonPos = PosName.find(':');
    size_t atPos = PosName.find('@');
    string Chrom = PosName.substr(0, colonPos);
    char *pEnd;
    int refCoord = static_cast<int>(strtol(
            PosName.substr(colonPos + 1, atPos - colonPos + 1).c_str(), &pEnd, 10)); // coordinate of variant site
    //size_t slashPos = PosName.find("/");
    //string ref = PosName.substr(atPos + 1, slashPos - atPos - 1); //ref allele
    int realCoord(0);
    int RefRealStart(0), RefRealEnd(0);
    if (PosName[PosName.size() - 1] == 'L') {
        realCoord = refCoord - opt->flank_long_len + pos - 1; //real coordinate of current reads on reference
        RefRealStart = refCoord - opt->flank_long_len;
        RefRealEnd = refCoord + opt->flank_long_len;
    } else {
        realCoord = refCoord - opt->flank_len + pos - 1; //real coordinate of current reads on reference
        RefRealStart = refCoord - opt->flank_len;
        RefRealEnd = refCoord + opt->flank_len;
    }
    //realCoord += 100; //
    int tmpCycle(0), tmpCycleVcfTable(0), left_to_right_coord(0), tmp_left_to_right_coord(0);//coord on reads
    char sign[2] =
            {1, -1};
    if (p->strand != 0) {
        tmpCycle = p->len - 1;
    }

    unsigned int tmp_index = 0;
    //cerr << p->name << "\t" << PosName<<"\t"<<pos<<"\t"<<Chrom<< "\t" << realCoord << "\t" << MD <<"\t";
    //cigar string
    if (p->cigar) {
        int n_cigar = p->n_cigar;
        bwa_cigar_t *cigar = p->cigar;
        for (int k = 0; k < n_cigar; ++k) {

            int cl = __cigar_len(cigar[k]);
            int cop = "MIDS"[__cigar_op(cigar[k])];
            switch (cop) {
                case 'M':
                    //cerr<<"M";
                    //cerr<< realCoord<<"-"<< realCoord + cl - 1;
                    /*************************variant site****************************************/

                    if (VcfTable.find(Chrom) != VcfTable.end()) {
//						tmpCycleVcfTable = tmpCycle;
//						tmp_left_to_right_coord = left_to_right_coord;
//						for (uint32_t i = realCoord; i != realCoord + cl - 1 + 1;
//							++i, tmpCycleVcfTable += 1 * sign[p->strand], ++tmp_left_to_right_coord)
//						{
//							if (VcfTable[Chrom].find(i) != VcfTable[Chrom].end())// actual snp site
//							{
//								tmp_index = VcfTable[Chrom][i];
//								SeqVec[tmp_index] += seq[tmp_left_to_right_coord];
//								QualVec[tmp_index] += qual[tmp_left_to_right_coord];
//								CycleVec[tmp_index].push_back(tmpCycleVcfTable);
//								MaqVec[tmp_index].push_back(p->mapQ + 33);
//								StrandVec[tmp_index].push_back(p->strand);
//							}
//						}
                        UpdateInfoVecAtMarker(tmpCycleVcfTable, tmpCycle,
                                              tmp_left_to_right_coord, left_to_right_coord,
                                              realCoord, cl, sign, p->strand, Chrom, tmp_index, seq,
                                              qual, p->mapQ);
                    }

                    /*****************************************************************************/
                    if (PositionTable.find(Chrom) != PositionTable.end()) //chrom exists
                        for (int i = realCoord; i != realCoord + cl - 1 + 1;
                             ++i, tmpCycle += 1 * sign[p->strand], ++left_to_right_coord) {
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
                                    == dbSNPTable[Chrom].end())    //not in dbsnp table
                                {
                                    StatVecDistUpdate(qual, left_to_right_coord,
                                                      tmp_index, RefSeq, seq, tmpCycle);
                                }
                            } else //coord not exists
                            {
                                //cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
                                AddBaseInfoToNewCoord(Chrom, i, qual,
                                                      left_to_right_coord, RefSeq, seq, tmpCycle);
                            }
                        }
                    else // chrom not exists
                    {
                        for (int i = realCoord; i != realCoord + cl - 1 + 1;
                             ++i, tmpCycle += 1 * sign[p->strand], ++left_to_right_coord) {
                            if (i < RefRealStart)
                                continue;
                            if (i > RefRealEnd)
                                break;
                            //cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
                            AddBaseInfoToNewCoord(Chrom, i, qual,
                                                  left_to_right_coord, RefSeq, seq, tmpCycle);
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
            }
        }            //end for
    } else {
        /*************************variant site****************************************/

        if (VcfTable.find(Chrom) != VcfTable.end()) {
//				tmpCycleVcfTable = tmpCycle;
//				tmp_left_to_right_coord = left_to_right_coord;
//				for (uint32_t i = realCoord; i != realCoord + p->len - 1 + 1;
//					++i, tmpCycleVcfTable += 1 * sign[p->strand], ++tmp_left_to_right_coord)
//					if (VcfTable[Chrom].find(i) != VcfTable[Chrom].end())// actual snp site
//					{
//					tmp_index = VcfTable[Chrom][i];
//					SeqVec[tmp_index] += seq[tmp_left_to_right_coord];
//					QualVec[tmp_index] += qual[tmp_left_to_right_coord];
//					CycleVec[tmp_index].push_back(tmpCycleVcfTable);
//					MaqVec[tmp_index].push_back(p->mapQ + 33);
//					StrandVec[tmp_index].push_back(p->strand);
//					}
            UpdateInfoVecAtMarker(tmpCycleVcfTable, tmpCycle,
                                  tmp_left_to_right_coord, left_to_right_coord,
                                  realCoord, p->len, sign, p->strand, Chrom, tmp_index, seq,
                                  qual, p->mapQ);

        }
        /************************************************************************************/
        if (PositionTable.find(Chrom) != PositionTable.end()) //chrom exists
            for (int i = realCoord; i != realCoord + p->len - 1 + 1;
                 ++i, tmpCycle += 1 * sign[p->strand], ++left_to_right_coord) {
                if (i < RefRealStart)
                    continue;
                if (i > RefRealEnd)
                    break;
                //std::pair<unordered_map<int, bool>::iterator, bool> ret =
                //		PositionTable[Chrom].insert(make_pair(i, 1));
                if (PositionTable[Chrom].find(i) !=
                    PositionTable[Chrom].end())    //!ret.second) //if insert failed, this coord exists
                {
                    //cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
                    tmp_index = PositionTable[Chrom][i];
                    DepthVec[tmp_index]++;
                    if (dbSNPTable[Chrom].find(i) == dbSNPTable[Chrom].end())//not in dbsnp table
                    {
                        StatVecDistUpdate(qual, left_to_right_coord, tmp_index,
                                          RefSeq, seq, tmpCycle);
                    }
                } else //coord not exists
                {
                    //cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
                    AddBaseInfoToNewCoord(Chrom, i, qual, left_to_right_coord,
                                          RefSeq, seq, tmpCycle);
                }
            }
        else // chrom not exists
        {
            for (int i = realCoord; i != realCoord + p->len - 1 + 1;
                 ++i, tmpCycle += 1 * sign[p->strand], ++left_to_right_coord) {
                if (i < RefRealStart)
                    continue;
                if (i > RefRealEnd)
                    break;
                //cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
                AddBaseInfoToNewCoord(Chrom, i, qual, left_to_right_coord,
                                      RefSeq, seq, tmpCycle);
            }
        }
    }

    return true;
}

bool StatCollector::addSingleAlignment(SamRecord &p, const gap_opt_t *opt) //
{

    if (p.getFlag() & SAM_FSU) {
        return false;
    }

    string seq(p.getSequence()), qual(p.getQuality()), newSeq, newQual;

    //ConstructFakeSeqQual(seq,qual,p->n_cigar,p->cigar,newSeq,newQual);
    string RefSeq(seq), MD;
    if (p.getStringTag("MD") != nullptr) {
        MD = string(p.getStringTag("MD")->c_str());
    } else
        return false;

    updateRefseqByMD(RefSeq, MD);

    string PosName = p.getReferenceName();

    int pos = p.get1BasedPosition();

    size_t colonPos = PosName.find(':');
    size_t atPos = PosName.find('@');
    string Chrom = PosName.substr(0, colonPos);

//    int refCoord = atoi(
//            PosName.substr(colonPos + 1, atPos - colonPos + 1).c_str()); // coordinate of variant site
    char *pEnd;
    int refCoord = static_cast<int>(strtol(
            PosName.substr(colonPos + 1, atPos - colonPos + 1).c_str(), &pEnd, 10));
    int realCoord(0);
    int RefRealStart(0), RefRealEnd(0);

    if (strchr(p.getReferenceName(), '@') != nullptr)//bam file is generated by FASTQuick itself
    {
        if (PosName[PosName.size() - 1] == 'L') {
            realCoord = refCoord - opt->flank_long_len + pos - 1; //real coordinate of current reads on reference
            RefRealStart = refCoord - opt->flank_long_len;
            RefRealEnd = refCoord + opt->flank_long_len;
        } else {
            realCoord = refCoord - opt->flank_len + pos - 1; //real coordinate of current reads on reference
            RefRealStart = refCoord - opt->flank_len;
            RefRealEnd = refCoord + opt->flank_len;
        }
    } else {
        if (VcfTable.find(Chrom) != VcfTable.end()) {//need to locate which artificial ref seq this reads aligned to
            auto bound = VcfTable[Chrom].lower_bound(pos);

            if (bound != VcfTable[Chrom].end())//right element exists
            {
                auto upper(bound);
                int overlap_right = realCoord + p.getReadLength() - 1 - (upper->first - opt->flank_len) + 1;
                if (bound != VcfTable[Chrom].begin())//both side exists
                {
                    auto lower(--bound);
                    int overlap_left = lower->first + opt->flank_len - realCoord + 1;
                    realCoord = pos;
                    if (overlap_left > overlap_right) {
                        RefRealStart = lower->first - opt->flank_len;
                        RefRealEnd = lower->first + opt->flank_len;
                    } else {
                        RefRealStart = upper->first - opt->flank_len;
                        RefRealEnd = upper->first + opt->flank_len;
                    }
                } else//only right exists
                {
                    RefRealStart = upper->first - opt->flank_len;
                    RefRealEnd = upper->first + opt->flank_len;
                }
            } else if (bound != VcfTable[Chrom].begin())//only left element exists
            {
                auto lower(--bound);
//                int overlap_left = lower->first + opt->flank_len - realCoord + 1;
                realCoord = pos;
                RefRealStart = lower->first - opt->flank_len;
                RefRealEnd = lower->first + opt->flank_len;
            } else//none neighbours, not possible
                return false;
        } else
            return false;
    }


    int tmpCycle(0), tmpCycleVcfTable(0), left_to_right_coord(0), tmp_left_to_right_coord(0);
    char sign[2] =
            {1, -1};

    bool strand = (p.getFlag() & SAM_FSR);
    u_char mapQ = p.getMapQuality();
    if (strand != 0) {
        tmpCycle = p.getReadLength() - 1;
    }

    unsigned int tmp_index = 0;

    if (std::string(p.getCigar()).find_first_of("SID") != std::string::npos)//DONE: change not mapped to only M
    {

        string cigar = p.getCigar();
        int n_cigar = static_cast<int>(cigar.size());
        int digitNow(0), digitLast(0);
        //for (int k = 0; k < n_cigar; ++k)
        while (digitNow != n_cigar) {
            while (isdigit(cigar[digitNow]))
                digitNow++;
            int cl = atoi(
                    cigar.substr(digitLast, digitNow - digitLast).c_str()); //digitNow is "SMID"
            digitLast = digitNow + 1;
            char cop = cigar[digitNow];
            digitNow++;
            //std::cerr<<cigar<<":\t"<<cop<<"\t"<<cl<<std::endl;
            switch (cop) {
                case 'M':
                    /*************************variant site****************************************/

                    if (VcfTable.find(Chrom) != VcfTable.end()) {
                        UpdateInfoVecAtMarker(tmpCycleVcfTable, tmpCycle,
                                              tmp_left_to_right_coord, left_to_right_coord,
                                              realCoord, cl, sign, strand, Chrom, tmp_index, seq,
                                              qual, mapQ);
                    }
                    /*****************************************************************************/
                    if (PositionTable.find(Chrom) != PositionTable.end()) //chrom exists
                        for (int i = realCoord; i != realCoord + cl - 1 + 1;
                             ++i, tmpCycle += 1 * sign[strand], ++left_to_right_coord) {
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
                                    == dbSNPTable[Chrom].end())    //not in dbsnp table
                                {
                                    StatVecDistUpdate(qual, left_to_right_coord,
                                                      tmp_index, RefSeq, seq, tmpCycle);
                                }
                            } else //coord not exists
                            {
                                AddBaseInfoToNewCoord(Chrom, i, qual,
                                                      left_to_right_coord, RefSeq, seq, tmpCycle);
                            }
                        }
                    else // chrom not exists
                    {
                        for (int i = realCoord; i != realCoord + cl - 1 + 1;
                             ++i, tmpCycle += 1 * sign[strand], ++left_to_right_coord) {
                            if (i < RefRealStart)
                                continue;
                            if (i > RefRealEnd)
                                break;
                            AddBaseInfoToNewCoord(Chrom, i, qual,
                                                  left_to_right_coord, RefSeq, seq, tmpCycle);
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
            }
        }            //end for
    } else {
        /*************************variant site****************************************/

        if (VcfTable.find(Chrom) != VcfTable.end()) {
//				tmpCycleVcfTable = tmpCycle;
//				tmp_left_to_right_coord = left_to_right_coord;
//				for (uint32_t i = realCoord;
//					i != realCoord + p.getReadLength() - 1 + 1;
//					++i, tmpCycleVcfTable += 1 * sign[strand], ++tmp_left_to_right_coord)
//					if (VcfTable[Chrom].find(i) != VcfTable[Chrom].end())// actual snp site
//					{
//					tmp_index = VcfTable[Chrom][i];
//					SeqVec[tmp_index] += seq[tmp_left_to_right_coord];
//					QualVec[tmp_index] += qual[tmp_left_to_right_coord];
//					CycleVec[tmp_index].push_back(tmpCycleVcfTable);
//					MaqVec[tmp_index].push_back(mapQ + 33);
//					StrandVec[tmp_index].push_back(strand);
//					}
            UpdateInfoVecAtMarker(tmpCycleVcfTable, tmpCycle,
                                  tmp_left_to_right_coord, left_to_right_coord,
                                  realCoord, p.getReadLength(), sign, strand, Chrom, tmp_index, seq,
                                  qual, mapQ);
        }
        /************************************************************************************/
        if (PositionTable.find(Chrom) != PositionTable.end()) //chrom exists
            for (int i = realCoord;
                 i != realCoord + p.getReadLength() - 1 + 1;
                 ++i, tmpCycle += 1 * sign[strand], ++left_to_right_coord) {
                if (i < RefRealStart)
                    continue;
                if (i > RefRealEnd)
                    break;
                //std::pair<unordered_map<int, bool>::iterator, bool> ret =
                //		PositionTable[Chrom].insert(make_pair(i, 1));
                if (PositionTable[Chrom].find(i) !=
                    PositionTable[Chrom].end())    //!ret.second) //if insert failed, this coord exists
                {
                    //cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
                    tmp_index = PositionTable[Chrom][i];
                    DepthVec[tmp_index]++;
                    if (dbSNPTable[Chrom].find(i) == dbSNPTable[Chrom].end())//not in dbsnp table
                    {
                        StatVecDistUpdate(qual, left_to_right_coord, tmp_index,
                                          RefSeq, seq, tmpCycle);
                    }
                } else //coord not exists
                {
                    //cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
                    AddBaseInfoToNewCoord(Chrom, i, qual, left_to_right_coord,
                                          RefSeq, seq, tmpCycle);
                }
            }
        else // chrom not exists
        {
            for (int i = realCoord;
                 i != realCoord + p.getReadLength() - 1 + 1;
                 ++i, tmpCycle += 1 * sign[strand], ++left_to_right_coord) {
                if (i < RefRealStart)
                    continue;
                if (i > RefRealEnd)
                    break;
                //cerr<<tmpCycle<<"\t"<<(int)sign[p->strand]<<"\t"<<(int)p->strand<<endl;
                AddBaseInfoToNewCoord(Chrom, i, qual, left_to_right_coord,
                                      RefSeq, seq, tmpCycle);
            }
        }
    }

    return true;
}

int StatCollector::IsDuplicated(const bntseq_t *bns, const bwa_seq_t *p,
                                const bwa_seq_t *q, const gap_opt_t *opt, int type,
                                ofstream &fout) {//in this setting p is fastq1 and q is fastq2
    //todo: remember to also update IsDuplicated overloaded version
    int maxInsert(-1), maxInsert2(-1);
    int readLength(0), seqid_p(-1), seqid_q(-1);
    int flag1 = 0;
    int flag2 = 0;
    int threshQual = 20;
    string status;

    int cl1(0), op1(0);//p left
    int cl2(0), op2(0);//p right
    int cl3(0), op3(0);//q left
    int cl4(0), op4(0);//q right

    if (p) {
        flag1 = p->extra_flag;
        if (p->type == BWA_TYPE_NO_MATCH) {
            flag1 |= SAM_FSU;
        }
        if (p->strand) flag1 |= SAM_FSR;
    }
    if (q) {
        flag2 = q->extra_flag;
        if (q->type == BWA_TYPE_NO_MATCH) {
            flag2 |= SAM_FSU;
        }
        if (q->strand) flag2 |= SAM_FSR;
    }

    if (type == 1) //q is aligned only, read2
    {
        readLength = static_cast<int>(pos_end(q) -
                                      q->pos); //length of read including cigar for bns tracking, length of the reference in the alignment
        bns_coor_pac2real(bns, q->pos, readLength, &seqid_q);
        if (q->mapQ >= threshQual) {
            status = "RevOnly";


            if (q->cigar) {
                op3 = "MIDS"[__cigar_op(q->cigar[0])];//left most cigar
                if (op3 == 'S') {
                    cl3 = __cigar_len(q->cigar[0]);
                }

                op4 = "MIDS"[__cigar_op(q->cigar[q->n_cigar - 1])];//right most cigar
                if (op4 == 'S') {
                    cl4 = __cigar_len(q->cigar[q->n_cigar - 1]);
                }
            }

            if (q->strand)//reverse
            {
                if (bns->anns[seqid_q].offset + bns->anns[seqid_q].len >=
                    (q->pos - cl3) + q->len)//soft clip stays within Flank Region
                    maxInsert2 = static_cast<int>((q->pos - cl3) + q->len - bns->anns[seqid_q].offset);
                else//has soft clip and not complete mapped
                    return 2;
            } else//forward
            {
                if ((q->pos - cl3) >= bns->anns[seqid_q].offset)//soft clip stays within Flank Region
                    maxInsert = static_cast<int>(bns->anns[seqid_q].offset + bns->anns[seqid_q].len - (q->pos - cl3));
                else
                    return 2;

                status = "FwdOnly";
            }

            fout << q->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
                 << "\t" << "*" << "\t" << "*" << "\t" << flag1 << "\t" << 0 << "\t" << "*"
                 << "\t" << bns->anns[seqid_q].name << "\t" << q->pos - bns->anns[seqid_q].offset + 1 << "\t" << flag2
                 << "\t" << q->len << "\t" << cigar_output(q->n_cigar, q->cigar, q->len)
                 << "\t" << status
                 << endl;
            return 0;
        } else {
            fout << q->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
                 << "\t" << "*" << "\t" << "*" << "\t" << flag1 << "\t" << 0 << "\t" << "*"
                 << "\t" << bns->anns[seqid_q].name << "\t" << q->pos - bns->anns[seqid_q].offset + 1 << "\t" << flag2
                 << "\t" << q->len << "\t" << cigar_output(q->n_cigar, q->cigar, q->len)
                 << "\t" << "LowQual"
                 << endl;
            return 2; //low quality
        }
    } else if (type == 3) //p is aligned only, read1
    {
        readLength = static_cast<int>(pos_end(p) -
                                      p->pos); //length of read including cigar for bns tracking, length of the reference in the alignment
        bns_coor_pac2real(bns, p->pos, readLength, &seqid_p);

        if (p->mapQ >= threshQual) {
            status = "FwdOnly";

            if (p->cigar) {
                op1 = "MIDS"[__cigar_op(p->cigar[0])];//left most cigar
                if (op1 == 'S') {
                    cl1 = __cigar_len(p->cigar[0]);
                }

                op2 = "MIDS"[__cigar_op(p->cigar[p->n_cigar - 1])];//right most cigar
                if (op2 == 'S') {
                    cl2 = __cigar_len(p->cigar[p->n_cigar - 1]);
                }
            }

            if (p->strand)//reverse
            {
                if (bns->anns[seqid_p].offset + bns->anns[seqid_p].len >=
                    (p->pos - cl1) + p->len)//soft clip stays within Flank Region
                    maxInsert2 = static_cast<int>((p->pos - cl1) + p->len - bns->anns[seqid_p].offset);
                else
                    return 2;

                status = "RevOnly";
            } else {
                if ((p->pos - cl1) >= bns->anns[seqid_p].offset)//soft clip stays within Flank Region
                    maxInsert = static_cast<int>(bns->anns[seqid_p].offset + bns->anns[seqid_p].len - (p->pos - cl1));
                else
                    return 2;
            }
            fout << p->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
                 << "\t" << bns->anns[seqid_p].name << "\t" << p->pos - bns->anns[seqid_p].offset + 1 << "\t" << flag1
                 << "\t" << p->len << "\t" << cigar_output(p->n_cigar, p->cigar, p->len)
                 << "\t" << "*" << "\t" << "*" << "\t" << flag2 << "\t" << 0 << "\t" << "*"
                 << "\t" << status
                 << endl;
            return 0;
        } else {
            fout << p->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
                 << "\t" << bns->anns[seqid_p].name << "\t" << p->pos - bns->anns[seqid_p].offset + 1 << "\t" << flag1
                 << "\t" << p->len << "\t" << cigar_output(p->n_cigar, p->cigar, p->len)
                 << "\t" << "*" << "\t" << "*" << "\t" << flag2 << "\t" << 0 << "\t" << "*"
                 << "\t" << "LowQual"
                 << endl;
            return 2;
        }
    } else if (type == 2) //both aligned
    {
        readLength = static_cast<int>(pos_end(p) - p->pos); //length of read
        bns_coor_pac2real(bns, p->pos, readLength, &seqid_p);
        readLength = static_cast<int>(pos_end(q) - q->pos); //length of read
        bns_coor_pac2real(bns, q->pos, readLength, &seqid_q);
        /*deal with cigar*/

        if (p->cigar) {
            op1 = "MIDS"[__cigar_op(p->cigar[0])];//left most cigar
            if (op1 == 'S') {
                cl1 = __cigar_len(p->cigar[0]);
            }

            op2 = "MIDS"[__cigar_op(p->cigar[p->n_cigar - 1])];//right most cigar
            if (op2 == 'S') {
                cl2 = __cigar_len(p->cigar[p->n_cigar - 1]);
            }
        }

        if (q->cigar) {
            op3 = "MIDS"[__cigar_op(q->cigar[0])];//left most cigar
            if (op3 == 'S') {
                cl3 = __cigar_len(q->cigar[0]);
            }

            op4 = "MIDS"[__cigar_op(q->cigar[q->n_cigar - 1])];//right most cigar
            if (op4 == 'S') {
                cl4 = __cigar_len(q->cigar[q->n_cigar - 1]);
            }
        }
        /*end deal with cigar*/
        if (!(p->strand) && q->strand && p->pos < q->pos)//FR
        {

            if ((p->pos - cl1) >= bns->anns[seqid_p].offset)//soft clip stays within Flank Region
                maxInsert = static_cast<int>(bns->anns[seqid_p].offset + bns->anns[seqid_p].len - (p->pos - cl1));
            else
                maxInsert = -1;

            if (bns->anns[seqid_q].offset + bns->anns[seqid_q].len >=
                (q->pos - cl3) + q->len)//soft clip stays within Flank Region
                maxInsert2 = static_cast<int>((q->pos - cl3) + q->len - bns->anns[seqid_q].offset);
            else
                maxInsert2 = -1;

        } else if (!(q->strand) && p->strand && q->pos < p->pos)//FR but rotated
        {

            if ((q->pos - cl3) >= bns->anns[seqid_q].offset)//soft clip stays within Flank Region
                maxInsert = static_cast<int>(bns->anns[seqid_q].offset + bns->anns[seqid_q].len - (q->pos - cl3));
            else
                maxInsert = -1;

            if (bns->anns[seqid_p].offset + bns->anns[seqid_p].len >=
                (p->pos - cl1) + p->len)//soft clip stays within Flank Region
                maxInsert2 = static_cast<int>((p->pos - cl1) + p->len - bns->anns[seqid_p].offset);
            else
                maxInsert2 = -1;
        } else//other than FR, treat as DiffChrom
        {
            fout << p->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
                 << "\t" << bns->anns[seqid_p].name << "\t" << p->pos - bns->anns[seqid_p].offset + 1 << "\t" << flag1
                 << "\t" << p->len << "\t" << cigar_output(p->n_cigar, p->cigar, p->len)
                 << "\t" << bns->anns[seqid_q].name << "\t" << q->pos - bns->anns[seqid_q].offset + 1 << "\t" << flag2
                 << "\t" << q->len << "\t" << cigar_output(q->n_cigar, q->cigar, q->len)
                 << "\tNotPair"
                 << endl;
            return 0;
        }
    }
//from now on only type 2 paired reads remained
    if (maxInsert >= INSERT_SIZE_LIMIT)
        maxInsert = INSERT_SIZE_LIMIT - 1;
    if (maxInsert2 >= INSERT_SIZE_LIMIT)
        maxInsert2 = INSERT_SIZE_LIMIT - 1;
//    MaxInsertSizeDist[maxInsert]++;
//    MaxInsertSizeDist[maxInsert2]++;
    if (seqid_p != seqid_q && seqid_p != -1 && seqid_q != -1) {
        //int ActualInsert(-1);
        //	ActualInsert=0;
        InsertSizeDist[0]++;
        //cerr<<"Duplicate function exit from diff chrom"<<endl;
        fout << p->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
             << "\t" << bns->anns[seqid_p].name << "\t" << p->pos - bns->anns[seqid_p].offset + 1 << "\t" << flag1
             << "\t" << p->len << "\t" << cigar_output(p->n_cigar, p->cigar, p->len)
             << "\t" << bns->anns[seqid_q].name << "\t" << q->pos - bns->anns[seqid_q].offset + 1 << "\t" << flag2
             << "\t" << q->len << "\t" << cigar_output(q->n_cigar, q->cigar, q->len)
             << "\tNotPair"
             << endl;
        return 0;
    }
    if (p->mapQ >= threshQual && q->mapQ >= threshQual) {
        int ActualInsert(-1), start(0), end(0);//[start, end)

        status = "PartialPair";
        if (!(p->strand) && q->strand && p->pos < q->pos)//FR
        {
            start = p->pos - cl1;
            end = q->pos - cl3 + q->len;
            ActualInsert = end - start;
        } else if (!(q->strand) && p->strand && q->pos < p->pos)//FR but rotated
        {
            start = q->pos - cl3;
            end = p->pos - cl1 + p->len;
            ActualInsert = end - start;
        }

        if (maxInsert != -1 and maxInsert2 != -1) status = "PropPair";

        InsertSizeDist[ActualInsert]++;
        fout << p->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << ActualInsert
             << "\t" << bns->anns[seqid_p].name << "\t" << p->pos - bns->anns[seqid_p].offset + 1 << "\t" << flag1
             << "\t" << p->len << "\t" << cigar_output(p->n_cigar, p->cigar, p->len)
             << "\t" << bns->anns[seqid_q].name << "\t" << q->pos - bns->anns[seqid_q].offset + 1 << "\t" << flag2
             << "\t" << q->len << "\t" << cigar_output(q->n_cigar, q->cigar, q->len)
             << "\t" << status
             << endl;

        char start_end[256];
        sprintf(start_end, "%d:%d", start, end);
        pair<unordered_map<string, bool>::iterator, bool> iter =
                duplicateTable.insert(make_pair(string(start_end), true));
        if (!iter.second) //insert failed, duplicated
        {
            //	cerr<<"Duplicate function exit from duplicate"<<endl;
            NumPCRDup++;
            return 1;
        }
    } else {
        //cerr<<"exit from LowQualf"<<endl;
        fout << p->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
             << "\t" << bns->anns[seqid_p].name << "\t" << p->pos - bns->anns[seqid_p].offset + 1 << "\t" << flag1
             << "\t" << p->len << "\t" << cigar_output(p->n_cigar, p->cigar, p->len)
             << "\t" << bns->anns[seqid_q].name << "\t" << q->pos - bns->anns[seqid_q].offset + 1 << "\t" << flag2
             << "\t" << q->len << "\t" << cigar_output(q->n_cigar, q->cigar, q->len)
             << "\tLowQual"
             << endl;
        return 2; //low quality
    }
    return 0;
}

//overload of function IsDuplicated for direct bam reading
int StatCollector::IsDuplicated(SamFileHeader &SFH, SamRecord &p, SamRecord &q,
                                const gap_opt_t *opt, int type, ofstream &fout) {
    int MaxInsert(-1), MaxInsert2(-1);

    if (type == 1) //q is aligned only
    {
        if (q.getFlag() & SAM_FR1) //if it's read one
            MaxInsert2 = atoi(SFH.getSQTagValue("LN", q.getReferenceName()))
                         - q.get1BasedPosition() + 1;
        else
            MaxInsert = q.get1BasedPosition() + q.getReadLength();

        int cl2;
        cl2 = atoi(q.getCigar());
        if (q.getCigar()[to_string(cl2).size()] != 'S')//if the first character is 'S'
        {
            cl2 = 0;
        } else {
            if (q.getFlag() & SAM_FSR)//reversed
                MaxInsert2 -= cl2;
            else
                MaxInsert2 += cl2;
        }


        fout << q.getReadName() << "\t" << -1 << "\t" << MaxInsert2 << "\t" << -1
             << "\t" << "*" << "\t" << "*" << "\t" << p.getFlag() << "\t" << p.getReadLength() << "\t" << p.getCigar()
             << "\t" << q.getReferenceName() << "\t" << q.get1BasedPosition() << "\t" << q.getFlag() << "\t"
             << q.getReadLength() << "\t" << q.getCigar()
             << "\tRevOnly" << endl;
        return 0;
    } else if (type == 3) //p is aligned only
    {
        if (p.getFlag() & SAM_FR1) //if it's read one
            MaxInsert = atoi(SFH.getSQTagValue("LN", p.getReferenceName()))
                        - p.get1BasedPosition() + 1;
        else
            MaxInsert = p.get1BasedPosition() + p.getReadLength();

        int cl2;
        cl2 = atoi(p.getCigar());
        if (p.getCigar()[to_string(cl2).size()] != 'S')//if the first character is 'S'
        {
            cl2 = 0;
        } else {
            if (p.getFlag() & SAM_FSR)//reversed
                MaxInsert -= cl2;
            else
                MaxInsert += cl2;
        }
        fout << p.getReadName() << "\t" << MaxInsert << "\t" << -1 << "\t" << -1
             << "\t" << p.getReferenceName() << "\t" << p.get1BasedPosition() << "\t" << p.getFlag() << "\t"
             << p.getReadLength() << "\t" << p.getCigar()
             << "\t*" << "\t" << "*" << "\t" << q.getFlag() << "\t" << q.getReadLength() << "\t" << q.getCigar()
             << "\tFwdOnly" << endl;
        return 0;
    } else if (type == 2) // both aligned
    {
        if (p.get1BasedPosition() < q.get1BasedPosition()) {
            //j = pos_end(q) - q->pos; //length of read
            //bns_coor_pac2real(bns, q->pos, j, &seqid_q);
            MaxInsert = p.get1BasedPosition() + p.getReadLength();
            //j = pos_end(p) - p->pos; //length of read
            //bns_coor_pac2real(bns, p->pos, j, &seqid_p);
            MaxInsert2 = atoi(SFH.getSQTagValue("LN", q.getReferenceName()))
                         - q.get1BasedPosition() + 1;
        } else {
            MaxInsert = q.get1BasedPosition() + q.getReadLength();
            //j = pos_end(p) - p->pos; //length of read
            //bns_coor_pac2real(bns, p->pos, j, &seqid_p);
            MaxInsert2 = atoi(SFH.getSQTagValue("LN", p.getReferenceName()))
                         - p.get1BasedPosition() + 1;
        }
    } else {
        warning("Alignment status fatal error!\n");
        exit(1);
    }

//	if (MaxInsert2 > MaxInsert)
//		MaxInsert = MaxInsert2;
    if (MaxInsert >= INSERT_SIZE_LIMIT)
        MaxInsert = INSERT_SIZE_LIMIT - 1;
    if (MaxInsert2 >= INSERT_SIZE_LIMIT)
        MaxInsert2 = INSERT_SIZE_LIMIT - 1;
//    MaxInsertSizeDist[MaxInsert]++;
    if (strcmp(p.getReferenceName(), q.getReferenceName()) != 0) {

        InsertSizeDist[0]++;
        //cerr<<"Duplicate function exit from diff chrom"<<endl;
        fout << p.getReadName() << "\t" << MaxInsert << "\t" << MaxInsert2 << "\t" << -1
             << "\t" << p.getReferenceName() << "\t" << p.get1BasedPosition() << "\t" << p.getFlag() << "\t"
             << p.getReadLength() << "\t" << p.getCigar()
             << "\t" << q.getReferenceName() << "\t" << q.get1BasedPosition() << "\t" << q.getFlag() << "\t"
             << q.getReadLength() << "\t" << q.getCigar()
             << "\tDiffChrom" << endl;
        return 0;
    }
    if (p.getMapQuality() >= 20 && q.getMapQuality() >= 20) {
        int ActualInsert(-1), start(0), end(0);

        if (p.get1BasedPosition() < q.get1BasedPosition()) {
            start = p.get1BasedPosition();
            end = q.get1BasedPosition() + q.getReadLength();
            ActualInsert = end - start;
        } else {
            start = q.get1BasedPosition();
            end = p.get1BasedPosition() + p.getReadLength();
            ActualInsert = end - start;
        }
        if (ActualInsert < INSERT_SIZE_LIMIT) {

            string cigar1 = p.getCigar();
            string cigar2 = q.getCigar();
            if (cigar1.find('S') == cigar1.npos && cigar2.find('S') == cigar1.npos)//no softclip
            {
                InsertSizeDist[ActualInsert]++;
                fout << p.getReadName() << "\t" << MaxInsert << "\t" << MaxInsert2 << "\t" << ActualInsert
                     << "\t" << p.getReferenceName() << "\t" << p.get1BasedPosition() << "\t" << p.getFlag() << "\t"
                     << p.getReadLength() << "\t" << p.getCigar()
                     << "\t" << q.getReferenceName() << "\t" << q.get1BasedPosition() << "\t" << q.getFlag() << "\t"
                     << q.getReadLength() << "\t" << q.getCigar()
                     << "\tPropPair" << endl;
            } else {
                int cl1, cl2;
                cl1 = atoi(p.getCigar());
                if (p.getCigar()[to_string(cl1).size()] != 'S')//if the first character is 'S'
                    cl1 = 0;
                cl2 = atoi(q.getCigar());
                if (q.getCigar()[to_string(cl2).size()] != 'S')//if the first character is 'S'
                    cl2 = 0;

                if (p.get1BasedPosition() < q.get1BasedPosition()) {
                    ActualInsert = ActualInsert + cl1 - cl2;
                    MaxInsert = MaxInsert + cl1;
                    MaxInsert2 = MaxInsert2 - cl2;
                } else {
                    ActualInsert = ActualInsert - cl1 + cl2;
                    MaxInsert = MaxInsert - cl1;
                    MaxInsert2 = MaxInsert2 + cl2;
                }
                InsertSizeDist[ActualInsert]++;
                fout << p.getReadName() << "\t" << MaxInsert << "\t" << MaxInsert2 << "\t" << ActualInsert
                     << "\t" << p.getReferenceName() << "\t" << p.get1BasedPosition() << "\t" << p.getFlag() << "\t"
                     << p.getReadLength() << "\t" << p.getCigar()
                     << "\t" << q.getReferenceName() << "\t" << q.get1BasedPosition() << "\t" << q.getFlag() << "\t"
                     << q.getReadLength() << "\t" << q.getCigar()
                     << "\tPartialPair" << endl;
            }
            char start_end[256];
            sprintf(start_end, "%d:%d", start, end);
            pair<unordered_map<string, bool>::iterator, bool> iter =
                    duplicateTable.insert(make_pair(string(start_end), true));
            if (!iter.second) //insert failed, duplicated
            {
                //	cerr<<"Duplicate function exit from duplicate"<<endl;
                NumPCRDup++;
                return 1;
            }
        } else {
            //cerr<<"exit from Inf"<<endl;
            fout << p.getReadName() << "\t" << MaxInsert << "\t" << MaxInsert2 << "\t" << -1
                 << "\t" << p.getReferenceName() << "\t" << p.get1BasedPosition() << "\t" << p.getFlag() << "\t"
                 << p.getReadLength() << "\t" << p.getCigar()
                 << "\t" << q.getReferenceName() << "\t" << q.get1BasedPosition() << "\t" << q.getFlag() << "\t"
                 << q.getReadLength() << "\t" << q.getCigar()
                 << "\tAbnormal" << endl;
        }
    } else {
        int ActualInsert(-1), start(0), end(0);
        if (p.get1BasedPosition() < q.get1BasedPosition()) {
            start = p.get1BasedPosition();
            end = q.get1BasedPosition() + q.getReadLength();
            ActualInsert = end - start;
        } else {
            start = q.get1BasedPosition();
            end = p.get1BasedPosition() + p.getReadLength();
            ActualInsert = end - start;
        }
        //cerr<<"exit from LowQualf"<<endl;
        fout << p.getReadName() << "\t" << MaxInsert << "\t" << MaxInsert2 << "\t" << -1
             << "\t" << p.getReferenceName() << "\t" << p.get1BasedPosition() << "\t" << p.getFlag() << "\t"
             << p.getReadLength() << "\t" << p.getCigar()
             << "\t" << q.getReferenceName() << "\t" << q.get1BasedPosition() << "\t" << q.getFlag() << "\t"
             << q.getReadLength() << "\t" << q.getCigar()
             << "\tLowQual" << endl;
        return 2; //low quality
    }
    return 0;
}

//return value: 0 failed, 1 add pair success, 2 add single success by using add pair interface
int StatCollector::addAlignment(const bntseq_t *bns, bwa_seq_t *p, bwa_seq_t *q,
                                const gap_opt_t *opt, ofstream &fout, int &total_add) {
    int seqid(0), seqid2(0), j(0), j2(0);
    /*checking if reads bridges two reference contigs*/
    if (p and p->type != BWA_TYPE_NO_MATCH) {
        j = pos_end(p) - p->pos; //length of read
        bns_coor_pac2real(bns, p->pos, j, &seqid);
        if (p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len) {
            p->type = BWA_TYPE_NO_MATCH; //this alignment bridges two adjacent reference sequences
        }
    }
    if (q and q->type != BWA_TYPE_NO_MATCH) {
        j2 = pos_end(q) - q->pos; //length of read
        bns_coor_pac2real(bns, q->pos, j2, &seqid2);
        if (q->pos + j - bns->anns[seqid2].offset > bns->anns[seqid2].len) {
            q->type = BWA_TYPE_NO_MATCH; //this alignment bridges two adjacent reference sequences
        }
    }
    /*done checking if reads bridges two reference contigs*/
    if (p == 0 || p->type == BWA_TYPE_NO_MATCH) {
        if (q != 0 && addSingleAlignment(bns, q, opt)) //adding single via pair interface
        {
            string qname(bns->anns[seqid2].name);
            if (string(qname).find("Y") != string::npos
                || string(qname).find("X") != string::npos) {
                if (!isPartialAlign(q)) {
                    contigStatusTable[qname].addNumOverlappedReads();
                    contigStatusTable[qname].addNumFullyIncludedReads();
                } else
                    contigStatusTable[qname].addNumOverlappedReads();
            }
            //TO-DO: adding IsDup function here to generate single end MAX insersize info
            IsDuplicated(bns, p, q, opt, 1, fout);//only q
            total_add++;
            return 2;
        }
        return 0;
    }

    //until now p is aligned
    string pname(bns->anns[seqid].name);
    if (q == 0 || q->type == BWA_TYPE_NO_MATCH) {
        if (addSingleAlignment(bns, p, opt)) //adding single via pair interface
        {
            if (string(pname).find('Y') != string::npos
                || string(pname).find('X') != string::npos) {
                if (!isPartialAlign(p)) {
                    contigStatusTable[pname].addNumOverlappedReads();
                    contigStatusTable[pname].addNumFullyIncludedReads();
                } else
                    contigStatusTable[pname].addNumOverlappedReads();
            }
            //TO-DO:adding IsDup function here to generate single end MAX insersize info
            IsDuplicated(bns, p, q, opt, 3, fout);//only p
            total_add++;
            return 2;
        }
        return 0;
    }

    string qname(bns->anns[seqid2].name);
    //until now both reads are aligned
    if (isPartialAlign(p)) //p is partially aligned
    {
        if (string(qname).find("Y") != string::npos
            || string(qname).find("X") != string::npos) {
            if (isPartialAlign(q)) {
                contigStatusTable[qname].addNumOverlappedReads();
            } else //q is perfectly aligned
            {
                contigStatusTable[qname].addNumOverlappedReads();
                contigStatusTable[qname].addNumFullyIncludedReads();
            }
            if (pname == qname) {
                contigStatusTable[qname].addNumPairOverlappedReads();
            }
            contigStatusTable[pname].addNumOverlappedReads();
        }

        if (IsDuplicated(bns, p, q, opt, 2, fout) != 1 || opt->cal_dup) {
            if (addSingleAlignment(bns, p, opt)) {
                if (addSingleAlignment(bns, q, opt)) {
                    total_add += 2;
                    return 1;
                } else {
                    total_add += 1;
                    return 2;
                }
            } else {
                if (addSingleAlignment(bns, q, opt)) {
                    total_add += 1;
                    return 2;
                } else
                    return 0;
            }
        }
        return 0;
    } else //p is perfectly aligned
    {
        // p and q are both aligned
        if (string(qname).find('Y') != string::npos
            || string(qname).find('X') != string::npos) {
            if (isPartialAlign(q)) //q is partially aligned
            {
                contigStatusTable[qname].addNumOverlappedReads();
                if (pname == qname) {
                    contigStatusTable[qname].addNumPairOverlappedReads();
                }
            } else //q is perfectly aligned
            {
                contigStatusTable[qname].addNumOverlappedReads();
                contigStatusTable[qname].addNumFullyIncludedReads();
                if (pname == qname) {
                    contigStatusTable[qname].addNumPairOverlappedReads();
                    contigStatusTable[qname].addNumFullyIncludedPairedReads();
                }
            }
            contigStatusTable[pname].addNumOverlappedReads();
            contigStatusTable[pname].addNumFullyIncludedReads();
        }
        if (IsDuplicated(bns, p, q, opt, 2, fout) != 1 || opt->cal_dup) {
            if (addSingleAlignment(bns, p, opt)) {
                if (addSingleAlignment(bns, q, opt)) {
                    total_add += 2;
                    return 1;
                } else {
                    total_add += 1;
                    return 2;
                }
            } else {
                if (addSingleAlignment(bns, q, opt)) {
                    total_add += 1;
                    return 2;
                } else
                    return 0;
            }
        } else {
            return 0;
        }
    }
    //cerr << "currently added reads " << total_add << endl;
    return 1;
}

//overload of addAlignment function for direct bam reading
//return value: 0 failed, 1 add pair success, 2 add single success by using add pair interface
int StatCollector::addAlignment(SamFileHeader &SFH, SamRecord *p,
                                SamRecord *q, const gap_opt_t *opt, ofstream &fout, int &total_add) {
    SamValidator Validator;
    SamValidationErrors VErrors;
    Validator.isValid(SFH, *p, VErrors);
    if (VErrors.numErrors() > 0)
        fprintf(stderr, "%s", VErrors.getNextError()->getMessage());
    //CycleDist[p->getReadLength()]++;
    //CycleDist[q->getReadLength()]++;
    //	if (p->type == BWA_TYPE_NO_MATCH)
    if (p == 0 || (p->getFlag() & SAM_FSU)) {
        if (q == 0 || (q->getFlag() & SAM_FSU)) //both end are not mapped
            return 0;
        else if (isPartialAlign(*q))  //q is partially mapped, p is not
        {
            string qname(q->getReferenceName());
            if (string(qname).find('Y') != string::npos
                || string(qname).find('X') != string::npos) {
                contigStatusTable[qname].addNumOverlappedReads();
            }
            if (addSingleAlignment(*q, opt)) //adding single via pair interface
            {
                //TO-DO:adding IsDup function here to generate single end MAX insersize info
                IsDuplicated(SFH, *p, *q, opt, 1, fout);// 1 is q
                total_add++;
                return 2;
            }
            return 0;
        } else // q is perfectly aligned, p is not
        {
            string qname(q->getReferenceName());
            if (string(qname).find('Y') != string::npos
                || string(qname).find('X') != string::npos) {
                contigStatusTable[qname].addNumOverlappedReads();
                contigStatusTable[qname].addNumFullyIncludedReads();
            }
            if (addSingleAlignment(*q, opt)) //adding single via pair interface
            {
                //TO-DO:adding IsDup function here to generate single end MAX insersize info
                IsDuplicated(SFH, *p, *q, opt, 1, fout);// 1 is q
                total_add++;
                return 2;
            } else
                return 0;
        }
    } else if (isPartialAlign(*p)) //p is partially aligned
    {
        if (q == 0 || (q->getFlag() & SAM_FSU)) //q is not aligned
        {
            string pname(p->getReferenceName());
            if (string(pname).find('Y') != string::npos
                || string(pname).find('X') != string::npos) {
                contigStatusTable[pname].addNumOverlappedReads();
            }
            if (addSingleAlignment(*p, opt)) //adding single via pair interface
            {
                //TO-DO:adding IsDup function here to generate single end MAX insersize info
                IsDuplicated(SFH, *p, *q, opt, 3, fout);// 3 is p
                total_add++;
                return 2;
            }
            return 0;
        } else // p is partially aligned q is aligned too
        {
            string pname(p->getReferenceName()), qname(q->getReferenceName());
            if (string(qname).find('Y') != string::npos
                || string(qname).find('X') != string::npos) {
                if (isPartialAlign(*q)) {
                    contigStatusTable[qname].addNumOverlappedReads();
                } else //q is perfectly aligned
                {
                    contigStatusTable[qname].addNumOverlappedReads();
                    contigStatusTable[qname].addNumFullyIncludedReads();
                }
                if (pname == qname) {
                    contigStatusTable[qname].addNumPairOverlappedReads();
                }
                contigStatusTable[pname].addNumOverlappedReads();
            }

            if (IsDuplicated(SFH, *p, *q, opt, 2, fout) != 1 || opt->cal_dup) {
                if (addSingleAlignment(*p, opt)) {
                    if (addSingleAlignment(*q, opt)) {
                        total_add += 2;
                        return 1;
                    } else {
                        total_add += 1;
                        return 2;
                    }
                } else {
                    if (addSingleAlignment(*q, opt)) {
                        total_add += 1;
                        return 2;
                    } else
                        return 0;
                }
            }
            return 0;
        }
    } else //p is perfectly aligned
    {
        if (q == 0 || (q->getFlag() & SAM_FSU)) // p is perfectly aligned, q is not
        {
            string pname(p->getReferenceName());
            if (string(pname).find('Y') != string::npos
                || string(pname).find('X') != string::npos) {
                contigStatusTable[pname].addNumOverlappedReads();
                contigStatusTable[pname].addNumFullyIncludedReads();
            }
            if (addSingleAlignment(*p, opt)) //adding single via pair interface
            {
                //TO-DO:adding IsDup function here to generate single end MAX insersize info
                IsDuplicated(SFH, *p, *q, opt, 3, fout);// 3 is p
                total_add++;
                return 2;
            } else
                return 0;
        } else // p and q are both aligned
        {
            string pname(p->getReferenceName()), qname(q->getReferenceName());
            if (string(qname).find('Y') != string::npos
                || string(qname).find('X') != string::npos) {
                if (isPartialAlign(*q)) //q is partially aligned
                {
                    contigStatusTable[qname].addNumOverlappedReads();
                    if (pname == qname) {
                        contigStatusTable[qname].addNumPairOverlappedReads();
                    }
                } else //q is perfectly aligned
                {
                    contigStatusTable[qname].addNumOverlappedReads();
                    contigStatusTable[qname].addNumFullyIncludedReads();
                    if (pname == qname) {
                        contigStatusTable[qname].addNumPairOverlappedReads();
                        contigStatusTable[qname].addNumFullyIncludedPairedReads();
                    }
                }
                contigStatusTable[pname].addNumOverlappedReads();
                contigStatusTable[pname].addNumFullyIncludedReads();
            }

            if (IsDuplicated(SFH, *p, *q, opt, 2, fout) != 1 || opt->cal_dup) {
                if (addSingleAlignment(*p, opt)) {
                    if (addSingleAlignment(*q, opt)) {
                        total_add += 2;
                        return 1;
                    } else {
                        total_add += 1;
                        return 2;
                    }
                } else {
                    if (addSingleAlignment(*q, opt)) {
                        total_add += 1;
                        return 2;
                    } else
                        return 0;
                }
            } else
                return 0;
        }
    }
    //cerr << "currently added reads " << total_add << endl;
    return 1;
}

int StatCollector::ReadAlignmentFromBam(const gap_opt_t *opt,
                                        const char *BamFile, std::ofstream &fout, int &total_add) {
    SamFileHeader SFH;
    SamFile SFIO;
    if (!SFIO.OpenForRead(BamFile, &SFH)) {
        cerr << SFIO.GetStatusMessage() << endl;
        warning("Reading Bam Header Failed!\n");
        exit(1);
    }
    FileStatCollector FSC(BamFile);
    notice("Start reading Bam file...\n");
    unordered_map<string, SamRecord *> pairBuffer;
    while (1) {
        SamRecord *SR = new SamRecord;
        //SFIO.ReadRecord( SFH, *SR);
        if (!SFIO.ReadRecord(SFH, *SR)) {
            cerr << SFIO.GetStatusMessage() << endl;
            notice("End Bam File Reading...\n");
            delete SR;
            break;
        }
        FSC.NumRead++;
        FSC.NumBase += SR->getReadLength();
        string readName;
        if (SR->getReadName()[SR->getReadNameLength() - 2] == '\\')
            readName = string(SR->getReadName()).substr(0,
                                                        SR->getReadNameLength() - 3);
        else
            readName = string(SR->getReadName());
        if (!(SR->getFlag() & SAM_FPP)) // read is not mapped in pair
        {
            addAlignment(SFH, SR, 0, opt, fout, total_add);
            delete SR;
            continue;
        } else if (pairBuffer.find(readName) == pairBuffer.end()) // mate not existed
        {
            pairBuffer.insert(make_pair(readName, SR));
        } else {
            SamRecord *SRtmp = pairBuffer[readName];
            addAlignment(SFH, SR, SRtmp, opt, fout, total_add);
            delete SR;
            delete SRtmp;
            pairBuffer.erase(readName);
        }
    } //end of while
    addFSC(FSC);
    for (unordered_map<string, SamRecord *>::iterator iter = pairBuffer.begin();
         iter != pairBuffer.end(); ++iter) {
        addAlignment(SFH, iter->second, 0, opt, fout, total_add);
        delete iter->second;
    }
    return 0;
}

int StatCollector::restoreVcfSites(const string &VcfPath, const gap_opt_t *opt) {

    VcfHeader header;
    VcfFileReader reader;
    string SelectedSite = VcfPath + ".SelectedSite.vcf.gz";
    if (!reader.open(SelectedSite.c_str(), header)) {
        warning("File open failed: %s\n", SelectedSite.c_str());
    }
    string GCpath = VcfPath + ".gc";
    ifstream FGC(GCpath, ios_base::binary);
    //int num_so_far(0);
    while (!reader.isEOF()) {
        VcfRecord *VcfLine = new VcfRecord;

        reader.readRecord(*VcfLine);
        VcfRecVec.push_back(VcfLine);

        string chr(VcfLine->getChromStr());
        int pos = VcfLine->get1BasedPosition();
        VcfTable[chr][pos] = VcfRecVec.size() - 1;
        _GCstruct GCstruct;
        GCstruct.read(FGC);
        int tmp_pos = pos - (GCstruct.len - 1) / 2;
        for (uint32_t i = 0; i != GCstruct.len; ++i) {
            GC[chr][tmp_pos + i] = GCstruct.GC[i];
        }
        //cerr<<chr<<"\t"<<pos<<"\t"<<num_so_far<<endl;
        //num_so_far++;

        SeqVec.push_back(string(""));
        QualVec.push_back(string(""));
        CycleVec.push_back(vector<int>(0));
        MaqVec.push_back(vector<unsigned char>(0));
        StrandVec.push_back(vector<bool>(0));
        //VcfTable[Chrom][i]=vcf_index;
        //vcf_index++;

    }
    reader.close();
    string BedFile = VcfPath + ".dpSNP.subset.vcf";
    if (!reader.open(BedFile.c_str(), header)) {
        notice("Open %s failed!\n", BedFile.c_str());
        exit(1);
    }
    while (!reader.isEOF()) {
        VcfRecord VcfLine;
        reader.readRecord(VcfLine);
        string chr(VcfLine.getChromStr());
        int pos = VcfLine.get1BasedPosition();
        dbSNPTable[chr][pos] = 1;
    }
    reader.close();
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

int StatCollector::releaseVcfSites() {

    for (int j = 0; j < VcfRecVec.size(); ++j) {
        delete VcfRecVec[j];
    }
    VcfRecVec.clear();
    VcfTable.clear();
    GC.clear();
    SeqVec.clear();
    QualVec.clear();
    CycleVec.clear();
    MaqVec.clear();
    StrandVec.clear();
    dbSNPTable.clear();
    return 0;
}

int StatCollector::getDepthDist(const string &outputPath, const gap_opt_t *opt) {

    ofstream fout(outputPath + ".DepthDist");
    int max_XorYmarker(0);
    if (opt->num_variant_short >= 100000)
        max_XorYmarker = 3000;
    else if (opt->num_variant_short >= 10000)
        max_XorYmarker = 300;
    else
        max_XorYmarker = 100;

    for (size_t i = 1; i != DepthDist.size(); ++i) {
        NumPositionCovered += DepthDist[i];
        if (i >= 2)
            NumPositionCovered2 += DepthDist[i];
        if (i >= 5)
            NumPositionCovered5 += DepthDist[i];
        if (i >= 10)
            NumPositionCovered10 += DepthDist[i];
    }
    total_region_size = (opt->flank_len  * 2 + 1) * opt->num_variant_short
                        + (opt->flank_long_len  * 2 + 1) * opt->num_variant_long
                        + 501 * max_XorYmarker;
    fout << 0 << "\t" << total_region_size - NumPositionCovered << endl;
    DepthDist[0] = total_region_size - NumPositionCovered;
    for (uint32_t i = 1; i != DepthDist.size(); ++i) {
        fout << i << "\t" << DepthDist[i] << endl;
    }
    fout.close();
    return 0;
}

int StatCollector::getGCDist(const string &outputPath,
                             const vector<int> &PosNum) {
    ofstream fout(outputPath + ".GCDist");
    double MeanDepth = NumBaseMapped / (double) NumPositionCovered;
    for (uint32_t i = 0; i != GCDist.size(); ++i) {
        fout << i << "\t" << GCDist[i] << "\t" << PosNum[i] << "\t";
        if (PosNum[i] == 0) {
            fout << 0;
        } else {
            fout << (double(GCDist[i]) / PosNum[i]) / MeanDepth;
        }
        fout << endl;
    }
    fout.close();
    return 0;
}

int StatCollector::getEmpRepDist(const string &outputPath) {
    ofstream fout(outputPath + ".EmpRepDist");
    for (uint32_t i = 0; i != EmpRepDist.size(); ++i) {
        fout << i << "\t" << (misEmpRepDist[i]) << "\t" << (EmpRepDist[i])
             << "\t"
             << PHRED((double) (misEmpRepDist[i] + 1) / (EmpRepDist[i] + 1))
             << endl;
    }
    fout.close();
    return 0;
}

int StatCollector::getEmpCycleDist(const string &outputPath) {
    ofstream fout(outputPath + ".EmpCycleDist");
    for (uint32_t i = 0; i != EmpCycleDist.size(); ++i) {
        fout << i << "\t" << misEmpCycleDist[i] << "\t" << EmpCycleDist[i]
             << "\t"
             << PHRED(
                        (double) (misEmpCycleDist[i] + 1)
                        / (EmpCycleDist[i] + 2)) << "\t" << CycleDist[i] << endl;
    }
    fout.close();
    return 0;
}

int StatCollector::getInsertSizeDist(const std::string &outputPath) {

    std::string InFileName(outputPath + ".InsertSizeTable");
    std::string OutFileName(outputPath + ".AdjustedInsertSizeDist");

    std::vector<double> f1;
    std::vector<double> f2;

    InsertSizeEstimator Estimator;
    Estimator.InputInsertSizeTable(InFileName, "FwdOnly");
    f1 = Estimator.UpdateWeight();
    Estimator.ReInit();
    Estimator.InputInsertSizeTable(InFileName, "RevOnly");
    f2 = Estimator.UpdateWeight();

    ofstream foutAdjust(OutFileName);
    for (int i = 0; i < f1.size(); ++i) {
        f1[i] = (f1[i] + f2[i]);
        foutAdjust << i << "\t" << f1[i] << endl;
    }
    foutAdjust.close();

    ofstream foutRaw(outputPath + ".RawInsertSizeDist");
    for (uint32_t i = 0; i != InsertSizeDist.size(); ++i) {
        foutRaw << i << "\t" << InsertSizeDist[i] << endl;
    }
    foutRaw.close();
    return 0;
}

int StatCollector::getSexChromInfo(const string &outputPath) {
    ofstream fout(outputPath + ".SexChromInfo");
    for (unordered_map<string, ContigStatus>::iterator iter =
            contigStatusTable.begin(); iter != contigStatusTable.end(); ++iter) {
        fout << iter->first << "\t" << iter->second.getNumOverlappedReads()
             << "\t" << iter->second.getNumFullyIncludedReads() << "\t"
             << iter->second.getNumPairOverlappedReads() << "\t"
             << iter->second.getNumFullyIncludedPairedReads() << endl;
    }
    fout.close();
    return 0;
}

int StatCollector::processCore(const string &statPrefix, const gap_opt_t *opt) {
    vector<int> PosNum(256, 0);
    //int total(0);
    for (auto i = PositionTable.begin();
         i != PositionTable.end(); ++i) //each chr
    {
        for (auto j =
                i->second.begin(); j != i->second.end(); ++j) //each site
        {
            /************DepthDist**************************************************************/
            {
                NumBaseMapped += DepthVec[j->second];
                if (DepthVec[j->second] > 1023)
                    DepthDist[1023]++;
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
    outputPileup(statPrefix, opt);
    SummaryOutput(statPrefix, opt);
    getGenoLikelihood(statPrefix);
    return 0;
}

int StatCollector::outputPileup(const string &outputPath, const gap_opt_t *opt) {
    ofstream fout(outputPath + ".Pileup");
    int qualoffset = 33;
    if (opt->mode | BWA_MODE_IL13) qualoffset = 64;
    for (sort_map::iterator i = VcfTable.begin(); i != VcfTable.end(); ++i) //each chr
    {
        for (std::map<int, unsigned int>::iterator j = i->second.begin();
             j != i->second.end(); ++j) //each site
        {
            if (SeqVec[j->second].size() == 0)
                continue;
            fout << i->first << "\t" << j->first << "\t.\t"
                 << StrandVec[j->second].size() << "\t";
            for (uint32_t k = 0; k != StrandVec[j->second].size(); ++k) {
                if (StrandVec[j->second][k])
                    fout << (char) toupper(SeqVec[j->second][k]);
                else
                    fout << (char) tolower(SeqVec[j->second][k]);
            }
            fout << "\t";
            for (uint32_t k = 0; k != QualVec[j->second].size(); ++k)
                fout << char(QualVec[j->second][k] + qualoffset);
            fout << "\t";
            for (uint32_t k = 0; k != MaqVec[j->second].size(); ++k) {
                fout << MaqVec[j->second][k];
            }
            fout << "\t";
            for (uint32_t k = 0; k != CycleVec[j->second].size(); ++k) {
                fout << CycleVec[j->second][k];
                if (k != CycleVec[j->second].size() - 1)
                    fout << ",";
            }
            fout << endl;
        }
    }
    fout.close();
    return 0;
}

vector<double> calLikelihood(const string &seq, const string &qual, const char &maj, const char &min)//maj:ref, min:alt
{
    double GL0(0), GL1(0), GL2(0);
    {
        for (uint32_t i = 0; i != seq.size(); ++i) {
            double seq_error = REV_PHRED(qual[i]);
            //fprintf(stderr,"Debug (qual:%d\tseq_error:%f)",qual[i],seq_error);
            if (seq[i] == maj) {
                GL0 += log10(1 - seq_error);
                GL1 += log10(0.5 - seq_error / 3);
                GL2 += log10(seq_error / 3);
            } else if (seq[i] == min) {
                GL0 += log10(seq_error / 3);
                GL1 += log10(0.5 - seq_error / 3);
                GL2 += log10(1 - seq_error);
            } else {
                GL0 += log10(2 * seq_error / 3);
                GL1 += log10(2 * seq_error / 3);
                GL2 += log10(2 * seq_error / 3);
            }
        }
        //fprintf(stderr, "\n");
    }
    vector<double> tmp(3, 0);
    tmp[0] = GL0 * (-10);
    tmp[1] = GL1 * (-10);
    tmp[2] = GL2 * (-10);
    double minimal(tmp[0]);
    if (tmp[1] < minimal) {
        minimal = tmp[1];
    }
    if (tmp[2] < minimal) {
        minimal = tmp[2];
    }
    tmp[0] -= minimal;
    tmp[1] -= minimal;
    tmp[2] -= minimal;
    return tmp;
}

int StatCollector::getGenoLikelihood(const string &outputPath) {
    ofstream fout(outputPath + ".Likelihood");
    size_t numAllele[4] =
            {0};
//    char majAllele, minAllele;
    int maxIndex = 0;
    //for (unsort_map::iterator i = PositionTable.begin();
    //	i != PositionTable.end(); ++i) //each chr
    //{
    //	for (std::unordered_map<int, unsigned int>::iterator j =
    //		i->second.begin(); j != i->second.end(); ++j) //each site
    //	{
    for (sort_map::iterator i = VcfTable.begin(); i != VcfTable.end(); ++i) {
        for (map<int, unsigned int>::iterator j = i->second.begin(); j != i->second.end(); ++j) {
            uint32_t markerIndex = j->second;
            VcfRecord *VcfLine = VcfRecVec[markerIndex];
            string RefStr(VcfLine->getRefStr()), AltStr(VcfLine->getAltStr());

            countAllele(numAllele, SeqVec[markerIndex]);
            maxIndex = findMaxAllele(numAllele, 4);
//            majAllele = "ACGT"[maxIndex];
            numAllele[maxIndex] = 0;
            maxIndex = findMaxAllele(numAllele, 4);
//            minAllele = "ACGT"[maxIndex];
            vector<double> tmpGL = calLikelihood(SeqVec[markerIndex],
                                                 QualVec[markerIndex],                //majAllele, minAllele);
                                                 RefStr[0], AltStr[0]);
            fout << i->first << "\t" << j->first << "\t"
                 << tmpGL[0] << "\t" << tmpGL[1] << "\t" << tmpGL[2];
            fout << endl;
        }
    }
    fout.close();
    return 0;
}

int StatCollector::addFSC(FileStatCollector a) {
    FSCVec.push_back(a);
    return 0;
}

int StatCollector::getGenomeSize(std::string RefPath) {
    ifstream fin(RefPath + ".fai");
    if (!fin.is_open()) {
        cerr << "Open file:" << RefPath + ".fai" << " failed!" << endl;
        exit(1);
    }
    string line, dummy, length;
    while (getline(fin, line)) {
        stringstream ss(line);
        ss >> dummy;
        ss >> length;
        ref_genome_size += atoi(length.c_str());
    }
    fin.close();
    return 0;
}

/*convert lambda function into regular functions*/
double StatCollector::Q20AvgDepth() {
    long long tmp(0);
    for (size_t i = 0; i != Q20DepthVec.size(); ++i) tmp += Q20DepthVec[i];
    return double(tmp) / NumPositionCovered;
}

double StatCollector::Q30AvgDepth() {
    long long tmp(0);
    for (size_t i = 0; i != Q30DepthVec.size(); ++i) tmp += Q30DepthVec[i];
    return double(tmp) / NumPositionCovered;
}

size_t StatCollector::MIS500() {
    long long tmp(0), total(0);
    for (size_t i = 500; i != InsertSizeDist.size(); ++i) total += InsertSizeDist[i];
    for (size_t i = 500; i != InsertSizeDist.size(); ++i) {
        tmp += InsertSizeDist[i];
        if (tmp > total / 2) return i;
    }
    return 0;
}

size_t StatCollector::MIS300() {
    long long tmp(0), total(0);
    for (size_t i = 300; i != InsertSizeDist.size(); ++i) total += InsertSizeDist[i];
    for (size_t i = 300; i != InsertSizeDist.size(); ++i) {
        tmp += InsertSizeDist[i];
        if (tmp > total / 2) return i;
    }
    return 0;
}

double StatCollector::Q20BaseFraction() {
    long long tmp(0);
    for (size_t i = 0; i != Q20DepthVec.size(); ++i) tmp += Q20DepthVec[i];
    return NumBaseMapped == 0 ? 0 : double(tmp) / NumBaseMapped;
}

double StatCollector::Q30BaseFraction() {
    long long tmp(0);
    for (size_t i = 0; i != Q30DepthVec.size(); ++i) tmp += Q30DepthVec[i];
    return NumBaseMapped == 0 ? 0 : double(tmp) / NumBaseMapped;
}

int StatCollector::SummaryOutput(const string &outputPath,const gap_opt_t *opt) {
    ofstream fout(outputPath + ".Summary");
    /*	int max_XorYmarker(0);
        if (opt->num_variant_short >= 100000)
        max_XorYmarker = 3000;
        else if (opt->num_variant_short >= 10000)
        max_XorYmarker = 300;
        else
        max_XorYmarker = 100;
        total_region_size = ((opt->flank_len-100) * 2 + 1) * opt->num_variant_short
        + ((opt->flank_long_len-100) * 2 + 1) * opt->num_variant_long
        + (501-100) * max_XorYmarker;*/
    long total_base(0);
    long total_reads(0);
    //output content for file HG00858.summary
    fout << "FILE 1|FILE 2|# Reads|Average Length" << endl;
    fout << "---|---|---|---" << endl;
    for (size_t i = 0; i != FSCVec.size(); ++i) {
        fout << FSCVec[i].FileName1 << "|" << FSCVec[i].FileName2 << "|"
             << FSCVec[i].NumRead << "|";
        fout << ((FSCVec[i].NumRead == 0) ? 0 : (FSCVec[i].NumBase / FSCVec[i].NumRead)) << endl;
        total_base += FSCVec[i].NumBase;
        total_reads += FSCVec[i].NumRead;
    }
    fout << "All" << "|-|" << total_reads << "|";
    fout << ((total_reads == 0) ? 0 : total_base / total_reads)
         << endl;
    fout << endl;
    fout << "Expected Read Depth : " << (double) total_base / ref_genome_size
         << " [" << total_base << "/" << ref_genome_size << "]" << endl;
    /*auto AvgDepth =
        [&]()->double
    {	long long tmp(0); for (size_t i = 0; i != DepthDist.size(); ++i) tmp += i*DepthDist[i]; return double(tmp) / total_region_size; };*/
    fout << "Estimated AvgDepth : ";
    fout << ((NumPositionCovered == 0) ? 0 : NumBaseMapped / (double) NumPositionCovered) << endl;
    fout << "Estimated percentage of accessible genome covered : "
         << (1 - (double) DepthDist[0] / total_region_size) * 100 << "/100"
         << endl;
    /*auto Q20AvgDepth =
        [&]()->double
    {	long long tmp(0); for (size_t i = 0; i != Q20DepthVec.size(); ++i) tmp += Q20DepthVec[i]; return double(tmp) / total_region_size; };*/
    fout << "Estimated AvgDepth for Q20 bases : " << Q20AvgDepth() << endl;
    /*auto Q30AvgDepth =
        [&]()->double
    {	long long tmp(0); for (size_t i = 0; i != Q30DepthVec.size(); ++i) tmp += Q30DepthVec[i]; return double(tmp) / total_region_size; };*/
    fout << "Estimated AvgDepth for Q30 bases : " << Q30AvgDepth() << endl;
    /*auto MIS500 =
        [&]()->double
    {	long long tmp(0), total(0);
        for (size_t i = 500; i != InsertSizeDist.size(); ++i) total += InsertSizeDist[i];
        for (size_t i = 500; i != InsertSizeDist.size(); ++i)
        {
            tmp += InsertSizeDist[i]; if (tmp > total / 2) return i;
        }
        return 0;
    };*/
    fout << "Median Insert Size(>=500bp) : " << MIS500() << endl;
    /*auto MIS300 =
        [&]()->double
    {	long long tmp(0), total(0);
        for (size_t i = 300; i != InsertSizeDist.size(); ++i) total += InsertSizeDist[i];
        for (size_t i = 300; i != InsertSizeDist.size(); ++i)
        {
        tmp += InsertSizeDist[i]; if (tmp > total / 2) return i;
        }
        return 0;
    };*/
    fout << "Median Insert Size(>=300bp) : " << MIS300() << endl;
    //output for fraction figure
    /*auto Q20BaseFraction = [&]()->double
    {	long long tmp(0); for (size_t i = 0; i != Q20DepthVec.size(); ++i) tmp += Q20DepthVec[i]; return NumBaseMapped==0?0:double(tmp) / NumBaseMapped; };
    auto Q30BaseFraction = [&]()->double
    {	long long tmp(0); for (size_t i = 0; i != Q30DepthVec.size(); ++i) tmp += Q30DepthVec[i]; return NumBaseMapped == 0 ? 0 : double(tmp) / NumBaseMapped; };*/
    double DP1fraction = NumPositionCovered / (double) total_region_size;
    double DP2fraction = NumPositionCovered2 / (double) total_region_size;
    double DP5fraction = NumPositionCovered5 / (double) total_region_size;
    double DP10fraction = NumPositionCovered10 / (double) total_region_size;
    fout << "Q20 Base Fraction :" << Q20BaseFraction() << endl;
    fout << "Q30 Base Fraction :" << Q30BaseFraction() << endl;
    fout << "Depth 1 or above position fraction :" << DP1fraction << endl;
    fout << "Depth 2 or above position fraction :" << DP2fraction << endl;
    fout << "Depth 5 or above position fraction :" << DP5fraction << endl;
    fout << "Depth 10 or above position fraction :" << DP10fraction << endl;

    return 0;
}

StatCollector::~StatCollector() {
    releaseVcfSites();
}


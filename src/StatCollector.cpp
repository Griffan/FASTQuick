#include "StatCollector.h"

#include "../libbwa/bwase.h"
#include "../libbwa/bwtaln.h"

#include "../misc/bam/SamFile.h"
#include "../misc/bam/SamRecord.h"
#include "../misc/bam/SamValidation.h"
#include "../misc/vcf/VcfFileReader.h"

#include "InsertSizeEstimator.h"
#include <algorithm>
#include <iostream>

using namespace std;

// extern string Prefix;
extern void notice(const char *, ...);

extern void warning(const char *, ...);

extern void error(const char *, ...);

const uint16_t __cigar_table[16] = {
    2 /*D*/, 0, 0, 0, 0, 1 /*I*/, 0, 0, 0, 0 /*M*/, 0, 0, 0, 0, 0, 3 /*S*/
};

const char sign[2] = {1, -1};

inline int StatCollector::IsPartialAlign(const bwa_seq_t *q) {
  for (int k = 0; k < q->n_cigar; ++k) {

    // int cl = __cigar_len(q->cigar[k]);
    int cop = "MIDS"[__cigar_op(q->cigar[k])];
    switch (cop) {
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

inline int StatCollector::IsPartialAlign(SamRecord &q) {
  if (std::string(q.getCigar()).find('S') != std::string::npos)
    return 1;
  else
    return 0;
}

#define __cigar_convert(__cigar) ((__cigar_table[__cigar - 68]))

static std::string Cigar2String(int n_cigar, const bwa_cigar_t *cigar,
                                int len) {
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

static std::vector<bwa_cigar_t> String2Cigar(const std::string &cigarString) {

  std::vector<bwa_cigar_t> cigar;
  int len = static_cast<int>(cigarString.size());
  int digitNow(0), digitLast(0);
  while (digitNow != len) {
    while (isdigit(cigarString[digitNow]))
      digitNow++;
    uint16_t cl = atoi(cigarString.substr(digitLast, digitNow - digitLast)
                           .c_str()); // digitNow is "SMID"
    digitLast = digitNow + 1;
    char cop = cigarString[digitNow];
    digitNow++;
    cigar.push_back(__cigar_create(__cigar_convert(cop), cl));
  }
  return cigar;
}

std::string
StatCollector::RecoverRefseqByMDandCigar(const std::string &readSeq,
                                         std::string MD,
                                         const std::string &cigarString) {

  std::vector<bwa_cigar_t> cigar = String2Cigar(cigarString);
  return RecoverRefseqByMDandCigar(readSeq, MD, &cigar[0], cigar.size());
}

std::string StatCollector::RecoverRefseqByMDandCigar(const std::string &readSeq,
                                                     std::string MD,
                                                     const bwa_cigar_t *cigar,
                                                     int n_cigar) {
  std::transform(MD.begin(), MD.end(), MD.begin(), ::toupper);
  if (MD.find_first_of("ATCGN") == std::string::npos &&
      atoi(MD.c_str()) == readSeq.size()) // means perfect match
    return readSeq;

  // only collect M pieces
  std::string refSeq;
  if (n_cigar > 0) {
    int relativeCoordOnRead = 0;
    for (int k = 0; k < n_cigar; ++k) {
      int cl = __cigar_len(cigar[k]); // cigar length
      int cop = "MIDS"[__cigar_op(cigar[k])];
      switch (cop) {
      case 'M':
        refSeq += readSeq.substr(relativeCoordOnRead, cl);
        relativeCoordOnRead += cl;
        break;
      case 'S': // ignore soft clip
        relativeCoordOnRead += cl;
        break;
      case 'D': // ref has base, read has no base;
        break;
      case 'I': // ref has no base, read has base
        relativeCoordOnRead += cl;
        break;
      default:
        warning("Unhandled cigar_op %d:%d\n", cop, cl);
      }
    } // end for
  } else {
    refSeq = readSeq;
  }

  // begin of MD modify
  int last(0), total_len(0) /*pos on seq*/;
  // remember we are making reference sequence
  for (uint32_t i = 0; i != MD.size(); ++i)
    if (isdigit(MD[i]))
      continue;
    else if (MD[i] == '^') // sign for deletion
    {
      // len from last to current deletion, '^' not included
      int len = atoi(MD.substr(last, i - last).c_str());
      total_len += len;
      int start_on_read = total_len; // 1 based
      i++;
      std::string tmp;
      while (!isdigit(MD[i])) // deletion sequence
      {
        tmp += MD[i];
        i++;
        total_len++;
      }
      std::string left = refSeq.substr(0, start_on_read);
      std::string right =
          refSeq.substr(start_on_read, refSeq.length() - start_on_read + 1);
      refSeq = left + tmp + right;
      last = i;
    } else { // SNP
      // len from last to current mismatch, include mismatch
      int len = atoi(MD.substr(last, i - last).c_str()) + 1;
      total_len += len;
      refSeq[total_len - 1] = MD[i];
      last = i + 1;
    }

  //  // assembly seq pieces
  //  if (n_cigar > 0) {
  //    int relativeCoordOnBackbone = 0;
  //    int relativeCoordOnRead = 0;
  //    std::string backboneSeq = refSeq;
  //    refSeq = "";
  //    for (int k = 0; k < n_cigar; ++k) {
  //      int cl = __cigar_len(cigar[k]); // cigar length
  //      int cop = "MIDS"[__cigar_op(cigar[k])];
  //      switch (cop) {
  //      case 'M':
  //        refSeq += backboneSeq.substr(relativeCoordOnBackbone, cl);
  //        relativeCoordOnRead += cl;
  //        relativeCoordOnBackbone += cl;
  //        break;
  //      case 'S':     // ignore soft clip
  //        if (k == 0) // left clip
  //          refSeq = readSeq.substr(relativeCoordOnRead, cl);
  //        else // right
  //          refSeq += readSeq.substr(relativeCoordOnRead, cl);
  //        relativeCoordOnRead += cl;
  //        break;
  //      case 'D': // ref has base, read has no base;
  //        refSeq += backboneSeq.substr(relativeCoordOnBackbone, cl);
  //        relativeCoordOnBackbone += cl;
  //        break;
  //      case 'I': // ref has no base, read has base
  //        relativeCoordOnRead += cl;
  //        break;
  //      default:
  //        warning("Unhandled cigar_op %d:%d\n", cop, cl);
  //      }
  //    } // end for
  //  }
  return refSeq;
}

static int FindMaxAllele(const size_t *a, size_t len) {
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

static int CountAllele(size_t *a, const string &seq) {
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
#define REV_PHRED(x) pow(10.0, (x / (-10.0)))
#define PHRED(x) (-10) * log10(x)

StatCollector::StatCollector() {
  // cerr << "NOTE:Using default initializer..." << endl;
  PositionTable.clear();
  VcfRecVec.clear();
  duplicateTable.clear();
  VariantProxyTable.clear();
  GC.clear();
  index = 0;
  total_region_size = 0;
  ref_genome_size = 0;
  NumXorYMarker = 0;
  NumShortMarker = 0;
  NumLongMarker = 0;
  NumPCRDup = 0;
  NumBaseMapped = 0;
  NumPositionCovered = 0;              // position with depth larger than 0
  NumPositionCovered2 = 0;             // larger than 1
  NumPositionCovered5 = 0;             // larger than 4
  NumPositionCovered10 = 0;            // larger than 9
  DepthDist = vector<size_t>(1024, 0); // up to depth 1024X
  CycleDist = vector<size_t>(512, 0);
  GCDist = vector<size_t>(101, 0);
  PosNum = vector<size_t>(101, 0);
  GCDist = vector<size_t>(256, 0);
  EmpRepDist = vector<size_t>(256, 0);
  misEmpRepDist = vector<size_t>(256, 0);
  EmpCycleDist = vector<size_t>(256, 0);
  misEmpCycleDist = vector<size_t>(256, 0);
  InsertSizeDist = vector<size_t>(INSERT_SIZE_LIMIT, 0);
  //    MaxInsertSizeDist = vector<size_t>(INSERT_SIZE_LIMIT, 0);
}

StatCollector::StatCollector(const std::string &OutFile) {
  PositionTable.clear();
  VcfRecVec.clear();
  duplicateTable.clear();
  VariantProxyTable.clear();
  GC.clear();
  index = 0;
  total_region_size = 0;
  ref_genome_size = 0;
  NumXorYMarker = 0;
  NumShortMarker = 0;
  NumLongMarker = 0;
  NumPCRDup = 0;
  NumBaseMapped = 0;
  NumPositionCovered = 0;              // position with depth larger than 0
  NumPositionCovered2 = 0;             // larger than 1
  NumPositionCovered5 = 0;             // larger than 4
  NumPositionCovered10 = 0;            // larger than 9
  DepthDist = vector<size_t>(1024, 0); // up to depth 1024X
  CycleDist = vector<size_t>(512, 0);
  GCDist = vector<size_t>(101, 0);
  PosNum = vector<size_t>(101, 0);
  EmpRepDist = vector<size_t>(256, 0);
  misEmpRepDist = vector<size_t>(256, 0);
  EmpCycleDist = vector<size_t>(256, 0);
  misEmpCycleDist = vector<size_t>(256, 0);
  InsertSizeDist = vector<size_t>(INSERT_SIZE_LIMIT, 0);
  //    MaxInsertSizeDist = vector<size_t>(INSERT_SIZE_LIMIT, 0);
}

void StatCollector::StatVecDistUpdate(const std::string &chrom, int i,
                                      int tmpCycle, char readBase,
                                      char baseQual, char refBase) {
  EmpRepDist[baseQual]++;
  EmpCycleDist[tmpCycle]++;
  if (readBase != 'N' and refBase != readBase and refBase != 'N' and
      dbSNPTable[chrom].find(i) == dbSNPTable[chrom].end()) {
    //    std::cerr << "pos\t" << chrom << "\t" << i << "\t" << refBase << "\t"
    //              << readBase << "\t" << baseQual << std::endl;
    misEmpRepDist[baseQual]++;
    misEmpCycleDist[tmpCycle]++;
  }
}

void StatCollector::AddBaseInfoToNewCoord(const std::string &chrom, int i,
                                          int tmpCycle, char readBase,
                                          char baseQual, char refBase) {
  DepthVec.push_back(0);
  Q20DepthVec.push_back(0);
  Q30DepthVec.push_back(0);
  DepthVec[index]++;
  if (baseQual >= 20) {
    Q20DepthVec[index]++;
    if (baseQual >= 30) {
      Q30DepthVec[index]++;
    }
  }
  StatVecDistUpdate(chrom, i, tmpCycle, readBase, baseQual, refBase);
  PositionTable[chrom][i] = index++;
}

void StatCollector::UpdateInfoVecAtMarker(
    int tmpCycle, int absoluteSite, int cl, bool strand, const string &chrom,
    const string &seq, const string &qual, u_char mapQ,
    int relativeCoordOnRead) // only update info at marker site, for pileup
                             // output
{
  if (VcfTable.find(chrom) != VcfTable.end()) {
    for (int i = absoluteSite; i != absoluteSite + cl - 1 + 1;
         ++i, tmpCycle += 1 * sign[strand], ++relativeCoordOnRead) {
      if (VcfTable[chrom].find(i) !=
          VcfTable[chrom].end()) // actual marker site
      {
        int tmpIndex = VcfTable[chrom][i];
        SeqVec[tmpIndex] += seq[relativeCoordOnRead];
        QualVec[tmpIndex] += qual[relativeCoordOnRead];
        CycleVec[tmpIndex].push_back(tmpCycle);
        MaqVec[tmpIndex].push_back(mapQ + 33);
        StrandVec[tmpIndex].push_back(strand);
      }
    }
  }
}

int StatCollector::AddMatchBaseInfo(const string &seq, const string &qual,
                                    const string &refSeq, const string &chr,
                                    int readRealStart, bool strand, u_char mapQ,
                                    int matchLen, int tmpCycle,
                                    int relativeCoordOnRead,
                                    int relativeCoordOnRef) {
  std::string Chrom(chr);
  std::transform(Chrom.begin(), Chrom.end(), Chrom.begin(), ::toupper);
  if (Chrom.find("CHR") != std::string::npos) {
    Chrom = Chrom.substr(3); // strip chr
  }

  UpdateInfoVecAtMarker(tmpCycle, readRealStart, matchLen, strand, Chrom, seq,
                        qual, mapQ, relativeCoordOnRead);

  /*****************************************************************************/
  return UpdateInfoVecAtRegularSite(seq, qual, refSeq, Chrom, readRealStart,
                                    strand, matchLen, tmpCycle,
                                    relativeCoordOnRead, relativeCoordOnRef);
}
#define FLANK_EDGE 0.65f

int StatCollector::UpdateInfoVecAtRegularSite(
    const string &seq, const string &qual, const string &refSeq,
    const string &chr, int readRealStart, bool strand, int matchLen,
    int tmpCycle, int relativeCoordOnRead, int relativeCoordOnRef) {
  char refBase = 0;
  char readBase = 0;
  char baseQual = 0;
  int total_effective_len = 0;
  if (PositionTable.find(chr) != PositionTable.end()) // chrom exists
    for (int i = readRealStart; i != readRealStart + matchLen - 1 + 1; ++i,
             tmpCycle += 1 * sign[strand], ++relativeCoordOnRead,
             ++relativeCoordOnRef) {

      if (not flankRegion.IsOverlapped(chr, i))
        continue;

      refBase = refSeq[relativeCoordOnRef];
      readBase = seq[relativeCoordOnRead];
      baseQual = qual[relativeCoordOnRead];

      if (PositionTable[chr].find(i) !=
          PositionTable[chr].end()) { // site already exists
        int tmpIndex = PositionTable[chr][i];
        DepthVec[tmpIndex]++;
        if (baseQual >= 20) {
          Q20DepthVec[tmpIndex]++;
          if (baseQual >= 30) {
            Q30DepthVec[tmpIndex]++;
          }
        }
        total_effective_len++;
        StatVecDistUpdate(chr, i, tmpCycle, readBase, baseQual, refBase);
      } else // coord not exists
      {
        total_effective_len++;
        AddBaseInfoToNewCoord(chr, i, tmpCycle, readBase, baseQual, refBase);
      }
    }
  else // chrom not exists
  {
    for (int i = readRealStart; i != readRealStart + matchLen - 1 + 1; ++i,
             tmpCycle += 1 * sign[strand], ++relativeCoordOnRead,
             ++relativeCoordOnRef) {

      if (not flankRegion.IsOverlapped(chr, i))
        continue;

      refBase = refSeq[relativeCoordOnRef];
      readBase = seq[relativeCoordOnRead];
      baseQual = qual[relativeCoordOnRead];
      total_effective_len++;
      AddBaseInfoToNewCoord(chr, i, tmpCycle, readBase, baseQual, refBase);
    }
  }
  return total_effective_len;
}

bool StatCollector::AddSingleAlignment(const bntseq_t *bns, bwa_seq_t *p,
                                       const gap_opt_t *opt) {
  // DEBUG related variables
  int total_effective_len(0);
  // DEBUG related variables end
  int threshQual = 20;
  int seqid(0);
  int j(0);

  if (p->type == BWA_TYPE_NO_MATCH or p->mapQ < threshQual) {
    return false;
  }

  j = static_cast<int>(pos_end(p) - p->pos); // length of read
  bns_coor_pac2real(bns, p->pos, j, &seqid);

  string seq, qual;

  if (p->strand == 0)
    for (j = 0; j != p->full_len; ++j) {
      seq += "ACGTN"[(int)(p)->seq[j]];
      qual += (char)(p->qual[j] - 33); // change from ASCII to Phred
    }
  else
    for (j = 0; j != p->full_len; ++j) {
      seq += "TGCAN"[(int)(p)->seq[p->full_len - 1 - j]];
      qual += (char)(p->qual[p->full_len - 1 - j] - 33); // 0 based
    }

  string chrName = string(bns->anns[seqid].name);

  int pos = static_cast<int>(p->pos - bns->anns[seqid].offset + 1);

  int readRealStart(0);
  int refRealStart(0), refRealEnd(0);
  size_t atPos;
  string chrom;
  char *pEnd;
  int refCoord; // absolute coordinate of variant
  size_t colonPos = chrName.find(':');
  if (colonPos != std::string::npos) // FASTQuick own alignment
  {
    atPos = chrName.find('@');
    chrom = chrName.substr(0, colonPos);
    refCoord = static_cast<int>(strtol(
        chrName.substr(colonPos + 1, atPos - colonPos + 1).c_str(), &pEnd,
        10)); // absolute coordinate of variant
    if (chrName[chrName.size() - 1] == 'L') {
      // absolute coordinate of refSeq start
      refRealStart = refCoord - opt->flank_long_len;
      refRealEnd = refCoord + opt->flank_long_len;
      // absolute coordinate of current reads on reference
      readRealStart = refRealStart + pos - 1;
    } else {
      refRealStart = refCoord - opt->flank_len;
      refRealEnd = refCoord + opt->flank_len;
      readRealStart = refRealStart + pos - 1;
    }
  } else // External alignment
  {
    readRealStart = pos;
    int readRealEnd = readRealStart + p->len;
    chrom = chrName;
    int overlap_right(0), overlap_left(0);
    int left_flank_len(0), right_flank_len(0) /*flank length of right element*/;
    // need to locate which artificial ref seq this reads aligned to
    if (VcfTable.find(chrom) != VcfTable.end()) {
      auto right = VcfTable[chrom].lower_bound(readRealStart);
      if (right != VcfTable[chrom].begin() and right != VcfTable[chrom].end())
      // meaning both right and left elements are valid
      {
        auto left(--right);
        if (std::string(VcfRecVec[right->second]->getIDStr()).back() ==
            'L') // long region
          right_flank_len = opt->flank_long_len;
        else
          right_flank_len = opt->flank_len;

        if (std::string(VcfRecVec[left->second]->getIDStr()).back() ==
            'L') // long region
          left_flank_len = opt->flank_long_len;
        else
          left_flank_len = opt->flank_len;

        overlap_right = readRealEnd - (right->first - right_flank_len) + 1;
        overlap_left = left->first + left_flank_len - readRealStart + 1;

        if (overlap_left > 0 and overlap_right > 0) { // both overlap exist
          if (overlap_left > overlap_right) {
            refRealStart = left->first - left_flank_len;
            refRealEnd = left->first + left_flank_len;
          } else {
            refRealStart = right->first - right_flank_len;
            refRealEnd = right->first + right_flank_len;
          }
        } else if (overlap_left > 0) { // only left overlap exists
          refRealStart = left->first - left_flank_len;
          refRealEnd = left->first + left_flank_len;

        } else if (overlap_right > 0) { // only right overlap exists
          refRealStart = right->first - right_flank_len;
          refRealEnd = right->first + right_flank_len;
        } else // not overlapped adjacent elements
          return false;
      } else if (right != VcfTable[chrom].begin()) // left element is valid, but
                                                   // not right element isn't
      {
        auto left(--right);
        if (std::string(VcfRecVec[left->second]->getIDStr()).back() ==
            'L') // long region
          left_flank_len = opt->flank_long_len;
        else
          left_flank_len = opt->flank_len;
        overlap_left = left->first + left_flank_len - readRealStart + 1;
        if (overlap_left > 0) { // left overlap exists
          refRealStart = left->first - left_flank_len;
          refRealEnd = left->first + left_flank_len;

        } else // not overlapped adjacent elements
          return false;
      } else if (right != VcfTable[chrom]
                              .end()) // only this element exists, which is rare
      {
        if (std::string(VcfRecVec[right->second]->getIDStr()).back() ==
            'L') // long region
          right_flank_len = opt->flank_long_len;
        else
          right_flank_len = opt->flank_len;
        overlap_right = readRealEnd - (right->first - right_flank_len) + 1;
        if (overlap_right > 0) { // right overlap exists
          refRealStart = right->first - right_flank_len;
          refRealEnd = right->first + right_flank_len;
        } else
          return false;
      } else // meaning empty set, impossible
      {
        warning(
            "read:%s cannot find adjacent variants on chrom:%s, impossible!",
            p->name, chrom.c_str());
        exit(EXIT_FAILURE);
      }
    } else
      return false;
  }

  string MD(p->md);
  string refSeq = RecoverRefseqByMDandCigar(seq, MD, p->cigar, p->n_cigar);

  int absoluteSite = readRealStart;
  int tmpCycle(0) /*sequencing cycle*/,
      relativeCoordOnRead(0) /*relative coord on reads*/, relativeCoordOnRef(0);
  if (p->strand != 0) {
    tmpCycle = p->full_len - 1;
  }

  if (p->cigar) {
    int n_cigar = p->n_cigar;
    bwa_cigar_t *cigar = p->cigar;
    for (int k = 0; k < n_cigar; ++k) {
      int cl = __cigar_len(cigar[k]); // cigar length
      char cop = "MIDS"[__cigar_op(cigar[k])];
      switch (cop) {
      case 'M':
        /*matching sites*/
        total_effective_len += AddMatchBaseInfo(
            seq, qual, refSeq, chrom, absoluteSite, p->strand, p->mapQ, cl,
            tmpCycle, relativeCoordOnRead, relativeCoordOnRef);
        absoluteSite += cl;
        tmpCycle += cl * sign[p->strand];
        relativeCoordOnRead += cl;
        relativeCoordOnRef += cl;
        break;
      case 'S': // ignore soft clip
        tmpCycle += cl * sign[p->strand];
        relativeCoordOnRead += cl;
        //        relativeCoordOnRef += cl;
        break;
      case 'D': // ref has base, read has no base;
        absoluteSite += cl;
        relativeCoordOnRef += cl;
        break;
      case 'I': // ref has no base, read has base;
        tmpCycle += cl * sign[p->strand];
        relativeCoordOnRead += cl;
        break;
      default:
        warning("Unhandled cigar_op %d:%d\n", cop, cl);
      }
    } // end for
  } else {
    /*matching sites*/
    total_effective_len += AddMatchBaseInfo(
        seq, qual, refSeq, chrom, absoluteSite, p->strand, p->mapQ, p->len,
        tmpCycle, relativeCoordOnRead, relativeCoordOnRef);
  }
  //    std::cerr<<p->name<<" effective len:"<<total_effective_len<<std::endl;
  return true;
}
// in this setting p is fastq1 and q is fastq2
int StatCollector::IsDuplicated(const bntseq_t *bns, const bwa_seq_t *p,
                                const bwa_seq_t *q, AlignStatus type,
                                ofstream &fout) {

  int maxInsert(-1), maxInsert2(-1);
  int readLength(0), seqid_p(-1), seqid_q(-1);
  int flag1 = 0;
  int flag2 = 0;
  int threshQual = 0;
  string status;

  int cl1(0), op1(0); // p left
  int cl2(0), op2(0); // p right
  int cl3(0), op3(0); // q left
  int cl4(0), op4(0); // q right

  if (p) {
    flag1 = p->extra_flag;
    if (p->type == BWA_TYPE_NO_MATCH) {
      flag1 |= SAM_FSU;
    }
    if (p->strand)
      flag1 |= SAM_FSR;
  }
  if (q) {
    flag2 = q->extra_flag;
    if (q->type == BWA_TYPE_NO_MATCH) {
      flag2 |= SAM_FSU;
    }
    if (q->strand)
      flag2 |= SAM_FSR;
  }

  if (type == AlignStatus::Second) // q is aligned only, read2
  {
    // length of read including cigar for bns tracking, length of
    // the reference in the alignment
    readLength = static_cast<int>(pos_end(q) - q->pos);
    bns_coor_pac2real(bns, q->pos, readLength, &seqid_q);

    if (q->mapQ > threshQual) {
      status = "RevOnly";
      if (q->cigar) {
        op3 = "MIDS"[__cigar_op(q->cigar[0])]; // left most cigar
        if (op3 == 'S') {
          cl3 = __cigar_len(q->cigar[0]);
        }

        op4 = "MIDS"[__cigar_op(q->cigar[q->n_cigar - 1])]; // right most cigar
        if (op4 == 'S') {
          cl4 = __cigar_len(q->cigar[q->n_cigar - 1]);
        }
      }

      if (q->strand) // reverse
      {
        if (bns->anns[seqid_q].offset + bns->anns[seqid_q].len >=
            (q->pos - cl3) + q->len) // soft clip stays within Flank Region
          maxInsert2 = static_cast<int>((q->pos - cl3) + q->len -
                                        bns->anns[seqid_q].offset);
        else // has soft clip and not complete mapped
          return 2;
      } else // forward
      {
        if ((q->pos - cl3) >=
            bns->anns[seqid_q].offset) // soft clip stays within Flank Region
          maxInsert = static_cast<int>(bns->anns[seqid_q].offset +
                                       bns->anns[seqid_q].len - (q->pos - cl3));
        else
          return 2;
        status = "FwdOnly";
      }
      fout << q->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
           << "\t"
           << "*"
           << "\t"
           << "*"
           << "\t" << flag1 << "\t" << 0 << "\t"
           << "*"
           << "\t" << bns->anns[seqid_q].name << "\t"
           << q->pos - bns->anns[seqid_q].offset + 1 << "\t" << flag2 << "\t"
           << q->len << "\t" << Cigar2String(q->n_cigar, q->cigar, q->len)
           << "\t" << status << endl;
      return 0;
    } else {
      fout << q->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
           << "\t"
           << "*"
           << "\t"
           << "*"
           << "\t" << flag1 << "\t" << 0 << "\t"
           << "*"
           << "\t" << bns->anns[seqid_q].name << "\t"
           << q->pos - bns->anns[seqid_q].offset + 1 << "\t" << flag2 << "\t"
           << q->len << "\t" << Cigar2String(q->n_cigar, q->cigar, q->len)
           << "\t"
           << "LowQual" << endl;
      return 2; // low quality
    }
  } else if (type == AlignStatus::First) // p is aligned only, read1
  {
    readLength = static_cast<int>(
        pos_end(p) -
        p->pos); // length of read including cigar for bns tracking, length of
                 // the reference in the alignment
    bns_coor_pac2real(bns, p->pos, readLength, &seqid_p);

    if (p->mapQ > threshQual) {
      status = "FwdOnly";

      if (p->cigar) {
        op1 = "MIDS"[__cigar_op(p->cigar[0])]; // left most cigar
        if (op1 == 'S') {
          cl1 = __cigar_len(p->cigar[0]);
        }

        op2 = "MIDS"[__cigar_op(p->cigar[p->n_cigar - 1])]; // right most cigar
        if (op2 == 'S') {
          cl2 = __cigar_len(p->cigar[p->n_cigar - 1]);
        }
      }

      if (p->strand) // reverse
      {
        if (bns->anns[seqid_p].offset + bns->anns[seqid_p].len >=
            (p->pos - cl1) + p->len) // soft clip stays within Flank Region
          maxInsert2 = static_cast<int>((p->pos - cl1) + p->len -
                                        bns->anns[seqid_p].offset);
        else
          return 2;

        status = "RevOnly";
      } else {
        if ((p->pos - cl1) >=
            bns->anns[seqid_p].offset) // soft clip stays within Flank Region
          maxInsert = static_cast<int>(bns->anns[seqid_p].offset +
                                       bns->anns[seqid_p].len - (p->pos - cl1));
        else
          return 2;
      }
      fout << p->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
           << "\t" << bns->anns[seqid_p].name << "\t"
           << p->pos - bns->anns[seqid_p].offset + 1 << "\t" << flag1 << "\t"
           << p->len << "\t" << Cigar2String(p->n_cigar, p->cigar, p->len)
           << "\t"
           << "*"
           << "\t"
           << "*"
           << "\t" << flag2 << "\t" << 0 << "\t"
           << "*"
           << "\t" << status << endl;
      return 0;
    } else {
      fout << p->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
           << "\t" << bns->anns[seqid_p].name << "\t"
           << p->pos - bns->anns[seqid_p].offset + 1 << "\t" << flag1 << "\t"
           << p->len << "\t" << Cigar2String(p->n_cigar, p->cigar, p->len)
           << "\t"
           << "*"
           << "\t"
           << "*"
           << "\t" << flag2 << "\t" << 0 << "\t"
           << "*"
           << "\t"
           << "LowQual" << endl;
      return 2;
    }
  } else if (type == AlignStatus::Both) // both aligned
  {
    readLength = static_cast<int>(pos_end(p) - p->pos); // length of read
    bns_coor_pac2real(bns, p->pos, readLength, &seqid_p);
    readLength = static_cast<int>(pos_end(q) - q->pos); // length of read
    bns_coor_pac2real(bns, q->pos, readLength, &seqid_q);
    /*deal with cigar*/

    if (p->cigar) {
      op1 = "MIDS"[__cigar_op(p->cigar[0])]; // left most cigar
      if (op1 == 'S') {
        cl1 = __cigar_len(p->cigar[0]);
      }

      op2 = "MIDS"[__cigar_op(p->cigar[p->n_cigar - 1])]; // right most cigar
      if (op2 == 'S') {
        cl2 = __cigar_len(p->cigar[p->n_cigar - 1]);
      }
    }

    if (q->cigar) {
      op3 = "MIDS"[__cigar_op(q->cigar[0])]; // left most cigar
      if (op3 == 'S') {
        cl3 = __cigar_len(q->cigar[0]);
      }

      op4 = "MIDS"[__cigar_op(q->cigar[q->n_cigar - 1])]; // right most cigar
      if (op4 == 'S') {
        cl4 = __cigar_len(q->cigar[q->n_cigar - 1]);
      }
    }
    /*end deal with cigar*/
    if (!(p->strand) && q->strand && p->pos < q->pos) // FR
    {
      if ((p->pos - cl1) >=
          bns->anns[seqid_p].offset) // soft clip stays within Flank Region
        maxInsert = static_cast<int>(bns->anns[seqid_p].offset +
                                     bns->anns[seqid_p].len - (p->pos - cl1));
      else
        maxInsert = -1;

      if (bns->anns[seqid_q].offset + bns->anns[seqid_q].len >=
          (q->pos - cl3) + q->len) // soft clip stays within Flank Region
        maxInsert2 = static_cast<int>((q->pos - cl3) + q->len -
                                      bns->anns[seqid_q].offset);
      else
        maxInsert2 = -1;

    } else if (!(q->strand) && p->strand && q->pos < p->pos) // FR but rotated
    {

      if ((q->pos - cl3) >=
          bns->anns[seqid_q].offset) // soft clip stays within Flank Region
        maxInsert = static_cast<int>(bns->anns[seqid_q].offset +
                                     bns->anns[seqid_q].len - (q->pos - cl3));
      else
        maxInsert = -1;

      if (bns->anns[seqid_p].offset + bns->anns[seqid_p].len >=
          (p->pos - cl1) + p->len) // soft clip stays within Flank Region
        maxInsert2 = static_cast<int>((p->pos - cl1) + p->len -
                                      bns->anns[seqid_p].offset);
      else
        maxInsert2 = -1;
    } else // other than FR, treat as DiffChrom
    {
      fout << p->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
           << "\t" << bns->anns[seqid_p].name << "\t"
           << p->pos - bns->anns[seqid_p].offset + 1 << "\t" << flag1 << "\t"
           << p->len << "\t" << Cigar2String(p->n_cigar, p->cigar, p->len)
           << "\t" << bns->anns[seqid_q].name << "\t"
           << q->pos - bns->anns[seqid_q].offset + 1 << "\t" << flag2 << "\t"
           << q->len << "\t" << Cigar2String(q->n_cigar, q->cigar, q->len)
           << "\tNotPair" << endl;
      return 0;
    }
  }
  // from now on only type 2 paired reads remained
  if (maxInsert >= INSERT_SIZE_LIMIT)
    maxInsert = INSERT_SIZE_LIMIT - 1;
  if (maxInsert2 >= INSERT_SIZE_LIMIT)
    maxInsert2 = INSERT_SIZE_LIMIT - 1;

  if (seqid_p != seqid_q && seqid_p != -1 && seqid_q != -1) {
    // int ActualInsert(-1);
    //	ActualInsert=0;
    InsertSizeDist[0]++;
    // cerr<<"Duplicate function exit from diff chrom"<<endl;
    fout << p->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
         << "\t" << bns->anns[seqid_p].name << "\t"
         << p->pos - bns->anns[seqid_p].offset + 1 << "\t" << flag1 << "\t"
         << p->len << "\t" << Cigar2String(p->n_cigar, p->cigar, p->len) << "\t"
         << bns->anns[seqid_q].name << "\t"
         << q->pos - bns->anns[seqid_q].offset + 1 << "\t" << flag2 << "\t"
         << q->len << "\t" << Cigar2String(q->n_cigar, q->cigar, q->len)
         << "\tNotPair" << endl;
    return 0;
  }

  if (p->mapQ > threshQual && q->mapQ > threshQual) {

    int ActualInsert(-1), start(0), end(0); //[start, end)
    status = "PartialPair";
    if (!(p->strand) && q->strand && p->pos < q->pos) // FR
    {
      start = p->pos - cl1;
      end = q->pos - cl3 + q->len;
      ActualInsert = end - start;
    } else if (!(q->strand) && p->strand && q->pos < p->pos) // FR but rotated
    {
      start = q->pos - cl3;
      end = p->pos - cl1 + p->len;
      ActualInsert = end - start;
    }

    if (maxInsert != -1 and maxInsert2 != -1)
      status = "PropPair";

    InsertSizeDist[ActualInsert]++;
    fout << p->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t"
         << ActualInsert << "\t" << bns->anns[seqid_p].name << "\t"
         << p->pos - bns->anns[seqid_p].offset + 1 << "\t" << flag1 << "\t"
         << p->len << "\t" << Cigar2String(p->n_cigar, p->cigar, p->len)
         << "\t" << bns->anns[seqid_q].name << "\t"
         << q->pos - bns->anns[seqid_q].offset + 1 << "\t" << flag2 << "\t"
         << q->len << "\t" << Cigar2String(q->n_cigar, q->cigar, q->len)
         << "\t" << status << endl;

    char start_end[256];
    sprintf(start_end, "%d:%d:%d", seqid_p, start, end);
    pair<unordered_map<string, bool>::iterator, bool> iter =
        duplicateTable.insert(make_pair(string(start_end), true));
    if (!iter.second) // insert failed, duplicated
    {
      NumPCRDup += 2;
//since we don't mark PCR duplicate
//      return 1;
    }
  } else {
    // cerr<<"exit from LowQualf"<<endl;
    fout << p->name << "\t" << maxInsert << "\t" << maxInsert2 << "\t" << -1
         << "\t" << bns->anns[seqid_p].name << "\t"
         << p->pos - bns->anns[seqid_p].offset + 1 << "\t" << flag1 << "\t"
         << p->len << "\t" << Cigar2String(p->n_cigar, p->cigar, p->len)
         << "\t" << bns->anns[seqid_q].name << "\t"
         << q->pos - bns->anns[seqid_q].offset + 1 << "\t" << flag2 << "\t"
         << q->len << "\t" << Cigar2String(q->n_cigar, q->cigar, q->len)
         << "\tLowQual" << endl;
    return 2; // low quality
  }
  return 0;
}

int StatCollector::AddAlignment(const bntseq_t *bns, bwa_seq_t *p, bwa_seq_t *q,
                                const gap_opt_t *opt, ofstream &fout,
                                long long &total_add_failed) {
  int seqid(0), seqid2(0), j(0), j2(0);

  /*checking if reads bridges two reference contigs*/
  if (p and p->type != BWA_TYPE_NO_MATCH) {
    j = pos_end(p) - p->pos; // length of read
    bns_coor_pac2real(bns, p->pos, j, &seqid);
    if (p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len) {
      p->type = BWA_TYPE_NO_MATCH; // this alignment bridges two adjacent
                                   // reference sequences
    }
  }
  if (q and q->type != BWA_TYPE_NO_MATCH) {
    j2 = pos_end(q) - q->pos; // length of read
    bns_coor_pac2real(bns, q->pos, j2, &seqid2);
    if (q->pos + j2 - bns->anns[seqid2].offset > bns->anns[seqid2].len) {
      q->type = BWA_TYPE_NO_MATCH; // this alignment bridges two adjacent
                                   // reference sequences
    }
  }
  /*done checking if reads bridges two reference contigs*/
  if (p == 0 || p->type == BWA_TYPE_NO_MATCH) {
    if (q != 0 && AddSingleAlignment(bns, q, opt)) {
      string qname(bns->anns[seqid2].name);
      if (string(qname).find("Y") != string::npos ||
          string(qname).find("X") != string::npos) {
        if (!IsPartialAlign(q)) {
          contigStatusTable[qname].addNumOverlappedReads();
          contigStatusTable[qname].addNumFullyIncludedReads();
        } else
          contigStatusTable[qname].addNumOverlappedReads();
      }
      // TO-DO: adding IsDup function here to generate single end MAX insersize
      // info
      IsDuplicated(bns, p, q, AlignStatus::Second, fout); // only q
      total_add_failed += 1;
      return 1;
    }
    total_add_failed += 2;
    return 0;
  }

  // until now p is aligned
  string pname(bns->anns[seqid].name);
  if (q == 0 || q->type == BWA_TYPE_NO_MATCH) {
    if (AddSingleAlignment(bns, p, opt)) {
      if (string(pname).find('Y') != string::npos ||
          string(pname).find('X') != string::npos) {
        if (!IsPartialAlign(p)) {
          contigStatusTable[pname].addNumOverlappedReads();
          contigStatusTable[pname].addNumFullyIncludedReads();
        } else
          contigStatusTable[pname].addNumOverlappedReads();
      }
      // TO-DO:adding IsDup function here to generate single end MAX insersize
      // info
      IsDuplicated(bns, p, q, AlignStatus::First, fout); // only p
      total_add_failed += 1;
      return 1;
    }
    total_add_failed += 2;
    return 0;
  }

  string qname(bns->anns[seqid2].name);
  // until now both reads are aligned
  if (IsPartialAlign(p)) // p is partially aligned
  {
    if (string(qname).find("Y") != string::npos ||
        string(qname).find("X") != string::npos) {
      if (IsPartialAlign(q)) {
        contigStatusTable[qname].addNumOverlappedReads();
      } else // q is perfectly aligned
      {
        contigStatusTable[qname].addNumOverlappedReads();
        contigStatusTable[qname].addNumFullyIncludedReads();
      }
      if (pname == qname) {
        contigStatusTable[qname].addNumPairOverlappedReads();
      }
      contigStatusTable[pname].addNumOverlappedReads();
    }

    if (IsDuplicated(bns, p, q, AlignStatus::Both, fout) != 1 ||
        opt->cal_dup) // test IsDuplicated first to get insert size info
    {
      if (AddSingleAlignment(bns, p, opt)) {
        if (AddSingleAlignment(bns, q, opt)) {
          return 2;
        } else {
          total_add_failed += 1;
          return 1;
        }
      } else {
        if (AddSingleAlignment(bns, q, opt)) {
          total_add_failed += 1;
          return 1;
        } else {
          total_add_failed += 2;
          return 0;
        }
      }
    }
    total_add_failed += 2;
    return 0;
  } else // p is perfectly aligned
  {
    // p and q are both aligned
    if (string(qname).find('Y') != string::npos ||
        string(qname).find('X') != string::npos) {
      if (IsPartialAlign(q)) // q is partially aligned
      {
        contigStatusTable[qname].addNumOverlappedReads();
        if (pname == qname) {
          contigStatusTable[qname].addNumPairOverlappedReads();
        }
      } else // q is perfectly aligned
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
    if (IsDuplicated(bns, p, q, AlignStatus::Both, fout) != 1 || opt->cal_dup) {
      if (AddSingleAlignment(bns, p, opt)) {
        if (AddSingleAlignment(bns, q, opt)) {
          return 2;
        } else {
          total_add_failed += 1;
          return 1;
        }
      } else {
        if (AddSingleAlignment(bns, q, opt)) {
          total_add_failed += 1;
          return 1;
        } else {
          total_add_failed += 2;
          return 0;
        }
      }
    } else {
      total_add_failed += 2;
      return 0;
    }
  }
  //  cerr << "currently added reads " << total_add_failed << endl;
  return 0;
}

// overload of AddAlignment function for direct bam reading
// return value: 0 failed, 1 add single success, 2 add pair success by using add
// pair interface
bool StatCollector::AddSingleAlignment(SamRecord &p, const gap_opt_t *opt) {
  // TODO:unify both versions
  if (p.getFlag() & SAM_FSU or p.getMapQuality() == 0) {
    return false;
  }

  string seq(p.getSequence()), qual(p.getQuality());
  std::transform(qual.begin(), qual.end(), qual.begin(),
                 [](unsigned char c) -> unsigned char { return c - 33; });

  string posName = p.getReferenceName();

  int pos = p.get1BasedPosition();

  string chrom = p.getReferenceName();

  string refSeq, MD;
  if (p.getStringTag("MD") != nullptr) {
    MD = string(p.getStringTag("MD")->c_str());
  } else
    return false;

  std::vector<bwa_cigar_t> cigar = String2Cigar(p.getCigar());
  refSeq = RecoverRefseqByMDandCigar(seq, MD, &cigar[0], cigar.size());

  int readRealStart(0);
  int refRealStart(0), refRealEnd(0);

  if (chrom.find('@') !=
      std::string::npos) // bam file is generated by FASTQuick itself
  {
    size_t colonPos = posName.find(':');
    size_t atPos = posName.find('@');
    chrom = posName.substr(0, colonPos);

    //    int refCoord = atoi(
    //            posName.substr(colonPos + 1, atPos - colonPos + 1).c_str());
    //            // coordinate of variant site
    char *pEnd;
    int refCoord = static_cast<int>(strtol(
        posName.substr(colonPos + 1, atPos - colonPos + 1).c_str(), &pEnd, 10));

    if (posName[posName.size() - 1] == 'L') {
      readRealStart = refCoord - opt->flank_long_len + pos -
                      1; // real coordinate of current reads on reference
      refRealStart = refCoord - opt->flank_long_len;
      refRealEnd = refCoord + opt->flank_long_len;
    } else {
      readRealStart = refCoord - opt->flank_len + pos -
                      1; // real coordinate of current reads on reference
      refRealStart = refCoord - opt->flank_len;
      refRealEnd = refCoord + opt->flank_len;
    }
  } else {
    readRealStart = pos;
    int readRealEnd = readRealStart + p.getReadLength() - 1 -
                      (__cigar_op(cigar[0]) == 'S' ? __cigar_len(cigar[0]) : 0);
    int overlap_right(0), overlap_left(0);
    int left_flank_len(0), right_flank_len(0) /*flank length of right element*/;
    if (VcfTable.find(chrom) !=
        VcfTable.end()) { // need to locate which artificial ref seq this reads
      // aligned to
      auto right = VcfTable[chrom].lower_bound(readRealStart);
      if (right != VcfTable[chrom].begin() and right != VcfTable[chrom].end())
      // meaning both right and left elements are valid
      {
        auto left(--right);
        if (std::string(VcfRecVec[right->second]->getIDStr()).back() ==
            'L') // long region
          right_flank_len = opt->flank_long_len;
        else
          right_flank_len = opt->flank_len;

        if (std::string(VcfRecVec[left->second]->getIDStr()).back() ==
            'L') // long region
          left_flank_len = opt->flank_long_len;
        else
          left_flank_len = opt->flank_len;

        overlap_right = readRealEnd - (right->first - right_flank_len) + 1;
        overlap_left = left->first + left_flank_len - readRealStart + 1;

        if (overlap_left > 0 and overlap_right > 0) { // both overlap exist
          if (overlap_left > overlap_right) {
            refRealStart = left->first - left_flank_len;
            refRealEnd = left->first + left_flank_len;
          } else {
            refRealStart = right->first - right_flank_len;
            refRealEnd = right->first + right_flank_len;
          }
        } else if (overlap_left > 0) { // only left overlap exists
          refRealStart = left->first - left_flank_len;
          refRealEnd = left->first + left_flank_len;

        } else if (overlap_right > 0) { // only right overlap exists
          refRealStart = right->first - right_flank_len;
          refRealEnd = right->first + right_flank_len;
        } else // not overlapped adjacent elements
          return false;
      } else if (right != VcfTable[chrom].begin()) // left element is valid, but
                                                   // not right element isn't
      {
        auto left(--right);
        if (std::string(VcfRecVec[left->second]->getIDStr()).back() ==
            'L') // long region
          left_flank_len = opt->flank_long_len;
        else
          left_flank_len = opt->flank_len;
        overlap_left = left->first + left_flank_len - readRealStart + 1;
        if (overlap_left > 0) { // left overlap exists
          refRealStart = left->first - left_flank_len;
          refRealEnd = left->first + left_flank_len;

        } else // not overlapped adjacent elements
          return false;
      } else if (right != VcfTable[chrom]
                              .end()) // only this element exists, which is rare
      {
        if (std::string(VcfRecVec[right->second]->getIDStr()).back() ==
            'L') // long region
          right_flank_len = opt->flank_long_len;
        else
          right_flank_len = opt->flank_len;
        overlap_right = readRealEnd - (right->first - right_flank_len) + 1;
        if (overlap_right > 0) { // right overlap exists
          refRealStart = right->first - right_flank_len;
          refRealEnd = right->first + right_flank_len;
        } else
          return false;
      } else // meaning empty set, impossible
      {
        warning(
            "read:%s cannot find adjacent variants on chrom:%s, impossible!",
            p.getReadName(), chrom.c_str());
        exit(EXIT_FAILURE);
      }
    } else
      return false;
  }

  int absoluteSite = readRealStart;
  int tmpCycle(0), relativeCoordOnRead(0) /*relative coord on reads*/,
      relativeCoordOnRef(0);

  bool strand = (p.getFlag() & SAM_FSR);
  u_char mapQ = p.getMapQuality();
  if (strand != 0) {
    tmpCycle = p.getReadLength() - 1;
  }

  if (cigar.size() > 1) // more than M only
  {
    //        return true;
    //        if("MIDS"[__cigar_op(cigar[0])]=='S' /*or
    //        "MIDS"[__cigar_op(cigar[cigar.size()-1])]=='S'*/) return true;
    for (int i = 0; i < cigar.size(); ++i) {
      int cl = __cigar_len(cigar[i]);
      char cop = "MIDS"[__cigar_op(cigar[i])];
      switch (cop) {
      case 'M':
        /*************************variant
         * site****************************************/
        {
          AddMatchBaseInfo(seq, qual, refSeq, chrom, absoluteSite, strand, mapQ,
                           cl, tmpCycle, relativeCoordOnRead,
                           relativeCoordOnRef);
        }
        absoluteSite += cl;
        tmpCycle += cl * sign[strand];
        relativeCoordOnRead += cl;
        relativeCoordOnRef += cl;
        break;
      case 'S':
        tmpCycle += cl * sign[strand];
        relativeCoordOnRead += cl;
        //        relativeCoordOnRef += cl;
        break;
      case 'D':
        absoluteSite += cl;
        relativeCoordOnRef += cl;
        break;
      case 'I':
        tmpCycle += cl * sign[strand];
        relativeCoordOnRead += cl;
        break;
      default:
        warning("Unhandled cigar_op %d:%d\n", cop, cl);
      }
    } // end for

  } else {
    /*************************variant
     * site****************************************/
    AddMatchBaseInfo(seq, qual, refSeq, chrom, absoluteSite, strand, mapQ,
                     p.getReadLength(), tmpCycle, relativeCoordOnRead,
                     relativeCoordOnRef);
  }

  return true;
}

int StatCollector::IsDuplicated(SamFileHeader &SFH, SamRecord &p, SamRecord &q,
                                const gap_opt_t *opt, AlignStatus type,
                                ofstream &fout) {
  int MaxInsert(-1), MaxInsert2(-1);

  if (type == AlignStatus::Second) // q is aligned only
  {
    if (q.getFlag() & SAM_FR1) // if it's read one
      MaxInsert2 = atoi(SFH.getSQTagValue("LN", q.getReferenceName())) -
                   q.get1BasedPosition() + 1;
    else
      MaxInsert = q.get1BasedPosition() + q.getReadLength();

    int cl2;
    cl2 = atoi(q.getCigar());
    if (q.getCigar()[to_string(cl2).size()] !=
        'S') // if the first character is 'S'
    {
      cl2 = 0;
    } else {
      if (q.getFlag() & SAM_FSR) // reversed
        MaxInsert2 -= cl2;
      else
        MaxInsert2 += cl2;
    }

    fout << q.getReadName() << "\t" << -1 << "\t" << MaxInsert2 << "\t" << -1
         << "\t"
         << "*"
         << "\t"
         << "*"
         << "\t"
         << "*"
         << "\t"
         << "*"
         << "\t"
         << "*"
         << "\t" << q.getReferenceName() << "\t" << q.get1BasedPosition()
         << "\t" << q.getFlag() << "\t" << q.getReadLength() << "\t"
         << q.getCigar() << "\tRevOnly" << endl;
    return 0;
  } else if (type == AlignStatus::First) // p is aligned only
  {
    if (p.getFlag() & SAM_FR1) // if it's read one
      MaxInsert = atoi(SFH.getSQTagValue("LN", p.getReferenceName())) -
                  p.get1BasedPosition() + 1;
    else
      MaxInsert = p.get1BasedPosition() + p.getReadLength();

    int cl2;
    cl2 = atoi(p.getCigar());
    if (p.getCigar()[to_string(cl2).size()] !=
        'S') // if the first character is 'S'
    {
      cl2 = 0;
    } else {
      if (p.getFlag() & SAM_FSR) // reversed
        MaxInsert -= cl2;
      else
        MaxInsert += cl2;
    }
    fout << p.getReadName() << "\t" << MaxInsert << "\t" << -1 << "\t" << -1
         << "\t" << p.getReferenceName() << "\t" << p.get1BasedPosition()
         << "\t" << p.getFlag() << "\t" << p.getReadLength() << "\t"
         << p.getCigar() << "\t*"
         << "\t"
         << "*"
         << "\t"
         << "*"
         << "\t"
         << "*"
         << "\t"
         << "*"
         << "\tFwdOnly" << endl;
    return 0;
  } else if (type == AlignStatus::Both) // both aligned
  {
    if (p.get1BasedPosition() < q.get1BasedPosition()) {
      // j = pos_end(q) - q->pos; //length of read
      // bns_coor_pac2real(bns, q->pos, j, &seqid_q);
      MaxInsert = p.get1BasedPosition() + p.getReadLength();
      // j = pos_end(p) - p->pos; //length of read
      // bns_coor_pac2real(bns, p->pos, j, &seqid_p);
      MaxInsert2 = atoi(SFH.getSQTagValue("LN", q.getReferenceName())) -
                   q.get1BasedPosition() + 1;
    } else {
      MaxInsert = q.get1BasedPosition() + q.getReadLength();
      // j = pos_end(p) - p->pos; //length of read
      // bns_coor_pac2real(bns, p->pos, j, &seqid_p);
      MaxInsert2 = atoi(SFH.getSQTagValue("LN", p.getReferenceName())) -
                   p.get1BasedPosition() + 1;
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
    // cerr<<"Duplicate function exit from diff chrom"<<endl;
    fout << p.getReadName() << "\t" << MaxInsert << "\t" << MaxInsert2 << "\t"
         << -1 << "\t" << p.getReferenceName() << "\t" << p.get1BasedPosition()
         << "\t" << p.getFlag() << "\t" << p.getReadLength() << "\t"
         << p.getCigar() << "\t" << q.getReferenceName() << "\t"
         << q.get1BasedPosition() << "\t" << q.getFlag() << "\t"
         << q.getReadLength() << "\t" << q.getCigar() << "\tDiffChrom" << endl;
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
      if (cigar1.find('S') == cigar1.npos &&
          cigar2.find('S') == cigar1.npos) // no softclip
      {
        InsertSizeDist[ActualInsert]++;
        fout << p.getReadName() << "\t" << MaxInsert << "\t" << MaxInsert2
             << "\t" << ActualInsert << "\t" << p.getReferenceName() << "\t"
             << p.get1BasedPosition() << "\t" << p.getFlag() << "\t"
             << p.getReadLength() << "\t" << p.getCigar() << "\t"
             << q.getReferenceName() << "\t" << q.get1BasedPosition() << "\t"
             << q.getFlag() << "\t" << q.getReadLength() << "\t" << q.getCigar()
             << "\tPropPair" << endl;
      } else {
        int cl1, cl2;
        cl1 = atoi(p.getCigar());
        if (p.getCigar()[to_string(cl1).size()] !=
            'S') // if the first character is 'S'
          cl1 = 0;
        cl2 = atoi(q.getCigar());
        if (q.getCigar()[to_string(cl2).size()] !=
            'S') // if the first character is 'S'
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
        fout << p.getReadName() << "\t" << MaxInsert << "\t" << MaxInsert2
             << "\t" << ActualInsert << "\t" << p.getReferenceName() << "\t"
             << p.get1BasedPosition() << "\t" << p.getFlag() << "\t"
             << p.getReadLength() << "\t" << p.getCigar() << "\t"
             << q.getReferenceName() << "\t" << q.get1BasedPosition() << "\t"
             << q.getFlag() << "\t" << q.getReadLength() << "\t" << q.getCigar()
             << "\tPartialPair" << endl;
      }
      char start_end[256];
      sprintf(start_end, "%d:%d", start, end);
      pair<unordered_map<string, bool>::iterator, bool> iter =
          duplicateTable.insert(make_pair(string(start_end), true));
      if (!iter.second) // insert failed, duplicated
      {
        //	cerr<<"Duplicate function exit from duplicate"<<endl;
        NumPCRDup++;
        return 1;
      }
    } else {
      // cerr<<"exit from Inf"<<endl;
      fout << p.getReadName() << "\t" << MaxInsert << "\t" << MaxInsert2 << "\t"
           << -1 << "\t" << p.getReferenceName() << "\t"
           << p.get1BasedPosition() << "\t" << p.getFlag() << "\t"
           << p.getReadLength() << "\t" << p.getCigar() << "\t"
           << q.getReferenceName() << "\t" << q.get1BasedPosition() << "\t"
           << q.getFlag() << "\t" << q.getReadLength() << "\t" << q.getCigar()
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
    // cerr<<"exit from LowQualf"<<endl;
    fout << p.getReadName() << "\t" << MaxInsert << "\t" << MaxInsert2 << "\t"
         << -1 << "\t" << p.getReferenceName() << "\t" << p.get1BasedPosition()
         << "\t" << p.getFlag() << "\t" << p.getReadLength() << "\t"
         << p.getCigar() << "\t" << q.getReferenceName() << "\t"
         << q.get1BasedPosition() << "\t" << q.getFlag() << "\t"
         << q.getReadLength() << "\t" << q.getCigar() << "\tLowQual" << endl;
    return 2; // low quality
  }
  return 0;
}

int StatCollector::AddAlignment(SamFileHeader &SFH, SamRecord *p, SamRecord *q,
                                const gap_opt_t *opt, ofstream &fout,
                                long long &total_add_failed) {
  SamValidator Validator;
  SamValidationErrors VErrors;
  Validator.isValid(SFH, *p, VErrors);
  if (VErrors.numErrors() > 0)
    fprintf(stderr, "%s", VErrors.getNextError()->getMessage());

  if (p == 0 || (p->getFlag() & SAM_FSU)) { // p is not aligned
    if (q == 0 || (q->getFlag() & SAM_FSU)) // both end are not mapped
      return 0;
    else if (IsPartialAlign(*q)) // q is partially mapped, p is not
    {
      string qname(q->getReferenceName());
      if (string(qname).find('Y') != string::npos ||
          string(qname).find('X') != string::npos) {
        contigStatusTable[qname].addNumOverlappedReads();
      }
      if (AddSingleAlignment(*q, opt)) {
        // TO-DO:adding IsDup function here to generate single end MAX insersize
        // info
        IsDuplicated(SFH, *p, *q, opt, AlignStatus::Second, fout); // 1 is q
        return 1;
      }
      return 0;
    } else // q is perfectly aligned, p is not
    {
      string qname(q->getReferenceName());
      if (string(qname).find('Y') != string::npos ||
          string(qname).find('X') != string::npos) {
        contigStatusTable[qname].addNumOverlappedReads();
        contigStatusTable[qname].addNumFullyIncludedReads();
      }
      if (AddSingleAlignment(*q, opt)) {
        // TO-DO:adding IsDup function here to generate single end MAX insersize
        // info
        IsDuplicated(SFH, *p, *q, opt, AlignStatus::Second, fout); // 1 is q
        return 1;
      } else
        return 0;
    }
  } else if (IsPartialAlign(*p)) // p is partially aligned
  {
    if (q == 0 || (q->getFlag() & SAM_FSU)) // q is not aligned
    {
      string pname(p->getReferenceName());
      if (string(pname).find('Y') != string::npos ||
          string(pname).find('X') != string::npos) {
        contigStatusTable[pname].addNumOverlappedReads();
      }
      if (AddSingleAlignment(*p, opt)) {

        IsDuplicated(SFH, *p, *q, opt, AlignStatus::First, fout); // 3 is p
        return 1;
      }
      return 0;
    } else // p is partially aligned q is aligned too
    {
      string pname(p->getReferenceName()), qname(q->getReferenceName());
      if (string(qname).find('Y') != string::npos ||
          string(qname).find('X') != string::npos) {
        if (IsPartialAlign(*q)) {
          contigStatusTable[qname].addNumOverlappedReads();
        } else // q is perfectly aligned
        {
          contigStatusTable[qname].addNumOverlappedReads();
          contigStatusTable[qname].addNumFullyIncludedReads();
        }
        if (pname == qname) {
          contigStatusTable[qname].addNumPairOverlappedReads();
        }
        contigStatusTable[pname].addNumOverlappedReads();
      }

      if (IsDuplicated(SFH, *p, *q, opt, AlignStatus::Both, fout) != 1 ||
          opt->cal_dup) {
        if (AddSingleAlignment(*p, opt)) {
          if (AddSingleAlignment(*q, opt)) {
            return 2;
          } else {
            return 1;
          }
        } else {
          if (AddSingleAlignment(*q, opt)) {
            return 1;
          } else
            return 0;
        }
      }
      total_add_failed += 2;
      return 0;
    }
  } else // p is perfectly aligned
  {
    if (q == 0 || (q->getFlag() & SAM_FSU)) // p is perfectly aligned, q is not
    {
      string pname(p->getReferenceName());
      if (string(pname).find('Y') != string::npos ||
          string(pname).find('X') != string::npos) {
        contigStatusTable[pname].addNumOverlappedReads();
        contigStatusTable[pname].addNumFullyIncludedReads();
      }
      if (AddSingleAlignment(*p, opt)) {

        IsDuplicated(SFH, *p, *q, opt, AlignStatus::First, fout); // 3 is p
        return 1;
      } else
        return 0;
    } else // p and q are both aligned
    {
      string pname(p->getReferenceName()), qname(q->getReferenceName());
      if (string(qname).find('Y') != string::npos ||
          string(qname).find('X') != string::npos) {
        if (IsPartialAlign(*q)) // q is partially aligned
        {
          contigStatusTable[qname].addNumOverlappedReads();
          if (pname == qname) {
            contigStatusTable[qname].addNumPairOverlappedReads();
          }
        } else // q is perfectly aligned
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

      if (IsDuplicated(SFH, *p, *q, opt, AlignStatus::Both, fout) != 1 ||
          opt->cal_dup) {
        if (AddSingleAlignment(*p, opt)) {
          if (AddSingleAlignment(*q, opt)) {
            return 2;
          } else {
            return 1;
          }
        } else {
          if (AddSingleAlignment(*q, opt)) {
            return 1;
          } else
            return 0;
        }
      } else {
        total_add_failed += 2;
        return 0;
      }
    }
  }
  // cerr << "currently added reads " << total_add << endl;
  return 1;
}

int StatCollector::ReadAlignmentFromBam(const gap_opt_t *opt,
                                        const char *BamFile,
                                        std::ofstream &fout, int &total_add) {
  long total_dup(0);
  SamFileHeader SFH;
  SamFile SFIO;
  if (!SFIO.OpenForRead(BamFile, &SFH)) {
    cerr << SFIO.GetStatusMessage() << endl;
    warning("Reading Bam Header Failed!\n");
    exit(1);
  }
  FileStatCollector FSC(BamFile);
  notice("Start reading Bam file...");
  unordered_map<string, SamRecord *> pairBuffer;
  while (1) {
    SamRecord *SR = new SamRecord;
    // SFIO.ReadRecord( SFH, *SR);
    if (!SFIO.ReadRecord(SFH, *SR)) {
      cerr << SFIO.GetStatusMessage() << endl;
      notice("End Bam File Reading...");
      delete SR;
      break;
    }
    FSC.NumRead++;
    FSC.NumBase += SR->getReadLength();
    string readName;
    if (SR->getReadName()[SR->getReadNameLength() - 2] == '\\')
      readName =
          string(SR->getReadName()).substr(0, SR->getReadNameLength() - 3);
    else
      readName = string(SR->getReadName());
    if (!(SR->getFlag() & SAM_FPP)) // read is not mapped in pair
    {
      FSC.TotalRetained += AddAlignment(SFH, SR, 0, opt, fout, FSC.TotalMAPQ);
      delete SR;
      continue;
    } else if (pairBuffer.find(readName) ==
               pairBuffer.end()) // mate not existed
    {
      pairBuffer.insert(make_pair(readName, SR));
    } else {
      SamRecord *SRtmp = pairBuffer[readName];
      FSC.TotalRetained +=
          AddAlignment(SFH, SR, SRtmp, opt, fout, FSC.TotalMAPQ);
      delete SR;
      delete SRtmp;
      pairBuffer.erase(readName);
    }
  } // end of while
  AddFSC(FSC);
  for (unordered_map<string, SamRecord *>::iterator iter = pairBuffer.begin();
       iter != pairBuffer.end(); ++iter) {
    FSC.TotalRetained +=
        AddAlignment(SFH, iter->second, 0, opt, fout, FSC.TotalMAPQ);
    delete iter->second;
  }
  return 0;
}

int StatCollector::RestoreVcfSites(const std::string &RefPath,
                                   const gap_opt_t *opt) {

  VcfHeader header;
  VcfFileReader reader;
  string SelectedSite = RefPath + ".SelectedSite.vcf";
  if (!reader.open(SelectedSite.c_str(), header)) {
    warning("File open failed: %s\n", SelectedSite.c_str());
  }
  string GCpath = RefPath + ".gc";
  ifstream FGC(GCpath, ios_base::binary);
  // int num_so_far(0);
  while (!reader.isEOF()) {
    VcfRecord *VcfLine = new VcfRecord;

    reader.readRecord(*VcfLine);
    VcfRecVec.push_back(VcfLine);

    string chr(VcfLine->getChromStr());
    std::transform(chr.begin(), chr.end(), chr.begin(), ::toupper);
    if (chr.find("CHR") != std::string::npos) {
      chr = chr.substr(3); // strip chr
    }

    int pos = VcfLine->get1BasedPosition();
    VcfTable[chr][pos] = VcfRecVec.size() - 1;
    _GCstruct GCstruct;
    GCstruct.read(FGC);
    int tmp_pos = pos - (GCstruct.len - 1) / 2;
    for (uint32_t i = 0; i != GCstruct.len; ++i) {
      GC[chr][tmp_pos + i] = GCstruct.GC[i];
    }
    // cerr<<chr<<"\t"<<pos<<"\t"<<num_so_far<<endl;
    // num_so_far++;
    if (chr == "X" or chr == "Y") {
      NumXorYMarker++;
      flankRegion.AddRegion(chr,
                            pos - opt->flank_len + opt->read_len * FLANK_EDGE,
                            pos + opt->flank_len - opt->read_len * FLANK_EDGE);
      //      flankRegion.AddRegion(chr,
      //                            pos - opt->flank_len,
      //                            pos + opt->flank_len);
    } else if (std::string(VcfLine->getIDStr()).back() == 'L') {
      NumLongMarker++;
      flankRegion.AddRegion(
          chr, pos - opt->flank_long_len + opt->read_len * FLANK_EDGE,
          pos + opt->flank_long_len - opt->read_len * FLANK_EDGE);
      //      flankRegion.AddRegion(
      //          chr, pos - opt->flank_long_len,
      //          pos + opt->flank_long_len );
    } else {
      NumShortMarker++;
      flankRegion.AddRegion(chr,
                            pos - opt->flank_len + opt->read_len * FLANK_EDGE,
                            pos + opt->flank_len - opt->read_len * FLANK_EDGE);
      //      flankRegion.AddRegion(chr,
      //                            pos - opt->flank_len,
      //                            pos + opt->flank_len);
    }

    SeqVec.emplace_back("");
    QualVec.emplace_back("");
    CycleVec.emplace_back(0);
    MaqVec.emplace_back(0);
    StrandVec.emplace_back(0);
  }
  reader.close();

  std::string dbSNPFile = RefPath + ".dbSNP.subset.vcf";

  if (!reader.open(dbSNPFile.c_str(), header)) {
    error("Open %s failed!", dbSNPFile.c_str());
  }
  while (!reader.isEOF()) {
    VcfRecord VcfLine;
    reader.readRecord(VcfLine);
    string chr(VcfLine.getChromStr());
    std::transform(chr.begin(), chr.end(), chr.begin(), ::toupper);
    if (chr.find("CHR") != std::string::npos) {
      chr = chr.substr(3); // strip chr
    }
    int pos = VcfLine.get1BasedPosition();
    dbSNPTable[chr][pos] = 1;
  }
  reader.close();
  for (auto kv : dbSNPTable)
    notice("Input dbSNP %d sites on chromosome %s", kv.second.size(),
           kv.first.c_str());
  FGC.close();
  return 0;
}

int StatCollector::ReleaseVcfSites() {

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

int StatCollector::GetDepthDist(const string &outputPath,
                                const gap_opt_t *opt) {

  for (auto i = PositionTable.begin(); i != PositionTable.end();
       ++i) // each chr
  {
    for (auto j = i->second.begin(); j != i->second.end(); ++j) // each site
    {
      /************DepthDist**************************************************************/
      assert(DepthVec[j->second] != 0);
      NumBaseMapped += DepthVec[j->second];
      if (DepthVec[j->second] > 1023)
        DepthDist[1023]++;
      else
        DepthDist[DepthVec[j->second]]++;
      GCDist[GC[i->first][j->first]] += DepthVec[j->second];
      PosNum[GC[i->first][j->first]]++;
    }
  }

  //    int max_XorYmarker(0);
  //    if (opt->num_variant_short >= 100000)
  //        max_XorYmarker = 3000;
  //    else if (opt->num_variant_short >= 10000)
  //        max_XorYmarker = 300;
  //    else
  //        max_XorYmarker = 100;

  for (size_t i = 1; i != DepthDist.size(); ++i) {
    NumPositionCovered += DepthDist[i];
    if (i >= 2)
      NumPositionCovered2 += DepthDist[i];
    if (i >= 5)
      NumPositionCovered5 += DepthDist[i];
    if (i >= 10)
      NumPositionCovered10 += DepthDist[i];
  }

  if (targetRegion.IsEmpty()) {
    total_region_size =
        ((opt->flank_len - opt->read_len * FLANK_EDGE) * 2 + 1) *
            NumShortMarker +
        ((opt->flank_long_len - opt->read_len * FLANK_EDGE) * 2 + 1) *
            NumLongMarker +
        ((opt->flank_len - opt->read_len * FLANK_EDGE) * 2 + 1) *
            NumXorYMarker; // X and Y each
  } else {
    total_region_size = flankRegion.Size();
  }
  ofstream fout(outputPath + ".DepthDist");
  DepthDist[0] = total_region_size - NumPositionCovered;
  for (uint32_t i = 0; i != DepthDist.size(); ++i) {
    fout << i << "\t" << DepthDist[i] << endl;
  }
  fout.close();
  return 0;
}

int StatCollector::GetGCDist(const string &outputPath) {
  ofstream fout(outputPath + ".GCDist");

  double MeanDepth =
      NumBaseMapped / (double)NumPositionCovered; /*i.e. coverage in qplot*/
  for (uint32_t i = 0; i != 101; ++i) {
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

int StatCollector::GetEmpRepDist(const string &outputPath) {
  ofstream fout(outputPath + ".EmpRepDist");
  for (uint32_t i = 0; i != EmpRepDist.size(); ++i) {
    fout << i << "\t" << (misEmpRepDist[i]) << "\t" << (EmpRepDist[i]) << "\t"
         << (EmpRepDist[i] == 0
                 ? 0
                 : PHRED((double)(misEmpRepDist[i] + 1) / (EmpRepDist[i] + 2)))
         << endl;
    //        if(misEmpRepDist[i]!=0)
    //            prevQual = PHRED((double) (misEmpRepDist[i] + 1e-6) /
    //            (EmpRepDist[i] + 1e-6));
  }
  fout.close();
  return 0;
}

int StatCollector::GetEmpCycleDist(const string &outputPath) {
  ofstream fout(outputPath + ".EmpCycleDist");
  double prevQual = 0;
  for (uint32_t i = 0; i != EmpCycleDist.size(); ++i) {
    fout << i + 1 << "\t" << misEmpCycleDist[i] << "\t" << EmpCycleDist[i]
         << "\t"
         << (misEmpCycleDist[i] == 0
                 ? prevQual
                 : PHRED((double)(misEmpCycleDist[i] + 1e-6) /
                         (EmpCycleDist[i] + 1e-6)))
         << "\t" << CycleDist[i] << endl;
    if (misEmpCycleDist[i] != 0)
      prevQual =
          PHRED((double)(misEmpCycleDist[i] + 1e-6) / (EmpCycleDist[i] + 1e-6));
  }
  fout.close();
  return 0;
}

int StatCollector::GetInsertSizeDist(const std::string &outputPath) {

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

int StatCollector::GetSexChromInfo(const string &outputPath) {
  ofstream fout(outputPath + ".SexChromInfo");
  for (unordered_map<string, ContigStatus>::iterator iter =
           contigStatusTable.begin();
       iter != contigStatusTable.end(); ++iter) {
    fout << iter->first << "\t" << iter->second.getNumOverlappedReads() << "\t"
         << iter->second.getNumFullyIncludedReads() << "\t"
         << iter->second.getNumPairOverlappedReads() << "\t"
         << iter->second.getNumFullyIncludedPairedReads() << endl;
  }
  fout.close();
  return 0;
}

int StatCollector::ProcessCore(const string &statPrefix, const gap_opt_t *opt) {

  GetDepthDist(statPrefix, opt);
  GetGCDist(statPrefix);
  GetEmpRepDist(statPrefix);
  GetEmpCycleDist(statPrefix);

  GetInsertSizeDist(statPrefix);

  GetSexChromInfo(statPrefix);
  GetPileup(statPrefix, opt);
  SummaryOutput(statPrefix);
  GetGenoLikelihood(statPrefix);
  return 0;
}

int StatCollector::GetPileup(const string &outputPath, const gap_opt_t *opt) {
  ofstream fout(outputPath + ".Pileup");
  int qualoffset = 33;
  if (opt->mode & BWA_MODE_IL13)
    qualoffset = 64;
  for (sort_map::iterator i = VcfTable.begin(); i != VcfTable.end();
       ++i) // each chr
  {
    for (auto j = i->second.begin(); j != i->second.end(); ++j) // each site
    {
      if (SeqVec[j->second].empty())
        continue;
      fout << i->first << "\t" << j->first << "\t.\t"
           << StrandVec[j->second].size() << "\t";
      for (uint32_t k = 0; k != StrandVec[j->second].size(); ++k) {
        if (StrandVec[j->second][k])
          fout << (char)toupper(SeqVec[j->second][k]);
        else
          fout << (char)tolower(SeqVec[j->second][k]);
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

  //    ofstream Dfout(outputPath + ".DebugPileup");
  //    int Dqualoffset = 33;
  //    if (opt->mode | BWA_MODE_IL13) Dqualoffset = 64;
  //    faidx_t *DEBUGREFFAI =
  //    fai_load("/Users/fanzhang/Downloads/FASTQuick_test/hs37d5.fa"); char
  //    DEBUGREGION[1024]; int dummy;
  //
  //    for (auto i = PositionTable.begin(); i != PositionTable.end(); ++i)
  //    //each chr
  //    {
  //        for (auto j = i->second.begin(); j != i->second.end(); ++j) //each
  //        site
  //        {
  //            if (DebugSeqVec[j->second].size() == 0)
  //                continue;
  //
  //            sprintf(DEBUGREGION, "%s:%d-%d", i->first.c_str(), j->first,
  //            j->first);
  //
  //            std::string DEBUGTEST(fai_fetch(DEBUGREFFAI, DEBUGREGION,
  //            &dummy));
  //
  //
  //            Dfout << i->first << "\t" << j->first << "\t" << DEBUGTEST <<
  //            "\t"
  //                  << DebugSeqVec[j->second].size() << "\t";
  //            for (uint32_t k = 0; k != DebugSeqVec[j->second].size(); ++k) {
  //                Dfout << (char) toupper(DebugSeqVec[j->second][k]);
  //            }
  //            Dfout << "\t";
  //            for (uint32_t k = 0; k != DebugQualVec[j->second].size(); ++k)
  //                Dfout << char(DebugQualVec[j->second][k] + Dqualoffset);
  //            Dfout << "\t";
  //            for (uint32_t k = 0; k != DebugCycleVec[j->second].size(); ++k)
  //            {
  //                Dfout << DebugCycleVec[j->second][k];
  //                if (k != DebugCycleVec[j->second].size() - 1)
  //                    Dfout << ",";
  //            }
  //            Dfout << endl;
  //        }
  //    }
  //    Dfout.close();
  return 0;
}

static vector<double> CalLikelihood(const string &seq, const string &qual,
                                    const char &maj,
                                    const char &min) // maj:ref, min:alt
{
  double GL0(0), GL1(0), GL2(0);
  {
    for (uint32_t i = 0; i != seq.size(); ++i) {
      double seq_error = REV_PHRED(qual[i]);
      // fprintf(stderr,"Debug (qual:%d\tseq_error:%f)",qual[i],seq_error);
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
    // fprintf(stderr, "\n");
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

int StatCollector::GetGenoLikelihood(const string &outputPath) {
  ofstream fout(outputPath + ".Likelihood");
  size_t numAllele[4] = {0};
  //    char majAllele, minAllele;
  int maxIndex = 0;
  // for (unsort_map::iterator i = PositionTable.begin();
  //	i != PositionTable.end(); ++i) //each chr
  //{
  //	for (std::unordered_map<int, unsigned int>::iterator j =
  //		i->second.begin(); j != i->second.end(); ++j) //each site
  //	{
  for (sort_map::iterator i = VcfTable.begin(); i != VcfTable.end(); ++i) {
    for (map<int, unsigned int>::iterator j = i->second.begin();
         j != i->second.end(); ++j) {
      uint32_t markerIndex = j->second;
      VcfRecord *VcfLine = VcfRecVec[markerIndex];
      string RefStr(VcfLine->getRefStr()), AltStr(VcfLine->getAltStr());

      CountAllele(numAllele, SeqVec[markerIndex]);
      maxIndex = FindMaxAllele(numAllele, 4);
      //            majAllele = "ACGT"[maxIndex];
      numAllele[maxIndex] = 0;
      maxIndex = FindMaxAllele(numAllele, 4);
      //            minAllele = "ACGT"[maxIndex];
      vector<double> tmpGL =
          CalLikelihood(SeqVec[markerIndex],
                        QualVec[markerIndex], // majAllele, minAllele);
                        RefStr[0], AltStr[0]);
      fout << i->first << "\t" << j->first << "\t" << tmpGL[0] << "\t"
           << tmpGL[1] << "\t" << tmpGL[2];
      fout << endl;
    }
  }
  fout.close();
  return 0;
}

int StatCollector::AddFSC(FileStatCollector a) {
  FSCVec.push_back(a);
  return 0;
}

int StatCollector::SetGenomeSize(long total_size, long total_N_size) {
  ref_genome_size = total_size;
  ref_N_size = total_N_size;
  return 0;
}

int StatCollector::SetTargetRegion(const std::string &targetRegionPath) {
  targetRegion.ReadRegionList(targetRegionPath);
  flankRegion.InnerJoin(targetRegion);
  return 0;
}

/*convert lambda function into regular functions*/
double StatCollector::Q20AvgDepth() {
  long long tmp(0);
  for (size_t i = 0; i != Q20DepthVec.size(); ++i)
    tmp += Q20DepthVec[i];
  return double(tmp) / NumPositionCovered;
}

double StatCollector::Q30AvgDepth() {
  long long tmp(0);
  for (size_t i = 0; i != Q30DepthVec.size(); ++i)
    tmp += Q30DepthVec[i];
  return double(tmp) / NumPositionCovered;
}

size_t StatCollector::MIS500() {
  long long tmp(0), total(0);
  for (size_t i = 500; i != InsertSizeDist.size(); ++i)
    total += InsertSizeDist[i];
  for (size_t i = 500; i != InsertSizeDist.size(); ++i) {
    tmp += InsertSizeDist[i];
    if (tmp > total / 2)
      return i;
  }
  return 0;
}

size_t StatCollector::MIS300() {
  long long tmp(0), total(0);
  for (size_t i = 300; i != InsertSizeDist.size(); ++i)
    total += InsertSizeDist[i];
  for (size_t i = 300; i != InsertSizeDist.size(); ++i) {
    tmp += InsertSizeDist[i];
    if (tmp > total / 2)
      return i;
  }
  return 0;
}

double StatCollector::Q20BaseFraction() {
  long long tmp(0);
  for (size_t i = 0; i != Q20DepthVec.size(); ++i)
    tmp += Q20DepthVec[i];
  return NumBaseMapped == 0 ? 0 : double(tmp) / NumBaseMapped;
}

double StatCollector::Q30BaseFraction() {
  long long tmp(0);
  for (size_t i = 0; i != Q30DepthVec.size(); ++i)
    tmp += Q30DepthVec[i];
  return NumBaseMapped == 0 ? 0 : double(tmp) / NumBaseMapped;
}

int StatCollector::SummaryOutput(const string &outputPath) {
  ofstream fastqOut(outputPath + ".FASTQ.csv");
  fastqOut << "FileIndex,PairEnd1,PairEnd2" << endl;
  for (size_t i = 0; i != FSCVec.size(); ++i) {
    size_t last_slash_idx = FSCVec[i].FileName1.find_last_of("\\/");
    if (std::string::npos != last_slash_idx) {
      FSCVec[i].FileName1.erase(0, last_slash_idx + 1);
    }

    last_slash_idx = FSCVec[i].FileName2.find_last_of("\\/");
    if (std::string::npos != last_slash_idx) {
      FSCVec[i].FileName2.erase(0, last_slash_idx + 1);
    }

    fastqOut << i + 1 << "," << FSCVec[i].FileName1 << ","
             << FSCVec[i].FileName2 << "\n";
  }
  fastqOut.close();
  ofstream fcsv(outputPath + ".Sequence.csv");
  long long total_base(0);
  long long total_reads(0);
  long long total_retained(0);
  long long total_unmapped(0);
  long long total_low_mapQ(0);
  fcsv << "FileIndex,NumOfBases,NumOfReads,NumOfUmappedReads,NumOfLowMAPQReads,"
          "NumOfQCPassReads,ReadLength"
       << endl;
  for (size_t i = 0; i != FSCVec.size(); ++i) {
    fcsv << i + 1 << "," << FSCVec[i].NumBase << "," << FSCVec[i].NumRead << ","
         << FSCVec[i].BwaUnmapped << "," << FSCVec[i].TotalMAPQ << ","
         << FSCVec[i].TotalRetained << ",";
    fcsv << ((FSCVec[i].NumRead == 0) ? 0
                                      : (FSCVec[i].NumBase / FSCVec[i].NumRead))
         << endl;
    total_base += FSCVec[i].NumBase;
    total_reads += FSCVec[i].NumRead;
    total_retained += FSCVec[i].TotalRetained;
    total_unmapped += FSCVec[i].BwaUnmapped;
    total_low_mapQ += FSCVec[i].TotalMAPQ;
  }
  /*note that the average read length will be rounded to integer */
  double avgReadLen = std::floor(
      0.5 + ((total_reads == 0) ? 0 : ((double)total_base / total_reads)));
  fcsv << "Total," << total_base << "," << total_reads << "," << total_unmapped
       << "," << total_low_mapQ << "," << total_retained << ",";
  fcsv << avgReadLen << endl;
  fcsv.close();

  ofstream fout(outputPath + ".Summary");
  fout << "Statistics : "
       << "Value\n";

  /*Estimate total number of mapped reads from report_genome_size region only*/

  double total_mapped_reads = (double)NumBaseMapped / avgReadLen *
                              (ref_genome_size / total_region_size);
  fout << "Estimated Read Mapping Rate : "
       << total_mapped_reads / total_reads
       << "\n";
//       << "[" << total_mapped_reads << "/" << total_reads << "]\n";
  fout << "Estimated Read PCR Duplication Rate : "
       << NumPCRDup / ((double)NumBaseMapped / avgReadLen)
            << "\n";
//       << "[" << NumPCRDup << "/" << (double)NumBaseMapped / avgReadLen << "]\n";

  fout << "Whole Genome Coverage : " << (double)total_base / ref_genome_size
       << "[" << total_base << "/" << ref_genome_size << "]\n";

  auto report_genome_size = targetRegion.IsEmpty()
                                ? (ref_genome_size - ref_N_size)
                                : targetRegion.Size();
  fout << "Expected Read Depth : " << (double)total_base / report_genome_size
       << "[" << total_base << "/" << report_genome_size << "]\n";
  /*auto AvgDepth =
      [&]()->double
  {	long long tmp(0); for (size_t i = 0; i != DepthDist.size(); ++i) tmp +=
  i*DepthDist[i]; return double(tmp) / total_region_size; };*/
  fout << "Estimated Read Depth : ";
  fout << ((NumPositionCovered == 0)
               ? 0
               : NumBaseMapped / (double)total_region_size)
       << "[" << NumBaseMapped << "/" << total_region_size << "]\n";
  fout << "Estimated Percentage of Accessible Genome Covered : "
       << (1. - (double)DepthDist[0] / total_region_size) * 100 << "%\n";
  // output for fraction figure
  /*auto Q20BaseFraction = [&]()->double
  {	long long tmp(0); for (size_t i = 0; i != Q20DepthVec.size(); ++i) tmp
  += Q20DepthVec[i]; return NumBaseMapped==0?0:double(tmp) / NumBaseMapped; };
  auto Q30BaseFraction = [&]()->double
  {	long long tmp(0); for (size_t i = 0; i != Q30DepthVec.size(); ++i) tmp
  += Q30DepthVec[i]; return NumBaseMapped == 0 ? 0 : double(tmp) /
  NumBaseMapped; };*/
  fout << "Total Accessible Genome Size : " << total_region_size << endl;
  double DP1fraction = NumPositionCovered / (double)total_region_size;
  double DP2fraction = NumPositionCovered2 / (double)total_region_size;
  double DP5fraction = NumPositionCovered5 / (double)total_region_size;
  double DP10fraction = NumPositionCovered10 / (double)total_region_size;

  fout << "Depth 1 or above position fraction : " << DP1fraction << endl;
  fout << "Depth 2 or above position fraction : " << DP2fraction << endl;
  fout << "Depth 5 or above position fraction : " << DP5fraction << endl;
  fout << "Depth 10 or above position fraction : " << DP10fraction << endl;
  /*auto Q20AvgDepth =
      [&]()->double
  {	long long tmp(0); for (size_t i = 0; i != Q20DepthVec.size(); ++i) tmp
  += Q20DepthVec[i]; return double(tmp) / total_region_size; };*/
  fout << "Q20 Base Fraction : " << Q20BaseFraction() << endl;
  fout << "Q30 Base Fraction : " << Q30BaseFraction() << endl;
  fout << "Estimated AvgDepth for Q20 bases : " << Q20AvgDepth() << endl;
  /*auto Q30AvgDepth =
      [&]()->double
  {	long long tmp(0); for (size_t i = 0; i != Q30DepthVec.size(); ++i) tmp
  += Q30DepthVec[i]; return double(tmp) / total_region_size; };*/
  fout << "Estimated AvgDepth for Q30 bases : " << Q30AvgDepth() << endl;
  /*auto MIS500 =
      [&]()->double
  {	long long tmp(0), total(0);
      for (size_t i = 500; i != InsertSizeDist.size(); ++i) total +=
  InsertSizeDist[i]; for (size_t i = 500; i != InsertSizeDist.size(); ++i)
      {
          tmp += InsertSizeDist[i]; if (tmp > total / 2) return i;
      }
      return 0;
  };*/

  fout << "Median Insert Size(>=500bp) : " << MIS500() << endl;
  /*auto MIS300 =
      [&]()->double
  {	long long tmp(0), total(0);
      for (size_t i = 300; i != InsertSizeDist.size(); ++i) total +=
  InsertSizeDist[i]; for (size_t i = 300; i != InsertSizeDist.size(); ++i)
      {
      tmp += InsertSizeDist[i]; if (tmp > total / 2) return i;
      }
      return 0;
  };*/
  fout << "Median Insert Size(>=300bp) : " << MIS300() << endl;
  //  fout << "Number of Long Region Variant : " << NumLongMarker << endl;
  //  fout << "Number of Short Region Variant : " << NumShortMarker << endl;
  //  fout << "Number of XorY Variant : " << NumXorYMarker << endl;

  return 0;
}

StatCollector::~StatCollector() { ReleaseVcfSites(); }

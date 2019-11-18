#include "RefBuilder.h"
#include "../libbwa/bwtaln.h"
#include "Error.h"
#include "TargetRegion.h"
#include "Utility.h"
#include "VcfFileReader.h"
#include "VcfHeader.h"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <random>
#include <sstream>

#define DEBUG 0
#define MIN_AF 0.01
#define CALLABLE_RATE 0.95f

// std::default_random_engine generator;
// std::uniform_real_distribution<double> distribution(0.0, 1.0);

static void CalculateGC(const int flank_len, const faidx_t *seq,
                        const std::string &Chrom, const int Position,
                        std::ofstream &FGC) { // Calculate GC content
  _GCstruct GCstruct(flank_len * 2 + 1);
  char region[1024];
  int dummy;
  for (int i = Position - flank_len, t = 0; i != Position + flank_len + 1;
       ++i, ++t) {
    sprintf(region, "%s:%d-%d", Chrom.c_str(), i - 50, i + 49);
    std::string Window(fai_fetch(seq, region, &dummy));
    int total = 0;
    for (unsigned int j = 0; j != Window.size(); ++j) {
      if (Window[j] == 'G' || Window[j] == 'C' || Window[j] == 'g' ||
          Window[j] == 'c')
        total++;
    }
    GCstruct.GC[t] = total;
  }
  GCstruct.write(FGC);
}

bool RefBuilder::VariantCheck(std::string &chr, int position,
                              VcfRecord *VcfLine, int chrFlag) {
  int flank_len = 0;
  if (not IsChromInWhiteList(chr)) // strip chr from e.g. chr11
  {
    warning("%s is not in chromosome white list:1-22,X,Y", chr.c_str());
    return true;
  }

  if (VcfLine->getNumRefBases() != 1 || strlen(VcfLine->getAltStr()) != 1 ||
      VcfLine->getNumAlts() != 1) // filtering indel sites
  {
    warning("%s:%d is around INDEL or MULTI_ALLELIC", chr.c_str(), position);
    return true;
  }

  std::string::size_type sz; // alias of size_t
  const std::string *AFstr = VcfLine->getInfo().getString("AF");
  if (AFstr != nullptr) {
    double AF = std::stod(*AFstr, &sz);
    if (AF < MIN_AF or AF > 1 - MIN_AF) {
      warning("%s:%d is a rare variant AF:%g", chr.c_str(), position, AF);
      return true;
    }
  } else {
    warning("%s:%d does not have AF field", chr.c_str(), position);
    return true;
  }

  if (chrFlag == 1)
    flank_len = flank_long_len;
  else
    flank_len = flank_short_len;

  // ensure no overlapping regions
  if (regionWhiteList.empty()) { // whole genome mode
    if (VcfTable.find(chr) != VcfTable.end()) {
      auto low_iter = VcfTable[chr].upper_bound(position);
      if (low_iter != VcfTable[chr].begin()) {
        auto up_iter = low_iter--;

        int adjacentPos = VcfVec[low_iter->second]->get1BasedPosition();
        int adjacentFlankLen = GetFlankLen(low_iter->second);

        if (abs(position - adjacentPos) < adjacentFlankLen + flank_len) {
          warning("%s:%d is too close to other variants", chr.c_str(),
                  position);
          return true; // overlap with left region
        }

        adjacentPos = VcfVec[up_iter->second]->get1BasedPosition();
        adjacentFlankLen = GetFlankLen(up_iter->second);

        if (up_iter != VcfTable[chr].end() and
            abs(position - adjacentPos) < adjacentFlankLen + flank_len) {
          warning("%s:%d is too close to other variants", chr.c_str(),
                  position);
          return true; // overlap with right region
        }
      }
    }
  }
  //
  if (MaskPath != "Empty") {
    std::string suffix = MaskPath.substr(MaskPath.size() - 3, 3);
    if (suffix == "bed" or suffix == "BED" or suffix == "Bed") {
      if (IsInCallableRegion(chr, position - flank_len, position + flank_len))
        return true;
    } else {
      int dummy_t;
      char region[1024];
      sprintf(region, "%s:%d-%d", chr.c_str(), position - flank_len,
              position + flank_len);
      std::string MaskSeq(fai_fetch(FastaMask, region, &dummy_t));
      double n = std::count(MaskSeq.begin(), MaskSeq.end(), 'P');
      if (n < CALLABLE_RATE * MaskSeq.size()) {
        return true;
      }
    }
  }
  return false;
}

int RefBuilder::GetFlankLen(int index) {
  int adjacentFlankLen = 0;
  if (std::string(VcfVec[index]->getIDStr()).back() == 'L')
    adjacentFlankLen = flank_long_len;
  else
    adjacentFlankLen = flank_short_len;
  return adjacentFlankLen;
}

bool RefBuilder::Skip(std::string &chr, int position, VcfRecord *VcfLine,
                      int chrFlag) {
  int flank_len = 0;
  if (not IsChromInWhiteList(chr)) // strip chr from e.g. chr11
    return true;

  if (VcfLine->getNumRefBases() != 1 || strlen(VcfLine->getAltStr()) != 1 ||
      VcfLine->getNumAlts() != 1) // filtering indel sites
    return true;

  std::string::size_type sz; // alias of size_t
  const std::string *AFstr = VcfLine->getInfo().getString("AF");
  if (AFstr != nullptr) {
    double AF = std::stod(*AFstr, &sz);
    if (AF < MIN_AF or AF > 1 - MIN_AF)
      return true;
  } else {
    warning("%s:%d has no AF tag in INFO field", chr.c_str(), position);
    return true;
  }

  if (chrFlag == 1)
    flank_len = flank_long_len;
  else
    flank_len = flank_short_len;

  // ensure no overlapping regions
  if (1 || regionWhiteList.empty()) { // whole genome mode
    if (VcfTable.find(chr) != VcfTable.end()) {
      auto low_iter = VcfTable[chr].upper_bound(position);
      if (low_iter != VcfTable[chr].begin()) {
        auto up_iter = low_iter--;

        int adjacentPos = VcfVec[low_iter->second]->get1BasedPosition();
        int adjacentFlankLen = GetFlankLen(low_iter->second);

        //        std::cerr<<"out of "<<VcfTable[chr].size()<<"\tcurrent
        //        marker:"<<position<<"("<<flank_len<<")"<<"\tleft:"
        //        <<adjacentPos<<"("<<adjacentFlankLen<<")";

        if (abs(position - adjacentPos) < adjacentFlankLen + flank_len) {
          //          std::cerr<<"\toverlapped"<<std::endl;
          return true; // overlap with left region
        }
        //        std::cerr<<"\tnon-overlapped"<<std::endl;
        if (up_iter != VcfTable[chr].end()) {
          adjacentPos = VcfVec[up_iter->second]->get1BasedPosition();
          adjacentFlankLen = GetFlankLen(up_iter->second);

          //          std::cerr << "out of " << VcfTable[chr].size()
          //                    << "\tcurrent marker:" << position << "(" <<
          //                    flank_len
          //                    << ")"
          //                    << "\tright:" << adjacentPos << "(" <<
          //                    adjacentFlankLen
          //                    << ")";

          if (up_iter != VcfTable[chr].end() and
              abs(position - adjacentPos) < adjacentFlankLen + flank_len) {
            //            std::cerr << "\toverlapped" << std::endl;
            return true; // overlap with right region
          }
          //          std::cerr << "\tnon-overlapped" << std::endl;
        }
      } else {
        int adjacentPos = VcfVec[low_iter->second]->get1BasedPosition();
        int adjacentFlankLen = GetFlankLen(low_iter->second);
        return abs(position - adjacentPos + 1) < adjacentFlankLen + flank_len;
      }
    }
  }

  // ensure no low complexity regions
  if (MaskPath != "Empty") {
    std::string suffix = MaskPath.substr(MaskPath.size() - 3, 3);
    if (suffix == "bed" or suffix == "BED" or suffix == "Bed") {
      if (not IsInCallableRegion(chr, position - flank_len,
                                 position + flank_len)) {
        std::cerr << "[DEBUG]Position:" << chr << "\t" << position
                  << " is in repetitive region." << std::endl;
        return true;
      }
    } else {
      int dummy_t;
      char region[1024];
      sprintf(region, "%s:%d-%d", chr.c_str(), position - flank_len,
              position + flank_len);
      std::string MaskSeq(fai_fetch(FastaMask, region, &dummy_t));
      double n = std::count(MaskSeq.begin(), MaskSeq.end(), 'P');
      if (n < CALLABLE_RATE * MaskSeq.size()) {
        std::cerr << "[DEBUG]Position:" << chr << "\t" << position
                  << " is in repetitive region." << std::endl;
        return true;
      }
    }
  }
  return false;
}

bool RefBuilder::IsChromInWhiteList(std::string &Chrom) {
  if (Chrom.find("chr") != std::string::npos or
      Chrom.find("CHR") != std::string::npos) {
    Chrom = Chrom.substr(3); // strip chr
  }
  return (chromWhiteList.find(Chrom) != chromWhiteList.end());
}

inline static int OverlapLen(int a, int c, int b, int d) {
  if (a <= d && c >= b) {
    // overlap
    return abs(std::min(c, d) - std::max(a, b));
  } else {
    return 0; // no overlap
  }
}

bool RefBuilder::IsInCallableRegion(std::string Chrom, int start, int end) {
  if (Chrom.find("chr") != std::string::npos or
      Chrom.find("CHR") != std::string::npos) {
    Chrom = Chrom.substr(3); // strip chr
  }
  if (repeatRegionList.find(Chrom) == repeatRegionList.end())
    return false;
  else {
    int length = end - start + 1;
    int overlap = 0;
    auto lower_iter = repeatRegionList[Chrom].lower_bound(start);
    if (lower_iter != repeatRegionList[Chrom].begin()) {
      lower_iter--;
      while (lower_iter->first <= end) {
        overlap +=
            OverlapLen(start, end, lower_iter->first, lower_iter->second);
        ++lower_iter;
      }
    } else // the first region
    {
      while (lower_iter->first <= end) {
        overlap +=
            OverlapLen(start, end, lower_iter->first, lower_iter->second);
        ++lower_iter;
      }
    }
    return (length * CALLABLE_RATE <= overlap);
  }
}

RefBuilder::RefBuilder(const std::string &Vcf, const std::string &Ref,
                       const std::string &New, const std::string &DBsnp,
                       const std::string &Mask, const int short_len,
                       const int long_len, const int short_num,
                       const int long_num)
    : VcfPath(Vcf), RefPath(Ref), NewRef(New), dbSNP(DBsnp), MaskPath(Mask),
      flank_short_len(short_len), flank_long_len(long_len),
      num_variant_short(short_num), num_variant_long(long_num) {
  if (num_variant_short >= 100000)
    maxXorYmarker = 3000;
  else if (num_variant_short >= 10000)
    maxXorYmarker = 300;
  else
    maxXorYmarker = 100;

  if (MaskPath != "Empty") {
    std::string suffix = MaskPath.substr(MaskPath.size() - 3, 3);
    if (suffix == "bed" or suffix == "BED" or suffix == "Bed") {
      std::ifstream fin(MaskPath);
      if (!fin.is_open()) {
        std::cerr << "Open file: " << MaskPath << " failed!" << std::endl;
        exit(EXIT_FAILURE);
      }
      std::string line;
      std::string chr;
      int start(0), end(0);
      while (std::getline(fin, line)) {
        std::stringstream ss(line);
        ss >> chr >> start >> end;
        if (chr.find("chr") != std::string::npos or
            chr.find("CHR") != std::string::npos)
          chr = chr.substr(3);
        if (repeatRegionList.find(chr) != repeatRegionList.end()) {
          if (repeatRegionList[chr].find(start) !=
              repeatRegionList[chr].end()) {
            if (repeatRegionList[chr][start] < end)
              repeatRegionList[chr][start] = end;
          } else {
            repeatRegionList[chr][start] = end;
          }
        } else {
          repeatRegionList[chr] = std::map<int, int>();
          repeatRegionList[chr][start] = end;
        }
      }
      fin.close();
      notice("Loading Mask Bed file done!\n");
    } else if (suffix == ".fa" or suffix == "sta" or suffix == ".gz") {
      FastaMask = fai_load(MaskPath.c_str());
      notice("Loading Mask fai file done!\n");
    } else {
      warning("Unknow file type for %s, fasta or bed file is required",
              MaskPath.c_str());
    }
  }

  chromWhiteList = {"1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",
                    "9",  "10", "11", "12", "13", "14", "15", "16",
                    "17", "18", "19", "20", "21", "22", "X",  "Y"};
}

bool RefBuilder::IsMaxNumMarker(
    const std::string &Chrom, int &chrFlag,
    bool isLong) // use this function after checking whiteList
{
  if (Chrom == "X") { /*0:short;1:long;2:Y;3:X*/
    if (nXMarker >= maxXorYmarker)
      return true;
    chrFlag = 3;
  } else if (Chrom == "Y") {
    if (nYMarker >= maxXorYmarker)
      return true;
    chrFlag = 2;
  } else { // must be from chr1-chr22
    if (isLong) {
      chrFlag = 1;
    } else {
      //      float ratio =
      //          num_variant_short / (float)(num_variant_short +
      //          num_variant_long);
      //      double sample = distribution(generator);
      if (/*sample > ratio and*/ nLongMarker < num_variant_long) {
        chrFlag = 1;
      } else if (                                  /*sample <= ratio and*/
                 nShortMarker < num_variant_short) // must be from chr1-chr22
      {
        chrFlag = 0;
      }
      if (nLongMarker >= num_variant_long and
          nShortMarker >= num_variant_short) {
        return true;
      }
    }
  }
  return false;
}

void RefBuilder::IncreaseNumMarker(int chrFlag) {
  switch (chrFlag) { /*0:short;1:long;2:Y;3:X*/
    ;
  case 0:
    nShortMarker++;
    break;
  case 1:
    nLongMarker++;
    break;
  case 2:
    nYMarker++;
    break;
  case 3:
    nXMarker++;
    break;
  default:
    error("Unexpected chromosome flag!");
  }
}

int RefBuilder::SelectMarker(const std::string &RegionList) {

  TargetRegion GI;
  if (RegionList != "Empty")
    GI.ReadRegionList(RegionList);
  notice("Start to select marker set...\n");
  std::string SelectedSite = NewRef + ".SelectedSite.vcf";
  InputFile FoutHapMapSelectedSite(SelectedSite.c_str(), "w");
  std::string BedPath = NewRef + ".bed";
  std::ofstream BedFile(BedPath);

  // read in vcf, hapmap3 sites
  VcfHeader header;
  VcfFileReader reader;
  reader.open(VcfPath.c_str(), header);
  header.write(&FoutHapMapSelectedSite);

  std::string Chrom;
  int Position;

  // begin select exome variants
  while (!reader.isEOF()) {
    int chrFlag(-1) /*0:short;1:long;2:Y;3:X*/;
    VcfRecord *VcfLine = new VcfRecord;
    reader.readRecord(*VcfLine);
    Position = VcfLine->get1BasedPosition();
    Chrom = VcfLine->getChromStr();

    if (IsMaxNumMarker(Chrom, chrFlag, false) ||
        (RegionList != "Empty" and not GI.IsOverlapped(Chrom, Position - 1))) {
      delete VcfLine;
      continue;
    }

    if (Skip(Chrom, Position, VcfLine, chrFlag)) {
      delete VcfLine;
      continue;
    }

    if (chrFlag == 1)
      VcfLine->setID((std::string(VcfLine->getIDStr()) + "$E|L").c_str());
    else
      VcfLine->setID((std::string(VcfLine->getIDStr()) + "$E").c_str());

    VcfTable[Chrom][Position] =
        nShortMarker + nLongMarker + nXMarker + nYMarker;
    VcfVec.push_back(VcfLine);

    IncreaseNumMarker(chrFlag);
  }
  reader.close();
  reader.open(VcfPath.c_str(), header);

  // begin select
  while (!reader.isEOF()) {
    int chrFlag(-1) /*0:short;1:long;2:Y;3:X*/;
    VcfRecord *VcfLine = new VcfRecord;
    reader.readRecord(*VcfLine);
    Position = VcfLine->get1BasedPosition();
    Chrom = VcfLine->getChromStr();

    if (IsMaxNumMarker(Chrom, chrFlag, false) ||
        (RegionList != "Empty" and GI.IsOverlapped(Chrom, Position - 1))) {
      delete VcfLine;
      continue;
    }

    if (Skip(Chrom, Position, VcfLine, chrFlag)) {
      delete VcfLine;
      continue;
    }

    if (chrFlag == 1)
      VcfLine->setID((std::string(VcfLine->getIDStr()) + "|L").c_str());

    VcfTable[Chrom][Position] =
        nShortMarker + nLongMarker + nXMarker + nYMarker;
    VcfVec.push_back(VcfLine);

    IncreaseNumMarker(chrFlag);
  }
  reader.close();

  if (nShortMarker + nLongMarker <
      num_variant_long + num_variant_short) // not enough essential markers
  {
    warning("there are insufficient candidate markers(%d/%d) in %s",
            nShortMarker + nLongMarker, num_variant_long + num_variant_short,
            VcfPath.c_str());
  }
  // output vcf and bed files
  int flank_len = 0;
  for (const auto &kv : VcfTable) {
    for (auto pq : kv.second) {
      if (!VcfVec[pq.second]->write(&FoutHapMapSelectedSite, 1)) {
        warning("Writing retained sites failed!\n");
        exit(EXIT_FAILURE);
      }
      flank_len = GetFlankLen(pq.second);
      BedFile << kv.first << "\t" << pq.first - flank_len << "\t"
              << pq.first + flank_len << std::endl;
      //            delete VcfVec[pq.second];
    }
  }
  BedFile.close();
  FoutHapMapSelectedSite.ifclose();

  reader.open(dbSNP.c_str(), header);
  reader.close();

  // create header for dbSNP subset vcf file
  InputFile FoutdbSNPSelectedSite(
      std::string(NewRef + ".dbSNP.subset.vcf").c_str(), "w");
  header.write(&FoutdbSNPSelectedSite);
  FoutdbSNPSelectedSite.ifclose();

  char cmdline[2048];
  // subset dbsnp
  sprintf(cmdline,
          "sort -k1,1 -k2,2n %s|tabix  -R  - %s  >> %s.dbSNP.subset.vcf",
          BedPath.c_str(), dbSNP.c_str(), NewRef.c_str());
  int ret = system(cmdline);
  if (ret != 0) {
    warning("Building dbsnp subset.vcf failed!\n");
    exit(EXIT_FAILURE);
  }
  return 0;
}

int RefBuilder::InputPredefinedMarker(const std::string &predefinedVcf) {
  notice("Start to load predefined marker set...");
  std::string SelectedSite = NewRef + ".SelectedSite.vcf";
  InputFile FoutHapMapSelectedSite(SelectedSite.c_str(), "w");
  std::string BedPath = NewRef + ".bed";
  std::ofstream BedFile(BedPath);

  // read in predefined vcf
  VcfHeader header;
  VcfFileReader reader;
  reader.open(predefinedVcf.c_str(), header);
  header.write(&FoutHapMapSelectedSite);

  std::string Chrom;
  int Position;

  // begin select
  while (!reader.isEOF()) {
    int chrFlag(-1) /*0:short;1:long;2:Y;3:X*/;
    VcfRecord *VcfLine = new VcfRecord;
    reader.readRecord(*VcfLine);
    bool isLong = (std::string(VcfLine->getIDStr()).back() == 'L');
    Position = VcfLine->get1BasedPosition();
    Chrom = VcfLine->getChromStr();

    if (IsMaxNumMarker(Chrom, chrFlag, isLong)) {
      delete VcfLine;
      continue;
    }

    if (VariantCheck(Chrom, Position, VcfLine, chrFlag)) {
      // warning("%s:%d is a low quality marker. Consider to filter
      // it.",VcfLine->getChromStr(),VcfLine->get1BasedPosition()); continue;
    }

    if (chrFlag == 1 and not isLong)
      VcfLine->setID((std::string(VcfLine->getIDStr()) + "|L").c_str());

    VcfTable[Chrom][Position] =
        nShortMarker + nLongMarker + nXMarker + nYMarker;
    VcfVec.push_back(VcfLine);

    IncreaseNumMarker(chrFlag);
  }

  if (nShortMarker + nLongMarker < num_variant_long + num_variant_short) {
    warning("Insufficient candidate markers %d/%d in %s.",
            nShortMarker + nLongMarker, num_variant_long + num_variant_short,
            predefinedVcf.c_str());

  } else
    notice("%s contains sufficient markers, consider filtered any markers that "
           "triggered warning above.",
           predefinedVcf.c_str());
  // output vcf and bed files
  int flank_len = 0;
  for (auto kv : VcfTable) {
    for (auto pq : kv.second) {
      if (!VcfVec[pq.second]->write(&FoutHapMapSelectedSite, 1)) {
        warning("Writing retained sites failed!\n");
        exit(EXIT_FAILURE);
      }
      if (std::string(VcfVec[pq.second]->getIDStr()).back() == 'L')
        flank_len = flank_long_len;
      else
        flank_len = flank_short_len;
      BedFile << kv.first << "\t" << pq.first - flank_len << "\t"
              << pq.first + flank_len << std::endl;
      //            delete VcfVec[pq.second];
    }
  }
  BedFile.close();
  FoutHapMapSelectedSite.ifclose();

  reader.open(dbSNP.c_str(), header);
  reader.close();

  // create header for dbSNP subset vcf file
  InputFile FoutdbSNPSelectedSite(
      std::string(NewRef + ".dbSNP.subset.vcf").c_str(), "w");
  header.write(&FoutdbSNPSelectedSite);
  FoutdbSNPSelectedSite.ifclose();

  char cmdline[2048];
  // subset dbsnp
  sprintf(cmdline,
          "sort -k1,1 -k2,2n %s|tabix  -R  - %s  >> %s.dbSNP.subset.vcf",
          BedPath.c_str(), dbSNP.c_str(), NewRef.c_str());
  int ret = system(cmdline);
  if (ret != 0) {
    warning("Building dbsnp subset.vcf failed!\n");
    exit(EXIT_FAILURE);
  }
  return 0;
}

void RefBuilder::SubstrRef(const faidx_t *seq, VcfRecord *VcfLine,
                           std::ofstream &FGC, std::ofstream &FaOut) {
  int flank_len;
  int dummy;
  char region[1024];
  char newChrName[1024];

  if (std::string(VcfLine->getIDStr()).find('L') !=
      std::string::npos) // Long region
  {
    flank_len = flank_long_len;
    sprintf(newChrName, ">%s:%d@%s/%s|L", VcfLine->getChromStr(),
            VcfLine->get1BasedPosition(), VcfLine->getRefStr(),
            VcfLine->getAltStr());
  } else {
    flank_len = flank_short_len;
    sprintf(newChrName, ">%s:%d@%s/%s", VcfLine->getChromStr(),
            VcfLine->get1BasedPosition(), VcfLine->getRefStr(),
            VcfLine->getAltStr());
  }
  sprintf(region, "%s:%d-%d", VcfLine->getChromStr(),
          VcfLine->get1BasedPosition() - flank_len,
          VcfLine->get1BasedPosition() + flank_len);
  std::string FetchedSeq(fai_fetch(seq, region, &dummy));

  CalculateGC(flank_len, seq, VcfLine->getChromStr(),
              VcfLine->get1BasedPosition(), FGC);

  //    SeqVec.push_back(FetchedSeq.substr(0, flank_len) +
  //    string(VcfLine.getRefStr()) +
  //                     FetchedSeq.substr(flank_len + 1, flank_len));
  //    sprintf(region, "%s:%d@%s/%s", VcfLine.getChromStr(),
  //    VcfLine.get1BasedPosition(), VcfLine.getRefStr(), VcfLine.getAltStr());
  //    RefTableIndex.insert(make_pair(string(region), nseqs));
  //    nseqs++;
  FaOut << newChrName << std::endl;
  FaOut << FetchedSeq.substr(0, flank_len) + std::string(VcfLine->getRefStr()) +
               FetchedSeq.substr(flank_len + 1, flank_len)
        << std::endl;
}

int RefBuilder::PrepareRefSeq() {
  faidx_t *seq = fai_load(RefPath.c_str());
  if (seq == 0)
    error("Please check if %s file exists!", RefPath.c_str());
  notice("Loading Ref fai file done!\n");
  std::string GCpath = NewRef + ".gc";
  std::ofstream FGC(GCpath, std::ios_base::binary);
  std::ofstream FaOut(NewRef);
  for (auto &kv : VcfTable) {
    for (auto &kv2 : kv.second) // each marker in VcfTable
    {
      SubstrRef(seq, VcfVec[kv2.second], FGC, FaOut);
      delete VcfVec[kv2.second];
    }
  }
  FaOut.close();
  FGC.close();
  return 0;
}
RefBuilder::~RefBuilder() {
  //    SeqVec.clear();
  //    RefTableIndex.clear();
}

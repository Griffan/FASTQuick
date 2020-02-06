#include "RefBuilder.h"

#include <algorithm>
#include <cctype>
#include <sstream>

#include "../libbwa/bwtaln.h"
#include "Error.h"
#include "TargetRegion.h"
#include "Utility.h"
#include "VcfFileReader.h"
#include "VcfHeader.h"

#define DEBUG 0
#define MIN_AF 0.01
#define CALLABLE_RATE 0.995f

static std::string ExtractSeq(const faidx_t *faidx, const std::string &chrom,
                              int beg, int end) {
  char region[128], chrRegion[128];
  int dummy;
  sprintf(region, "%s:%d-%d", chrom.c_str(), beg, end);
  char *seq = fai_fetch(faidx, region, &dummy);
  if (!seq) {
    sprintf(chrRegion, "chr%s:%d-%d", chrom.c_str(), beg, end);
    seq = fai_fetch(faidx, chrRegion, &dummy);
    if (!seq) {
      error("Cannot find %s from the reference file!\n", region);
    }
  }
  std::string ret(seq);
  free(seq);
  return ret;
}

static void CalculateGC(const int flank_len, const faidx_t *seq,
                        const std::string &Chrom, const int Position,
                        std::ofstream &FGC) { // Calculate GC content
  _GCstruct GCstruct(flank_len * 2 + 1);
  for (int i = Position - flank_len, t = 0; i != Position + flank_len + 1;
       ++i, ++t) {
    std::string Window(ExtractSeq(seq, Chrom.c_str(), i - 50, i + 49));
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
      std::string MaskSeq(ExtractSeq(
          FastaMask, chr.c_str(), position - flank_len, position + flank_len));
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

  int flank_len = 0;
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
        if (abs(position - adjacentPos) < adjacentFlankLen + flank_len) {
          return true; // overlap with left region
        }
        if (up_iter != VcfTable[chr].end()) {
          adjacentPos = VcfVec[up_iter->second]->get1BasedPosition();
          adjacentFlankLen = GetFlankLen(up_iter->second);
          if (up_iter != VcfTable[chr].end() and
              abs(position - adjacentPos) < adjacentFlankLen + flank_len) {
            return true; // overlap with right region
          }
        }
      } else { // first item
        int adjacentPos = VcfVec[low_iter->second]->get1BasedPosition();
        int adjacentFlankLen = GetFlankLen(low_iter->second);
        if (abs(position - adjacentPos + 1) < adjacentFlankLen + flank_len)
          return true;
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
    } else { // fasta file
      std::string MaskSeq(ExtractSeq(
          FastaMask, chr.c_str(), position - flank_len, position + flank_len));
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

bool RefBuilder::IsChromInWhiteList(const std::string &Chrom) {
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

  if (callableRegionList.find(Chrom) == callableRegionList.end())
    return false;
  else {
    int length = end - start + 1;
    int overlap = 0;
    auto lower_iter = callableRegionList[Chrom].lower_bound(start);
    if (lower_iter != callableRegionList[Chrom].begin()) {
      lower_iter--;
      while (lower_iter != callableRegionList[Chrom].end() and
             lower_iter->first <= end) {
        overlap +=
            OverlapLen(start, end, lower_iter->first, lower_iter->second);
        ++lower_iter;
      }
    } else // the first region
    {
      while (lower_iter != callableRegionList[Chrom].end() and
             lower_iter->first <= end) {
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
        if (callableRegionList.find(chr) != callableRegionList.end()) {
          if (callableRegionList[chr].find(start) !=
              callableRegionList[chr].end()) {
            if (callableRegionList[chr][start] < end)
              callableRegionList[chr][start] = end;
          } else {
            callableRegionList[chr][start] = end;
          }
        } else {
          callableRegionList[chr] = std::map<int, int>();
          callableRegionList[chr][start] = end;
        }
      }
      fin.close();
      notice("Loading Mask Bed file done!");
    } else if (suffix == ".fa" or suffix == "sta" or suffix == ".gz") {
      FastaMask = fai_load(MaskPath.c_str());
      if (!FastaMask)
        notice("Loading Mask fai file done!");
      else
        error("Check if %s exists and Mask fai chromosome name is consistent "
              "with reference fai!",
              MaskPath.c_str());
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
    const std::string &Chrom, int &chrFlag, bool isForcedLong,
    bool isForcedShort) // use this function after checking whiteList
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
    if (nLongMarker >= num_variant_long and nShortMarker >= num_variant_short) {
      return true;
    }
    if (isForcedLong) {
      if (isForcedShort)
        error("isForcedLong and isForcedShort cannot be set at the same time");
      else {
        chrFlag = 1;
      }
    } else if (isForcedShort) {
      chrFlag = 0;
    } else {
      if (nLongMarker < num_variant_long) {
        chrFlag = 1;
      } else if (nShortMarker < num_variant_short) // must be from chr1-chr22
      {
        chrFlag = 0;
      }
    }
  }
  return false;
}

void RefBuilder::IncreaseNumMarker(int chrFlag) {
  switch (chrFlag) { /*0:short;1:long;2:Y;3:X*/
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

/*
 * The marker selection priority is as follows:
 * 1. markers in target region with long flank region
 * 2. markers in target region with short flank region
 * 3. markers not in target region with long flank region
 * 4. markers not in target region with short flank region
 */
int RefBuilder::SelectMarker(const std::string &RegionList) {
  notice("Start to select markers...");
  std::string SelectedSite = NewRef + ".SelectedSite.vcf";

  // read in vcf, hapmap3 sites
  std::string Chrom;
  int Position(0);
  int chrFlag(-1);
  bool isForcedShort(false);
  VcfHeader header;
  VcfFileReader reader;
  reader.open(VcfPath.c_str(), header);
  // begin select variants in target regions
  if (RegionList != "Empty") {
    notice("Start to select markers from target regions...");
    TargetRegion targetRegion;
    targetRegion.ReadRegionList(RegionList);
    while (!reader.isEOF()) {
      chrFlag = -1 /*0:short;1:long;2:Y;3:X*/;
      isForcedShort = false;
      VcfRecord *VcfLine = new VcfRecord;
      reader.readRecord(*VcfLine);
      Position = VcfLine->get1BasedPosition();

      Chrom = VcfLine->getChromStr();
      std::transform(Chrom.begin(), Chrom.end(), Chrom.begin(), ::toupper);
      if (Chrom.find("chr") != std::string::npos or
          Chrom.find("CHR") != std::string::npos) {
        Chrom = Chrom.substr(3); // strip chr
      }
      notice("Now test if %s:%d is selected...", Chrom.c_str(), Position);

    RESCUE:
      if (IsMaxNumMarker(Chrom, chrFlag, false, isForcedShort)) {
        delete VcfLine;
        continue;
      }

      if (not targetRegion.IsOverlapped(Chrom, Position - 1)) {
        delete VcfLine;
        continue;
      }

      if (Skip(Chrom, Position, VcfLine, chrFlag)) {
        if (not isForcedShort) {
          isForcedShort = true;
          goto RESCUE;
        }
        delete VcfLine;
        continue;
      }

      if (chrFlag == 1)
        VcfLine->setID((std::string(VcfLine->getIDStr()) + "$E|L").c_str());
      else
        VcfLine->setID((std::string(VcfLine->getIDStr()) + "$E").c_str());
      notice("Now confirm if %s:%d is selected as %s...", Chrom.c_str(),
             Position, VcfLine->getIDStr());

      VcfTable[Chrom][Position] =
          nShortMarker + nLongMarker + nXMarker + nYMarker;
      VcfVec.push_back(VcfLine);

      IncreaseNumMarker(chrFlag);
    }
    reader.close();
    reader.open(VcfPath.c_str(), header);
  }

  // begin select markers outside of target region
  while (!reader.isEOF()) {
    chrFlag = -1 /*0:short;1:long;2:Y;3:X*/;
    VcfRecord *VcfLine = new VcfRecord;
    reader.readRecord(*VcfLine);
    Position = VcfLine->get1BasedPosition();

    Chrom = VcfLine->getChromStr();
    std::transform(Chrom.begin(), Chrom.end(), Chrom.begin(), ::toupper);
    if (Chrom.find("chr") != std::string::npos or
        Chrom.find("CHR") != std::string::npos) {
      Chrom = Chrom.substr(3); // strip chr
    }

    if (IsMaxNumMarker(Chrom, chrFlag, false, false)) {
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
  InputFile FoutHapMapSelectedSite(SelectedSite.c_str(), "w");
  std::string BedPath = NewRef + ".bed";
  std::ofstream BedFile(BedPath);
  header.write(&FoutHapMapSelectedSite);
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

    std::transform(Chrom.begin(), Chrom.end(), Chrom.begin(), ::toupper);
    if (Chrom.find("chr") != std::string::npos or
        Chrom.find("CHR") != std::string::npos) {
      Chrom = Chrom.substr(3); // strip chr
    }

    if (IsMaxNumMarker(Chrom, chrFlag, isLong, false)) {
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

  std::string FetchedSeq(ExtractSeq(seq, VcfLine->getChromStr(),
                                    VcfLine->get1BasedPosition() - flank_len,
                                    VcfLine->get1BasedPosition() + flank_len));

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
  notice("Loading Ref fai file done!");
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
RefBuilder::~RefBuilder() {}

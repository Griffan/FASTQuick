//
// Created by Fan Zhang on 11/14/19.
//

#include "TargetRegion.h"
#include "Error.h"
#include <algorithm>
#include <fstream>
#include <sstream>

TargetRegion::TargetRegion(const std::string &regionListPath) {
  notice("Read target region list...");
  std::ifstream fin(regionListPath);
  if (!fin.is_open())
    error("Region list bed file:%s open failed!", regionListPath.c_str());
  std::string line;
  while (std::getline(fin, line)) {
    std::string chr;
    int start(0), end(0);
    std::stringstream ss(line);
    ss >> chr >> start >> end;
    if (chr.find("chr") != std::string::npos)
      chr = chr.substr(3);
    if (regionWhiteList.find(chr) != regionWhiteList.end()) {
      if (regionWhiteList[chr].find(start) != regionWhiteList[chr].end()) {
        if (regionWhiteList[chr][start] < end)
          regionWhiteList[chr][start] = end;
      } else {
        regionWhiteList[chr][start] = end;
      }
    } else {
      regionWhiteList[chr] = std::map<int, int>();
      regionWhiteList[chr][start] = end;
    }
  }
  fin.close();
}

int TargetRegion::ReadRegionList(const std::string &regionListPath) {
  notice("Read target region list...");
  std::ifstream fin(regionListPath);
  if (!fin.is_open())
    error("Region list bed file:%s open failed!", regionListPath.c_str());
  std::string line;
  while (std::getline(fin, line)) {
    std::string chr;
    int start(0), end(0);
    std::stringstream ss(line);
    ss >> chr >> start >> end;
    if (chr.find("chr") != std::string::npos)
      chr = chr.substr(3);
    if (regionWhiteList.find(chr) != regionWhiteList.end()) {
      if (regionWhiteList[chr].find(start) != regionWhiteList[chr].end()) {
        if (regionWhiteList[chr][start] < end)
          regionWhiteList[chr][start] = end;
      } else {
        regionWhiteList[chr][start] = end;
      }
    } else {
      regionWhiteList[chr] = std::map<int, int>();
      regionWhiteList[chr][start] = end;
    }
  }
  fin.close();
  return 0;
}

bool TargetRegion::IsOverlapped(std::string Chrom, int start) {

  std::transform(Chrom.begin(), Chrom.end(), Chrom.begin(), ::toupper);
  if (Chrom.find("chr") != std::string::npos or
      Chrom.find("CHR") != std::string::npos) {
    Chrom = Chrom.substr(3); // strip chr
  }
  if (regionWhiteList.find(Chrom) == regionWhiteList.end())
    return false;
  else {
    auto lower_iter = regionWhiteList[Chrom].lower_bound(start);
    if (lower_iter != regionWhiteList[Chrom].begin())
      lower_iter--;
    // std::cerr<<"lower_bound:(size:"<<regionWhiteList[Chrom].size()<<")"<<Chrom<<"
    // "<<lower_iter->first<<" "<<lower_iter->second<<"
    // input:"<<start<<std::endl;
    if (lower_iter->first <= (start) and lower_iter->second > (start))
      return true;
    else
      return (++lower_iter)->first <= (start) and lower_iter->second > (start);
  }
}

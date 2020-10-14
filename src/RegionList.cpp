//
// Created by Fan Zhang on 11/14/19.
//

#include "RegionList.h"
#include "Error.h"
#include <algorithm>
#include <fstream>
#include <sstream>

RegionList::RegionList(const std::string &regionListPath) {
  ReadRegionList(regionListPath);
}

int RegionList::ReadRegionList(const std::string &regionListPath) {
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
    std::transform(chr.begin(), chr.end(), chr.begin(), ::toupper);
    if (chr.find("CHR") != std::string::npos) {
      chr = chr.substr(3); // strip chr
    }
    if (regionList.find(chr) != regionList.end()) {
      if (regionList[chr].find(start) != regionList[chr].end()) {
        if (regionList[chr][start] < end)//update
          regionList[chr][start] = end;
      } else {
        regionList[chr][start] = end;
      }
    } else {
      regionList[chr][start] = end;
    }
  }
  fin.close();
  return 0;
}

bool RegionList::IsOverlapped(const std::string & Chrom, int start) {

  if (regionList.empty())
    return false;
  else if (regionList.find(Chrom) == regionList.end())
    return false;
  else {
    auto lower_iter = regionList[Chrom].lower_bound(start);
    if (lower_iter != regionList[Chrom].begin())
      lower_iter--;
    if (lower_iter->first <= start and lower_iter->second >= start)
      return true;
    lower_iter++;
    if (lower_iter->first <= start and lower_iter->second >= start)
      return true;
  }
  return false;
}

bool RegionList::AddRegion(const std::string& Chrom, int start, int end)
{
  std::string chr(Chrom);
  std::transform(chr.begin(), chr.end(), chr.begin(), ::toupper);
  if (chr.find("CHR") != std::string::npos) {
    chr = chr.substr(3); // strip chr
  }
  regionList[chr][start] = end;
  return true;
}

bool RegionList::Collapse()//union of internal regions
{
  std::map<std::string, std::map<int, int> > tmpList;

  for(auto kv : regionList)
  {
    int beg1,end1,beg2,end2;
    auto holder = kv.second.begin();
    for(auto iter = kv.second.begin(); iter != kv.second.end(); ++iter)
    {
      beg1 = holder->first;
      end1 = holder->second;
      beg2 = iter->first;
      end2 = iter->second;
      if(end1 >= end2)//holder includes iter
        continue;
      else if(end1 < beg2) {//holder and iter do not overlap
        tmpList[kv.first][beg1] = end1;
        holder = iter;
      }
      else// partial overlap
      {
        tmpList[kv.first][beg1] = end2;
        holder->second = end2;
      }
    }
    tmpList[kv.first][holder->first] = holder->second;//last holder, regardless if dup or not
  }

  regionList = tmpList;
  uint64_t len = 0;
  for (const auto & kv : regionList) {
    for (const auto & kv2 : kv.second) {
      if(kv2.second < kv2.first) error("abnormal region %s:%d-%d",kv.first.c_str(), kv2.first, kv2.second);
      len += (kv2.second - kv2.first)+1;
    }
  }
  length = len;
  return true;
}

bool RegionList::Join(RegionList & b, bool isUnion = false)
{
  if(isUnion) {
    for (auto &kv : b.regionList) {
        for (auto kv2 : kv.second) {
          AddRegion(kv.first, kv2.first, kv2.second);
        }
    }
  }
  else//isIntersection
  {
    Collapse();
    b.Collapse();

    RegionList tmpList;
    for (auto &kv : b.regionList) {
      if(regionList.find(kv.first) != regionList.end()) {//chr shown in List A
        auto thisIter = regionList[kv.first].begin();
        auto iter = kv.second.begin();
        int beg1,end1,beg2,end2;
        while(thisIter != regionList[kv.first].end() && iter!=kv.second.end())
        {
          beg1 = thisIter->first;
          end1 = thisIter->second;
          beg2 = iter->first;
          end2 = iter->second;
          if(beg1 <= beg2)
          {
            if(end1 > end2) {// [1,4] and [2,3]
              tmpList.AddRegion(kv.first, beg2, end2);
              iter++;
            }
            else if (end1 > beg2) {//[1,3] and [2,4]
              tmpList.AddRegion(kv.first, beg2, end1);
              thisIter++;
            }
            else{ //[1,2] and [3,4]
              thisIter++;
            }
          }
          else//beg1 >beg2
          {
            if(end1 <= end2) {//  [2,3] and [1,4]
              tmpList.AddRegion(kv.first, beg1, end1);
              thisIter++;
            }
            else if (end1 > beg2 and beg1 < end2) {//[2,4] and [1,3]
              tmpList.AddRegion(kv.first, beg1, end2);
              iter++;
            }
            else{ //[3,4] and [1,2]
              iter++;
            }
          }
        }
      }
    }
    regionList = tmpList.regionList;
  }
  //collapse
  return Collapse();
}


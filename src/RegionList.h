//
// Created by Fan Zhang on 11/14/19.
//

#ifndef FASTQUICK_REGIONLIST_H
#define FASTQUICK_REGIONLIST_H

#include <map>
#include <string>
#include <iostream>

class RegionList {
public:
  RegionList() = default;
  explicit RegionList(const std::string &bedFile);
  int ReadRegionList(const std::string &regionListPath, bool autoCollapse = true);
  bool IsOverlapped(const std::string& Chrom, int start);
  //caller ensure b is collapsed before pass in
  bool Join(const RegionList &b, bool isUnion);
  inline bool InnerJoin(const RegionList & b)
  {
    return Join(b, false);
  }
  inline bool OuterJoin(const RegionList & b)
  {
    return Join(b, true);
  }
  bool AddRegion(const std::string& Chrom, int start, int end);
  bool Collapse();
  bool IsEmpty(){return regionList.empty();}
  void PrintRegion(const std::string & label)
  {
    std::cerr<<"Print RegionList:\n";
    for(auto kv : regionList)
    {
      for(auto kv2 : kv.second)
      {
        std::cerr<<label<<"\t"<<kv.first<<"\t"<<kv2.first<<"\t"<<kv2.second<<std::endl;
      }
    }
  }

  bool operator==(const RegionList& b)
  {
    return regionList.size() == b.regionList.size()
           && std::equal(regionList.begin(), regionList.end(),
                         b.regionList.begin());
  }

  void Clear()
  {
    regionList.clear();
  }

  uint64_t Size()
  {
    return length;
  }
private:
  std::map<std::string, std::map<int, int>> regionList;
  uint64_t length = 0;
};

#endif // FASTQUICK_REGIONLIST_H

//
// Created by Fan Zhang on 11/14/19.
//

#ifndef FASTQUICK_TARGETREGION_H
#define FASTQUICK_TARGETREGION_H

#include <map>
#include <string>
class TargetRegion {
public:
  bool isEmpty = false;
  TargetRegion() = default;
  TargetRegion(const std::string &bedFile);
  int ReadRegionList(const std::string &regionListPath);
  bool IsOverlapped(const std::string& Chrom, int start);

private:
  std::map<std::string, std::map<int, int>> regionWhiteList;
};

#endif // FASTQUICK_TARGETREGION_H

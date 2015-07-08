//
// Created by Fan Zhang on 7/3/15.
//

#ifndef INSERTSIZEESTIMATOR_INSERTSIZEESTIMATOR_H
#define INSERTSIZEESTIMATOR_INSERTSIZEESTIMATOR_H

#include <iosfwd>
#include <string>
#include <vector>
#include <map>

class InsertSizeRecord {
public:
    std::string ReadName;
    int ObservedInsertSize;
    int MaxInsertSize;
    double Weight;
    InsertSizeRecord(std::string &R,int Observed,int Max)
    {
        ReadName=R;
        ObservedInsertSize=Observed;
        MaxInsertSize=Max;
        Weight=1.;
    }

};
class InsertSizeEstimator {

public:
    std::map<int,double> Distribution;
    std::vector<InsertSizeRecord> ObsRecordVec;
    //std::vector<InsertSizeRecord> UnobsRecordVec;
    //std::unordered_map<int,double> DistWeight;
    bool init;
    InsertSizeEstimator():init(true) { }

    int InputInsertSizeTable(std::string FileName);
    int UpdateWeight();
    int Sort();
    int GetInsertDist(std::string OutFileName);
    int UpdateInsertDist();

};


#endif //INSERTSIZEESTIMATOR_INSERTSIZEESTIMATOR_H

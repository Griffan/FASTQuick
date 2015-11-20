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
    InsertSizeRecord(std::string &R,int Observed,int Max, double Wt)
    {
        ReadName=R;
        ObservedInsertSize=Observed;
        MaxInsertSize=Max;
        Weight=Wt;
    }
    InsertSizeRecord()
    {
        ReadName="Empty";
        ObservedInsertSize=-1;
        MaxInsertSize=-1;
        Weight=0;
    }
};


typedef std::vector<InsertSizeRecord>  RecordVec;
typedef std::vector<RecordVec> RecordVecDist;

typedef std::vector<double> Dist;
#define INSERT_LIMIT 4096
class InsertSizeEstimator {

public:
    int totalPair=0;
    std::map<int,double> Distribution;
    RecordVecDist ObsDistVec;
    RecordVecDist MisDistVec;
    Dist ObsDist;
    Dist MisDist;

    //std::vector<InsertSizeRecord> UnobsRecordVec;
    //std::unordered_map<int,double> DistWeight;
    bool init;
    InsertSizeEstimator():init(true) {
        RecordVec tmpVec(0,InsertSizeRecord());
        ObsDistVec=RecordVecDist(INSERT_LIMIT,tmpVec);
        MisDistVec=RecordVecDist(INSERT_LIMIT,tmpVec);
        ObsDist=Dist(INSERT_LIMIT,0.);
        MisDist=Dist(INSERT_LIMIT,0.);
    }

    int InputInsertSizeTable(std::string FileName);
    int UpdateWeight();
    int Sort();
    int GetInsertDist();
    int UpdateInsertDist();

};


#endif //INSERTSIZEESTIMATOR_INSERTSIZEESTIMATOR_H

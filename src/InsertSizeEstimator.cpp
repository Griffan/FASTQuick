//
// Created by Fan Zhang on 7/3/15.
//

#include <fstream>
#include <iostream>
#include <sstream>
#include "InsertSizeEstimator.h"


int InsertSizeEstimator::InputInsertSizeTable(std::string FileName) {
    std::ifstream fin(FileName);
    if(!fin.is_open()) std::cerr<<"File "<<FileName<<" Open Failed!"<<std::endl;
    std::string line;
    while(getline(fin,line))
    {
        std::stringstream ss(line);
        std::string ReadName,Chr1,Chr2,Status;
        int Max,Obs,Pos1,Pos2;
        ss>>ReadName;
        ss>>Max;
        ss>>Obs;
        ss>>Chr1;
        ss>>Pos1;
        ss>>Chr2;
        ss>>Pos2;
        ss>>Status;
        if(Status =="DiffChrom") continue;
        else if(Status != "PropPair")//Obs is NA
        {
            ObsRecordVec.push_back(InsertSizeRecord(ReadName,-1,Max));
        }
        else
        ObsRecordVec.push_back(InsertSizeRecord(ReadName,Obs,Max));
    }
    return 0;
}

int InsertSizeEstimator::UpdateWeight() {
    for (int i = ObsRecordVec.size()-1; i !=0; --i) {
        if(ObsRecordVec[i].ObservedInsertSize != -1) continue;//skip when
        double sum(0.0);
        for (int j = i; j < ObsRecordVec.size(); ++j) {
            if ( (ObsRecordVec[j].ObservedInsertSize != -1) && (ObsRecordVec[j].ObservedInsertSize > ObsRecordVec[i].MaxInsertSize) )
            sum+=ObsRecordVec[j].Weight;
        }
        for (int k = i; k < ObsRecordVec.size(); ++k) {
            if ( (ObsRecordVec[k].ObservedInsertSize != -1) && (ObsRecordVec[k].ObservedInsertSize > ObsRecordVec[i].MaxInsertSize) )
            ObsRecordVec[k].Weight+=1.0/sum;
        }
    }
    return 0;
}

int InsertSizeEstimator::Sort() {
    std::sort(ObsRecordVec.begin(), ObsRecordVec.end(),[&](InsertSizeRecord a,InsertSizeRecord b){ return a.MaxInsertSize < b.MaxInsertSize;});
    return 0;
}

int InsertSizeEstimator::GetInsertDist(std::string FileName) {

    std::ofstream fout(FileName);
    if(!fout.is_open()) std::cerr<<"File "<<FileName<<" Open Failed!"<<std::endl;
    for(auto kv:Distribution)
    {
        fout<<kv.first<<"\t"<<kv.second<<std::endl;
    }
//    for(auto kv:ObsRecordVec)
//    {
//        std::cout<<kv.ReadName<<"\t"<<kv.ObservedInsertSize<<"\t"<<kv.MaxInsertSize<<"\t"<<kv.Weight<<std::endl;
//    }
    return 0;
}

int InsertSizeEstimator::UpdateInsertDist() {
    //if(init)
        for(auto kv:ObsRecordVec)
        {
            if(kv.ObservedInsertSize!=-1)
            {
                if (Distribution.find(kv.ObservedInsertSize) != Distribution.end())
                {
                    Distribution[kv.ObservedInsertSize] += kv.Weight;
                }
                else
                {
                    Distribution[kv.ObservedInsertSize]=kv.Weight;
                }
            }
        }
   /* else
        for(auto kv:ObsRecordVec)
        {
            Distribution[kv.ObservedInsertSize]=Distribution[kv.ObservedInsertSize]*DistWeight[kv.ObservedInsertSize];
        }*/

    return 0;
}

//
// Created by Fan Zhang on 7/3/15.
//

#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include "InsertSizeEstimator.h"

using namespace std;
void split( vector<string> & theStringVector,  /* Altered/returned value */
       const  string  & theString,
       const  string  & theDelimiter)
{
    assert( (theDelimiter.size() > 0)); // My own ASSERT macro.
    theStringVector.clear();
    size_t  start = 0, end = 0;

    while ( end != string::npos)
    {
        end = theString.find( theDelimiter, start);

        // If at end, use length=maxLength.  Else use length=end-start.
        theStringVector.push_back( theString.substr( start,
                                                     (end == string::npos) ? string::npos : end - start));

        // If at end, use start=maxSize.  Else use start=end+delimiter.
        start = (   ( end > (string::npos - theDelimiter.size()) )
                    ?  string::npos  :  end + theDelimiter.size());
    }
}

int InsertSizeEstimator::InputInsertSizeTable(const std::string& FileName) {
    std::ifstream fin(FileName);
    if(!fin.is_open()) std::cerr<<"File "<<FileName<<" Open Failed!"<<std::endl;
    std::string line;
    std::istringstream ss;
    vector<string> stringVector;
    while(getline(fin,line))
    {
//        ss.str(line);
        std::string ReadName,Chr1,Chr2,Status;
        int Max,Max2,Obs,Pos1,Pos2;
//        ss>>ReadName;
//        ss>>Max;
//        ss>>Max2;
//        ss>>Obs;
//        //getline(ss,Chr1,'\t');
//        ss>>Chr1;
//        ss>>Pos1;
//        ss>>Chr1;
//        ss>>Chr1;
//        ss>>Chr1;
//        ss>>Chr2;
//        ss>>Pos2;
//        ss>>Chr2;
//        ss>>Chr2;
//        ss>>Chr2;
//        ss>>Status;
//        ss.clear();
        split(stringVector,line,string("\t"));
        ReadName=stringVector[0];
        Obs=atoi(stringVector[3].c_str());
        Pos2=atoi(stringVector[10].c_str());
        Chr2=stringVector[9];
        Max=atoi(stringVector[1].c_str());
        Max2=atoi(stringVector[2].c_str());
        Status=stringVector[14];
        //cerr<<"come here and "<<ReadName<<"\t"<<Obs<<"\t"<<Max<<"\t"<<Max2<<"\t"<<Status<<endl;
        //continue;
        if(Max>=INSERT_LIMIT) Max=INSERT_LIMIT-1;
        if(Max2>=INSERT_LIMIT) Max2=INSERT_LIMIT-1;
        if(Obs>=INSERT_LIMIT) Obs=INSERT_LIMIT-1;
        if(Status =="Abnormal" or Status == "DiffChrom") continue;
        else if(Status == "FwdOnly")//Obs is NA
        {
            MisDistVec[Max].push_back(InsertSizeRecord(ReadName, -1, Max,1));
            MisDist[Max]++;
        }
        else if(Status == "RevOnly")
        {
            MisDistVec[Max2].push_back(InsertSizeRecord(ReadName, -1, Max2,1));
            MisDist[Max2]++;
        }
        else if(Status == "PropPair")//||Status == "PartialPair")
        {
//            MisDistVec[Max].push_back(InsertSizeRecord(ReadName, -1, Max,0.5));
//            MisDist[Max]++;
//            MisDistVec[Max2].push_back(InsertSizeRecord(ReadName, -1, Max2,0.5));
//            MisDist[Max2]++;
            ObsDistVec[Obs].push_back(InsertSizeRecord(ReadName, Obs, Max,1));
            ObsDist[Obs]++;
           // cerr<<"ObsDist["<<Obs<<"]:"<<ObsDist[Obs]<<endl;
        }
        else {//pair end info available
            //continue;
            MisDistVec[Max].push_back(InsertSizeRecord(ReadName, -1, Max,0.5));
            MisDist[Max]+=0.5;
            MisDistVec[Max2].push_back(InsertSizeRecord(ReadName, -1, Max2,0.5));
            MisDist[Max2]+=0.5;
        }
        totalPair++;
        //ObsRecordVec.push_back(InsertSizeRecord(ReadName,Obs,Max));
    }
    return 0;
}

int InsertSizeEstimator::UpdateWeight(const std::string & outputPath) {
	ofstream fout(outputPath);
    vector<double> F(1000,0.),f(1000,0.);
    vector<double> G(F),g(f);
    for (int k = 0; k <1000 ; ++k) {
        double m(0),n(0);
        m=MisDist[k];
        n=ObsDist[k];
        if(k != 0) {
            f[k] = n / (1 - G[k - 1]) * 1 / double(totalPair);
            F[k] = F[k - 1] + f[k];
        }
        else {
            f[k] = n / double(totalPair);
            F[k]=f[k];
        }

        if(k != 0) {
            g[k] = m / (1 - F[k]) * 1 / double(totalPair);
            G[k] = G[k - 1] + g[k];
        }
        else {
            g[k] = m / double(totalPair);
            G[k]=g[k];
        }

        fout << k << "\t" << f[k]<< std::endl;
    }

    fout.close();
    return 0;
}

int InsertSizeEstimator::Sort() {
//    std::sort(ObsRecordVec.begin(), ObsRecordVec.end(),[&](InsertSizeRecord a,InsertSizeRecord b){ return a.MaxInsertSize < b.MaxInsertSize;});
    return 0;
}

int InsertSizeEstimator::GetInsertDist() {

    for(auto kv:Distribution)
    {
        std::cerr<<kv.first<<"\t"<<kv.second<<std::endl;
    }
//    for(auto kv:ObsRecordVec)
//    {
//        std::cout<<kv.ReadName<<"\t"<<kv.ObservedInsertSize<<"\t"<<kv.MaxInsertSize<<"\t"<<kv.Weight<<std::endl;
//    }
    return 0;
}

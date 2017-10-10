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
#define SAM_FPD   1 // paired
#define SAM_FPP   2 // properly paired
#define SAM_FSU   4 // self-unmapped
#define SAM_FMU   8 // mate-unmapped
#define SAM_FSR  16 // self on the reverse strand
#define SAM_FMR  32 // mate on the reverse strand
#define SAM_FR1  64 // this is read one
#define SAM_FR2 128 // this is read two
#define SAM_FSC 256 // secondary alignment

int InsertSizeEstimator::InputInsertSizeTable(std::string FileName, std::string Orientation) {
    std::ifstream fin(FileName);
    if(!fin.is_open()) std::cerr<<"File "<<FileName<<" Open Failed!"<<std::endl;
    std::string line;
    std::istringstream ss;
    vector<string> stringVector;
    double count = 1.;
    int leftCount=0,rightCount=0;
    while(getline(fin,line))
    {
        std::string ReadName,Chr1,Chr2,Status,Cigar1,Cigar2;
        int Max(0),Max2(0),Obs(0),Pos1(0),Pos2(0),ReadLen1(0),ReadLen2(0),Flag1(0),Flag2(0);
        split(stringVector,line,string("\t"));
        ReadName=stringVector[0];
        Max=atoi(stringVector[1].c_str());
        Max2=atoi(stringVector[2].c_str());
        Obs=atoi(stringVector[3].c_str());
        Chr1=stringVector[4];
        Pos1=atoi(stringVector[5].c_str());
        Flag1=atoi(stringVector[6].c_str());
        ReadLen1=atoi(stringVector[7].c_str());
        Cigar1=stringVector[8];
        Chr2=stringVector[9];
        Pos2=atoi(stringVector[10].c_str());
        Flag2=atoi(stringVector[11].c_str());
        ReadLen2=atoi(stringVector[12].c_str());
        Cigar2=stringVector[13];
        Status=stringVector[14];
        if(Max>=INSERT_LIMIT) Max=INSERT_LIMIT-1;
        if(Max2>=INSERT_LIMIT) Max2=INSERT_LIMIT-1;
        if(Obs>=INSERT_LIMIT) Obs=INSERT_LIMIT-1;
        if (Status == "Abnormal" or Status == "LowQual" or Status == "NotPair" or
            Status == Orientation/*or Status == "FwdOnly" or Status == "RevOnly"*/)
            continue;
        else if(Status == "FwdOnly")//Obs is NA
        {
            MisDistVec[Max].push_back(InsertSizeRecord(ReadName, -1, Max,1));
            MisDist[Max]+=count;
        }
        else if(Status == "RevOnly")
        {
            MisDistVec[Max2].push_back(InsertSizeRecord(ReadName, -1, Max2,1));
            MisDist[Max2]+=count;
        }
        else if(Status == "PropPair")
        {
            ObsDistVec[Obs].push_back(InsertSizeRecord(ReadName, Obs, Max,1));
            ObsDist[Obs]+=count;
        }
        else if (Status == "PartialPair")
        {
            if (Cigar1.find('S') == Cigar1.npos && Cigar2.find('S') != Cigar2.npos)//read1 no S
            {
                if(Flag1&SAM_FSR) { // R mapped
                    MisDistVec[Max2].push_back(InsertSizeRecord(ReadName, -1, Max2, 1));
                    MisDist[Max2] += 1.;
                    rightCount++;
                } else // F mapped
                {
                    MisDistVec[Max].push_back(InsertSizeRecord(ReadName, -1, Max, 1));
                    MisDist[Max] += 1.;
                    leftCount++;
                }
            } else if (Cigar1.find('S') != Cigar1.npos && Cigar2.find('S') == Cigar2.npos)//read2 no S
            {
                if(Flag2&SAM_FSR) {// R mapped
                    MisDistVec[Max2].push_back(InsertSizeRecord(ReadName, -1, Max2, 1));
                    MisDist[Max2] += 1.;
                    rightCount++;
                } else// F mapped
                {
                    MisDistVec[Max].push_back(InsertSizeRecord(ReadName, -1, Max, 1));
                    MisDist[Max] += 1.;
                    leftCount++;
                }
            } else
            {
                continue;
            }
        }
        else if ( Status == "NotPair")
        {
            MisDistVec[Max].push_back(InsertSizeRecord(ReadName, -1, Max,0.5));
            MisDist[Max]+=0.5;
            MisDistVec[Max2].push_back(InsertSizeRecord(ReadName, -1, Max2,0.5));
            MisDist[Max2]+=0.5;
        }
        else {//pair end info available
            exit(EXIT_FAILURE);
            continue;
            MisDistVec[Max].push_back(InsertSizeRecord(ReadName, -1, Max,0.5));
            MisDist[Max]+=0.5;
            MisDistVec[Max2].push_back(InsertSizeRecord(ReadName, -1, Max2,0.5));
            MisDist[Max2]+=0.5;
        }
        totalPair++;
        //ObsRecordVec.push_back(InsertSizeRecord(ReadName,Obs,Max));
    }
//    cerr<<"leftCount:"<<leftCount<<"\trightCount:"<<rightCount<<std::endl;
    return 0;
}

vector<double> InsertSizeEstimator::UpdateWeight() {
    vector<double> F(2000, 0.), f(2000, 0.);
    vector<double> G(F),g(f);
    for (int k = 0; k < 2000; ++k) {
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
        //       cout << k << "\t" << f[k]<<"\tMaximalIS:"<<m<<"\tObsIS:"<<n<< endl;
    }

    return f;
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

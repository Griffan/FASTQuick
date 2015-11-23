/* The MIT License

   Copyright (c) 2009 Genome Research Ltd (GRL), 2010 Broad Institute

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Fan Zhang <fanzhang@umich.edu> */

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

    int InputInsertSizeTable(const std::string& FileName);
    int UpdateWeight(const std::string &FileName);
    int Sort();
    int GetInsertDist();
    int UpdateInsertDist();

};


#endif //INSERTSIZEESTIMATOR_INSERTSIZEESTIMATOR_H

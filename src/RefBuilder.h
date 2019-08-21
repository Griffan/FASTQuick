/*The MIT License (MIT)
Copyright (c) 2017 Fan Zhang, Hyun Min Kang
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */
/* Contact: Fan Zhang <fanzhang@umich.edu> */

#ifndef REFBUILDER_H_
#define REFBUILDER_H_

#include <faidx.h>
#include "Utility.h"
#include "../libbwa/bwtaln.h"
#include "VcfRecord.h"

class RefBuilder {
public:
//    std::vector<std::string> SeqVec;
//    std::unordered_map<std::string, uint32_t> RefTableIndex;//seq name -> index
    std::map<std::string, std::map<int, int> > VcfTable;
    std::vector<VcfRecord* > VcfVec;

    RefBuilder();

//    RefBuilder(const std::string &VcfPath, const std::string &RefPath, const std::string &NewRefPath,
//               const std::string &DBsnpPath, const std::string &MaskPath, const gap_opt_t *opt, bool reselect);

    RefBuilder(const std::string &Vcf, const std::string &Ref, const std::string &New,
                         const std::string &DBsnp, const std::string &Mask,
                         int short_len, int long_len,
                         int short_num, int long_num);

    virtual ~RefBuilder();

    bool VariantCheck(VcfRecord* VcfLine, int chrFlag);

    bool Skip(VcfRecord* VcfLine, int  chrFlag);

    int SelectMarker();

    int InputPredefinedMarker(const std::string & predefinedVcf);

    void SubstrRef(const faidx_t *seq, VcfRecord *VcfLine, std::ofstream &FGC, std::ofstream &FaOut);

    int PrepareRefSeq();

    bool IsChromInWhiteList(std::string &Chrom);

    bool IsMaxNumMarker(const std::string &Chrom, int &chrFlag,
                        bool isLong = false);//check if reach max number of marker, and also set marker type

    void IncreaseNumMarker(int chrFlag);

private:
    int nYMarker = 0, nXMarker = 0, nShortMarker = 0, nLongMarker =0;
    int maxXorYmarker =0;


    const std::string VcfPath;
    const std::string RefPath;
    const std::string NewRef;
    const std::string dbSNP;
    const std::string MaskPath;
    const faidx_t *FastaMask;

    int flank_short_len=0;
    int flank_long_len=0;
    int num_variant_short=0;
    int num_variant_long=0;

    std::set<std::string> chromWhiteList;

};

#endif /* REFBUILDER_H_ */

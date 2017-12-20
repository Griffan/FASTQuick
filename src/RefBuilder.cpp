#include <algorithm>
#include "RefBuilder.h"
#include "Utility.h"
#include "../misc/general/Error.h"
#include "VcfFileReader.h"
#include "VcfHeader.h"
#include "../libbwa/bwtaln.h"
#include <fstream>
#include <string>
#include <cctype>

using namespace std;
#define DEBUG 0

//extern string Prefix;
extern void notice(const char *, ...);

extern void warning(const char *, ...);

extern void error(const char *, ...);

static void CalculateGC(const int flank_len, const faidx_t *seq, const string &Chrom, const int Position,
                        ofstream &FGC)  {//Calculate GC content
    _GCstruct GCstruct(flank_len * 2 + 1);
    char region[1024];
    int dummy;
    for (int i = Position - flank_len, t = 0; i != Position + flank_len + 1; ++i, ++t) {
        sprintf(region, "%s:%d-%d", Chrom.c_str(), i - 50, i + 49);
        string Window(fai_fetch(seq, region, &dummy));
        int total = 0;
        for (unsigned int j = 0; j != Window.size(); ++j) {
            if (Window[j] == 'G' || Window[j] == 'C' || Window[j] == 'g' || Window[j] == 'c')
                total++;
        }
        GCstruct.GC[t] = total;
    }
    GCstruct.write(FGC);
}

bool RefBuilder::Skip(VcfRecord* VcfLine, int & chrFlag) {
    std::string chr = VcfLine->getChromStr();
    int position = VcfLine->get1BasedPosition();
    int flank_len =0;
    if(not IsChromInWhiteList(chr))//strip chr from e.g. chr11
        return true;

    if(IsMaxNumMarker(chr, chrFlag))
        return true;

    if (VcfLine->getNumRefBases() != 1||strlen(VcfLine->getAltStr()) !=1||VcfLine->getNumAlts() != 1)// filtering indel sites
        return true;

    std::string::size_type sz;     // alias of size_t
    const std::string* AFstr = VcfLine->getInfo().getString("AF");
    if(AFstr != NULL) {
        double AF = std::stod(*AFstr, &sz);
        if (AF < 0.05 or AF > 0.95) return true;
    }
    else return true;

    if(chrFlag == 1)
        flank_len = flank_long_len;
    else
        flank_len = flank_short_len;

    //ensure no overlapping regions
    if(VcfTable.find(chr)!=VcfTable.end())
    {
        auto low_iter = VcfTable[chr].upper_bound(position);
        if(low_iter != VcfTable[chr].begin() and low_iter != VcfTable[chr].end()) {
            auto up_iter = low_iter--;
            if (abs(position - VcfVec[low_iter->second]->get1BasedPosition()) < 2 * flank_len or
                abs(position - VcfVec[up_iter->second]->get1BasedPosition()) < 2 * flank_len)
                return true;
        }
    }
    //
    int dummy_t;
    char region[1024];
    sprintf(region, "%s:%d-%d", chr.c_str(), position - flank_len, position + flank_len);
    if (MaskPath != "Empty") {
        string MaskSeq(fai_fetch(FastaMask, region, &dummy_t));
        double n = std::count(MaskSeq.begin(), MaskSeq.end(), 'P');
        if (n < MaskSeq.size())
            return true;
    }
    return false;
}

bool RefBuilder::IsChromInWhiteList(std::string &Chrom) {
    std::transform(Chrom.begin(), Chrom.end(), Chrom.begin(), ::toupper);
    size_t tmpPos = Chrom.find("chr");
    if (tmpPos != std::string::npos) {
        Chrom = Chrom.substr(tmpPos + 3, Chrom.size() - 3);//strip chr
    }
    return (chromWhiteList.find(Chrom) != chromWhiteList.end());
}


RefBuilder::RefBuilder(const std::string &Vcf, const std::string &Ref, const std::string &New,
                       const std::string &DBsnp, const std::string &Mask,
                        const int short_len, const int long_len,
                       const int short_num, const int long_num):
        VcfPath(Vcf), RefPath(Ref), NewRef(New), dbSNP(DBsnp),MaskPath(Mask),
        flank_short_len(short_len),flank_long_len(long_len),num_variant_short(short_num),num_variant_long(long_num)
{
    if (num_variant_short >= 100000)
        maxXorYmarker = 3000;
    else if (num_variant_short >= 10000)
        maxXorYmarker = 300;
    else
        maxXorYmarker = 100;

    if (MaskPath != "Empty") {
        FastaMask = fai_load(MaskPath.c_str());
        notice("Loading Mask fai file done!\n");
    }

    chromWhiteList =
            {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
             "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
             "21", "22", "X", "Y"};

}

RefBuilder::RefBuilder() {
}

//RefBuilder::RefBuilder(const string &VcfPath, const string &RefPath, const string &NewRef, const string &DBsnpPath,
//                       const string &MaskPath, const gap_opt_t *opt,
//                       bool reselect = false) {
//
//    faidx_t *FastaMask = 0;
//    if (MaskPath != "Empty") {
//        FastaMask = fai_load(MaskPath.c_str());
//        notice("Loading Mask fai file done!\n");
//    }
//
//    string SelectedSite = NewRef + ".SelectedSite.vcf";
//    InputFile FoutHapMapSelectedSite(SelectedSite.c_str(), "w");
//    string BedPath = NewRef + ".bed";
//    ofstream BedFile(BedPath);
//
//    //read in vcf, hapmap3 sites
//    VcfHeader header;
//    VcfFileReader reader;
//    reader.open(VcfPath.c_str(), header);
//    header.write(&FoutHapMapSelectedSite);
//    char region[128];
//
//    unsigned int nmarker(0);
//    unsigned int totalMarker(0);
//
//    srand(time(NULL));
//
//    string Chrom;
//    int Position;
//
//    std::set<std::string> autoRegionWL =
//            {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
//             "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
//             "21", "22"};
//    while (!reader.isEOF()) {
//        if (nmarker >= opt->num_variant_short)// for short region
//        {
//            break;
//        }
//        VcfRecord* VcfLine=new VcfRecord;
//        reader.readRecord(*VcfLine);
//        if (VcfLine->getNumRefBases() != 1||strlen(VcfLine->getAltStr()) !=1||VcfLine->getNumAlts() != 1 )// filtering indel sites
//            continue;
//
//        std::string::size_type sz;     // alias of size_t
//        double AF = std::stod(*VcfLine->getInfo().getString("AF"),&sz);
//        if(AF<0.05 or AF >0.95) continue;
//
//        Chrom = VcfLine->getChromStr();
//        Position = VcfLine->get1BasedPosition();
//
//        sprintf(region, "%s:%d-%d", Chrom.c_str(), Position - opt->flank_len, Position + opt->flank_len);
//
//        if (Skip(Chrom, Position, region, MaskPath, FastaMask, autoRegionWL, opt->flank_len)) continue;
//
//        VcfTable[Chrom][Position]=totalMarker+nmarker;
//        VcfVec.push_back(VcfLine);
//
//        nmarker++;
//    }
//    totalMarker+=nmarker;
//    // below is for long ref
//    nmarker = 0;
//    while (!reader.isEOF()) {
//        if (nmarker >= opt->num_variant_long)// for long region
//        {
//            break;
//        }
//        VcfRecord* VcfLine=new VcfRecord;
//        reader.readRecord(*VcfLine);
//        if (VcfLine->getNumRefBases() != 1||strlen(VcfLine->getAltStr()) !=1||VcfLine->getNumAlts() != 1)// filtering indel sites
//            continue;
//
//        std::string::size_type sz;     // alias of size_t
//        double AF = std::stod(*VcfLine->getInfo().getString("AF"),&sz);
//        if(AF<0.05 or AF >0.95) continue;
//
//        Chrom = VcfLine->getChromStr();
//        Position = VcfLine->get1BasedPosition();
//        sprintf(region, "%s:%d-%d", Chrom.c_str(), Position - opt->flank_long_len, Position + opt->flank_long_len);
//
//        VcfLine->setID((std::string(VcfLine->getIDStr())+"|L").c_str());
//
//        if (Skip(Chrom, Position, region, MaskPath, FastaMask, autoRegionWL, opt->flank_long_len)) continue;
//        VcfTable[Chrom][Position]=totalMarker+nmarker;
//        VcfVec.push_back(VcfLine);
//        nmarker++;
//    }
//    totalMarker+=nmarker;
//    int max_XorYmarker(0);
//    if (opt->num_variant_short >= 100000)
//        max_XorYmarker = 3000;
//    else if (opt->num_variant_short >= 10000)
//        max_XorYmarker = 300;
//    else
//        max_XorYmarker = 100;
//    // choosing chrX and chrY
//    int Ynmarker = 0, Xnmarker = 0, chr_flag(-1);
//
//    while (!reader.isEOF()) {//combine two in case not enough X/Y markers
//        if (Ynmarker >= max_XorYmarker && Xnmarker >= max_XorYmarker)
//        {
//            break;
//        }
//        VcfRecord* VcfLine=new VcfRecord;
//        reader.readRecord(*VcfLine);
//        if (VcfLine->getNumRefBases() != 1||strlen(VcfLine->getAltStr()) !=1||VcfLine->getNumAlts() != 1)// filtering indel sites
//            continue;
//
//        std::string::size_type sz;     // alias of size_t
//        double AF = std::stod(*VcfLine->getInfo().getString("AF"),&sz);
//        if(AF<0.05 or AF >0.95) continue;
//
//        Chrom=VcfLine->getChromStr();
//        if (Chrom == "X" || Chrom == "x" || Chrom == "chrX" || Chrom == "chrx") {
//            if (Xnmarker >= max_XorYmarker)
//                continue;
//            chr_flag = 0;
//        } else if (Chrom == "Y" || Chrom == "y" || Chrom == "chrY" || Chrom == "chry") {
//            if (Ynmarker >= max_XorYmarker)
//                continue;
//            chr_flag = 1;
//        } else
//            continue;
//
//
//        Position = VcfLine->get1BasedPosition();
//        int dummy;
//        sprintf(region, "%s:%d-%d", Chrom.c_str(), Position - opt->flank_len, Position + opt->flank_len);
//
//        //ensure no overlapping regions
//        if(VcfTable.find(Chrom)!=VcfTable.end())
//        {
//            auto low_iter = VcfTable[Chrom].upper_bound(Position);
//            if(low_iter != VcfTable[Chrom].begin() and low_iter != VcfTable[Chrom].end()) {
//                auto up_iter = low_iter--;
//                if (abs(Position - VcfVec[low_iter->second]->get1BasedPosition()) < 2 * opt->flank_len or
//                    abs(Position - VcfVec[up_iter->second]->get1BasedPosition()) < 2 * opt->flank_len)
//                    continue;
//            }
//        }
//        if (MaskPath != "Empty") {
//            string MaskSeq(fai_fetch(FastaMask, region, &dummy));
//            size_t n = std::count(MaskSeq.begin(), MaskSeq.end(), 'P');
//            if (n < MaskSeq.size())
//                continue;
//        }
//
//        VcfTable[Chrom][Position]=totalMarker+Xnmarker+Ynmarker;
//        VcfVec.push_back(VcfLine);
//
//
//        if (chr_flag == 0)
//            Xnmarker++;
//        else
//            Ynmarker++;
//    }
//
//
//    //output vcf and bed files
//    int flank_len=0;
//    for(auto kv:VcfTable)
//    {
//        for (auto pq:kv.second) {
//            if (!VcfVec[pq.second]->write(&FoutHapMapSelectedSite, 1)) {
//                warning("Writing retained sites failed!\n");
//                exit(EXIT_FAILURE);
//            }
//            if(std::string(VcfVec[pq.second]->getIDStr()).back()=='L')
//                flank_len=opt->flank_long_len;
//            else
//                flank_len=opt->flank_len;
//            BedFile<<kv.first<<"\t"<<pq.first - flank_len<<"\t"<<pq.first + flank_len<<endl;
//            delete VcfVec[pq.second];
//        }
//    }
//    BedFile.close();
//    FoutHapMapSelectedSite.ifclose();
//
//    reader.open(DBsnpPath.c_str(), header);
//    reader.close();
//
//    InputFile FoutdbSNPSelectedSite(std::string(NewRef+".dbSNP.subset.vcf").c_str(), "w");
//    header.write(&FoutdbSNPSelectedSite);
//    FoutdbSNPSelectedSite.ifclose();
//
//    char cmdline[2048];
//    //subset dbsnp
//    sprintf(cmdline, "sort -k1,1 -k2,2n %s|tabix  -R  - %s  >> %s.dbSNP.subset.vcf", BedPath.c_str(),DBsnpPath.c_str(),
//             NewRef.c_str());
//    int ret = system(cmdline);
//    if (ret != 0) {
//        warning("Building dbsnp subset.vcf failed!\n");
//        exit(EXIT_FAILURE);
//    }
//}

bool RefBuilder::IsMaxNumMarker(const std::string& Chrom, int& chrFlag)//use this function after checking whiteList
{
    if (Chrom == "X") {/*0:short;1:long;2:Y;3:X*/;
        if (nXMarker >= maxXorYmarker)
            return true;
        chrFlag = 3;
    } else if (Chrom == "Y") {
        if (nYMarker >= maxXorYmarker)
            return true;
        chrFlag = 2;
    }else { //must be from chr1-chr22
        if(nLongMarker < num_variant_long)
        {
            chrFlag = 1;
        }
        else if (nShortMarker < num_variant_short)//must be from chr1-chr22
        {
            chrFlag = 0;
        } else
            return true;
    }
    return false;
}

void RefBuilder::IncreaseNumMarker(int chrFlag)
{
    switch(chrFlag){/*0:short;1:long;2:Y;3:X*/;
        case 0:
            nShortMarker++;
            break;
        case 1:
            nLongMarker++;
            break;
        case 2:
            nYMarker++;
            break;
        case 3:
            nXMarker++;
            break;
    }
}

int RefBuilder::SelectMarker()
{
    string SelectedSite = NewRef + ".SelectedSite.vcf";
    InputFile FoutHapMapSelectedSite(SelectedSite.c_str(), "w");
    string BedPath = NewRef + ".bed";
    ofstream BedFile(BedPath);

    //read in vcf, hapmap3 sites
    VcfHeader header;
    VcfFileReader reader;
    reader.open(VcfPath.c_str(), header);
    header.write(&FoutHapMapSelectedSite);

    string Chrom;
    int Position;

    //begin select
    while (!reader.isEOF()) {
        int chrFlag(-1)/*0:short;1:long;2:Y;3:X*/;
        VcfRecord *VcfLine = new VcfRecord;
        reader.readRecord(*VcfLine);

        if(Skip(VcfLine,chrFlag))
            continue;

        Position = VcfLine->get1BasedPosition();
        Chrom = VcfLine->getChromStr();

        if(chrFlag==1)
        VcfLine->setID((std::string(VcfLine->getIDStr())+"|L").c_str());


        VcfTable[Chrom][Position]=nShortMarker+nLongMarker+nXMarker+nYMarker;
        VcfVec.push_back(VcfLine);

        IncreaseNumMarker(chrFlag);
    }

    if(nShortMarker+nLongMarker+nXMarker < maxXorYmarker + num_variant_long + num_variant_short)//not enough essential markers
    {
        warning("there are insufficient candidate markers in %s",VcfPath.c_str());
    }
    //output vcf and bed files
    int flank_len=0;
    for(auto kv:VcfTable)
    {
        for (auto pq:kv.second) {
            if (!VcfVec[pq.second]->write(&FoutHapMapSelectedSite, 1)) {
                warning("Writing retained sites failed!\n");
                exit(EXIT_FAILURE);
            }
            if(std::string(VcfVec[pq.second]->getIDStr()).back()=='L')
                flank_len=flank_long_len;
            else
                flank_len=flank_short_len;
            BedFile<<kv.first<<"\t"<<pq.first - flank_len<<"\t"<<pq.first + flank_len<<endl;
//            delete VcfVec[pq.second];
        }
    }
    BedFile.close();
    FoutHapMapSelectedSite.ifclose();

    reader.open(dbSNP.c_str(), header);
    reader.close();

    //create header for dbSNP subset vcf file
    InputFile FoutdbSNPSelectedSite(std::string(NewRef+".dbSNP.subset.vcf").c_str(), "w");
    header.write(&FoutdbSNPSelectedSite);
    FoutdbSNPSelectedSite.ifclose();

    char cmdline[2048];
    //subset dbsnp
    sprintf(cmdline, "sort -k1,1 -k2,2n %s|tabix  -R  - %s  >> %s.dbSNP.subset.vcf",
            BedPath.c_str(),dbSNP.c_str(), NewRef.c_str());
    int ret = system(cmdline);
    if (ret != 0) {
        warning("Building dbsnp subset.vcf failed!\n");
        exit(EXIT_FAILURE);
    }
    return 0;
}

int RefBuilder::InputPredefinedMarker(const std::string & predefinedVcf)
{
    string SelectedSite = NewRef + ".SelectedSite.vcf";
    InputFile FoutHapMapSelectedSite(SelectedSite.c_str(), "w");
    string BedPath = NewRef + ".bed";
    ofstream BedFile(BedPath);

    //read in predefined vcf
    VcfHeader header;
    VcfFileReader reader;
    reader.open(predefinedVcf.c_str(), header);
    header.write(&FoutHapMapSelectedSite);

    string Chrom;
    int Position;

    //begin select
    while (!reader.isEOF()) {
        int chrFlag(-1)/*0:short;1:long;2:Y;3:X*/;
        VcfRecord* VcfLine= new VcfRecord;
        reader.readRecord(*VcfLine);

        if(Skip(VcfLine,chrFlag))
            continue;

        Position = VcfLine->get1BasedPosition();
        Chrom = VcfLine->getChromStr();

        if(chrFlag==1)
            VcfLine->setID((std::string(VcfLine->getIDStr())+"|L").c_str());


        VcfTable[Chrom][Position]=nShortMarker+nLongMarker+nXMarker+nYMarker;
        VcfVec.push_back(VcfLine);

        IncreaseNumMarker(chrFlag);
    }

    if(nShortMarker+nLongMarker+nXMarker< maxXorYmarker + num_variant_long + num_variant_short)
    {
        warning("there are insufficient candidate markers in %s",predefinedVcf.c_str());

    }
    //output vcf and bed files
    int flank_len=0;
    for(auto kv:VcfTable)
    {
        for (auto pq:kv.second) {
            if (!VcfVec[pq.second]->write(&FoutHapMapSelectedSite, 1)) {
                warning("Writing retained sites failed!\n");
                exit(EXIT_FAILURE);
            }
            if(std::string(VcfVec[pq.second]->getIDStr()).back()=='L')
                flank_len=flank_long_len;
            else
                flank_len=flank_short_len;
            BedFile<<kv.first<<"\t"<<pq.first - flank_len<<"\t"<<pq.first + flank_len<<endl;
//            delete VcfVec[pq.second];
        }
    }
    BedFile.close();
    FoutHapMapSelectedSite.ifclose();

    reader.open(dbSNP.c_str(), header);
    reader.close();

    //create header for dbSNP subset vcf file
    InputFile FoutdbSNPSelectedSite(std::string(NewRef+".dbSNP.subset.vcf").c_str(), "w");
    header.write(&FoutdbSNPSelectedSite);
    FoutdbSNPSelectedSite.ifclose();

    char cmdline[2048];
    //subset dbsnp
    sprintf(cmdline, "sort -k1,1 -k2,2n %s|tabix  -R  - %s  >> %s.dbSNP.subset.vcf", BedPath.c_str(),dbSNP.c_str(),
            NewRef.c_str());
    int ret = system(cmdline);
    if (ret != 0) {
        warning("Building dbsnp subset.vcf failed!\n");
        exit(EXIT_FAILURE);
    }
    return 0;
}

void
RefBuilder::SubstrRef(const faidx_t *seq, VcfRecord *VcfLine, std::ofstream &FGC, std::ofstream &FaOut) {
    int flank_len;
    int dummy;
    char region[1024];
    char newChrName[1024];

    if (std::string(VcfLine->getIDStr()).find('L')!=std::string::npos)//Long region
    {
        flank_len=flank_long_len;
        sprintf(newChrName, ">%s:%d@%s/%s|L", VcfLine->getChromStr(), VcfLine->get1BasedPosition(), VcfLine->getRefStr(), VcfLine->getAltStr());
    } else{
        flank_len=flank_short_len;
        sprintf(newChrName, ">%s:%d@%s/%s", VcfLine->getChromStr(), VcfLine->get1BasedPosition(), VcfLine->getRefStr(), VcfLine->getAltStr());
    }
    sprintf(region, "%s:%d-%d", VcfLine->getChromStr(), VcfLine->get1BasedPosition() - flank_len,
            VcfLine->get1BasedPosition() + flank_len);
    std::string FetchedSeq(fai_fetch(seq, region, &dummy));

    CalculateGC(flank_len, seq, VcfLine->getChromStr(), VcfLine->get1BasedPosition(), FGC);

//    SeqVec.push_back(FetchedSeq.substr(0, flank_len) + string(VcfLine.getRefStr()) +
//                     FetchedSeq.substr(flank_len + 1, flank_len));
//    sprintf(region, "%s:%d@%s/%s", VcfLine.getChromStr(), VcfLine.get1BasedPosition(), VcfLine.getRefStr(), VcfLine.getAltStr());
//    RefTableIndex.insert(make_pair(string(region), nseqs));
//    nseqs++;
    FaOut<<newChrName<<std::endl;
    FaOut << FetchedSeq.substr(0, flank_len) + std::string(VcfLine->getRefStr()) +
             FetchedSeq.substr(flank_len + 1, flank_len) << std::endl;
}
int RefBuilder::PrepareRefSeq()
{

    faidx_t *seq = fai_load(RefPath.c_str());
    notice("Loading Ref fai file done!\n");
    string GCpath = NewRef + ".gc";
    ofstream FGC(GCpath, ios_base::binary);
    ofstream FaOut(NewRef);
    for(auto& kv: VcfTable)
    {
        for(auto& kv2: kv.second)//each marker in VcfTable
        {
            SubstrRef(seq, VcfVec[kv2.second], FGC, FaOut);
            delete VcfVec[kv2.second];
        }
    }
    FaOut.close();
    FGC.close();
    return 0;
}
RefBuilder::~RefBuilder() {
//    SeqVec.clear();
//    RefTableIndex.clear();
}


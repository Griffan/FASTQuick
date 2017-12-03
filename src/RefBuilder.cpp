#include <algorithm>
#include "RefBuilder.h"
#include "Utility.h"
#include "../misc/general/Error.h"
#include "VcfFileReader.h"
#include "VcfHeader.h"
#include "../libbwa/bwtaln.h"
#include <fstream>
#include <string>

using namespace std;
#define DEBUG 0

//extern string Prefix;
extern void notice(const char *, ...);

extern void warning(const char *, ...);

extern void error(const char *, ...);

static void CalculateGC(const int flank_len, const faidx_t *seq, char *region, const string &Chrom, const int Position,
                        ofstream &FGC, int &dummy)  {//Calculate GC content
    _GCstruct GCstruct(flank_len * 2 + 1);
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

bool RefBuilder::Skip(const string &Chrom, const int Position, const string &last_chr, const int last_pos, char *region,
                 const string &MaskPath, const faidx_t *FastaMask, const std::set<std::string> &chromWhiteList, int flank_len) {
    std::string data = Chrom;
    std::transform(data.begin(), data.end(), data.begin(), ::tolower);
    size_t tmpPos = data.find("chr");
    if (tmpPos != std::string::npos) {
        data = data.substr(tmpPos + 3, data.size() - 3);//strip chr
    }
    if (chromWhiteList.find(data) == chromWhiteList.end())
        return true;
//    if (Chrom == last_chr && abs(Position - last_pos) < 2*flank_len)//ensure no overlapping regions
//        return true;
    //ensure no overlapping regions
    if(VcfTable.find(Chrom)!=VcfTable.end())
    {
        auto low_iter = VcfTable[Chrom].upper_bound(Position);
        if(low_iter != VcfTable[Chrom].begin() and low_iter != VcfTable[Chrom].end()) {
            auto up_iter = low_iter--;
            if (abs(Position - VcfVec[low_iter->second]->get1BasedPosition()) < 2 * flank_len or
                abs(Position - VcfVec[up_iter->second]->get1BasedPosition()) < 2 * flank_len)
                return true;
        }
    }
    //
    int dummy_t;
    if (MaskPath != "Empty") {
        string MaskSeq(fai_fetch(FastaMask, region, &dummy_t));
        double n = std::count(MaskSeq.begin(), MaskSeq.end(), 'P');
        if (n < MaskSeq.size())
            return true;
    }
    return false;
}


RefBuilder::RefBuilder() {
}

RefBuilder::RefBuilder(const string &VcfPath, const string &RefPath, const string &NewRef, const string &DBsnpPath,
                       const string &MaskPath, const gap_opt_t *opt,
                       bool reselect = false) {
    notice("Initialization of RefBwt...\n");
    //read in ref.fa and ref.fai
    //string RefFaiPath=RefPath+".fai";
    faidx_t *seq = fai_load(RefPath.c_str());
    notice("Loading Ref fai file done!\n");

    faidx_t *FastaMask = 0;
    if (MaskPath != "Empty") {
        FastaMask = fai_load(MaskPath.c_str());
        notice("Loading Mask fai file done!\n");
    }
    //read in vcf, hapmap3 sites
    VcfHeader header;
    VcfFileReader reader;

    string SelectedSite = NewRef + ".SelectedSite.vcf";
    InputFile FoutHapMapSelectedSite(SelectedSite.c_str(), "w");
    string GCpath = NewRef + ".gc";
    ofstream FGC(GCpath, ios_base::binary);
    string BedPath = NewRef + ".bed";
    ofstream BedFile(BedPath);
    //FGC.write((char*)opt->num_variant_short,sizeof(int));
    //_GCstruct * GCstruct=new _GCstruct [opt->num_variant_long+opt->num_variant_short];

    reader.open(VcfPath.c_str(), header);
    header.write(&FoutHapMapSelectedSite);
    char region[128];
    unsigned int nseqs(0);
    unsigned int nmarker(0);
    unsigned int totalMarker(0);
    int last_pos = 0;
    string last_chr;
    srand(time(NULL));

    int dummy;
//    VcfRecord VcfLine;
    string Chrom;
    int Position;
    /* generate secret number between 1 and 10: */
    //double iSelect = (rand() % 1000 + 1)/1000.0;
    //int num_so_far=0;
    std::set<std::string> autoRegionWL =
            {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
             "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
             "21", "22"};
    while (!reader.isEOF()) {
        if (nmarker >= opt->num_variant_short)// for short region
        {
            break;
        }
        VcfRecord* VcfLine=new VcfRecord;
        reader.readRecord(*VcfLine);
        if (VcfLine->getNumRefBases() != 1||strlen(VcfLine->getAltStr()) !=1||VcfLine->getNumAlts() != 1 )// filtering indel sites
            continue;

        std::string::size_type sz;     // alias of size_t
        double AF = std::stod(*VcfLine->getInfo().getString("AF"),&sz);
        if(AF<0.05 or AF >0.95) continue;

        Chrom = VcfLine->getChromStr();
        Position = VcfLine->get1BasedPosition();

//        std::cerr<<Chrom<<"\t"<<Position<<"\t"<<AF<<std::endl;

        sprintf(region, "%s:%d-%d", Chrom.c_str(), Position - opt->flank_len, Position + opt->flank_len);

        if (Skip(Chrom, Position, last_chr, last_pos, region, MaskPath, FastaMask, autoRegionWL, opt->flank_len)) continue;

//        if (!VcfLine.write(&FoutHapMapSelectedSite, true)) {
//            warning("Writing retained sites failed!\n");
//            exit(EXIT_FAILURE);
//        }
        VcfTable[Chrom][Position]=totalMarker+nmarker;
        VcfVec.push_back(VcfLine);

        string FetchedSeq(fai_fetch(seq, region, &dummy));

        CalculateGC(opt->flank_len, seq, region, Chrom, Position, FGC, dummy);

        if (DEBUG) {
            string RefAllele(VcfLine->getAlleles(0));
            string RefSeq = FetchedSeq.substr(0, opt->flank_len) + RefAllele +
                            FetchedSeq.substr(opt->flank_len + 1, opt->flank_len);
            if (RefSeq != FetchedSeq) {
                cerr << "WARNING:Coordinate problem!!!!!!" << endl << "Number of Alt:" << VcfLine->getNumAlts() << endl;
                //cerr<<"Region: "<<Coord<<endl;
                cerr << VcfLine->getAlleles(0) << endl;
                cerr << VcfLine->getAlleles(1) << endl;
                cerr << RefSeq << endl;
                cerr << FetchedSeq << endl;
                exit(EXIT_FAILURE);
            }
        }

        SeqVec.push_back(FetchedSeq.substr(0, opt->flank_len) + string(VcfLine->getRefStr()) +
                         FetchedSeq.substr(opt->flank_len + 1, opt->flank_len));
        sprintf(region, "%s:%d@%s/%s", Chrom.c_str(), Position, VcfLine->getRefStr(), VcfLine->getAltStr());
        RefTableIndex.insert(make_pair(string(region), nseqs));

//        sprintf(region, "%s\t%d\t%d", Chrom.c_str(), Position - opt->flank_len, Position + opt->flank_len);
//        BedFile << region << endl;

        nseqs++;
        nmarker++;
        last_pos = Position;
        last_chr = Chrom;
    }
    totalMarker+=nmarker;
    // below is for long ref
    nmarker = 0;
    while (!reader.isEOF()) {
        if (nmarker >= opt->num_variant_long)// for long region
        {
            break;
        }
        VcfRecord* VcfLine=new VcfRecord;
        reader.readRecord(*VcfLine);
        if (VcfLine->getNumRefBases() != 1||strlen(VcfLine->getAltStr()) !=1||VcfLine->getNumAlts() != 1)// filtering indel sites
            continue;

        std::string::size_type sz;     // alias of size_t
        double AF = std::stod(*VcfLine->getInfo().getString("AF"),&sz);
        if(AF<0.05 or AF >0.95) continue;

        Chrom = VcfLine->getChromStr();
        Position = VcfLine->get1BasedPosition();
        sprintf(region, "%s:%d-%d", Chrom.c_str(), Position - opt->flank_long_len, Position + opt->flank_long_len);

        VcfLine->setID((std::string(VcfLine->getIDStr())+"|L").c_str());

        if (Skip(Chrom, Position, last_chr, last_pos, region, MaskPath, FastaMask, autoRegionWL, opt->flank_long_len)) continue;
//
//        if (!VcfLine.write(&FoutHapMapSelectedSite, true)) {
//            warning("Writing retained sites failed!\n");
//            exit(EXIT_FAILURE);
//        }

        VcfTable[Chrom][Position]=totalMarker+nmarker;
        VcfVec.push_back(VcfLine);

        string FetchedSeq(fai_fetch(seq, region, &dummy));

        CalculateGC(opt->flank_long_len, seq, region, Chrom, Position, FGC, dummy);

        SeqVec.push_back(FetchedSeq.substr(0, opt->flank_long_len) + string(VcfLine->getRefStr()) +
                         FetchedSeq.substr(opt->flank_long_len + 1, opt->flank_long_len));

        sprintf(region, "%s:%d@%s/%s|L", Chrom.c_str(), Position, VcfLine->getRefStr(), VcfLine->getAltStr());
        RefTableIndex.insert(make_pair(string(region), nseqs));

//        sprintf(region, "%s\t%d\t%d", Chrom.c_str(), Position - opt->flank_long_len, Position + opt->flank_long_len);
//        BedFile << region << endl;

        nseqs++;
        nmarker++;
        last_pos = Position;
        last_chr = Chrom;

    }
    totalMarker+=nmarker;
    int max_XorYmarker(0);
    if (opt->num_variant_short >= 100000)
        max_XorYmarker = 3000;
    else if (opt->num_variant_short >= 10000)
        max_XorYmarker = 300;
    else
        max_XorYmarker = 100;
    // choosing chrX and chrY
    int Ynmarker = 0, Xnmarker = 0, chr_flag(-1);

    while (!reader.isEOF()) {//combine two in case not enough X/Y markers
        if (Ynmarker >= max_XorYmarker && Xnmarker >= max_XorYmarker)
        {
            break;
        }
        VcfRecord* VcfLine=new VcfRecord;
        reader.readRecord(*VcfLine);
        if (VcfLine->getNumRefBases() != 1||strlen(VcfLine->getAltStr()) !=1||VcfLine->getNumAlts() != 1)// filtering indel sites
            continue;

        std::string::size_type sz;     // alias of size_t
        double AF = std::stod(*VcfLine->getInfo().getString("AF"),&sz);
        if(AF<0.05 or AF >0.95) continue;

        Chrom=VcfLine->getChromStr();
        if (Chrom == "X" || Chrom == "x" || Chrom == "chrX" || Chrom == "chrx") {
            if (Xnmarker >= max_XorYmarker)
                continue;
            chr_flag = 0;
        } else if (Chrom == "Y" || Chrom == "y" || Chrom == "chrY" || Chrom == "chry") {
            if (Ynmarker >= max_XorYmarker)
                continue;
            chr_flag = 1;
        } else
            continue;


        Position = VcfLine->get1BasedPosition();
        int dummy;
        sprintf(region, "%s:%d-%d", Chrom.c_str(), Position - opt->flank_len, Position + opt->flank_len);
//        if (Chrom == last_chr && abs(Position - last_pos) < opt->flank_len * 2)
//            continue;
        //ensure no overlapping regions
        if(VcfTable.find(Chrom)!=VcfTable.end())
        {
            auto low_iter = VcfTable[Chrom].upper_bound(Position);
            if(low_iter != VcfTable[Chrom].begin() and low_iter != VcfTable[Chrom].end()) {
                auto up_iter = low_iter--;
                if (abs(Position - VcfVec[low_iter->second]->get1BasedPosition()) < 2 * opt->flank_len or
                    abs(Position - VcfVec[up_iter->second]->get1BasedPosition()) < 2 * opt->flank_len)
                    continue;
            }
        }
        if (MaskPath != "Empty") {
            string MaskSeq(fai_fetch(FastaMask, region, &dummy));
            size_t n = std::count(MaskSeq.begin(), MaskSeq.end(), 'P');
            if (n < MaskSeq.size())
                continue;
//            if (!VcfLine.write(&FoutHapMapSelectedSite, 1)) {
//                warning("Writing retained sites failed!\n");
//                exit(EXIT_FAILURE);
//            }
        }

        VcfTable[Chrom][Position]=totalMarker+Xnmarker+Ynmarker;
        VcfVec.push_back(VcfLine);

        string FetchedSeq(fai_fetch(seq, region, &dummy));

        CalculateGC(opt->flank_len, seq, region, Chrom, Position, FGC, dummy);

        SeqVec.push_back(FetchedSeq.substr(0, opt->flank_len) + string(VcfLine->getRefStr()) + FetchedSeq.substr(opt->flank_len + 1, opt->flank_len));
        sprintf(region, "%s:%d@%s/%s", Chrom.c_str(), Position, VcfLine->getRefStr(), VcfLine->getAltStr());
        RefTableIndex.insert(make_pair(string(region), nseqs));
        nseqs++;

//        sprintf(region, "%s\t%d\t%d", Chrom.c_str(), Position - opt->flank_len, Position + opt->flank_len);
//        BedFile << region << endl;

        if (chr_flag == 0)
            Xnmarker++;
        else
            Ynmarker++;

        last_pos = Position;
        last_chr = Chrom;

    }
    //cerr<<"the total gc sites:"<<num_so_far<<endl;
    //FGC.write((char*)GCstruct,(opt->num_variant_long*(4*opt->flank_len+1)+opt->num_variant_short*(2*opt->flank_len+1))*sizeof(_GCstruct));
    //int total=(opt->num_variant_long*(4*opt->flank_len+1)+opt->num_variant_short*(2*opt->flank_len+1))*sizeof(_GCstruct);
    //FGC.write((char*) &nmarker, sizeof( unsigned int));

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
                flank_len=opt->flank_long_len;
            else
                flank_len=opt->flank_len;
            BedFile<<kv.first<<"\t"<<pq.first - flank_len<<"\t"<<pq.first + flank_len<<endl;
            delete VcfVec[pq.second];
        }
    }
    FGC.close();
    BedFile.close();
    FoutHapMapSelectedSite.ifclose();


    reader.open(DBsnpPath.c_str(), header);
    reader.close();

    InputFile FoutdbSNPSelectedSite(std::string(NewRef+".dbSNP.subset.vcf").c_str(), "w");
    header.write(&FoutdbSNPSelectedSite);
    FoutdbSNPSelectedSite.ifclose();

    char cmdline[2048];
    //subset dbsnp
    sprintf(cmdline, "sort -k1,1 -k2,2n %s|tabix  -R  - %s  >> %s.dbSNP.subset.vcf", BedPath.c_str(),DBsnpPath.c_str(),
             NewRef.c_str());
    int ret = system(cmdline);
    if (ret != 0) {
        warning("Building dbsnp subset.vcf failed!\n");
        exit(EXIT_FAILURE);
    }
//
//    sprintf(cmdline,
//            "(grep ^# %s.SelectedSite.vcf && grep -v  ^# %s.SelectedSite.vcf|sort -k1,1 -k2,2n) >%s.SelectedSite.vcf.tmp",
//            NewRef.c_str(), NewRef.c_str(), NewRef.c_str());
//    if (system(cmdline) != 0) {
//        warning("Call command line:\n%s\nfailed!\n", cmdline);
//        exit(EXIT_FAILURE);
//    }
//    sprintf(cmdline, "mv %s.SelectedSite.vcf.tmp %s.SelectedSite.vcf", NewRef.c_str(), NewRef.c_str());
//    if (system(cmdline) != 0) {
//        warning("Call command line:\n%s\nfailed!\n", cmdline);
//        exit(EXIT_FAILURE);
//    }
//    sprintf(cmdline, "bgzip -f %s.SelectedSite.vcf", NewRef.c_str());
//    if (system(cmdline) != 0) {
//        warning("Call command line:\n%s\nfailed!\n", cmdline);
//        exit(EXIT_FAILURE);
//    }
//    sprintf(cmdline, "tabix -pvcf %s.SelectedSite.vcf.gz", NewRef.c_str());
//    if (system(cmdline) != 0) {
//        warning("Call command line:\n%s\nfailed!\n", cmdline);
//        exit(EXIT_FAILURE);
//    }
}

RefBuilder::~RefBuilder() {
    SeqVec.clear();
    RefTableIndex.clear();
}


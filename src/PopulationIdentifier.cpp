#include "PopulationIdentifier.h"
#include <iostream>
#include <algorithm>
#include "InputFile.h"

#ifdef ARMADILLO
#include <armadillo>
using namespace arma;

/*Functions for GenotypeMatrix*/
GenotypeMatrix::GenotypeMatrix(){}
GenotypeMatrix::GenotypeMatrix(const char* vcfFile, bool siteOnly, bool findBest, /*std::vector<std::string>& subsetInds,*/ double minAF, double minCallRate)
{
	// open a VCF file
	std::vector<std::string> subsetInds;
	tvcf.infoAF = (siteOnly || ((subsetInds.size() <= 1) && !findBest));

	int i = 0, m = 0, M = 0;
	double af, maf;
	int unit = 10000;

	std::set<std::string> idset;
	for (i = 0; i < (int)subsetInds.size(); ++i) {
		if (idset.find(subsetInds[i]) != idset.end()) {
			error("ERROR: Duplicate individual ID %s", subsetInds[i].c_str());
		}
		idset.insert(subsetInds[i]);
	}
	tvcf.ignoreEmpty = findBest;
	tvcf.load(vcfFile, NULL, "GT", NULL, true, idset);

	indids = tvcf.inds;

	// set bytesPerMarker attribute
	//error("indids.size() = %d",tvcf.nInds);

	if (siteOnly) {
		bytesPerMarker = 0;
	}
	else {
		bytesPerMarker = (tvcf.nInds + 3) / 4;
	}

	for (M = 0; tvcf.readMarkers(unit);) {
		M += tvcf.nMarkers;
		fprintf(stderr, "Processing %d markers across %d individuals..\n", M, tvcf.nInds);
		for (i = 0, m = 0; i < tvcf.nMarkers; ++i) { // for each marker
			af = tvcf.alleleFreq(i);
			maf = af > 0.5 ? 1 - af : af;
			if ((maf >= minAF) &&
				(tvcf.callRate(i) >= minCallRate)) {
				addMarker(tvcf.chroms[i].c_str(), tvcf.pos1s[i], tvcf.refs[i][0], tvcf.alts[i][0], af, tvcf.nInds);
				for (int j = 0; j < tvcf.nInds; ++j) {
					setGenotype(tvcf.genos[i*tvcf.nInds + j], j, m);
				}
				++m;
			}
		}
	}
}
int GenotypeMatrix::addMarker(const char* chrom, int position, char refBase, char altBase, double alleleFreq, int numIndividuals)
{
	chroms.push_back(chrom);
	positions.push_back(position);
	refBases.push_back(refBase);
	altBases.push_back(altBase);
	alleleFrequencies.push_back(alleleFreq);

	//notice("%s:%d\t%d",chrom,position,(int)genotypes.size());
	for (int i = 0; i < bytesPerMarker; ++i) {
		genotypes.push_back(std::vector<float>(numIndividuals, 0));
	}
	//notice("%s:%d\t%d",chrom,position,(int)genotypes.size());

	return (int)chroms.size();
}
void GenotypeMatrix::setGenotype(float geno, int indIndex, int markerIndex = -1)
{
	float ngeno = isnan(geno) ? 0 : geno;
	genotypes[markerIndex][indIndex] = ngeno;
}
float GenotypeMatrix::getGenotype(int indIndex, int markerIndex)
{
	return genotypes[markerIndex][indIndex];
}
void GenotypeMatrix::printVCF(std::string path)
{
	std::ofstream fout(path);
	if (!fout.is_open()) { warning("Open file %s failed!\n", path.c_str()); }
	for (size_t i = 0; i != tvcf.headers.size(); ++i)
	{
		fout << tvcf.headers[i] << std::endl;
	}

	for (size_t i = 0; i != tvcf.chroms.size(); ++i)
	{
		for (size_t j = 0; j != tvcf.markers.size(); ++j)
		{
			fout << tvcf.chroms[i] << "\t" << tvcf.pos1s[j] << "\t" << tvcf.markers[j] << "\t" << tvcf.refs[j] << "\t" << tvcf.alts[j] << "\tPASS\tAF=" << tvcf.AFs[j] << std::endl;
		}
	}
	fout.close();
}
GenotypeMatrix::~GenotypeMatrix(){}
#endif
/*Functions for fullLLKFunc*/
PopulationIdentifier::fullLLKFunc::fullLLKFunc(){}
PopulationIdentifier::fullLLKFunc::~fullLLKFunc(){}
/*Functions for PopulationIdentifier*/
PopulationIdentifier::PopulationIdentifier()
{
}
#ifdef ARMADILLO
PopulationIdentifier::PopulationIdentifier(const std::string& VCF, PopArgs *p, const std::string & GLpath) :pArgs(p), GenoMatrix(VCF.c_str(), pArgs->bSiteOnly, pArgs->bFindBest/*subsetInds*/, pArgs->minAF, pArgs->minCallRate)
{//using selected sites from union
	/*PopArgs* pArgs =new PopArgs;*/
	NumMarker = GenoMatrix.tvcf.nMarkers;
	NumIndividual = GenoMatrix.tvcf.nInds;
	//fullLLKFunc fn;
	std::vector<PCtype> tmpPC(2, 0);
	UD = std::vector<std::vector<PCtype> >(NumMarker, tmpPC);
	PC = std::vector<std::vector<PCtype> >(NumIndividual, tmpPC);
	ReadMatrixGL(GLpath);
	FormatMarkerIntersection();
	ImputeMissing();
	RunSVD();
	fn = PopulationIdentifier::fullLLKFunc(this);
}
int PopulationIdentifier::ImputeMissing()
{

	auto mean = [&](std::vector<float>& vec){float sum = 0; int num(0); std::for_each(vec.begin(), vec.end(), [&](float a){if (!isnan(a)) { sum += a; num++; }}); return sum / num; };
	for (size_t i = 0; i != NumMarker; ++i)
	{
		means[i] = mean(GenoMatrix.genotypes[i]);
		for (size_t j = 0; j != NumIndividual; ++j)
		{
			if (isnan(GenoMatrix.genotypes[i][j]))
				GenoMatrix.genotypes[i][j] = means[i];
		}
	}
	return 0;
}
int PopulationIdentifier::RunSVD()
{
	fmat GenoM;// (GenoMatrix.numMarker, GenoMatrix.numIndividual);
	for (size_t i = 0; i != NumMarker; ++i)
	{
		for (size_t j = 0; j != NumIndividual; ++j)
		{
			GenoM << GenoMatrix.genotypes[i][j];
		}
		GenoM << endr;
	}
	fmat tU;
	fvec ts;
	fmat tV;
	svd(tU, ts, tV, GenoM);
	fmat tUD = tU*ts;
	for (size_t i = 0; i != NumMarker; ++i)
	{
		UD[i][0] = tUD(i, 0);
		UD[i][1] = tUD(i, 1);
	}
	for (size_t j = 0; j != NumIndividual; ++j)
	{
		PC[j][0] = tV(j, 0);
		PC[j][1] = tV(j, 1);
	}
	return 0;
}
#endif

#include <fstream>
#include <sstream>

std::vector<std::string> bedContainer;
int PopulationIdentifier::writeVcfFile(const std::string& path)
{
	std::ofstream fout(path);
	if (!fout.is_open()) { abort(); }
	fout
		<< "##fileformat=VCFv4.0\n"
		<< "##FILTER=<ID=NOT_POLY_IN_1000G,Description=\"Alternate allele count = 0\">\n"
		<< "##FILTER=<ID=badAssayMapping,Description=\"The mapping information for the SNP assay is internally inconsistent in the chip metadata\">\n"
		<< "##FILTER=<ID=dup,Description=\"Duplicate assay at same position with worse Gentrain Score\">\n"
		<< "##FILTER=<ID=id10,Description=\"Within 10 bp of an known indel\">\n"
		<< "##FILTER=<ID=id20,Description=\"Within 20 bp of an known indel\">\n"
		<< "##FILTER=<ID=id5,Description=\"Within 5 bp of an known indel\">\n"
		<< "##FILTER=<ID=id50,Description=\"Within 50 bp of an known indel\">\n"
		<< "##FILTER=<ID=refN,Description=\"Reference base is N. Assay is designed for 2 alt alleles\">\n"
		<< "##FORMAT=<ID=GC,Number=.,Type=Float,Description=\"Gencall Score\">\n"
		<< "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
		<< "##INFO=<ID=CR,Number=.,Type=Float,Description=\"SNP Callrate\">\n"
		<< "##INFO=<ID=GentrainScore,Number=.,Type=Float,Description=\"Gentrain Score\">\n"
		<< "##INFO=<ID=HW,Number=.,Type=Float,Description=\"Hardy-Weinberg Equilibrium\">\n"
		<< "##reference=human_g1k_v37.fasta\n"
		<< "##source=infiniumFinalReportConverterV1.0\n"
		<< "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
	
	for (size_t i = 0; i < bedContainer.size(); i++)
	{
		std::stringstream ss(bedContainer[i]);
		std::string chr;
		int pos;
		char ref, alt;
		ss >> chr;
		fout << chr << "\t";
		ss >> pos;
		fout << pos+1 << "\t"<<"SNPID"<<pos+1<<"\t";
		ss >> pos;
		ss >> ref;
		fout << ref<<"\t";
		ss >> alt;
		fout << alt << "\t"<<".\tPASS\tAF=";
		fout << this->AFs[i] << "\n";


	}
	return 0;
}
PopulationIdentifier::PopulationIdentifier(const std::string& UDpath, const std::string &PCpath, const std::string &Mean, const std::string & GLpath, const std::string &Bed)
{
	ReadChooseBed(Bed);
	NumMarker = ReadMatrixUD(UDpath);
	NumIndividual = ReadMatrixPC(PCpath);
	
	//ReadMatrixGL(GLpath);
	ReadPileup(GLpath);
	FormatMarkerIntersection();
	ReadMean(Mean);
	fn = PopulationIdentifier::fullLLKFunc(this);
}
int PopulationIdentifier::ReadMatrixUD(const std::string &path)
{
	std::ifstream fin(path);
	std::string line;
	uint32_t index(0);
	std::vector<PCtype> tmpUD(2, 0);
	if (!fin.is_open()) { warning("Open file %s failed!\n", path.c_str()); }
	while (std::getline(fin, line))
	{
		std::stringstream ss(line);
		//std::string chr;
		//int pos;
		//ss >> chr >> pos;
		ss >> tmpUD[0] >> tmpUD[1];
		UD.push_back(tmpUD);
	}
	fin.close();
	return UD.size();
}
int PopulationIdentifier::ReadChooseBed(const std::string &path)
{
	std::ifstream fin(path);
	std::string line,chr;
	uint32_t index(0),pos(0);
	char ref(0), alt(0);

	if (!fin.is_open()) { warning("Open file %s failed!\n", path.c_str()); }
	while (std::getline(fin, line))
	{
		index++;
		std::stringstream ss(line);
		//std::string chr;
		//int pos;
		ss >> chr >> pos>>pos;
		ss >> ref >> alt;
	
		bedContainer.push_back(line);
		ChooseBed[chr][pos] = std::make_pair(ref,alt);
	}
	fin.close();
	return UD.size();
}
int PopulationIdentifier::ReadMatrixPC(const std::string &path)
{
	std::ifstream fin(path);
	std::string line;
	uint32_t index(0);
	std::vector<PCtype> tmpPC(2, 0);
	if (!fin.is_open()) { warning("Open file %s failed!\n", path.c_str()); }
	while (std::getline(fin, line))
	{
		std::stringstream ss(line);
		std::string dummy;
		//std::string chr;
		//int pos;
		//ss >> chr >> pos;
		//MarkerIndex[chr][pos] = index;
		ss >> dummy >> tmpPC[0] >> tmpPC[1];
		PC.push_back(tmpPC);
	}
	fin.close();
	return PC.size();
}
int PopulationIdentifier::ReadMean(const std::string &path)
{
	std::ifstream fin(path);
	std::string line;
	uint32_t index(0),pos(0);
	double mu(0);
	std::string snpName,chr;
	if (!fin.is_open()) { warning("Open file %s failed!\n", path.c_str()); }
	while (std::getline(fin, line))
	{
		std::stringstream ss(line);
		ss >> snpName;
		chr=snpName.substr(0,snpName.find(':',0));
		pos = atoi(snpName.substr(snpName.find(':', 0)+1,snpName.find('_',0)).c_str());
		ss >> mu;
		//std::cerr << chr << "\t" << pos << std::endl;
		PosVec.push_back(make_pair(chr, pos));
		means[index]=mu;
		index++;
	}
	fin.close();
	return means.size();
}
int PopulationIdentifier::CheckMarkerSetConsistency()
{
	return UD.size() == GL.size();
}
int PopulationIdentifier::FormatMarkerIntersection()
{
	//if (CheckMarkerSetConsistency())
	{
		means = std::vector<PCtype>(NumMarker, 0.5);
		AFs = std::vector<double>(NumMarker, 0);
	}
	//else//implement intersection behavior if needed
	//{
	//	std::cerr << "[Waring] - Marker Sets are not consistent UD:"<<UD.size()<<"\tGL:"<<GL.size() << std::endl;
	//	exit(EXIT_FAILURE);
	//}
	return 0;
}
int PopulationIdentifier::ReadMatrixGL(const std::string& path)//Reading GL information for individual need to be classified
{
	std::ifstream fin(path);
	std::string line;
	uint32_t index(0);
	std::vector<double> tmpGL(3, 0);
	if (!fin.is_open()) { warning("Open file %s failed!\n", path.c_str()); }
	while (std::getline(fin, line))
	{
		std::stringstream ss(line);
		std::string chr;
		int pos;
		ss >> chr >> pos;
		if (ChooseBed.size()!=0)
		{
			if (ChooseBed.find(chr) == ChooseBed.end()) continue;
			else if (ChooseBed[chr].find(pos) == ChooseBed[chr].end()) continue;
		}
		MarkerIndex[chr][pos] = index;
		ss >> tmpGL[0] >> tmpGL[1] >> tmpGL[2];
		//cout << chr <<"\t"<< pos << "\t"<<tmpGL[0] << "\t"<<tmpGL[1] <<"\t"<< tmpGL[2] << endl;
		GL.push_back(tmpGL);
		index++;
	}
	fin.close();
	return 0;
}

#define REV_PHRED(x)	pow(10.0,x/(-10.0))
static std::vector<double>  calLikelihood(const std::string & seq, const std::string & qual,char ref,char alt)//maj:ref, min:alt
{	

	double lik(0), GL0(0), GL1(0), GL2(0);
	{
		for (uint32_t i = 0; i != seq.size(); ++i)
		{
			double seq_error = REV_PHRED(qual[i]);
			//fprintf(stderr,"Debug (qual:%d\tseq_error:%f)\n",qual[i],seq_error);
			if (seq[i] == '.'||seq[i]==',')
			{
				GL0 += log10(1 - seq_error);
				GL1 += log10(0.5 - seq_error / 3);
				GL2 += log10(seq_error / 3);
			}
			else if (seq[i] == alt||seq[i]==toupper(alt))
			{
				GL0 += log10(seq_error / 3);
				GL1 += log10(0.5 - seq_error / 3);
				GL2 += log10(1 - seq_error);
			}
			else
			{
				GL0 += log10(2 * seq_error / 3);
				GL1 += log10(2 * seq_error / 3);
				GL2 += log10(2 * seq_error / 3);
			}
		}
		//fprintf(stderr, "\n");
	}
	std::vector<double> tmp(3, 0);
	tmp[0] = GL0*(-10); tmp[1] = GL1*(-10); tmp[2] = GL2*(-10);
	double minimal(tmp[0]);
	if (tmp[1] < minimal)
	{
		minimal = tmp[1];
	}
	if (tmp[2] < minimal)
	{
		minimal = tmp[2];
	}
	tmp[0] -= minimal; tmp[1] -= minimal; tmp[2] -= minimal;
	return tmp;
}
int PopulationIdentifier::ReadPileup(const std::string& path)
{
	//std::ifstream fin(path);
	InputFile fin(path.c_str(), "r", InputFile::BGZF);
	std::string line;
	uint32_t index(0);

	std::vector<double> tmpGL(3, 0);
	char ref, alt;
	if (!fin.isOpen()) { warning("Open file %s failed!\n", path.c_str()); }

	int sephore = -1;
	while ((line="",fin.readLine(line)!=-1))
	{
		if (line[0] == '#') continue;
		std::stringstream ss(line);
		std::string chr,seq,qual;
		int pos;
		ss >> chr >> pos;
		if (ChooseBed[chr].find(pos) == ChooseBed[chr].end()||MarkerIndex[chr].find(pos)!=MarkerIndex[chr].end()){continue; }
		else
		{
			ref = ChooseBed[chr][pos].first;
			alt = ChooseBed[chr][pos].second;
		}
		MarkerIndex[chr][pos] = index;
		ss >> seq;
		ss >> seq;
		ss >> seq;
		ss >> qual;
		ss >> qual;
		ss >> qual;
		//std::cout << "calculating likelihood for: " << chr << "\t" << pos << "\t" << ref << "\t" << alt << std::endl;
		tmpGL = calLikelihood(seq,qual,ref,alt);

		ss >> tmpGL[0] >> tmpGL[1] >> tmpGL[2];
		//std::cout << chr <<"\t"<< pos << "\t"<<tmpGL[0] << "\t"<<tmpGL[1] <<"\t"<< tmpGL[2] << std::endl;
		GL.push_back(tmpGL);
		index++;
	}
	
	{//for marker existed in UD but not in pileup
		tmpGL[0] = tmpGL[1] = tmpGL[2] = 0.0;
		GL.push_back(tmpGL);
	}
	fin.ifclose();
	return 0;
}

int PopulationIdentifier::OptimizeLLK()
{
	//std::cerr << "Now the label is:1" << std::endl;
	AmoebaMinimizer myMinimizer;
	//std::cerr << "Now the label is:2" << std::endl;
	Vector startingPoint("TestPoint",2);
	//std::cerr << "Now the label is:3" << std::endl;
	startingPoint[0] = PC[0][0];  
	startingPoint[1] = PC[0][1];   
	//startingPoint.label = "startPoint";
	//std::cerr << "Now the label is:" << startingPoint.label << std::endl;
	myMinimizer.func = &fn;
	//std::cerr << "Before minimizing:1" << startingPoint.label << std::endl;
	myMinimizer.Reset(2);
	//std::cerr << "Before minimizing:2" << startingPoint.label << std::endl;
	myMinimizer.point = startingPoint;
	//std::cerr << "Before minimizing:3" << startingPoint.label << std::endl;
	myMinimizer.Minimize(1e-6);
	double optimalPC1 = myMinimizer.point[0];// fullLLKFunc::invLogit(myMinimizer.point[0]);
	double optimalPC2 = myMinimizer.point[1]; //fullLLKFunc::invLogit(myMinimizer.point[1]);
	std::cout << "PCs in OptimizaLLK():" << std::endl;
	std::cout << "PC1:" << optimalPC1 << "\tPC2:" << optimalPC2 << std::endl;

	return 0;
}
int PopulationIdentifier::RunMapping()
{
	return 0;
}
#ifdef ARMADILLO
int PopulationIdentifier::PrintVcf(const std::string & path)
{
	GenoMatrix.printVCF(path);
	std::cout << "Population PCs for this individuals:" << std::endl;
	std::cout << "PC1:" << fn.PC1 << "\tPC2:" << fn.PC2 << std::endl;
	return 0;
}
#endif
PopulationIdentifier::~PopulationIdentifier()
{
}


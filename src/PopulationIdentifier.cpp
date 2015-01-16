#include "PopulationIdentifier.h"
#include <fstream>
#include <algorithm>
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
/*Functions for fullLLKFunc*/
PopulationIdentifier::fullLLKFunc::fullLLKFunc(){}
PopulationIdentifier::fullLLKFunc::~fullLLKFunc(){}
/*Functions for PopulationIdentifier*/
PopulationIdentifier::PopulationIdentifier()
{
}
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

PopulationIdentifier::PopulationIdentifier(const std::string& UDpath, const std::string &PCpath, const std::string & GLpath)
{
	NumMarker = ReadMatrixUD(UDpath);
	NumIndividual = ReadMatrixPC(PCpath);
	ReadMatrixGL(GLpath);
	FormatMarkerIntersection();
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
int PopulationIdentifier::CheckMarkerSetConsistency()
{
	return UD.size() == GL.size();
}
int PopulationIdentifier::FormatMarkerIntersection()
{
	if (CheckMarkerSetConsistency())
	{
		means = std::vector<PCtype>(NumMarker, 0.5);
		AFs = std::vector<double>(NumMarker, 0);
	}
	else//implement intersection behavior if needed
	{
		std::cerr << "[Waring] - Marker Sets are not consistent" << std::endl;
		exit(EXIT_FAILURE);
	}
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
		MarkerIndex[chr][pos] = index;
		ss >> tmpGL[0] >> tmpGL[1] >> tmpGL[2];
		GL.push_back(tmpGL);
	}
	fin.close();
	return 0;
}

int PopulationIdentifier::OptimizeLLK()
{
	//fullLLKFunc myFunc(UD,PC, GL);
	AmoebaMinimizer myMinimizer;
	Vector startingPoint(2);
	startingPoint[0] = PC[0][0];  // start with fMix = 0.01
	startingPoint[1] = PC[0][1];       // pRefHet = 0.5

	myMinimizer.func = &fn;
	myMinimizer.Reset(2);
	myMinimizer.point = startingPoint;
	myMinimizer.Minimize(1e-3);
	double optimalPC1 = fullLLKFunc::invLogit(myMinimizer.point[0]);
	double optimalPC2 = fullLLKFunc::invLogit(myMinimizer.point[1]);
	std::cout << "PCs in OptimizaLLK():" << std::endl;
	std::cout << "PC1:" << optimalPC1 << "\tPC2:" << optimalPC2 << std::endl;
	return 0;
}
int PopulationIdentifier::RunMapping()
{
	return 0;
}
int PopulationIdentifier::PrintVcf(const std::string & path)
{
	GenoMatrix.printVCF(path);
	std::cout << "Population PCs for this individuals:" << std::endl;
	std::cout << "PC1:" << fn.PC1 << "\tPC2:" << fn.PC2 << std::endl;
	return 0;
}

PopulationIdentifier::~PopulationIdentifier()
{
}

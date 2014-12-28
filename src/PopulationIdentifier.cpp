#include "PopulationIdentifier.h"
#include <fstream>
#include <algorithm>
#include <armadillo>
using namespace arma;
/*Functions for GenotypeMatrix*/
GenotypeMatrix::GenotypeMatrix(){}
GenotypeMatrix::GenotypeMatrix(const char* vcfFile, bool siteOnly, bool findBest, std::vector<std::string>& subsetInds, double minAF, double minCallRate)
{
	// open a VCF file

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
				addMarker(tvcf.chroms[i].c_str(), tvcf.pos1s[i], tvcf.refs[i][0], tvcf.alts[i][0], af,tvcf.nInds);
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
	genotypes[markerIndex][indIndex]= ngeno;
}
float GenotypeMatrix::getGenotype(int indIndex, int markerIndex)
{
	return genotypes[markerIndex][indIndex];
}
double GenotypeMatrix::computeAlleleFrequency(int markerIndex)
{

}
void GenotypeMatrix::printVCF(std::string path)
{
	std::ofstream fout(path);
	if (!fout.is_open()) { notice("Open file %s failed!\n",path.c_str()); }
	for (size_t i; i != tvcf.headers.size(); ++i)
	{
		fout << tvcf.headers[i] << std::endl;
	}

	for (size_t i; i != tvcf.chroms.size(); ++i)
	{
		for (size_t j; j != tvcf.markers.size(); ++j)
		{
			fout << tvcf.chroms[i] << "\t" << tvcf.pos1s[j] << "\t" << tvcf.markers[j] <<"\t"<<tvcf.refs[j]<<"\t"<<tvcf.alts[j]<<"\tPASS\tAF="<<tvcf.AFs[j]<<std::endl;
		}
	}
}
GenotypeMatrix::~GenotypeMatrix(){}
/*Functions for fullLLKFunc*/
fullLLKFunc::fullLLKFunc(){}
fullLLKFunc::~fullLLKFunc(){}
/*Functions for PopulationIdentifier*/
PopulationIdentifier::PopulationIdentifier()
{
}
PopulationIdentifier::PopulationIdentifier(std::string& VCF) :fn(this)
{

}
int PopulationIdentifier::ReadingGL(const std::string& path)
{

}
int PopulationIdentifier::ImputeMissing()
{
	means = std::vector<float>(NumMarker, 0);//mean genotype for each marker
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
		PC[j][0] = tV(j,0);
		PC[j][1] = tV(j,1);
	}

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
	myMinimizer.Minimize(1e-6);
	double optimalPC1 = fullLLKFunc::invLogit(myMinimizer.point[0]);
	double optimalPC2 = fullLLKFunc::invLogit(myMinimizer.point[1]);

}
int PopulationIdentifier::RunMapping()
{

}
int PopulationIdentifier::PrintVcf()
{

}

PopulationIdentifier::~PopulationIdentifier()
{
}

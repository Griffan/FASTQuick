#ifndef POPULATIONIDENTIFIER_H_
#define POPULATIONIDENTIFIER_H_
#include <string>
#include "../libmpu/fVcf.h"
#include "../libmpu/MathVector.h"
#include "../libmpu/MathGenMin.h"
class GenotypeMatrix{
public:
	fVcf tvcf;
	//int nInds;
	std::vector<std::string> chroms;
	std::vector<int> positions;
	std::vector<double> alleleFrequencies;
	std::vector<char> refBases;
	std::vector<char> altBases;

	std::vector<std::string> indids;

	std::vector<std::vector<float> > genotypes;
	//uint32_t numIndividual;
	//uint32_t numMarker;
	int bytesPerMarker;
	GenotypeMatrix();
	GenotypeMatrix(const char* vcfFile, bool siteOnly, bool findBest, std::vector<std::string>& subsetInds, double minAF, double minCallRate);
	int addMarker(const char* chrom, int position, char refBase, char altBase, double alleleFreq,int numIndividuals);
	void setGenotype(float genotype, int indIndex, int markerIndex = -1);
	float getGenotype(int indIndex, int markerIndex);
	double computeAlleleFrequency(int markerIndex);
	void printVCF(std::string path);
	~GenotypeMatrix();
};
#define PCtype float
class fullLLKFunc : public VectorFunc {
public:

	double llk0;
	double llk1;
	PopulationIdentifier * ptr;
	double PC1, PC2;
	fullLLKFunc();
	~fullLLKFunc();
	inline static double fullLLKFunc::invLogit(double & x){double e = exp(x); return e / (1. + e); };
	inline double computeMixLLKs(double PC1, double PC2)
	{
		double min_af(0.5 / ptr->NumIndividual), max_af((ptr->NumIndividual - 0.5) / ptr->NumIndividual);
		double sumLLK(0),GF0(0),GF1(0),GF2(0);
		for (size_t i = 0; i != ptr->NumMarker; ++i)
		{
			ptr->AFs[i] = ((ptr->UD[i][0] * PC1 + ptr->UD[i][1] * PC2) + ptr->means[i]) / 2;
			if (ptr->AFs[i] < min_af) ptr->AFs[i] = min_af;
			if (ptr->AFs[i] > max_af) ptr->AFs[i] = max_af;
			GF0 = (1 - ptr->AFs[i])*(1 - ptr->AFs[i]);
			GF1 = 2*(ptr->AFs[i])*(1-ptr->AFs[i]);
			GF2 = (ptr->AFs[i])*(ptr->AFs[i]);
			sumLLK += log(ptr->GL[i][0] * GF0 + ptr->GL[i][1] * GF1 + ptr->GL[i][2] * GF2);
		}
	}
	fullLLKFunc(PopulationIdentifier* inPtr){
		ptr = inPtr;
		llk1 = llk0 = (0 - computeMixLLKs(static_cast<double>(ptr->PC[0][0]), static_cast<double>(ptr->PC[0][1])));
	}

	virtual double Evaluate(Vector& v) {
		if (v.Length() != 2)
			error("fullMixLLKFunc(): Input vector must be length of 3");

		double tmpPC1 = invLogit(v[0]);
		double tmpPC2 = invLogit(v[1]);

		double smLLK = 0 - computeMixLLKs(tmpPC1,tmpPC2);

		if (smLLK < llk1) {
			llk1 = smLLK;
			PC1 = tmpPC1;
			PC2 = tmpPC2;
		}

		return smLLK;
	}
};
class PopulationIdentifier
{
public:
	GenotypeMatrix GenoMatrix;
	int NumMarker;
	int NumIndividual;
	fullLLKFunc fn;

	std::vector<std::vector<PCtype> > UD;
	std::vector<std::vector<PCtype> > PC;
	std::vector<std::vector<PCtype> > GL;
	std::vector<float> means;
	std::vector<double> AFs;

	PopulationIdentifier();
	PopulationIdentifier(std::string& VCF);
	int ReadingGL(const std::string& path);
	int ImputeMissing();
	int RunSVD();
	int OptimizeLLK();
	int RunMapping();
	int PrintVcf();
	~PopulationIdentifier();
};
#endif POPULATIONIDENTIFIER_H_

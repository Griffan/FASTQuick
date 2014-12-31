#ifndef POPULATIONIDENTIFIER_H_
#define POPULATIONIDENTIFIER_H_
#include <string>
#include <unordered_map>
#include "fVcf.h"
#include "MathVector.h"
#include "MathGenMin.h"
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
	GenotypeMatrix(const char* vcfFile, bool siteOnly, bool findBest, double minAF, double minCallRate);
	int addMarker(const char* chrom, int position, char refBase, char altBase, double alleleFreq,int numIndividuals);
	void setGenotype(float genotype, int indIndex, int markerIndex);
	float getGenotype(int indIndex, int markerIndex);
	//double computeAlleleFrequency(int markerIndex);
	void printVCF(std::string path);
	~GenotypeMatrix();
};



class PopArgs {
public:
	std::string sSubsetInds;
	std::string sVcfFile;
	std::string sMpuFile;
	std::string sOutFile;
	std::string sSMID;
	std::string sScanPair;

	bool bSelfOnly;
	bool bSiteOnly;
	bool bFindBest;
	bool bIgnoreSex;

	//bool bIgnoreOverlapPair;
	bool bPrecise;
	bool bVerbose;
	bool bSilent;

	double genoError;
	double minAF;
	double minCallRate;
	double contamThres;

	double pRefRef;
	double pRefHet;
	double pRefAlt;

	bool bFreeNone;
	bool bFreeMixOnly;
	bool bFreeRefBiasOnly;
	bool bFreeFull;

	bool bChipNone;
	bool bChipMixOnly;
	bool bChipRefBiasOnly;
	bool bChipFull;

	int minMapQ;
	int maxDepth;
	int minQ;
	int maxQ;

	double grid;

	PopArgs() {
		bSelfOnly = false;
		bSiteOnly = false;
		bFindBest = false;
		bIgnoreSex = false;

		bPrecise = false;

		genoError = 0.001;
		minAF = 0.01;
		minCallRate = 0.50;
		contamThres = 0.02;

		pRefRef = 1;
		pRefHet = 0.5;
		pRefAlt = 0;

		bFreeNone = false;
		bFreeMixOnly = true;
		bFreeRefBiasOnly = false;
		bFreeFull = false;

		bChipNone = false;
		bChipMixOnly = true;
		bChipRefBiasOnly = false;
		bChipFull = false;

		bVerbose = false;

		minMapQ = 10;
		maxDepth = 20;
		minQ = 13;
		maxQ = 40;

		grid = 0.05;
	}
};

class PopulationIdentifier
{
public:

#define PCtype float
	class fullLLKFunc : public VectorFunc {
	public:

		double llk0;
		double llk1;
		PopulationIdentifier* ptr;
		double PC1, PC2;
		fullLLKFunc();
		~fullLLKFunc();
		inline static double invLogit(double & x){ double e = exp(x); return e / (1. + e); };
		inline double computeMixLLKs(double tPC1, double tPC2)
		{
			double min_af(0.5 / ptr->NumIndividual), max_af((ptr->NumIndividual - 0.5) / ptr->NumIndividual);
			double sumLLK(0), GF0(0), GF1(0), GF2(0);
			for (size_t i = 0; i != ptr->NumMarker; ++i)
			{
				ptr->AFs[i] = ((ptr->UD[i][0] * tPC1 + ptr->UD[i][1] * tPC2) + ptr->means[i]) / 2;
				if (ptr->AFs[i] < min_af) ptr->AFs[i] = min_af;
				if (ptr->AFs[i] > max_af) ptr->AFs[i] = max_af;
				GF0 = (1 - ptr->AFs[i])*(1 - ptr->AFs[i]);
				GF1 = 2 * (ptr->AFs[i])*(1 - ptr->AFs[i]);
				GF2 = (ptr->AFs[i])*(ptr->AFs[i]);
				sumLLK += log(ptr->GL[i][0] * GF0 + ptr->GL[i][1] * GF1 + ptr->GL[i][2] * GF2);
			}
			return sumLLK;
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

			double smLLK = 0 - computeMixLLKs(tmpPC1, tmpPC2);

			if (smLLK < llk1) {
				llk1 = smLLK;
				PC1 = tmpPC1;
				PC2 = tmpPC2;
			}

			return smLLK;
		}
	};


	PopArgs* pArgs;
	GenotypeMatrix GenoMatrix;
	int NumMarker;
	int NumIndividual;
	fullLLKFunc fn;

	std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t> > MarkerIndex;

	std::vector<std::vector<PCtype> > UD;
	std::vector<std::vector<PCtype> > PC;
	std::vector<std::vector<PCtype> > GL;
	std::vector<PCtype> means;
	std::vector<double> AFs;

	PopulationIdentifier();
	PopulationIdentifier(std::string& VCF,PopArgs* p);
	int ReadingGL(const std::string& path);
	int ImputeMissing();
	int RunSVD();
	int OptimizeLLK();
	int RunMapping();
	int PrintVcf(const std::string & path);
	~PopulationIdentifier();
};


#endif

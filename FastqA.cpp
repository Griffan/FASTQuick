//============================================================================
// Name        : FastqA.cpp
// Author      : FanZhang
// Version     :
// Copyright   : Copyright Reserved.
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>
#include "./libbwa/utils.h"
#include "BwtIndexer.h"
#include "BwtMapper.h"
#include "params.h"
//using namespace std;

#define USE_BWT 1
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.1"
#endif

int main(int argc, char ** argv)
{

	double t_real;
	t_real = realtime();
	int c, opte = -1;
	gap_opt_t *opt;
	pe_opt_t *popt;
	popt = bwa_init_pe_opt();
	opt = gap_init_opt();

	/*
	 * Parameters
	 *
	 */

	std::string RefPath, VcfPath, MaskPath("Empty"), Fastq_1("Empty"), Fastq_2(
			"Empty"), BamOut("Empty"), ReadGroup("default"),DepthDist, SitePileup,FaList("Empty");

	bool loggap(0), compread(0), nonstop(0), IL13(0);
	paramList pl;

	BEGIN_LONG_PARAMS(longParameters) LONG_PARAM_GROUP("Input/Output Files","Input/Output files for the program[Complete Path Recommended]")
	LONG_STRING_PARAM("vcf",&VcfPath,"Input VCF file name")
	LONG_STRING_PARAM("ref",&RefPath,"Reference FASTA file.")
	LONG_STRING_PARAM("mask",&MaskPath,"Repeat Mask FASTA file.")
	LONG_STRING_PARAM("fastq_1",&Fastq_1,"Pair end 1 fastq file.")
	LONG_STRING_PARAM("fastq_2",&Fastq_2,"Pair end 2 fastq file.[Leave empty if using single end]")
	LONG_STRING_PARAM("bam_out",&BamOut,"Output file prefix")
	LONG_STRING_PARAM("fq_list",&FaList,"Path of input fastq files, tab-delimited, one pair-end files per line(one file per line for single end) ")


	//LONG_STRING_PARAM("out",&outf,"Output file prefix")
	LONG_PARAM_GROUP("Parameters for Reference Sequence ", "Parameters being used to extract reference sequences.")
	// LONG_INT_PARAM("K",&Indexer.RollParam.kmer_size,"kmer size for Rolling Hash filtering")
	// LONG_INT_PARAM("T",&Indexer.RollParam.thresh,"threshold for filtering specific read")
	// LONG_INT_PARAM("S",&Indexer.RollParam.read_step_size,"step size for extracting kmer")
	LONG_INT_PARAM("var_long",&opt->num_variant_long,"number of variants with long flanking region")
	LONG_INT_PARAM("var_short",&opt->num_variant_short,"number of variants with short flanking region")
	LONG_INT_PARAM("flank_len",&opt->flank_len,"flanking region length around each marker")

	LONG_PARAM_GROUP("Parameters for Alignment ", "Parameters the are universal for both single end and pair end alignment.")
	LONG_DOUBLE_PARAM("n",&opt->fnr,"Max #diff (int) or missing prob under 0.02 error rate [float]")
	LONG_INT_PARAM("o",&opt->max_gapo,"maximum number or fraction of gap opens")
	LONG_INT_PARAM("e",&opte,"maximum number of gap extensions, -1 for disabling long gaps [-1]")
	LONG_INT_PARAM("i",&opt->indel_end_skip,"do not put an indel within INT bp towards the ends")
	LONG_INT_PARAM("d",&opt->max_del_occ,"maximum occurrences for extending a long deletion")
	LONG_INT_PARAM("l",&opt->seed_len,"seed length")
	LONG_INT_PARAM("k",&opt->max_seed_diff,"maximal seed difference")
	LONG_INT_PARAM("m",&opt->max_entries,"maximal stack entries")
	LONG_INT_PARAM("t",&opt->n_threads,"number of threads")
	// LONG_INT_PARAM("L",&opt->seed_len,"seed length")
	LONG_INT_PARAM("R",&opt->max_top2,"stop searching when there are >INT equally best hits")
	LONG_INT_PARAM("q",&opt->trim_qual,"quality threshold for read trimming down to 35dbp")
	LONG_STRING_PARAM("RG",&ReadGroup,"set ReadGroup name")
	//Roll Hash

	// LONG_PARAM("c",&compread,"seed length")
	LONG_PARAM("N",&nonstop,"non-iterative mode: search for all n-difference hits (slooow)")
	LONG_PARAM("I",&IL13,"the input is in the Illumina 1.3+ FASTQ-like format")
	LONG_PARAM("L",&loggap,"log-scaled gap penalty for long deletions")

	LONG_PARAM_GROUP("Parameters for PairEnd ", "Parameters specified for Pair end mapping .")
	LONG_INT_PARAM("max_isize",&popt->max_isize," maximum insert size")
	LONG_INT_PARAM("max_occ",&popt->max_occ,"maximum occurrences for one end ")
	LONG_PARAM("is_sw",&popt->is_sw,"disable Smith-Waterman for the unmapped mate")
	LONG_INT_PARAM("n_multi",&popt->n_multi,"maximum hits to output for paired reads")
	LONG_INT_PARAM("N_multi",&popt->N_multi,"maximum hits to output for discordant pairs")
	LONG_DOUBLE_PARAM("ap_prior",&popt->ap_prior,"prior of chimeric rate (lower bound) ")
	LONG_PARAM("force_isize",&popt->force_isize," disable insert size estimate")

	LONG_PARAM_GROUP("Parameters for Statistics ", "Parameters specified for statistics and summary .")
	LONG_PARAM("cal_dup",&opt->cal_dup,"enable the calculation of duplicated reads in depth calculation ")
	//LONG_STRING_PARAM("depth_dist",&DepthDist,"Output file for depth distribution ")
	//LONG_STRING_PARAM("site_pileup",&SitePileup,"Output file for Pileup information on specified sites ")

	END_LONG_PARAMS();

	pl.Add(new longParams("Available Options", longParameters));
	pl.Read(argc, argv);
	pl.Status();
	if(BamOut!="Empty")
	{
	opt->bam_name=strdup(BamOut.c_str());
	}
	opt->RG=strdup(ReadGroup.c_str());
	if (opte > 0)
	{
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}
	if (nonstop)
	{
		opt->mode |= BWA_MODE_NONSTOP;
		opt->max_top2 = 0x7fffffff;
	}
	if (IL13)
	{
		std::cerr << "Using Illumina 1.3 version quality system..." << std::endl;
		opt->mode |= BWA_MODE_IL13;
	}
	if (loggap)
	{
		opt->mode |= BWA_MODE_LOGGAP;
	}
	//Extract ref seqs from VCF and Ref

	struct stat sb;

	if (USE_BWT) // using bwt methods
	{
		//build ref index
		BwtIndexer Indexer;
		std::string BwtPath = RefPath + ".bwt";
		//string BwtPathR = RefPath + ".rbwt";
		if (stat(BwtPath.c_str(), &sb) != 0) //|| stat(BwtPathR.c_str(), &sb)!=0)
		{
			std::cerr << "Index file doesn't exist, building...\n";
			RefBuilder ArtiRef(VcfPath, RefPath, MaskPath, opt);
					//Indexer.longRefTable);
			Indexer.BuildIndex(ArtiRef, RefPath, opt);
			if(FaList!="Empty")
			{
				BwtMapper Mapper(Indexer, FaList, VcfPath, popt, opt);
			}
			else
			{
				BwtMapper Mapper(Indexer, Fastq_1, Fastq_2, VcfPath, popt, opt);
			}
		}
		else //load ref index
		{
			std::cerr << "Index file exists, loading...\n";
			Indexer.LoadIndex(RefPath);
			if(FaList!="Empty")
			{
				BwtMapper Mapper(Indexer, FaList, VcfPath, popt, opt);
			}
			else
			{
				BwtMapper Mapper(Indexer, Fastq_1, Fastq_2, VcfPath, popt, opt);
			}
		}

	}
	else // using hash method
	{
		//build ref index
		//main loop for mapping

		//output alignment file
	}

	fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
	fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__,
			realtime() - t_real, cputime());

	return 0;
}

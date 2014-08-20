//============================================================================
// Name        : FastqA.cpp
// Author      : FanZhang
// Version     :
// Copyright   : Copyright Reserved.
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <string>
#include <string.h>
#include <vector>
#include <unistd.h>
//#include "./libbwa/bwtaln.h"
#include "./libbwa/utils.h"
#include "BwtIndexer.h"
#include "BwtMapper.h"
#include "params.h"
using namespace std;

#define USE_BWT 1
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.1"
#endif


BwtIndexer Indexer;


int main(int argc, char ** argv)
{
	cerr << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	//parameters handling
	double t_real;
	t_real = realtime();
	int c, opte = -1;
	gap_opt_t *opt;
	pe_opt_t *popt;
	popt = bwa_init_pe_opt();
	opt = gap_init_opt();
	/*
	while ((c = getopt(argc, argv, "n:o:e:i:d:l:k:cLR:m:t:NM:O:E:P")) >= 0)
	{
		switch (c)
		{
		case 'n':
			if (strstr(optarg, "."))
				opt->fnr = atof(optarg), opt->max_diff = -1;
			else
				opt->max_diff = atoi(optarg), opt->fnr = -1.0;
			break;
		case 'o':
			opt->max_gapo = atoi(optarg);
			break;
		case 'e':
			opte = atoi(optarg);
			break;
		case 'M':
			opt->s_mm = atoi(optarg);
			break;
		case 'O':
			opt->s_gapo = atoi(optarg);
			break;
		case 'E':
			opt->s_gape = atoi(optarg);
			break;
		case 'd':
			opt->max_del_occ = atoi(optarg);
			break;
		case 'i':
			opt->indel_end_skip = atoi(optarg);
			break;
		case 'l':
			opt->seed_len = atoi(optarg);
			break;
		case 'k':
			opt->max_seed_diff = atoi(optarg);
			break;
		case 'm':
			opt->max_entries = atoi(optarg);
			break;
		case 't':
			opt->n_threads = atoi(optarg);
			break;
		case 'L':
			opt->mode |= BWA_MODE_LOGGAP;
			break;
		case 'R':
			opt->max_top2 = atoi(optarg);
			break;
		case 'c':
			opt->mode &= ~BWA_MODE_COMPREAD;
			break;
		case 'N':
			opt->mode |= BWA_MODE_NONSTOP;
			opt->max_top2 = 0x7fffffff;
			break;
		case 'P':
			opt->pair_end = true;
			break;
		default:
			return 1;
		}
	}
	if (opte > 0)
	{
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 2 > argc)
	{
		fprintf(stderr, "\n");
		fprintf(stderr,
				"Usage:   FastqA  [options] <reference.fa> <hapmap.vcf> <in.fq>\n\n");
		fprintf(stderr,
				"Options: -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n",
				BWA_AVG_ERR, opt->fnr);
		fprintf(stderr,
				"         -o INT    maximum number or fraction of gap opens [%d]\n",
				opt->max_gapo);
		fprintf(stderr,
				"         -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
		fprintf(stderr,
				"         -i INT    do not put an indel within INT bp towards the ends [%d]\n",
				opt->indel_end_skip);
		fprintf(stderr,
				"         -d INT    maximum occurrences for extending a long deletion [%d]\n",
				opt->max_del_occ);
		fprintf(stderr, "         -l INT    seed length [%d]\n", opt->seed_len);
		fprintf(stderr,
				"         -k INT    maximum differences in the seed [%d]\n",
				opt->max_seed_diff);
		fprintf(stderr,
				"         -m INT    maximum entries in the queue [%d]\n",
				opt->max_entries);
		fprintf(stderr, "         -t INT    number of threads [%d]\n",
				opt->n_threads);
		fprintf(stderr, "         -M INT    mismatch penalty [%d]\n",
				opt->s_mm);
		fprintf(stderr, "         -O INT    gap open penalty [%d]\n",
				opt->s_gapo);
		fprintf(stderr, "         -E INT    gap extension penalty [%d]\n",
				opt->s_gape);
		fprintf(stderr,
				"         -R INT    stop searching when there are >INT equally best hits [%d]\n",
				opt->max_top2);
		fprintf(stderr,
				"         -c        input sequences are in the color space\n");
		fprintf(stderr,
				"         -L        log-scaled gap penalty for long deletions\n");
		fprintf(stderr,
				"         -N        non-iterative mode: search for all n-difference hits (slooow)\n");
		fprintf(stderr, "\n");
		return 1;
	}*/


/*
 * Parameters
 *
 */

	string RefPath, VcfPath, MaskPath("Empty"), Fastq_1("Empty"), Fastq_2("Empty"), OutPrefix, DepthDist,SitePileup;

	bool loggap(0),compread(0),nonstop(0),IL13(0);
	paramList pl;

	  BEGIN_LONG_PARAMS(longParameters)
	    LONG_PARAM_GROUP("Input/Output Files","Input/Output files for the program[Complete Path Recommended]")
	    LONG_STRING_PARAM("vcf",&VcfPath,"Input VCF file name")
	    LONG_STRING_PARAM("ref",&RefPath,"Reference FASTA file.")
	    LONG_STRING_PARAM("mask",&MaskPath,"Repeat Mask FASTA file.")
	    LONG_STRING_PARAM("fastq_1",&Fastq_1,"Pair end 1 fastq file.")
	    LONG_STRING_PARAM("fastq_2",&Fastq_2,"Pair end 2 fastq file.[Leave empty if using single end]")
	    LONG_STRING_PARAM("out",&OutPrefix,"Output file prefix")
	    //LONG_STRING_PARAM("out",&outf,"Output file prefix")
	    LONG_PARAM_GROUP("Parameters for Hash ", "Parameters the are universal for construction of HashTable.")
	    LONG_INT_PARAM("K",&Indexer.RollParam.kmer_size,"kmer size for Rolling Hash filtering")
	    LONG_INT_PARAM("T",&Indexer.RollParam.thresh,"threshold for filtering specific read")
	    LONG_INT_PARAM("S",&Indexer.RollParam.read_step_size,"step size for extracting kmer")
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
	    LONG_STRING_PARAM("depth_dist",&DepthDist,"Output file for depth distribution ")
	    LONG_STRING_PARAM("site_pileup",&SitePileup,"Output file for Pileup information on specified sites ")


	  END_LONG_PARAMS();

	  pl.Add(new longParams("Available Options", longParameters));
	  pl.Read(argc, argv);
	  pl.Status();



		if (opte > 0) {
			opt->max_gape = opte;
			opt->mode &= ~BWA_MODE_GAPE;
		}
		if(nonstop)
		{
			opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0x7fffffff;
		}
		if(IL13)
		{
			opt->mode |= BWA_MODE_IL13;
		}
		if(loggap)
		{
			opt->mode |= BWA_MODE_LOGGAP;
		}
	//Extract ref seqs from VCF and Ref

	//string RefPath(argv[optind]), VcfPath(argv[optind + 1]), FastQ(argv[optind + 2]);



	//const char *fn_fa;
	struct stat sb;


	if (USE_BWT) // using bwt methods
	{
		//build ref index
		string BwtPath = RefPath + ".bwt";
		string BwtPathR = RefPath + ".rbwt";
		if (stat(BwtPath.c_str(), &sb) != 0) //|| stat(BwtPathR.c_str(), &sb)!=0)
		{
			cerr << "Index file doesn't exist, building...\n";
			RefBuilder ArtiRef(VcfPath, RefPath, MaskPath,opt);
			Indexer.BuildIndex(ArtiRef, RefPath, opt);
			BwtMapper Mapper(Indexer,Fastq_1, Fastq_2,  popt, opt);
		}
		else //load ref index
		{
			cerr << "Index file exists, loading...\n";
			Indexer.LoadIndex(RefPath);
			BwtMapper Mapper(Indexer,Fastq_1, Fastq_2, popt, opt);
		}

		//output alignment file
	}
	else // using hash method
	{
		//build ref index
		//main loop for mapping

		//output alignment file
	}

	//keep records for pileup-like vcf file
	//do statistic
	fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
	fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__,
			realtime() - t_real, cputime());

	return 0;
}

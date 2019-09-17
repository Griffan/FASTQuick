[![Build Status](https://travis-ci.org/Griffan/FASTQuick.png?branch=master)](https://travis-ci.org/Griffan/FASTQuick)
[![GitHub Downloads](https://img.shields.io/github/downloads/Griffan/FASTQuick/total.svg?style=flat)](https://github.com/Griffan/FASTQuick/releases)
### NAME
   FASTQuick, a Fastq file based **ultra-rapid** QC tool which skips full-genome alignment and incoporates **ancestry estimation**, **contamination estimation** and variety of QC analysis. 
   
   
### Tutorial
   For more detailed tutorial information please refer to wiki page:[https://github.com/Griffan/FASTQuick/wiki]
   
   
### CONTENTS

- [SYNOPSIS](#synopsis)
- [DESCRIPTION](#description)
- [DOWNLOAD AND INSTALL](#download-and-install)
- [COMMANDS AND OPTIONS](#commands-and-options)
- [EXAMPLES](#examples)
- [USEFUL TIPS](#useful-tips)
- [BUGS](#bugs)
- [AUTHOR](#author)
- [COPYRIGHT](#copyright)

### SYNOPSIS

Below is a simple example to demonstrate its usage, the reference version and dbSNP version are not limited to hg19.

```
FASTQuick index --siteVCF hapmap.test.vcf.gz --dbsnpVCF dbsnp.test.vcf.gz --ref hg19.fa --out_prefix reduced_ref_index

FASTQuick align  --index_prefix reduced_ref_index --fq_list NA12878.fq.list --out_prefix NA12878 

FASTQuick pop+con --BamFile NA12878.bam --SVDPrefix resource/hapmap_3.3.b37.dat --Reference hg19.fa
```
### DESCRIPTION
   FASTQuick is designed for fast quality control analysis of fastq files. It rapidly map reads to selected region and generate a variety of quality control statistics. In principal, you can choose any common genetic variants list for your dataset. 
### DOWNLOAD AND INSTALL
   git clone https://github.com/Griffan/FASTQuick.git
   
   mkdir build
   
   cd build
   
   cmake ..
```
In case any required libraries is missing, you may specify customized installing path by replacing "cmake .." with:

For libhts:
  - cmake -DHTS_INCLUDE_DIRS=/hts_absolute_path/include/  -DHTS_LIBRARIES=/hts_absolute_path/lib/libhts.a ..

For bzip2:
  - cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

For lzma:
  - cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..
```
   
   make
   
   make test

   Installation complete successfully if pass all  tests finished successfully.

####Notice that if you use docker to deploy, the minimal memory requirement is 4GB.

### EXAMPLES

You can find various example scripts in example directory as templates for your own usage.

For example:
 * the script **example.sh** is the template for one-stop analysis.
 * the script **example.index.sh** is the template for selection new marker set and indexing reference data structures.
 * the script **example.align.sh** is the template for primary analysis.
 * the script **example.pop+con.sh** is the template to estimate contamination level and genetic ancestry of the intended sample.
 * the script **example.predefine.marker.index.sh** is the template to use pre-defined marker set to build indices.
   

### COMMANDS AND OPTIONS

**index**

    FASTQuick index --siteVCF [hapmap site vcf] --dbsnpVCF [dbsnp site vcf]  --ref [reference fasta]  --flank_len [250] --var_short [9000] --flank_long_len [1000] --var_long [1000] --mask [repeat_mask.fasta] --out_prefix [reduced_ref_index]

Index database sequences, using known variant sites to anchor informative region.
```
OPTIONS
--siteVCF        [String] VCF file with candidate variant sites(if predefinedVCF not specified)
--predefinedVCF  [String] VCF file with predefined variant sites(if siteVCF not specified)
--regionList     [String] Bed file with target region list [Optional]")
--dbsnpVCF       [String] VCF file with dbsnp site 
--ref            [String] Fasta file with reference genome
--mask           [String] Fasta file with repetitive region mask 
--flank_len      [Int] Flanking region length of short-flanking-region variant
--var_short      [Int] Number of short-flanking-region variant
--flank_long_len [Int] Flanking region length of long-flanking-region variant
--var_long       [Int] Number of long-flanking-region variant
--out_prefix     [String] Prefix of all the output index files
```

In principal, you can choose any common genetic variants database for --siteVCF; or if you have a specific list of variants, you can specify --predefinedVCF to skip marker slelection stage. 

To further simplify this step, we also provided a bundle of resource files with pre-defined genetic variant list in $(FASTQUICK_HOME)/resource/ directory(site-only vcf files).

A example for repeat_mask.fasta: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/StrictMask/20140520.allAutosome.strict_mask.fasta.gz


**align**

    FASTQuick align --index_prefix [reduced_ref_index] --fq_list [sample’s fastq list file] [--bam_out] [--cal_dup] [--I] --t [2]  --out_prefix [NA12878] --frac_samp [1.0] 

Align short reads 70~300 bp to selected reference region to generate comprehensive quality control related statistics in very short time.

```    
OPTIONS 
--fq_list        [String] Path of fastq file list, format: [path of pair-end 1]\t[path of pair-end 2(optional)]
--bam_in         [String] Input reads in bam format
--index_prefix   [String] Input Prefix of all the index files
--out_prefix     [String] Prefix of variety of output files
--n              [Float] Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. [0.04]
--o              [Int] Maximum number of gap opens [1]
--e              [Int] Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps) [-1]
--i              [Int] Disallow an indel within INT bp towards the ends [5]
--d              [Int] Disallow a long deletion within INT bp towards the 3’-end [16]
--l              [Int] Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. For long reads, this option is typically ranged from 25 to 35 for ‘-k 2’. [32]
--k              [Int] Maximum edit distance in the seed [2]
--m              [Int] Maximum entries in the queue [2000000]
--t              [Int] Number of threads [1]
--R              [Int] Stop searching when there are >INT equally best hits [30]
--q              [Int] Quality threshold for read trimming down to 35bp [0]
--RG             [String] Read group name
--N              [Bool] Non-iterative mode: search for all n-difference hits
--NonI           [Bool] Input fastq quality is sanger format 
--L              [Bool] Log-scaled gap penalty for long deletions 
--max_isize      [Int] Maximal insert size for a read pair to be considered being mapped properly.[500] 
--max_occ        [Int] Maximum occurrences of a read for pairing. A read with more occurrences will be treated as a single-end read. Reducing this parameter helps faster pairing. [100000]
--is_sw          [Bool] Enable Smith-Waterman for the unmapped mate.[True]
--n_multi        [Int] Maximum hits to output for paired reads.[3]
--N_multi        [Int] Maximum hits to output for discordant pairs.[10]
--ap_prior       [Float] Prior of chimeric rate (lower bound)[1.0e-05]
--force_isize    [Bool] Disable insert size estimate.[False]
--cal_dup        [Bool] Calculate PCR duplicated reads in all the statistics.[False]
--frac_samp      [Float] Overall reads downsampling rate.[1]
```
This step utilize a optimized version of BWA to rapidly screen unmap reads without losing sensitivity to potential mismatches on reads.

Each line of fq_list contains full path of each fastq file. Pair-end fastq files are in the same line and delimated by tab; single-end fastq file occpies whole line.


**pop+con**

    FASTQuick pop+con --BamFile NA12878.bam --SVDPrefix resource/hapmap_3.3.b37.dat --Reference hg19.fa
Jointly estimate sample genetic ancestry and contamination rate(VerifyBAMID2)

```
OPTIONS
--SVDPrefix      [String] SVD related files prefix(normally shared by .UD, .mu and .bed files)[Required]
--BamFile        [String] Bam or Cram file for the sample[Required]
--Reference      [String] reference file[Required]
--Seed           [Int] Random number seed(default:12345)
--NumPC          [Int] Number of Principal Components used in estimation
--NumThread      [Int] Set number of threads in likelihood calculation[default:4]
--FixPC          [String] Specify known PC coordinates for the sample[format PC1:PC2:PC3...]
--FixAlpha       [Double] Specify known contamination level
--WithinAncestry [Bool] Enabling withinAncestry assume target sample and contamination source are from the same populations,[default:betweenAncestry] otherwise")
--KnownAF        [String] A Bed file that provide known allele frequency for each marker, similar behaviour with VerifyBamID 1.0
--Epsilon        [Double] Minimization procedure convergence threshold, usually a trade-off bettween accuracy and running time[default:1e-10]
--OutputPileup   [Bool] If output temp pileup file
--Verbose        [Bool] If print the progress of the method on the screen
/*To construct SVDPrefix auxillary files*/
--RefVCF         [String] Reference panel VCF with genotype information, for generation of .UD .mu .bed files[Optional]
``` 

This step is a wrapper of our previous tool VerifyBAMID2[https://github.com/Griffan/VerifyBamID].

### Generate Final Report

Usage:
```
Rscript $(FASTQUICK_HOME)/bin/RPlotScript.R <FASTQuick align out_prefix>  <SVDPrefix> <FASTQuickInstallDir>
```

Example:
```
Rscript $(FASTQUICK_HOME)/bin/RPlotScript.R regular_size_predefine_b37_10k_NWD315195 1000g.phase3.10k.b37.vcf.gz.dat ~/Downloads/FASTQuick/resource/
```

### More Details on Generating PC plot

The final report generated above will include the PC plot to indicate sample ancestry.

If you want to visualize customized background population information, the PC coordinates files(ending with .V) in ``$(FASTQUICK_HOME)/resource/`` might help you which
provides background PC points of 1000 Genomes Project samples(e.g. 1000g.100k.b38.vcf.gz.dat.V) or of Human Genome Diversity Project samples(e.g. hgdp.100k.b38.vcf.gz.dat.V)

You can use this script to generate PC plot with customized dataset as background points, for example:
```
sh $(FASTQUICK_HOME)/bin/run.plot.sh -i ./resource/hapmap_3.3.b37.dat.V -o ./resource/hapmap.test -r 1000g -g grey
```
You may run ``sh $(FASTQUICK_HOME)/bin/run.plot.sh -h`` for further help.

   
   
### USEFUL TIPS

#### Workflow
   The recommended flow is first indexing your reference and then align your fastq file to this reference, and then infer the population identity and the contamination level.
#### Resource file preparation.
   FASTQuick was released along with pre-selected variant sites for information collection, which could be found in resource directory. 
   
   If you want to use your own abitrary variant sites, you can also refer to https://github.com/Griffan/VerifyBamID for resource files preparation.
   
   or you can find generate_new_matrix.sh in $(FASTQUICK_HOME)/bin/ to update your own variant list and then you can update everything you need with the auxiliary tools in bin directory. 
   
   
### BUGS
   List known bugs.
### AUTHOR
Fan Zhang (email:fanzhang@umich.edu)
### COPYRIGHT
   The full FASTQuick package is distributed under MIT License.



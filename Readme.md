[![Build Status](https://travis-ci.org/Griffan/FASTQuick.png?branch=master)](https://travis-ci.org/Griffan/FASTQuick)
[![GitHub Downloads](https://img.shields.io/github/downloads/Griffan/FASTQuick/total.svg?style=flat)](https://github.com/Griffan/FASTQuick/releases)
### OVERVIEW
   FASTQuick is an **ultra-fast** QC tool for NGS sequencing fastq files. It generates a comprehensive list of QC statistics, including **ancestry estimation** and **contamination estimation**, at 50x faster turnaround time.
   
### CONTENTS

- [QUICK START](#quick-start)
- [SYNOPSIS](#synopsis)
- [DESCRIPTION](#description)
- [INSTALL](#install)
- [EXAMPLES](#examples)
- [RESOURCE FILES](#resource-files)
- [WIKI PAGE](#wiki-page)
- [COMMANDS AND OPTIONS](#commands-and-options)
- [RECOMMENDED WORKFLOW](#recommended-workflow)
- [BUGS](#bugs)
- [AUTHOR](#author)
- [COPYRIGHT](#copyright)

### QUICK START
First, to start using FASTQuick, **clone the repository** and refer to [INSTALL](#download-and-install) to install FASTQuick first.

To perform a **test run** FASTQuick with a very small-sized example to understand how the software tool works, see [EXAMPLES](#examples).

To run FASTQuick **with your own FASTQ files**, you need to **download** the [RESOURCE FILES](#resource-files) first. 

For simplicity, we prepared an all-in-one script to process the whole FASTQuick pipeline or choose any start point of the pipeline (All | AllButIndex | Index | Align | Contamination | Ancestry | Visualize) in one command line.

```
${FASTQuick_HOME}/bin/FASTQuick.sh 
--steps All \
--reference <hs37d5.fa> \
--dbSNP <dbsnp132_20101103.vcf.gz> \
--callableRegion <20141020.strict_mask.whole_genome.bed> \
--output <output.prefix> \
--fastqList <input.fq.list> \
--candidateVCF <candidate.variant.vcf.gz> \
[--targetRegion <targetRegion.bed>] 
```

**Note** that you only need to build indices once, hence "--steps AllButIndex" should be the preferred option once indices are ready.

##### RESOURCE FILES

You can download commonly used resource files from:

**reference genome**(**--reference**) 

[ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz](http://tinyurl.com/jvflzg3)

**dbSNP VCF**(**--dbSNP**) 

[ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp//technical/reference/dbsnp132_20101103.vcf.gz](http://tinyurl.com/sl2kgof)

**1000 strict masked region**(**--callableRegion**) 

[ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed](http://tinyurl.com/sjhb5nn)

All these resource files are in the version of build37/hg19, which should be sufficient for the purpose of QC. 

**Note** that if other reference genomes are needed, input of **--dbSNP** and **--callableRegion** are also required to be the same build version with the reference genome.


**_INPUT FILES_**

**--fastqList** expects tab-delimited format as follows:

```
read.group.A.read_1.fq.gz   read.group.A.read_2.fq.gz
read.group.A.single.end.fq.gz
read.group.B.read_1.fq.gz   read.group.B.read_2.fq.gz
read.group.C.read_1.fq.gz   read.group.C.raed_2.fq.gz
read.group.C.single.end.fq.gz
```

**--candidateVCF** expects a list of variants with VCF format. You can provide your own candidate variant list or simply use dbsnp132_20101103.vcf.gz listed in resource files.

**Optionally**, we can enable **_target region_** mode by specifying **--targetRegion** with a bed format file(which should be the same build version as reference genome)

**_OUTPUT FILES_**

Once the process finished, you'll find summary statistics in various files starting with the same prefix(provided by **--output**). You also will find a similar [FinalReport.html](https://www.dropbox.com/s/7fbtpq82zduk4la/FinalReport.html?dl=1) in your output directory(base directory of prefix provided by **--output**). 


### SYNOPSIS


```
FASTQuick index --siteVCF hapmap.test.vcf.gz --dbsnpVCF dbsnp.test.vcf.gz --ref ref.test.fa --out_prefix test_out_ref

FASTQuick align --fq_list fq.test.list --index_prefix test_out_ref --out_prefix test_out 

FASTQuick pop+con pop+con --DisableSanityCheck --BamFile test_out.sorted.bam --SVDPrefix resource/hapmap_3.3.b37.dat --Reference ref.test.fa --Output test_out

```
You can also directly run FASTQuick without using wrapping script. Below are simple examples to demonstrate its usage, the reference version and dbSNP version are not limited to hg19.


### DESCRIPTION

   FASTQuick is designed for fast quality control analysis of fastq files. It rapidly map reads to selected region and generate a variety of quality control statistics. In principal, you can choose any common genetic variants list for your data set. 
   
### INSTALL

To install FASTQuick, run the following series of commands.

```
git clone https://github.com/Griffan/FASTQuick.git   
mkdir build
cd build
cmake ..
make   
make test
```

Installation is complete if all tests finish successfully.

```
In case any required libraries is missing, you may specify customized installing path by replacing "cmake .." with:

For libhts:
  - cmake -DHTS_INCLUDE_DIRS=/hts_absolute_path/include/  -DHTS_LIBRARIES=/hts_absolute_path/lib/libhts.a ..

For bzip2:
  - cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

For lzma:
  - cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..
```


**Note** that if you use docker to deploy, the minimal memory requirement is 4GB.

### EXAMPLES

You can find example scripts for each single step in example directory as template for customized usage.

For example:
 * the script **example.sh** is the template for one-stop analysis.
 * the script **example.index.sh** is the template for selection new marker set and indexing reference data structures.
 * the script **example.align.sh** is the template for primary analysis.
 * the script **example.pop+con.sh** is the template to estimate contamination level and genetic ancestry of the intended sample.
 * the script **example.predefine.marker.index.sh** is the template to use pre-defined marker set to build indices.

### WIKI PAGE

You can always refer to our wiki page for more detailed introduction about the design and usage of FASTQuick. [https://github.com/Griffan/FASTQuick/wiki]

### COMMANDS AND OPTIONS

**index**

    FASTQuick index --siteVCF [hapmap site vcf] --dbsnpVCF [dbsnp site vcf]  --ref [reference fasta]  --flank_len [250] --var_short [9000] --flank_long_len [1000] --var_long [1000] --callableRegion [repeat_mask.fasta] --out_prefix [reduced_ref_index]

Index database sequences, using known variant sites to anchor informative region.

```
OPTIONS
--siteVCF        [String] VCF file with candidate variant sites(if predefinedVCF not specified)
--predefinedVCF  [String] VCF file with predefined variant sites(if siteVCF not specified)
--regionList     [String] Bed file with target region list
--dbsnpVCF       [String] VCF file with dbsnp site 
--ref            [String] Fasta file with reference genome
--callableRegion [String] Masked fasta or bed file to specify callable regions 
--flank_len      [Int] Flanking region length of short-flanking-region variant
--var_short      [Int] Number of short-flanking-region variant
--flank_long_len [Int] Flanking region length of long-flanking-region variant
--var_long       [Int] Number of long-flanking-region variant
--out_prefix     [String] Prefix of all the output index files
```

In principal, you can choose any common genetic variants database for --siteVCF; or if you have a specific list of variants, you can specify --predefinedVCF to skip marker selection stage. 

To further simplify this step, we also provided a bundle of resource files with pre-defined genetic variant list in $(FASTQUICK_HOME)/resource/

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

Each line of fq_list contains full path of each fastq file. Pair-end fastq files are in the same line and delimited by tab; single-end fastq file occupies the whole line.


**pop+con**

    FASTQuick pop+con --BamFile [sample bam file] --SVDPrefix [SVD resource files prefix] --Reference [reference genome fasta file]

Jointly estimate sample genetic ancestry and contamination rate(same method as in VerifyBAMID2)

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

### Generate Final Report

Usage:
```
Rscript $(FASTQUICK_HOME)/bin/RPlotScript.R <FASTQuick align out_prefix>  <SVDPrefix> <FASTQuickInstallDir>
```

Example:
```
Rscript $(FASTQUICK_HOME)/bin/RPlotScript.R  FASTQuick_align_out_prefix $(FASTQUICK_HOME)/resource/1000g.phase3.10k.b37.vcf.gz.dat $(FASTQUICK_HOME)
```
   
### Recommended Workflow
   The FASTQuick.sh script demonstrated the recommended flow to first indexing your reference and then align your fastq file to this reference, and then infer the population identity and the contamination level.
 
   
### BUGS
   List known bugs.
   
### AUTHOR
Fan Zhang (email:fanzhang@umich.edu)

### COPYRIGHT
   The full FASTQuick package is distributed under MIT License.



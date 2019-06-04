[![Build Status](https://travis-ci.org/Griffan/FASTQuick.png?branch=master)](https://travis-ci.org/Griffan/FASTQuick)
### NAME
   FASTQuick, a Fastq file based **ultra-rapid** QC tool which incorporates the alignment, population identification, contamination estimation and variety of QC analysis. 
   
   
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
```
FASTQuick index --siteVCF hapmap.test.vcf.gz --dbsnpVCF dbsnp.test.vcf.gz --ref hg19.fa --out_prefix reduced_ref_index

FASTQuick align  --index_prefix reduced_ref_index --fq_list NA12878.fq.list --out_prefix NA12878 

FASTQuick pop+con --BamFile NA12878.bam --SVDPrefix resource/hapmap_3.3.b37.dat --Reference hg19.fa

/*Below are deprecated but still available*/
FASTQuick pop --SVDPrefix resource/hapmap_3.3.b37.dat --Pileup NA12878.Pileup.gz

FASTQuick con --SVDPrefix resource/hapmap_3.3.b37.dat --Pileup NA12878.Pileup.gz —-Out test
```
### DESCRIPTION
   FASTQuick is designed for fast quality control analysis of fastq files. It rapidly map reads to selected region and generate a variety of quality control statistics.
### DOWNLOAD AND INSTALL
   git clone https://github.com/Griffan/FASTQuick.git
   
   mkdir build
   
   cd build
   
   cmake ..
   
   make
   
   (required libraries):<http://www.htslib.org/doc/tabix.html>
### COMMANDS AND OPTIONS

**index**

    FASTQuick index --siteVCF [hapmap site vcf] --dbsnpVCF [dbsnp site vcf]  --ref [reference fasta]  --flank_len [250] --var_short [9000] --flank_long_len [1000] --var_long [1000] --mask [repeat_mask.fasta] --out_prefix [reduced_ref_index]

Index database sequences, using known variant sites to anchor informative region.
```
OPTIONS
--siteVCF        [String] Path of VCF file with candidate variant sites(if predefinedVCF not specified)
--dbsnpVCF       [String] Path of VCF file with dbsnp site 
--ref            [String] Path of fasta file with reference genome
--mask           [String] Path of fasta file with repetitive region mask 
--predefinedVCF  [String] Path of VCF file with predefined marker set(if siteVCF not specified)
--flank_len      [Int] Flanking region length of short-flanking-region variant
--var_short      [Int] Number of short-flanking-region variant
--flank_long_len [Int] Flanking region length of long-flanking-region variant
--var_long       [Int] Number of long-flanking-region variant
--out_prefix     [String] Prefix of all the output index files
```
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
/*Below are deprecated but still available*/
--UDPath         [String] .UD matrix file from SVD result of genotype matrix[Required]
--MeanPath       [String] .mu matrix file of genotype matrix[Required]
--BedPath        [String] .Bed file for markers used in this analysis,format(chr\tpos-1\tpos\trefAllele\taltAllele)[Required]
```    
**pop**

    FASTQuick pop --SVDPrefix [resource/hapmap_3.3.b37.dat] --Pileup [NA12878.Pileup.gz]
Identify individual’s population identity, ancestry information. The geometric distance in plot represents how close the relatedness is.
```
OPTIONS
--SVDPrefix      [String] Specify the prefix used by SVD matrices. If you are using FASTQuick default marker set, you may find them in resource directory.[Required]
--GL             [String] Input genotype likelihood file generated from align step.[Required if no pileup file]
--Pileup         [String] Input pileup file generated from align[Required if no gl file]
```
**con**

    FASTQuick con --SVDPrefix [resource/hapmap_3.3.b37.dat] --Pileup [NA12878.Pileup.gz]
Estimate the probability that this sample is contaminated with other genomic material.
```
OPTIONS
--SVDPrefix      [String] Specify the prefix used by SVD matrices. If you are using FASTQuick default marker set, you may find them in resource directory.[Required if VCF file doesn't provide site allele frequency]
--VCF            [String] Specify VCF file that contains site allele frequency.[Required if SVD matrices don't exist]
--Pileup         [String] Specify pileup file for current individual, could be the one generated from align step.[Required]
--Out            [String] Specify the output prefix.[Required]
```
### Generate Final Report

```
Rscript ../bin/RPlotScript.R [FASTQuick align out_prefix]
```

### EXAMPLES
   Some examples of common usage.
   See wiki page tutorial.
   [https://github.com/Griffan/FASTQuick/wiki]
   
   
### USEFUL TIPS
   The recommended flow is first indexing your reference and then align your fastq file to this reference, and then infer the population identity or you infer the contamination level.
   
   FASTQuick was released along with pre-selected variant sites for information collection, which could be found in resource directory. If you want to use your own abitrary variant sites, you may look into the bin directory to use generate_new_matrix.sh to update your own variants set and then you can update everything you need with the auxiliary tools in bin directory.
### BUGS
   List known bugs.
### AUTHOR
Fan Zhang (email:fanzhang@umich.edu)
### COPYRIGHT
   The full FASTQuick package is distributed under MIT License.



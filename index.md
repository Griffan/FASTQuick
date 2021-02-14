###NAME
   FASTQuick, a Fastq file based Population identification and Contamination detection tool
###CONTENTS

###SYNOPSIS
```
FASTQuick index --vcf hapmap.vcf --dbsnp 00-All.vcf.gz  --ref hs37d5.fa --flank_len 250 --var_short 9000 --flank_long_len 1000 --var_long 1000 --mask 20141007.all.strict_mask.fasta 

FASTQuick align --ref hs37d5.fa --fq_list NA12878.fq.list --bam_out --cal_dup --flank_len 250 --var_short 9000 --flank_long_len 1000 --var_long 1000  --I --t 2  --prefix NA12878 --frac_samp 1.0 

FASTQuick pop --UD resource/hapmap.dat.UD --PC resource/hapmap.dat.V --mu resource/hapmap.dat.mu --gl NA12878.likelihood --bed resource/choose.bed.post.bed.allele.bed

FASTQuick con --prefix NA12878
```
###DESCRIPTION
   FASTQuick is short for fastq file based population identification and contamination detection tool. It is designed for fast quality control analysis of fastq files. It rapidly map reads to selected region and generate a variety of quality control statistics.
###COMMANDS AND OPTIONS

**index**	

    FASTQuick index --vcf [hapmap site vcf] --dbsnp [dbsnp site vcf]  --ref [reference fasta]  --flank_len [250] --var_short [9000] --flank_long_len [1000] --var_long [1000] --mask [repeat_mask.fasta] 

Index database sequences, using known variant sites to anchor informative region.

    OPTIONS
    --vcf	STR	path of input hapmap site vcf file 
    --dbsnp	STR	path of input dbsnp site vcf file
    --ref	STR	path of reference genome fasta file
    --mask	STR	path of repetitive region  mask fasta file
    --flank_len	INT	flanking region length of short-flanking-region variant
    --var_short	INT	number of short-flanking-region variant
    --flank_long_len	INT flanking region length of long-flanking-region variant
    --var_long	INT	number of long-flanking-region variant

**align**

    FASTQuick align --ref [reference genome fasta file] --fq_list [sample’s fastq list file] [--bam_out] [--cal_dup] --flank_len [250] --var_short [9000] --flank_long_len [1000] --var_long [1000]  [--I] --t [2]  --prefix [NA12878] --frac_samp [1.0]

Align short reads 70~300 bp to selected reference region to generate comprehensive quality control related statistics in very short time.
    
    OPTIONS
    --ref	STR	path of reference genome fasta file
    --fq_list	STR path of fastq file list, format: [path of pair-end 1]\t[path of pair-end]
    --fastq_1	STR path of pair end 1 fastq file
    --fastq_2	STR path of pair end 2 fastq file
    --bam_in	STR path of already aligned bam file
    --bam_out	BOOL	output is bam format or not
    --prefix	STR	prefix of variety of output file associated with specific sample
    --flank_len	INT	flanking region length of short-flanking-region variant
    --var_short	INT	number of short-flanking-region variant
    --flank_long_len	INT flanking region length of long-flanking-region variant
    --var_long	INT	number of long-flanking-region variant
    --n	float	Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. [0.04]
    --o	INT	Maximum number of gap opens [1]
    --e	INT	Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps) [-1]
    --i	INT	Disallow an indel within INT bp towards the ends [5]
    --d	INT	Disallow a long deletion within INT bp towards the 3’-end [16]
    --l	INT	Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. For long reads, this option is typically ranged from 25 to 35 for ‘-k 2’. [inf]
    --k	INT	Maximum edit distance in the seed [2]
    --m	INT	Maximum entries in the queue [2000000]
    --t	INT	Number of threads [1]
    --R	INT	Stop searching when there are >INT equally best hits [30]
    --q	INT	Quality threshold for read trimming down to 35bp [0]
    --RG	STR	Read group name
    --N	BOOL	Non-iterative mode: search for all n-difference hits
    --I	BOOL	The input is in the Illumina 1.3+ read format (quality equals ASCII-64)
    --L	BOOL	Log-scaled gap penalty for long deletions 
    --max_isize	INT	Maximal insert size for a read pair to be considered being mapped properly.[500] 
    --max_occ	INT	Maximum occurrences of a read for pairing. A read with more occurrences will be treated as a single-end read. Reducing this parameter helps faster pairing. [100000]
    --is_sw	BOOL	Enable Smith-Waterman for the unmapped mate.[True]
    --n_multi	INT	Maximum hits to output for paired reads.[3]
    --N_multi	INT	Maximum hits to output for discordant pairs.[10]
    --ap_prior	FLOAT	Prior of chimeric rate (lower bound)[1.0e-05]
    --force_isize	BOOL	Disable insert size estimate.[False]
    --cal_dup	BOOL	Calculate PCR duplicated reads in all the statistics.[False]
    --frac_samp	FLOAT	Overall reads downsampling rate.[1]
**pop**

    FASTQuick pop --UD [resource/hapmap.dat.UD] --PC [resource/hapmap.dat.V] --mu [resource/hapmap.dat.mu] --gl [prefix.likelihood] --bed [resource/choose.bed.post.bed.allele.bed]
Identify individual’s population identity, ancestry information. The geometric distance in plot represents how close the relatedness is.

    OPTIONS
    --UD	STR	Path of UD matrix in resource directory
    --PC	STR	Path of V matrix in resource directory
    --mu	STR	Path of mu matrix in resource directory
    --gl	STR	Path of output likelihood file in align step, with prefix specified in align step and suffix likelihood.
    --bed	STR	Bed format file that specified markers used in pop inference, also can be found in resource directory.
**con**

    FASTQuick con --prefix [NA12878]
Estimate the probability that this sample is contaminated with other genomic material.

    OPTIONS
    --prefix STR Specify the prefix used in previous steps, which will be used to retrieve all the information needed in this step.
###EXAMPLES
   Some examples of common usage.
   
    src/FASTQuick index --vcf new_sites/All.hapmap.omni.HDGP.recode.vcf --dbsnp dbSNP/b137/00-All.vcf.gz  --ref hs37d5.fa --flank_len 250 --var_short 9000 --flank_long_len 1000 --var_long 1000 --mask 20141007.all.strict_mask.fasta
    
    src/FASTQuick align --ref hs37d5.fa --fq_list NA12878.fq.list --bam_out --cal_dup --flank_len 250 --var_short 9000 --flank_long_len 1000 --var_long 1000  --I --t 2  --prefix NA12878 --frac_samp 1
###Useful Tips
   The recommended flow is first indexing your reference and then align your fastq file to this reference, and then infer the population identity or you infer the contamination level.
   
   FASTQuick was released along with pre-selected variant sites for information collection, which could be found in resource directory. If you want to use your own abitrary variant sites, you may look into the bin directory to use generate_new_matrix.sh to update your own variants set and then you can update everything you need with the auxilary tools in bin directory.
###BUGS
   List known bugs.
###AUTHOR
Fan Zhang (email:fanzhang@umich.edu)
###COPYRIGHT
   The full FASTQuick package is distributed under GPLv3.



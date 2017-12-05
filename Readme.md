### NAME
   FASTQuick, a Fastq file based **ultra-rapid** QC tool which incorporates the alignment, population identification, contamination estimation and variety of QC analysis. 
   
   
### Tutorial
   For more detailed tutorial information please refer to wiki page:[https://github.com/Griffan/FASTQuick/wiki]
   
   
### CONTENTS

- [SYNOPSIS](#synopsis)
- [DESCRIPTION](#description)
- [COMMANDS AND OPTIONS](#commands-and-options)
- [EXAMPLES](#examples)
- [USEFUL TIPS](#useful-tips)
- [BUGS](#bugs)
- [AUTHOR](#author)
- [COPYRIGHT](#copyright)

### SYNOPSIS
```
FASTQuick index --siteVCF hapmap.test.vcf.gz --dbsnpVCF dbsnp.test.vcf.gz --ref test.fa --out_index_prefix reduced_ref_index

FASTQuick align  --in_index_prefix reduced_ref_index --fq_list NA12878.fq.list --out_prefix NA12878 

FASTQuick pop --SVD_prefix resource/hapmap_3.3.b37.dat --pileup NA12878.Pileup.gz --BED resource/choose.bed

FASTQuick con --SVD_prefix resource/hapmap_3.3.b37.dat --pileup NA12878.Pileup.gz --BED resource/choose.bed —out test
```
### DESCRIPTION
   FASTQuick is designed for fast quality control analysis of fastq files. It rapidly map reads to selected region and generate a variety of quality control statistics.
### COMMANDS AND OPTIONS

**index**	

    FASTQuick index --siteVCF [hapmap site vcf] --dbsnpVCF [dbsnp site vcf]  --ref [reference fasta]  --flank_len [250] --var_short [9000] --flank_long_len [1000] --var_long [1000] --mask [repeat_mask.fasta] --out_index_prefix [reduced_ref_index]

Index database sequences, using known variant sites to anchor informative region.

    OPTIONS
    --siteVCF	STR	Path of selected Sites VCF file,e.g. hapmap vcf
    --dbsnpVCF	STR	Path of input dbsnp site vcf file
    --ref	STR	Path of reference genome fasta file
    --mask	STR	Path of repetitive region  mask fasta file
    --predefinedVCF STR path of predefined marker set vcf file
    --flank_len	INT	Flanking region length of short-flanking-region variant
    --var_short	INT	Number of short-flanking-region variant
    --flank_long_len	INT Flanking region length of long-flanking-region variant
    --var_long	INT	Number of long-flanking-region variant
    --out_index_prefix   STR   Prefix of all the output index files

**align**

    FASTQuick align --in_index_prefix [reduced_ref_index] --fq_list [sample’s fastq list file] [--bam_out] [--cal_dup] [--I] --t [2]  --out_prefix [NA12878] --frac_samp [1.0] 

Align short reads 70~300 bp to selected reference region to generate comprehensive quality control related statistics in very short time.
    
    OPTIONS
    
    --fq_list	STR path of fastq file list, format: [path of pair-end 1]\t[path of pair-end]
    --fastq_1	STR path of pair end 1 fastq file
    --fastq_2	STR path of pair end 2 fastq file
    --bam_in	STR path of already aligned bam file
    --sam_out Bool  If output sam file[default output bam file]
    --in_idx_prefix	STR	Input prefix of all the index files
    --out_prefix	STR	prefix of variety of output files
    --n	FLOAT	Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. [0.04]
    --o	INT	Maximum number of gap opens [1]
    --e	INT	Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps) [-1]
    --i	INT	Disallow an indel within INT bp towards the ends [5]
    --d	INT	Disallow a long deletion within INT bp towards the 3’-end [16]
    --l	INT	Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. For long reads, this option is typically ranged from 25 to 35 for ‘-k 2’. [32]
    --k	INT	Maximum edit distance in the seed [2]
    --m	INT	Maximum entries in the queue [2000000]
    --t	INT	Number of threads [1]
    --R	INT	Stop searching when there are >INT equally best hits [30]
    --q	INT	Quality threshold for read trimming down to 35bp [0]
    --RG	STR	Read group name
    --N	BOOL	Non-iterative mode: search for all n-difference hits
    --NonI  Bool  the input fastq quality is sanger format 
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

    FASTQuick pop --SVD_prefix [resource/hapmap_3.3.b37.dat] --pileup [NA12878.Pileup.gz] --BED [resource/choose.bed]
Identify individual’s population identity, ancestry information. The geometric distance in plot represents how close the relatedness is.

    OPTIONS
    --SVD_prefix	STR	Specify the prefix used by SVD matrices. If you are using FASTQuick default marker set, you may find them in resource directory.[Required]
    --gl	STR	Input genotype likelihood file generated from align step.[Required if no pileup file]
    --pileup STR Input pileup file generated from align[Required if no gl file]
    --BED	STR	Specify the matching BED format file that contains marker information. If you are using FASTQuick default marker set, you may find choose.bed file in resource directory.[Required]
**con**

    FASTQuick con --SVD_prefix [resource/hapmap_3.3.b37.dat] --pileup [NA12878.Pileup.gz] --BED [resource/choose.bed]
Estimate the probability that this sample is contaminated with other genomic material.

    OPTIONS
    --SVD_prefix  STR	Specify the prefix used by SVD matrices. If you are using FASTQuick default marker set, you may find them in resource directory.[Required if VCF file doesn't provide site allele frequency]
    --VCF   STR   Specify VCF file that contains site allele frequency.[Required if SVD matrices don't exist]
    --pileup   STR   Specify pileup file for current individual, could be the one generated from align step.[Required]
    --BED   STR   Specify the matching BED format file that contains marker information, which should match markers in SVD matrices.[Required]
    --out   STR   Specify the output prefix.[Required]
    --RG STR   set ReadGroup name
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



samtools sort -O BAM regular_size_predefine_b37_10k_output.bam >regular_size_predefine_b37_10k_output.bam.sort.bam
samtools index regular_size_predefine_b37_10k_output.bam.sort.bam
bin/FASTQuick pop+con --Reference hs37d5.fa --BamFile regular_size_predefine_b37_10k_output.bam.sort.bam --SVDPrefix ./resource/1000g.phase3.10k.b37.vcf.gz.dat --DisableSanityCheck --Output regular_size_predefine_b37_10k_output

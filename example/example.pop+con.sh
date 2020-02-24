samtools sort test_out.bam >test_out.sorted.bam
samtools index test_out.sorted.bam
../bin/FASTQuick pop+con --DisableSanityCheck --BamFile test_out.sorted.bam --SVDPrefix ../resource/hapmap_3.3.b37.dat --Reference ref.test.fa --Output test_out&>test_out.e.log

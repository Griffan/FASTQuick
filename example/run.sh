../bin/FASTQuick index --siteVCF hapmap.test.vcf.gz --dbsnpVCF dbsnp.test.vcf.gz --ref test.fa --out_index_prefix NA12878
../bin/FASTQuick align --fq_list test.fq.list --in_index_prefix NA12878 --out_prefix test_out
../bin/FASTQuick pop --SVD_prefix ../resource/hapmap_3.3.b37.dat --pileup fake.Pileup.gz --BED ../resource/choose.bed
../bin/FASTQuick con --SVD_prefix ../resource/hapmap_3.3.b37.dat --pileup fake.Pileup.gz --BED ../resource/choose.bed --out test
if [ $? -ne 0 ]	
then exit $?;
else echo "Test done successfully.";
fi

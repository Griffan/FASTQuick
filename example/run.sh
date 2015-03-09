../src/FASTQuick index --siteVCF hapmap.test.vcf.gz --dbsnpVCF dbsnp.test.vcf.gz --ref test.fa --out_idx_prefix NA12878
../src/FASTQuick align --fq_list test.fq.list --in_idx_prefix NA12878 --out_prefix test_out
../src/FASTQuick pop --UD ../resource/hapmap_3.3.b37.dat.UD --PC ../resource/hapmap_3.3.b37.dat.V --mu ../resource/hapmap_3.3.b37.dat.mu --gl fake.likelihood --bed ../resource/choose.bed.post.bed.allele.bed
../src/FASTQuick con --in_idx_prefix NA12878  --in_prefix fake 
if [ $? -ne 0 ]	
then exit $?;
	else echo "Test done successfully.";
	fi

../src/FASTQuick index --hapmap hapmap.test.vcf.gz --dbsnp dbsnp.test.vcf.gz --ref test.fa --index_prefix NA12878
../src/FASTQuick align --fq_list test.fq.list --index_prefix NA12878 --prefix test_out
if [ $? -ne 0 ]	
then exit $?;
	else echo "Test done successfully.";
	fi

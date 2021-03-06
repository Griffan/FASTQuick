# Validate tools exist on path
for tool in samtools; do
	if ! which $tool >/dev/null; then
	  echo "Error: unable to find $tool on \$PATH, please install before continue" 1>&2 ; exit 2; fi
	echo "Found $(which $tool)" 1>&2
done

../bin/FASTQuick index --siteVCF hapmap.test.vcf.gz --dbsnpVCF dbsnp.test.vcf.gz --ref ref.test.fa --out_prefix test_out_ref 1>test_out.o.log 2>test_out.e.log
if [ $? -ne 0 ]
then exit $?;
else echo "Test on index finished successfully.";
fi
../bin/FASTQuick align --fq_list fq.test.list --index_prefix test_out_ref --out_prefix test_out 1>test_out.o.log 2>test_out.e.log
if [ $? -ne 0 ]
then exit $?;
else echo "Test on align finished successfully.";
fi
samtools sort test_out.bam >test_out.sorted.bam
samtools index test_out.sorted.bam
../bin/FASTQuick pop+con --DisableSanityCheck --BamFile test_out.sorted.bam --SVDPrefix ../resource/hapmap_3.3.b37.dat --Reference ref.test.fa --Output test_out 1>test_out.o.log 2>test_out.e.log
if [ $? -ne 0 ]
then exit $?;
else
  echo "Test on pop+con finished successfully.";
  echo "All tests finished successfully.";
fi

../bin/FASTQuick index --siteVCF hapmap.test.vcf.gz --dbsnpVCF dbsnp.test.vcf.gz --ref test.fa --out_prefix test_ref &>test.e.log
if [ $? -ne 0 ]
then exit $?;
else echo "Test on index finished successfully.";
fi
../bin/FASTQuick align --fq_list test.fq.list --index_prefix test_ref --out_prefix test_out &>test.e.log
if [ $? -ne 0 ]
then exit $?;
else echo "Test on align finished successfully.";
fi
samtools index test_out.bam
../bin/FASTQuick pop+con --DisableSanityCheck --BamFile test_out.bam --SVDPrefix ../resource/hapmap_3.3.b37.dat --Reference test.fa &>test.e.log
if [ $? -ne 0 ]
then exit $?;
else echo "Test on pop+con finished successfully.";
fi
../bin/FASTQuick pop --SVDPrefix ../resource/hapmap_3.3.b37.dat --Pileup fake.Pileup.gz &>test.e.log
if [ $? -ne 0 ]
then exit $?;
else echo "Test on pop finished successfully.";
fi
../bin/FASTQuick con --SVDPrefix ../resource/hapmap_3.3.b37.dat --Pileup fake.Pileup.gz --Out test_out &>test.e.log
if [ $? -ne 0 ]
then exit $?;
else echo "Test on con finished successfully.\nAll tests finished successfully.";
fi

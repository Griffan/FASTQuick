#!/usr/bin/perl -w
use strict;
die "perl a03-extract-GLs-from-phase1.pl <resource directory> <1000g release 20110521 directory>\n" unless @ARGV==2;
my $resourceDIR=$ARGV[0];
my $glDIR=$ARGV[1];
#print @ARGV,"\n";
my %hpops = ();
my @ids = ();
#print "$resourceDIR/phase1_integrated_calls.20101123.ALL.panel\n";
open(IN,"$resourceDIR/phase1_integrated_calls.20101123.ALL.panel") || die "Cannot open file 1$!\n";
while(<IN>) {
    my ($id,$con,$pop,$pl) = split;
    $hpops{$id} = $pop;
    push(@ids,$id);
}
close IN;

open(OUT,">$resourceDIR/1kg.phase1.selected.GLs.dat") || die "Cannot open file 2$!\n";
for(my $i=0; $i < @ids; ++$i) {
    print OUT "\t" if ( $i > 0 );
    print OUT join("\t","$ids[$i]-$hpops{$ids[$i]}-0","$ids[$i]-$hpops{$ids[$i]}-1","$ids[$i]-$hpops{$ids[$i]}-2");
}
print OUT "\n";

open(IN,"<$resourceDIR/choose.bed.post.bed.allele.bed") || die "Cannot open file 3$!\n";
while(<IN>) {
    my ($chrom,$beg,$bp,$refB,$altB) = split;
#    print "tabix $glDIR/ALL.chr$chrom.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz $chrom:$bp-$bp |\n";
    open(IN2,"tabix $glDIR/ALL.chr$chrom.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz $chrom:$bp-$bp |") || die "Cannot open file\n";
    my $flag = 0;
    while(<IN2>) {
	my ($chr,$pos,$id,$ref,$alt,$qual,$filt,$info,$fmt,@F) = split;
	next unless ( ( length($ref) == 1 ) && ( length($alt) == 1 ) ); ## snps only
	next if ( $pos != $bp );

	print OUT "$chrom:$bp:$refB:$altB";
	for(my $i=0; $i < @F; ++$i) {
	    my ($gt,$ds,$gl) = split(/:/,$F[$i]);
	    my @gls = split(/,/,$gl);
	    print OUT "\t";
	    print OUT join("\t",@gls);
	}
	print OUT "\n";
	$flag = 1;
    };
    close IN2;

    if ( $flag == 0 ) {
	print OUT "$chrom:$bp:$refB:$altB";	
	for(my $i=0; $i < @ids; ++$i) {
	    print OUT "\t";
	    print OUT join("\t",0,0,0);
	}
	print OUT "\n";
    }
}
close IN;
close OUT;

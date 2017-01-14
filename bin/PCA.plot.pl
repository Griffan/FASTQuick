#!/usr/bin/perl -w

use strict;
die "Usage:
	perl PCA.plot.pl <SampleID> <PC1> <PC2> <resource dir>\n" unless @ARGV==4;
	my $sampleID=shift;
	my $PC1=shift;
	my $PC2=shift;
	my $resourceDir=shift;
my %hmpops = ();
open(IN,$resourceDir."/relationships_w_pops_041510.txt") || die "Cannot open file $!\n";
while(<IN>) {
    my ($fid,$iid,$dad,$mom,$sex,$pheno,$pop) = split;
    $pop = "MXL" if ( $pop eq "MEX" );
    $hmpops{$iid} = $pop;
}
close IN;

my %kgpops = ();
open(IN,$resourceDir."/phase1_integrated_calls.20101123.ALL.panel") || die "Cannot open file $!\n";
while(<IN>) {
    my ($id,$pop,$con,$pl) = split;
    $kgpops{$id} = [$pop,$con];
}
close IN;

my %hmpcs = ();
open(IN,$resourceDir."/hapmap_3.3.b37.dat.V") || die "Cannot open file $!\n";
while(<IN>) {
    my ($id,$pc1,$pc2) = split;
    $id = substr($id,0,7);
    $hmpcs{$id} = [$pc1,$pc2];
}
close IN;

my %kgpcs = ();
open(IN,$resourceDir."/1kg.phase1.selected.pc2.out") || die "Cannot open file $!\n";
while(<IN>) {
    my ($id,$pc1,$pc2) = split;
    $kgpcs{$id} = [$pc1,$pc2];
}
close IN;

open(OUT,">1kg.hapmap.pcs.dat") || die "Cannot open file $!\n";
print OUT join("\t","SET","ID","POPULATION","PC1","PC2","OVERLAP")."\n";
foreach my $id ( sort keys %kgpcs ) {
    print OUT join("\t","1000G.Likelihoods",$id,$kgpops{$id}->[0],@{$kgpcs{$id}},defined($hmpops{$id}) ? 1 : 0)."\n";
}
foreach my $id ( sort keys %hmpcs ) {
    print OUT join("\t","HapMap3.Genotypes",$id,$hmpops{$id},@{$hmpcs{$id}},defined($kgpops{$id}) ? 1 : 0)."\n";
}
print OUT join("\t","HapMap3.Genotypes","$sampleID","TBA",$PC1,$PC2,0)."\n";
print OUT join("\t","1000G.Likelihoods","$sampleID","TBA",$PC1,$PC2,0)."\n";

close OUT;

open(R,">1kg.hapmap.pcs.r") || die "Cannot open file\n";
print R "library(ggplot2)\n";
print R "df <- as.data.frame(read.table('1kg.hapmap.pcs.dat',header=TRUE))\n";
print R "pdf('1kg.hapmap.pcs.all.pdf',width=12,height=8)\n";
print R "ggplot(df,aes(x=PC1,y=PC2,colour=POPULATION)) + geom_point() + geom_text(data=df[df\$POPULATION==\"TBA\", ],aes(label=ID),hjust=0, vjust=0,color=\"blue\")+geom_point(data=df[df\$POPULATION==\"TBA\", ],size=4,color=\"blue\")+facet_grid(. ~ SET,scales=\"free\") + theme(legend.position='bottom')\n";
print R "dev.off()\n";
#print R "pdf('1kg.hapmap.pcs.overlap.pdf',width=8,height=6)\n";
#print R "ggplot(subset(df,OVERLAP==1),aes(PC1,PC2,colour=POPULATION)) + geom_point(size=1.5) + facet_grid(. ~ SET,scales=\"free\") + theme(legend.position='bottom')\n";
#print R "dev.off()\n";
close R;

print `Rscript 1kg.hapmap.pcs.r`;

#!/usr/bin/perl -w
use strict;
my $sampleID=shift;
my $pc1=shift;
my $pc2=shift;
my $popID="ToBeAssigned";

my %hmpops = ();
open(IN,"/net/fantasia/home/hmkang/data/GATK-resources/relationships_w_pops_041510.txt") || die "Cannot open file\n";
while(<IN>) {
    my ($fid,$iid,$dad,$mom,$sex,$pheno,$pop) = split;
    $pop = "MXL" if ( $pop eq "MEX" );
    $hmpops{$iid} = $pop;
}
close IN;

#my %kgpops = ();
#open(IN,"/net/1000g/1000g/release/20110521/phase1_integrated_calls.20101123.ALL.panel") || die "Cannot open file\n";
#while(<IN>) {
#    my ($id,$pop,$con,$pl) = split;
#    $kgpops{$id} = [$pop,$con];
#}
#close IN;

my %hmpcs = ();
open(IN,"hapmap_3.3.b37.selected.dat.bk.V") || die "Cannot open file\n";
while(<IN>) {
    my ($id,$pc1,$pc2) = split;
    $id = substr($id,0,7);
    $hmpcs{$id} = [$pc1,$pc2];
}
close IN;

#my %kgpcs = ();
#open(IN,"1kg.phase1.selected.pc2.out") || die "Cannot open file\n";
#while(<IN>) {
#    my ($id,$pc1,$pc2) = split;
#    $kgpcs{$id} = [$pc1,$pc2];
#}
#close IN;

open(OUT,">hapmap.pcs.dat") || die "Cannot open file\n";
print OUT join("\t","SET","ID","POPULATION","PC1","PC2")."\n";
#foreach my $id ( sort keys %kgpcs ) {
#    print OUT join("\t","1000G.Likelihoods",$id,$kgpops{$id}->[0],@{$kgpcs{$id}},defined($hmpops{$id}) ? 1 : 0)."\n";
#}
foreach my $id ( sort keys %hmpcs ) {
    print OUT join("\t","HapMap3.Genotypes",$id,$hmpops{$id},@{$hmpcs{$id}})."\n";
}
#print OUT join("\t","HapMap3.Genotypes",$sampleID,$popID,$pc1,$pc2)."\n";
print OUT join("\t","HapMap3.Genotypes","NA12878T","TBA",-0.0156387,0.0432906)."\n";
print OUT join("\t","HapMap3.Genotypes","HG00553T","TBA",-0.00626694,0.024657)."\n";
close OUT;
my $lastLineNum=(keys %hmpcs)+1;
open(R,">hapmap.pcs.r") || die "Cannot open file\n";
print R "library(ggplot2)\n";
print R "df <- as.data.frame(read.table('hapmap.pcs.dat',header=TRUE))\n";
print R "pdf('hapmap.pcs.all.pdf',width=12,height=8)\n";
print R "ggplot(df,aes(x=PC1,y=PC2,colour=POPULATION)) + geom_point(size=1.5) + facet_grid(. ~ SET,scales=\"free\") + geom_point(data=df[df\$POPULATION==\"CEU\",], aes(x=PC1, y=PC2), colour=\"blue\", size=2.5)+geom_point(data=df[".$lastLineNum.", ], aes(x=PC1, y=PC2), size=4.5)+geom_text(data=df[".$lastLineNum.", ],aes(label=ID),hjust=0, vjust=0)+theme(legend.position='bottom')\n";
print R "dev.off()\n";
#print R "pdf('1kg.hapmap.pcs.overlap.pdf',width=8,height=6)\n";
#print R "ggplot(subset(df,OVERLAP==1),aes(PC1,PC2,colour=POPULATION)) + geom_point(size=1.5) + facet_grid(. ~ SET) + theme(legend.position='bottom')\n";
#print R "dev.off()\n";
close R;

print `Rscript hapmap.pcs.r`;

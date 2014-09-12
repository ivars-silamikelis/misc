#!/usr/bin/perl;
use strict;
use warnings;
#input vcf file in stream
print "chr\tpos\tsnp_baf\tdel_baf\tdel_varw\tins_baf\tins_varw\n";
my $vcfname = shift;
open VCF, "<", $vcfname;
while(<VCF>){
	next if $_=~/^#/;
	chomp($_);
	#print STDERR $_,"\n";
	$_=~s/\|/\\\|/g;
	my @fields=split("\t",$_);
	system("echo $fields[0]\t$fields[1] > $vcfname.tmp.bed");
	my $pos=$fields[1];
	my $ref_base=$fields[3];
	#print $pos,"\n";
	system("samtools view -F 1024 -f ~/Data/Act/references/TCDC-AB0715/TCDC-AB0715.fasta ~/Data/Act/data/aln/a-2885-12.sorted.reheaded.marked.realigned.bam -L $vcfname.bed >$vcfname.tmp.sam");
	system("perl from_bam.pl $pos $ref_base $fields[0] $vcfname.tmp.sam");

}
close VCF;
exit;

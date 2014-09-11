#!/usr/bin/perl;
use strict;
use warnings;
while(<>){
	chomp($_);
	print STDERR $_,"\n";
	$_=~s/\|/\\\|/g;
	system("echo $_ > tmp.bed");
	my @fields=split("\t",$_);
	my $pos=$fields[1];
	#print $pos,"\n";
	system("samtools view -F 1024 -f ~/Data/Act/references/TCDC-AB0715/TCDC-AB0715.fasta ~/Data/Act/data/aln/a-2885-12.sorted.reheaded.marked.realigned.bam -L tmp.bed >tmp.sam");
	system("perl from_bam.pl $pos tmp.sam ");

}

#!/usr/bin/perl;
use strict;
use warnings;

my @header = ( );
my $vcf_name= shift;
my @chroms = (1..22,"X","Y","MT","M");
#print @chroms;
#exit
open my $vcf_handle,"<",$vcf_name;
while (<$vcf_handle>){
	if ($_=~/^#/){
		push (@header, $_);
		next;
	}
	my @fields = split("\t",$_);
	my $chr_field = $fields[0];
	for my $chr (@chroms){
		if ("chr".$chr eq $chr_field){
			print "writing to chrom$chr.vcf\n";
			open my $fh, ">>", "chrom".$chr.".vcf";
			print $fh $_;
			close $fh;
		}
	}
}
close $vcf_handle;

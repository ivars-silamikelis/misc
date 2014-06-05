#!/usr/bin/perl;
use strict;
use warnings;
my $fname=shift;
open FH ,"<", $fname or die "No such file\n";
$name=~s/\.cov$//;
open OUT,">", "$fname.repaired.cov";
while (<FH>){
	chomp($_);
	my @fields=split("\t",$_);
	unless ($fields[2]){
		$fields[2]=$fields[1];	
	}
	my $line=join("\t", @fields);
	print OUT $line,"\n";
}
close FH;
close OUT;



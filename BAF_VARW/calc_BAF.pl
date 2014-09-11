#!/usr/bin/perl;
use strict;
use warnings;

while(<>){
	my $ref_count = 0;
	my @fields = split("\t",$_);
	my $pl = $fields[4];
	$pl =~s/\$//g;
	$pl =~s/\^.//g;
	#print $pl,"\n";
	while ($pl=~/[.,]/g){
		$ref_count+=1;
	}
	unless (length($pl)==0){
		print $pl,"\t",$ref_count,"\t",length($pl),"\t",(length($pl)-$ref_count)/length($pl),"\n";
	}
}



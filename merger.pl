#!/usr/bin/perl;
use strict;
use warnings;
use Data::Dumper;
my %gene_read_count;
my $last_annot="dummy";
my $last_end=0;
my $last_start=0;
my $gene_start;
my %data;
my $is_gene_start=1;
my $current_annot;
my $total_reads_on_region=0;
my @unannot_data;
my %record_unannot;
print "sample\treference\tstart\tend\tlength\tread_count\tannotation\n";
while (<>){
	next if $_=~/^sample/;
	chomp($_);
	my @fields=split("\t",$_);
	my $sample=$fields[0];
	my $ref=$fields[1];
	my $start=$fields[2];
	my $end=$fields[3];
	my $len=$fields[4];
	my $av_cov=$fields[5];
	my $median_cov=$fields[6];
	my $read_count=$fields[7];
	my $total=$fields[8];
	my $annot=$fields[9];
	
	#if ($annot eq $last_annot){
	#	$total_reads_on_region+=$read_count;
	#}
	
	unless (!$annot){
		#print $annot,"\n";
		if ($annot eq $last_annot){
			$total_reads_on_region+=$read_count;
			$last_end=$end;
			
			next;
		} else {
			
			print $last_end,"\t",$last_end-$last_start,"\t",$total_reads_on_region,"\t",$last_annot,"\n";
			$total_reads_on_region=0;
			print $sample,"\t",$ref,"\t",$start,"\t";
			$total_reads_on_region+=$read_count;
			$last_start=$start;
			$last_annot=$annot;
			$last_end=$end;
			next;
		}
	} else {
			push (@unannot_data, [$start,$end,$sample,$ref, $read_count]);
	}
	#	}
	
	#$last_end=$end;
	#}
	
}
print $last_end,"\t",$last_end-$last_start,"\t",$total_reads_on_region,"\n";
$total_reads_on_region=0;
$record_unannot{"start"}=$unannot_data[0]->[0];
$record_unannot{"end"}=$unannot_data[0]->[1];
$record_unannot{"read_count"}=0; #$unannot_data[0]->[4];

for (my $j=0; $j<=$#unannot_data-1; $j++){
	if ($unannot_data[$j+1]->[0]-$unannot_data[$j]->[1]<=300 && $unannot_data[$j+1]->[0]-$unannot_data[$j]->[1]>0){
		$record_unannot{"end"}=$unannot_data[$j+1]->[1];
		$record_unannot{"read_count"}+=$unannot_data[$j+1]->[4];
		next;
	} else { 
		$record_unannot{"end"}=$unannot_data[$j]->[1];
		$record_unannot{"read_count"}+=$unannot_data[$j]->[4];
	}
	print $unannot_data[$j]->[2],"\t",$unannot_data[$j]->[3],"\t", $record_unannot{"start"},"\t",$record_unannot{"end"},"\t",$record_unannot{"end"}-$record_unannot{"start"},"\t",$record_unannot{"read_count"},"\n";
	
	$record_unannot{'start'}=$unannot_data[$j+1]->[0];
	$record_unannot{"read_count"}=0;
}
$record_unannot{"read_count"}+=$unannot_data[$#unannot_data]->[4];
print $unannot_data[$#unannot_data]->[2],"\t",$unannot_data[$#unannot_data]->[3],"\t",$record_unannot{'start'},"\t",$unannot_data[$#unannot_data]->[1],"\t",$unannot_data[$#unannot_data]->[1]-$record_unannot{'start'},"\t",$record_unannot{"read_count"},"\n";


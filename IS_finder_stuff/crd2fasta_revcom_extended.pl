#!/usr/bin/perl;
use strict;
use warnings;
use Bio::Perl;
use Data::Dumper;

my %coords;
my $suffix=shift;
my $alls=shift;
open FH, "<", $alls;
my @crds=glob "$suffix/*.crd";
my @coords_2try;
map{$_=~s/\.$suffix\.crd//} @crds;
map{$_=~s/$suffix\///} @crds;
map {print "$_\n"} @crds;
#exit;
my %data;
my $i=0;
while(<FH>){
	chomp($_);
	my @fields=split("\t",$_);
	my $start_other=$fields[0];
	my $start=$fields[1];
	my $end=$fields[2];
	my $end_other=$fields[3];
	my $name=$fields[5];
	push (@coords_2try,[$start_other,$start,$end,$end_other,$name]);
}
close FH;
#Ja starta koordinātes frekvence mazāka par beigu koordinātes frekvenci, tad izmanto beigu koordināti

for my $record (@coords_2try){
	my $start=$record->[1];
	my $end=$record->[2];
	my $start_other=$record->[0];
	my $end_other=$record->[3];
	#noskaidro kura koordināte biežāka un veido hash ar korektām koordinātēm
	my $correct_idx=correct_coord($record->[1],$record->[2],\@coords_2try);
	#print $correct_idx,"\t",$record->[$correct_idx],"\n";
	my $strand;
	if ($end_other-$start_other>=0){
		$strand=1;
	} else {
		$strand=-1;	
	}
	push(@{$coords{$record->[$correct_idx]}}, [$record->[4], $strand]);
}

my %sequences;
for my $key (sort {$a<=>$b} keys %coords){
	print $key,"\n";
	for my $sample (@crds){
		if (grep{$_->[0] eq $sample and $_->[1]==1} @{$coords{$key}}){		#pir TC - pur GA
			push(@{$sequences{$sample}},"A");								#A	A	TTT
		} elsif (grep{$_->[0] eq $sample and $_->[1]==-1} @{$coords{$key}}) {
			push(@{$sequences{$sample}},"G");								#A	G	TTC
		} else {
			push(@{$sequences{$sample}},"T");								#A	T	CTA
		}
	}
}
print Dumper (\%sequences);
my @aln_sq;
for my $key (keys %sequences){
	my $seq=new_sequence(join("",@{$sequences{$key}}),$key);
	push(@aln_sq,$seq);
}
write_sequence(">$suffix/all.nuc.fasta","fasta",@aln_sq);

sub correct_coord{
	my $start=shift;
	my $end=shift;
	my $rec_ref=shift;
	my $start_count=0;
	my $end_count=0;

	for my $record (@$rec_ref){
		$start_count+=1 if $record->[1]==$start;
		$end_count+=1 if $record->[2]==$end;
	}
	#print "Counts for $start and $end: ", $start_count,"\t",$end_count,"\n";
	if ($start_count>=$end_count){
		return 1;
	} else {
		return 2;
	}
}

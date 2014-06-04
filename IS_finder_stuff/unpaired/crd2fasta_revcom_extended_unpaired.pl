#!/usr/bin/perl;
use strict;
use warnings;
use Bio::Perl;
use Data::Dumper;
my @crds=glob "*.crd";
my %coords;
my @coords_2try;
map{$_=~s/\.CTRI-2.Russia.CP002992.crd//} @crds;
#print @crds;
my %data;
my $i=0;
while(<>){
	chomp($_);
	my @fields=split("\t",$_);
	my $end=$fields[0];
	my $coord=$fields[1];
	my $coord_other=$fields[2];
	my $name=$fields[3];
	#push(@{$coords{$start}},$name);
	#push(@{$coords{}{$end}});
	push (@coords_2try,[$coord, $coord_other,$name, $end]);
}
#Ja starta koordinātes frekvence mazāka par beigu koordinātes frekvenci, tad izmanto beigu koordināti

for my $record (@coords_2try){
	my $coord=$record->[0];
	my $coord_other=$record->[1];
	my $end=$record->[3];
	my $strand;
	if ($end eq "#start_unpaired"){
		if ($coord-$coord_other>=0){
			$strand=1;
		} else {
			$strand=-1;	
		}
	}
	if ($end eq "#end_unpaired"){
		if ($coord_other-$coord>=0){
			$strand=1;
		} else {
			$strand=-1;	
		}
	}
	push(@{$coords{$coord}}, [$record->[2], $strand]);
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
write_sequence(">all.nuc.fasta","fasta",@aln_sq);

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

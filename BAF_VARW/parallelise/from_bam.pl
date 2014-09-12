#!/usr/bin/perl;
use strict;
use warnings;
use Data::Dumper;
use Statistics::Descriptive;
#Calculate VARW
my @dels=( );
my @inserts= ( );
my $snp_pos_on_ref = shift;	#C
my $ref_base = shift;
my $chrom = shift;
my @SNPs;
my @alt_snps;
open BAM,"<",$ARGV[0];
my $depth=0;
while(<BAM>){
	$depth+=1;
	my $cigars_after = 0;
	my @fields=split("\t",$_);
	#print @fields;
	my $read_pos_on_ref = $fields[3];
	my $snp_pos_on_read = $snp_pos_on_ref - $read_pos_on_ref +1; # one based position if cigar has all matches
	#print $read_pos_on_ref,"\n";
	my $cigar = $fields[5];

	my $cigar_position=0;
	my %indels_so_far = ("I"=>0,"D"=>0);	#how many indels there are on read before snp is found
	my $flag_position = 0;
	my $flag_position2 = 0;
	while ($cigar =~/(\d+\w)/g){
		my $cig_frag = $1;
		#print $cig_frag,"\n";
		if ($cigar_position < $snp_pos_on_read){
			
			#print $cigar_position,"\t",$snp_pos_on_read,"\n";
			if ($cig_frag=~/M/){
				$cig_frag =~ s/M//g;
				$cigar_position+=$cig_frag;
     	   } elsif ($cig_frag=~/D/){
				$cig_frag =~ s/D//g;
				$cigar_position+=$cig_frag;
				$indels_so_far{"D"}+=$cig_frag;
			} elsif ($cig_frag=~/I/){
				$cig_frag =~ s/I//g;
				$indels_so_far{"I"}+=$cig_frag;
			} elsif ($cig_frag=~/S/){
				$cig_frag =~ s/S//g;
				$cigar_position+=$cig_frag;
			}
		} else {
			$cigars_after+=1;
		}
		#inserciju cigar
		if ($cigar_position == $snp_pos_on_read && $cigars_after==1 && ($cig_frag=~/I/) && $flag_position2 == 0 && ($cig_frag=~/I/)) {

			#print $cigar_position,"\t",$snp_pos_on_read,"\t",$cigar,"\t",$fields[0],"\t",$fields[1],"\n";
			#print $cig_frag,"\n";
			$indels_so_far{"cigar"}=$cig_frag; #nakamais cigar str
			$flag_position2=1;
			last;
		}
		#deleciju cigar
		if ($flag_position == 0 && ($cigar_position == $snp_pos_on_read)&& $cigars_after==1 && ($cig_frag=~/D/) ){
			#print $cigar_position,"\t",$snp_pos_on_read,"\t",$cigar,"\t",$fields[0],"\t",$fields[1],"\n";
			$indels_so_far{"cigar"}=$cig_frag; #esošais cigar str
			$flag_position=1;
			#print $cig_frag,"\n";

			last;
		}
		if ($indels_so_far{"cigar"}){
			#print $indels_so_far{"cigar"},"\n";
		}
	}

	#TODO  nepraktiski ilgi katru bam failu atsevišķi ģenerēt
	#apkopo informāciju par snp un indeliem
	$snp_pos_on_read+=$indels_so_far{"I"}-$indels_so_far{"D"};
	my $read = $fields[9];
	if ($indels_so_far{"cigar"} && ($indels_so_far{"cigar"}=~/D/)){
		$indels_so_far{"cigar"}=~s/D//;
	    push @dels, $indels_so_far{"cigar"}; 
	} elsif ($indels_so_far{"cigar"} && ($indels_so_far{"cigar"}=~/I/)){
		$indels_so_far{"cigar"}=~s/I//;
		push @inserts, $indels_so_far{"cigar"};
	} else {
		my $SNP = substr $read, ($snp_pos_on_read-1), 1;
		push @SNPs, $SNP;
	}
}
close BAM;
for my $snp (@SNPs){
	if (length $ref_base == 1 && length $snp == 1){
		if ($snp ne $ref_base){
			push @alt_snps,$snp;
		}
	}
}
#map {print $_,"\t",$ref_base,"\n"} @alt_snps;
#print "chr\tpos\tsnp_baf\tdel_baf\tdel_varw\tins_baf\tins_varw\n";
print $chrom,"\t",$snp_pos_on_ref;
if (scalar @alt_snps > 0){
	#my @snp_length = map {length $_} @alt_snps;
	my $stat = Statistics::Descriptive::Full->new();
	#$stat->add_data(\@snp_length);
	print "\t",((scalar @alt_snps)/$depth);
} else {
	print "\t0";
}
if (scalar @dels > 0){
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(\@dels);
	print "\t",((scalar @dels)/$depth),"\t",$stat->variance();
} else {
	print "\t0\t0";
}
if (scalar @inserts > 0) {
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(\@inserts);
	print "\t",(scalar @inserts)/$depth,"\t",$stat->variance(),"\n";
} else {
	print "\t0\t0\n";
}

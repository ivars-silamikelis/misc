#!/usr/bin/perl;
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Perl;
use Data::Dumper;
my $vcfname=shift;

my $gbname="H37Rv.gb";
my %data;
my $gb = Bio::SeqIO->new(-file=>$gbname);
my $gb_seq = $gb->next_seq;
for my $feat_obj ($gb_seq->get_SeqFeatures){
	if ($feat_obj->primary_tag eq "CDS"){
		if ($feat_obj->has_tag("gene")){
			my @gene_names=$feat_obj->get_tag_values("gene");
			my @notes=$feat_obj->get_tag_values("note") if $feat_obj->has_tag("note");
			my @locus_tags=$feat_obj->get_tag_values("locus_tag") if $feat_obj->has_tag("locus_tag");
			#my $locus_tag=join("",@locus_tags);
			my $note=join("",@notes);
			for my $gene_name (@locus_tags){
				$data{$gene_name}{"start"}=$feat_obj->start;
				$data{$gene_name}{"end"}=$feat_obj->end;
				$data{$gene_name}{"note"}=$note;
			#	$data{$gene_name}{"locus_tag"}=$locus_tag;
			}
			#print $feat_obj->start,"\t",$feat_obj->end,"\t",@gene_names,"\n";
		}
	}
}
open VCF, "<", $vcfname;
my @gene_list;
while (<VCF>){
	next if $_=~/^#/;
	my @fields=split("\t",$_);
	my $chrom=$fields[0];
	my $pos=$fields[1];
	my $ref=$fields[3];
	my $alt=$fields[4];
	push(@gene_list, annotate($pos, $ref, $alt, \%data));
}
close VCF;
my $lineage;
#print Dumper(@gene_list);
my $is_euro=1;
if (grep {$_->[1]==321 && $_->[2] eq "Rv0557" && $_->[4] eq "C"} @gene_list){
	$is_euro=0;
	$lineage="not_Euro-American";
} else {
	$lineage="Euro-American";
}
if ($is_euro==0){
	if (grep {$_->[1]==810 && $_->[2] eq "Rv0557" && $_->[4] eq "T"} @gene_list){
		$lineage="EAI";
	} elsif (grep {$_->[1]==221 && $_->[2] eq "Rv0557" && $_->[4] eq "T"} @gene_list){
		$lineage="M_africanum";
	} elsif (grep {$_->[1]==1050 && $_->[2] eq "Rv0557" && $_->[4] eq "A"} @gene_list){
		$lineage="M_caprae";
	} elsif (grep {$_->[1]==1066 && $_->[2] eq "Rv0557" && $_->[4] eq "T"} @gene_list){
		$lineage="M_canettii";
	} else {
		if (grep {$_->[1]==472 && $_->[2] eq "Rv0129c" && $_->[4] eq "A"} @gene_list){
			$lineage="Delhi/CAS";
		} elsif (grep {$_->[1]==735 && $_->[2] eq "Rv1009" && $_->[4] eq "T"} @gene_list){
			$lineage="M. microti, M. pinipedii";
		} elsif (grep {$_->[1]==1070 && $_->[2] eq "Rv1009" && $_->[4] eq "T"} @gene_list){
			$lineage="M. bovis";
		} elsif (grep {$_->[1]==284 && $_->[2] eq "Rv1811" && $_->[4] eq "T"} @gene_list){
			$lineage="M. africanum 1b";
		} elsif (grep {$_->[1]==320 && $_->[2] eq "Rv1811" && $_->[4] eq "T"} @gene_list){
			$lineage="M. africanum 2";
		} elsif (grep {$_->[1]==191 && $_->[2] eq "Rv2629" && $_->[4] eq "C"} @gene_list){
			$lineage="Beijing";
		} else {
			$lineage="not_defined";
		} 
	}
} elsif ($is_euro==1) {
	if (grep {$_->[1]==455 && $_->[2] eq "Rv0557" && $_->[4] eq "C"} @gene_list){
		$lineage="Haarlem";
	} elsif (grep {$_->[1]==457 && $_->[2] eq "Rv0557" && $_->[4] eq "G"} @gene_list){
		$lineage="S-type";
	} elsif (grep {$_->[1]==532 && $_->[2] eq "Rv0557" && $_->[4] eq "G"} @gene_list){
		$lineage="Cameroon";
	} else {
		if (grep {$_->[1]==309 && $_->[2] eq "Rv0129c" && $_->[4] eq "A"} @gene_list){
			$lineage="LAM";
		} elsif (grep {$_->[1]==1075 && $_->[2] eq "Rv1009" && $_->[4] eq "A"} @gene_list){
			$lineage="TUR";
		} elsif (grep {$_->[1]==12 && $_->[2] eq "Rv1811" && $_->[4] eq "A"} @gene_list){
			$lineage="Ural";
		} elsif (grep {$_->[1]==965 && $_->[2] eq "Rv2629" && $_->[4] eq "T"} @gene_list){
			$lineage="Ghana";
		} else {
			$lineage="not defined";
		}
	}

}
print $lineage,"\n";
#for my $anno (@gene_list){
#	for my $record (@$anno){
#		if ($record->[1]==321 && $record->[2] eq "Rv0557" && $record->[4] eq "C"){
#			$is_euro=0;
#			
#		}
#	}
#}
#if ($is_euro==0){
#	for my $anno (@gene_list){
#		for my $record (@$anno){
#			if ($record->[1]==321 && $record->[2] eq "Rv0557" && $record->[4] eq "C"){
#				$is_euro=0;
#				
#			}
#		}
#	}
#}



sub annotate {
	my $pos=shift;
	my $ref=shift;
	my $alt=shift;
	my $db=shift;
	#my $resdb=shift;
	my @gene_list_data;
	for my $gene (keys %$db){
		#print $gene,"\t", $db->{$gene}{"start"},"\n";
		if ($pos>=$db->{$gene}{"start"} && $pos<=$db->{$gene}{"end"}){
			my $pos_in_gene=$pos-$db->{$gene}{"start"}+1;
			#my $debug_ref= substr $mutated_dna, $pos_in_gene-1, 1, ",";#$alt;
			#print $pos,"\t",$pos_in_gene,"\t",$gene,"\t",$ref,"\t",$alt,"\n";#,"\t",translate_as_string($dna),"\n";	#$db->{$gene}{"note"},"\n";
			push (@gene_list_data, [$pos,$pos_in_gene,$gene,$ref,$alt]);
		}
	}
	return @gene_list_data;
}

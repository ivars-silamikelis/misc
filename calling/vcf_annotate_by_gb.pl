#!/usr/bin/perl;
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Perl;
my $res_db="resistance_db.txt";
my $false_regions="rep_ppe_reg.txt";
open REG, "<", $false_regions;
my @regions=<REG>;
close REG;
my %res_data;
my $gbname="H37Rv.gb";
open RES, "<", $res_db;
while (<RES>){
	next if $_=~/^ID/;
	my @fields=split("\t",$_);
	my $aachange=$fields[14];
	my $pos=$fields[13];
	my $gene=$fields[1];
	my $drug=$fields[3];
	my $gene_id=$fields[2];
	next if !$pos or $aachange!~/\// or $aachange=~/ERROR/i or $pos<0;
	#print aa_abrv($aachange,$pos),"\n";
	push @{$res_data{$gene}{"changes"}},aa_abrv($aachange,$pos);
	$res_data{$gene}{aa_abrv($aachange,$pos)}=$drug;
}
close RES;
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
	#exit;
#open VCF, "<", $vcfname;
while (<>){
	next if $_=~/^#/;
	chomp($_);
	my $resistance;#=[".",".","."];
	my @fields=split("\t",$_);
	my $chrom=$fields[0];
	my $pos=$fields[1];
	my $ref=$fields[3];
	my $alt=$fields[4];
	my @info=split(";",$fields[7]);
	my $eff= $info[$#info];
	next if is_in_region($pos,\@regions)==1;
	my @gene_list=annotate($fields[1],$fields[3],$fields[4],\%data);
	for my $gene (@gene_list){
		unless (!$res_data{$gene}){
			my @changes=grep {$eff=~/$_/} @{$res_data{$gene}{"changes"}};
			unless (!@changes){
				$resistance=[$gene, join("",@changes),$res_data{$gene}{$changes[0]}];
			}
		}
	}
	 if ($resistance){
		print "$pos\t$ref\t$alt\t@gene_list\t",$resistance->[0],"\t",$resistance->[1],"\t",$resistance->[2],"\t$eff\n";
	}
}
#close VCF;

#############################subroutines##########################################################################
sub annotate {
	my $pos=shift;
	my $ref=shift;
	my $alt=shift;
	my $db=shift;
	#my $resdb=shift;
	my @gene_list;
	for my $gene (keys %$db){
		#print $gene,"\t", $db->{$gene}{"start"},"\n";
		if ($pos>=$db->{$gene}{"start"} && $pos<=$db->{$gene}{"end"}){
			#my $pos_in_gene=$pos-$db->{$gene}{"start"}+1;
			#my $debug_ref= substr $mutated_dna, $pos_in_gene-1, 1, ",";#$alt;
			#print $pos,"\t",$gene,"\t",$ref,"\t",$debug_ref,"\t",$mutated_dna,"\n";#,"\t",translate_as_string($dna),"\n";	#$db->{$gene}{"note"},"\n";
			push (@gene_list, $gene);
		}
	}
	return @gene_list;
}

sub aa_abrv {
	my $aastr=shift;
	my $pos=shift;
	my %aminos=('Ala' => 'A',
		'Arg' => 'R',
		'Asn' => 'N',
		'Asp' => 'D',
		'Cys' => 'C',
		'Gln' => 'Q',
		'Glu' => 'E',
		'Gly' => 'G',
		'His' => 'H',
		'Ile' => 'I',
		'Leu' => 'L',
		'Lys' => 'K',
		'Met' => 'M',
		'Phe' => 'F',
		'Pro' => 'P',
		'Pyl' => 'O',
		'Sec' => 'U',
		'Ser' => 'S',
		'Thr' => 'T',
		'Trp' => 'W',
		'Tyr' => 'Y',
		'Val' => 'V',
		'Asx' => 'B',
		'Glx' => 'Z',
		'Xle' => 'J',
		'Xaa' => 'X',
		'STOP'=> '*'
	);
	my @aas=split("/",$aastr);
	my $abrv=join("",$aminos{$aas[0]},$pos,$aminos{$aas[1]});
	return $abrv;
}
sub is_in_region{
	my $pos=shift;
	my $regs=shift;
	for my $line (@$regs){
		my @fields=split("\t",$line);
		if ($pos>=$fields[0] && $pos<=$fields[1]){
			return 1;
		}
	}
	return 0;
}

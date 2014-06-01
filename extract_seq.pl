#!/usr/bin/perl;
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Perl;
my @gb_names=qw/NC_000962.gbk NC_009525.gbk NC_012943.gbk NC_016934.gbk NC_017524.gbk NC_002755.gbk NC_009565.gbk NC_017026.gbk NC_017528.gbk/;
open my $fh, "genes";
my @genes=<$fh>;
my %gene_seq;
my @gene;
close $fh;
#Katram failam izvelk gena sekvenci no genu saraksta un saglaba hashaa %gene_seq
for my $name (@gb_names){
  my $io_obj=Bio::SeqIO->new(-file=>"/home/ivars/Data/myco_tuber/references/MTuber_ref_diff_strains/$name");
  my $sq_obj=$io_obj->next_seq;
  for my $feat_obj ($sq_obj->get_SeqFeatures){
    if ($feat_obj->primary_tag eq "CDS"){
      if ($feat_obj->has_tag("gene")){
        #print $feat_obj->get_tag_values("gene"),"\n";
        @gene=$feat_obj->get_tag_values("gene");
        $gene_seq{$name}{$gene[0]}=$feat_obj->spliced_seq->seq if (grep{$_=~/$gene[0]/}@genes);
      }
    }
  }
}
#print "CTRI-2\tSN2-1\n";
for my $gene (keys %{$gene_seq{"NC_017524.gbk"}}){
  my $length=length($gene_seq{"NC_017524.gbk"}{$gene});  #blast query length
  my $rez=&blaster($gene_seq{"NC_017524.gbk"}{$gene},$gene); #blast rezult for SN2-1 contig
  while (my $hit=$rez->next_hit()){#parsin result
    while (my $hsp=$hit->next_hsp()){
      #print $hsp->query_string,"\n";
      #print $length,"\t",length($hsp->hit_string),"\n" if $length>length($hsp->hit_string);
      if (($length-20)<=length($hsp->hit_string)){
        #print $length,"\t",length($hsp->hit_string),"\n";
        $gene_seq{"SN2-1"}{$gene}=$hsp->hit_string;
      } else { 
        for my $name (keys %gene_seq){
          delete $gene_seq{$name}{$gene};
        }
      }
    }
  }
}
#print Dumper(\%gene_seq);
my @seq_objs;
for my $name (keys %gene_seq){
  my @frags;
  print $name,"\t",scalar keys %{$gene_seq{$name}},"\n";
  for my $gene (sort keys %{$gene_seq{$name}}){
    push @frags,$gene_seq{$name}{$gene};
  }
  my $sekv=join("",@frags);
  my $sq_obj=new_sequence($sekv,$name);
  push @seq_objs,$sq_obj;
}
write_sequence(">MTB_filogen3.fasta",'fasta',@seq_objs);





sub blaster {
  my $fac=Bio::Tools::Run::StandAloneBlastPlus->new(-db_name=>"/home/ivars/Blast/databases/mira/SN2-1/SN2-1_contigs");
  my $q_sekv=shift;
  my $q_name=shift;
  my $query=new_sequence($q_sekv,$q_name);
  my $result=$fac->blastn(-query=>$query);
  $fac->cleanup;
  return $result;
}

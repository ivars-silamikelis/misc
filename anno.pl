#!/usr/bin/perl;
use strict;
use warnings;
#use Bio::DB::SeqFeature::Store;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::SeqFeature::Generic;
use Bio::Perl;
#Anotē pēc viena genbank faila gēniem citu sekvenci un izveido tai savu genbank failu
#fasta fails kursh tiek anotets (tadai pashai jābūt blast database)
my $in_fasta=shift;#"/home/ivars/Data/Mtb/references/MTuber_ref_diff_strains/nederigas/AE000516.fasta";
#genu saraksts, kuri tiek pēc tam meklēti blast database
my $gene_gb="/home/ivars/Data/Act/references/MDR-TJ.gb";
#anotētā genbank faila nosaukums
my $out_gb=shift; #"/home/ivars/Data/Mtb/references/MTuber_ref_diff_strains/nederigas/AE000516_annotated.gb";
my $gene_name;
my $locus_tag;
my $seq_obj=Bio::SeqIO->new(-file=>$in_fasta,-format=>"fasta");
my $gb_obj= Bio::SeqIO->new(-file => ">$out_gb",
                           -format => 'GenBank');
my $gb_seq=$seq_obj->next_seq();
my $str = Bio::SeqIO->new(-file=>"$gene_gb", -format=>"GenBank");
my $seq=$str->next_seq();
#my @feats=$seq->get_SeqFeatures();
for my $feat_obj ($seq->get_SeqFeatures()){
  if ($feat_obj->primary_tag eq "CDS"){
    if ($feat_obj->has_tag('locus_tag')){
      #print $feat_obj->spliced_seq->seq,"\n";
      #last;
      #CDS features iemet blastāi un iegūst to koordinātes
      my @coords=&extr_cds_coords($feat_obj->spliced_seq->seq,$feat_obj->get_tag_values('locus_tag'));
      unless($coords[0]==0){
        for my $kref (@coords){
          my $feat=&create_tag($kref,$feat_obj->get_tag_values('locus_tag'));
          $gb_seq->add_SeqFeature($feat);
          print "Adding ", $feat_obj->get_tag_values('locus_tag'),"\n";
        }
      }
    }
  }
}
$gb_obj->write_seq($gb_seq);



sub create_tag {
my $kref=shift;
my $gene=shift;
my $feat=new Bio::SeqFeature::Generic(-start=>$kref->[0]->[0],
                                      -end=>$kref->[0]->[1],
                                      -strand=>$kref->[1],
                                      -primary_tag=>"CDS",
                                      -tag=>{locus_tag=>$gene}
                                     );
return $feat;
}


sub extr_cds_coords {
  my $fac=Bio::Tools::Run::StandAloneBlastPlus->new(-db_data=>$in_fasta,create=>1);
  my $seq=shift;
  my $gene=shift;
  my @ret_crd;
  my $query=new_sequence($seq, $gene);
  my $result=$fac->blastn(-query=>$query);
    my $hit = $result->next_hit;
    #print "\n",$hit,"\n";
    if($hit){
    my $length=0;
    # insert code here for hit processing
    while( my $hsp = $hit->next_hsp ) {
        if ($length-5<=length $hsp->query_string()){  #ja ieprieksheja alig
            $length=length $hsp->query_string();
            push(@ret_crd,[[$hsp->range('hit')],$hit->strand('hit')]);
        }
    }
    $fac->cleanup;
    return @ret_crd;
  } else {
    $fac->cleanup;
    return 0;
  }
}


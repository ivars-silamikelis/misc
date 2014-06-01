#!/usr/bin/perl;
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Perl;
$filename=shift;
my @seq=read_all_sequences($filename,"fasta")

sub blaster {
  my $query=shift;

  my $fac=Bio::Tools::Run::StandAloneBlastPlus->new(-db_name=>"/home/ivarz/Blast/databases/CTRI-2/CTRI-2");
  my $result=$fac->blastn(-query=>$query);
  return $result;
}

sub parse_result{
  my $rez =shift;
  #print $rez;
  while (my $hit=$rez->next_hit()){
    my $cov=0;
    while (my $hsp=$hit->next_hsp()){
      my @range=$hsp->range('hit');
      my $q_cov=abs(($range[1]-$range[0])/$window);
      if ($cov<= $q_cov) {  #rāda tikai labāko query coverage
        #print $hit->start
        $cov=$q_cov;

       
      }
    }
  }
}

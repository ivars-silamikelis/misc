#!/usr/bin/perl;
use strict;
use warnings;
use Bio::DB::Sam;
use Bio::Perl;
use Statistics::Descriptive;
use Bio::Tools::Run::StandAloneBlastPlus;
my $sam = Bio::DB::Sam->new(-bam  =>"/home/ivars/Data/myco_tuber/reads/SN2-2_CTRI-2.sorted.bam",
                             -fasta=>"/home/ivars/Data/myco_tuber/references/CTRI-2/CTRI-2.fasta"
                             );

my $header=$sam->header;
my @targets    = $sam->seq_ids;
my ($coverage) = $sam->features(-type =>"coverage",-seq_id => 'gi|344217818|gb|CP002992.1|');
my $idx=0;
my @cover=$coverage->coverage;

#my $stat1=Statistics::Descriptive::Full->new();
#my $stat2=Statistics::Descriptive::Full->new();
#my $ix=1806104+9;
#$stat1->add_data(@cover[$ix-20..$ix-10]);
#$stat2->add_data(@cover[$ix-9..$ix]);
#print $stat_s->median()-$stat_e->median(),"\t",$stat_s->sample_range()<=25 && $stat_e->sample_range()<=25);
for my $point (@cover){
  if ($idx=>40){
  my $stat_s=Statistics::Descriptive::Full->new();
   my $stat_e=Statistics::Descriptive::Full->new();
    $stat_s->add_data(@cover[$idx-40..$idx-20]);
    $stat_e->add_data(@cover[$idx-19..$idx]);
    if(abs($stat_s->mean()-$stat_e->mean())>35 && ($stat_s->sample_range()<=20 && $stat_e->sample_range()<=20)){
      print $idx-319,"-",$idx+281,"\n";
      
    }

  }
$idx++;
#if ($idx==$#cover){
#  last;
#  }
}


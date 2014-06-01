#!/usr/bin/perl;
use strict;
use warnings;
use Statistics::Descriptive;
my $stat = Statistics::Descriptive::Full->new();
my $genome_length=4400000;
my @seq_length;
my $num_reads=0;
while(<>){
  my $seq=<>;
  <>;
  <>;
  push(@seq_length,length($seq));
  $num_reads+=1;
}
$stat->add_data(@seq_length);
my $sum = $stat->sum();
my $av_cov=$sum/$genome_length;
print "Average Coverage: $av_cov X\n";

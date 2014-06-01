#!/usr/bin/perl;
use strict;
use warnings;
use Bio::Perl;
#use Statistics::Descriptive;
my $outname=shift;
while(<>){
  my $id=$_;
  chomp($id);
  $id=~s/^@//;
  $id=~s/ /_/g;
  my $sequence=<>;
  chomp($sequence);
  <>;
  <>;
  #my $stat=Statistics::Descriptive::Full->new();
  #my @fields=split("\t",$_);
  #my @qual_chr=split("",$fields[10]);
  #my @qual = map{ord($_)-33} @qual_chr;
  #$stat->add_data(@qual);

  my $seq=new_sequence($sequence, $id);
  write_sequence(">>$outname.fasta", "fasta", $seq);
  #  print "\@$fields[0]\n$fields[9]\n+\n$fields[10]\n";# if $stat->mean>=$mean;
  #   print $fields[10],"\n" if $stat->mean>=27;
  
  }

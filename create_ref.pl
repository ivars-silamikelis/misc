#!/usr/bin/perl;
use strict;
use warnings;
use Bio::Perl;
my $name=shift;
my @seq_objs=read_all_sequences($name,"fasta");
my @large_seqs;
my @bpoints;
my $i=0;
for my $seq_obj (@seq_objs){
  if (length($seq_obj->seq())>=5000){
    $i+=length($seq_obj->seq());
    my $sequence=$seq_obj->seq();
    $sequence=~tr/wsmkrybdhvn/acagaccaaaa/;
    print $i,"\n";
    push(@bpoints,$i);
    push(@large_seqs,$sequence);
  }
}
my $large_seq=join("",@large_seqs);
write_sequence(">act1.fasta","fasta",new_sequence($large_seq,"act1_denovo"));

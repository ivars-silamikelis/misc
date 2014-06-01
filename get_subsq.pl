#!/usr/bin/perl;
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Perl;
my $fasta=shift;
my $start=shift;
my $end=shift;
my $in  = Bio::SeqIO->new(-file => $fasta ,
                           -format => 'Fasta');
my $seq= $in->next_seq;
my $subseq=$seq->trunc($start,$end);
print ">",$subseq->display_id(),"\n";
print $subseq->seq(),"\n";
exit;


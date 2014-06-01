#!/usr/bin/perl;
use strict;
use warnings;
use Bio::Perl;
my $in_name=shift;
my $seq=read_sequence($in_name,'genbank');
 $in_name=~s/\.\w*?$/\.fasta/g;
write_sequence(">$in_name",'fasta',$seq);

#!/usr/bin/perl;
use strict;
use warnings;
use Bio::DB::Sam;
my $ref = "/home/ivars/Data/human/references/Ampliseq_hg19/hg19.fasta";
my $sam = Bio::DB::Sam->new(-bam  =>"/home/ivars/Data/human/AmpliSeq/aligned/SN2-3/IonXpress_001_R_2013_05_29_19_11_56_user_SN2-3-Inetas-Mody-1_Inetas-Mody-1.5.bam",
                             -fasta=>$ref,
			     );
$sam->header;
my $i=1;
my @targets    = $sam->seq_ids;
print @targets;
my ($coverage) = $sam->features(-type =>"coverage",-seq_id=>"chr1");#seq id

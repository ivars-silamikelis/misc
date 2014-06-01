#!/usr/bin/perl;
use strict;
use warnings;
use Bio::DB::Sam;
use Bio::Perl;
use Bio::Tools::Run::StandAloneBlastPlus;
my $sam = Bio::DB::Sam->new(-bam  =>"/home/ivars/Data/myco_tuber/reads/SN2-1_CTRI-2.sorted.bam",
                             -fasta=>"/home/ivars/Data/myco_tuber/reads/references/CTRI-2/CTRI-2.fasta"
			     );
 my $header=$sam->header;
 my @targets    = $sam->seq_ids;
 my @alignments = $sam->get_features_by_location(-seq_id => 'gi|344217818|gb|CP002992.1|',
                                                 -start  => 2635800,
                                                 -end    => 2636000);

#my $out=Bio::DB::Bam->open("CTRI-2_reads.bam","w") or die "could not create bam file";
#$out->header_write($header);
#output fastq formātā
for (@alignments){
  if (clasify($_->query->dna)==1){
    #tad ieraksta bam failā
    print "@",$_->query->name,"\n";
    print $_->query->dna,"\n+\n";
    print map {chr($_+33)} $_->qscore;
    print "\n";
    #$out->write1($_) or warn "could not write read";
  }
}
sub clasify {
  my $sekv=shift;
  my $fac=Bio::Tools::Run::StandAloneBlastPlus->new(-db_name=>"/home/ivars/Blast/databases/H37Rv_CTRI-2/H37Rv_CTRI-2");
  my $result=$fac->blastn(-query=>new_sequence($sekv));
  my %e_val=( );
  while (my $hit=$result->next_hit){
    if ($hit->name=~/344217818/){
     # print "CTRI-2\n";
     # print $hit->significance,"\n";
      $e_val{"CTRI-2"}=$hit->significance;
    }elsif($hit->name=~/444893469/){
     # print "H37Rv\n";
     # print $hit->significance,"\n";
      $e_val{"H37Rv"}=$hit->significance;
    }
  }
  $fac->cleanup;
  if ($e_val{"CTRI-2"} && $e_val{"H37Rv"}){
    if ($e_val{"CTRI-2"}<=$e_val{"H37Rv"}){   #ja reads ir tuvaks CTRI-2 tad atgriež 1
      return 1;
    } else {
    return 1;
    }
  } elsif ($e_val{"CTRI-2"}){
      return 1;
  } elsif ($e_val{"H37Rv"}){
      return 1;
  } else {
      return 0;
  }
#$fac->cleanup;
}


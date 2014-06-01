#!/usr/bin/perl;
use strict;
use warnings;
use Bio::Tools::Run::RemoteBlast;
use Bio::Perl;
my $filename=shift;
my @seq_objs=read_all_sequences($filename);
for my $quer (@seq_objs){
  print $quer->display_name(),"\t";
  &blast_seq_and_get_title($quer);
}
sub blast_seq_and_get_title {
  my $score=0;
  my $query=shift;
  my $prog='blastn';
  my $db = 'nr';
  my $e_val= '1e-10';

  my @params=('-prog' => $prog,
              '-data' => $db,
              '-expect' => $e_val,
              '-readmethod' => 'SearchIO');
  my $factory = Bio::Tools::Run::RemoteBlast->new(@params);
  my $r = $factory->submit_blast($query);
  while (my @rids = $factory->each_rid){
    foreach my $rid ( @rids ) {
      my $rc = $factory->retrieve_blast($rid);
      if( !ref($rc) ) {
        if( $rc < 0 ) {
          $factory->remove_rid($rid);
        }
        print STDERR ".";
        sleep 5;
      } else {
        print STDERR "\n";
        my $result=$rc->next_result();
        $factory->remove_rid($rid);
        while (my $hit=$result->next_hit){
          if ($score<=$hit->score){
            $score=$hit->score;
            print $hit->description,"\t";
            print "query length\t", length($query->seq),"\t";
            print "hit length\t", $hit->length,"\t", "score\t", $hit->score,"\n";
          }
          
        }
      }
    }
  }
}




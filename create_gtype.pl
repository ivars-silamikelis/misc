#!/usr/bin/perl;
use strict;
use warnings;
use Bio::Perl;
use Data::Dumper;
#problēma ar MNP, kuru no alternatīvajām alēlēm izmantot kā reprezentatīvo
my $head;
my %seq;
my @head_fields;
my $t3=0;
my $r=0;
my $outname=shift;
while(<>){
  $_=~s/\n//;
  if ($_=~/^#CHR/){
    $head=$_;
    $head=~s/^#||\n//g;
    @head_fields=split("\t",$head);
    #map {print "$_\n"} @head_fields;
  }
  unless ($_=~/^#/){
    my @fields=split("\t", $_);
    my $ref_al=$fields[3];
    my @alt_al=split(",",$fields[4]);
    my $i=0;
    for my $variant_idx (9..$#fields){
      my @gtype_fields=split(":",$fields[$variant_idx]);
      ++$i if $gtype_fields[0]=~/(\.\/\.)/;
    }
    next if $i >0; #izlaizh vietas kur vismaz vienam nav coverage
    for my $variant_idx (9..$#fields){
      my @gtype_fields=split(":",$fields[$variant_idx]);    
      
      my @AD=split(",",$gtype_fields[1]);
      if($gtype_fields[0]=~/0\/0/){
        push(@{$seq{$head_fields[$variant_idx]}},$fields[3]);
        print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[3],"\tref","\n";
      } 
      if ($gtype_fields[0]=~/1\/1/){
        print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[4],"\talt","\n";
        push(@{$seq{$head_fields[$variant_idx]}},$alt_al[0]);
      } 
      if ($gtype_fields[0]=~/2\/2/){
        print "MNP###### ",$head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[4],"\talt","\n";
        push(@{$seq{$head_fields[$variant_idx]}},$alt_al[1]);
      }
      if ($gtype_fields[0]=~/3\/3/){
        print "MNP###### ",$head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[4],"\talt","\n";
        push(@{$seq{$head_fields[$variant_idx]}},$alt_al[2]);
      }
      if ($gtype_fields[0]=~/0\/1/){
        
        #ko darīt ja ir heterozigotisks, kuru alēli likt?
        #1. degenerate bases
        #2. ņem mazhoro ja tā parsniedz kādu threshold
        if($AD[1]/( $AD[0]+$AD[1])>=0.5){
          push(@{$seq{$head_fields[$variant_idx]}},$alt_al[0]);
          print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[4],"alt","\n";
        } else {
          push(@{$seq{$head_fields[$variant_idx]}},$fields[3]);
          print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[3],"ref","\n";
        }
      }
      if ($gtype_fields[0]=~/0\/2/){

        if($AD[2]/( $AD[0]+$AD[2])>=0.5){
          push(@{$seq{$head_fields[$variant_idx]}},$alt_al[1]);
          print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[4],"alt","\n";
        } else {
          push(@{$seq{$head_fields[$variant_idx]}},$fields[3]);
          print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[3],"ref","\n";
        }       
      }
      if ($gtype_fields[0]=~/0\/3/){

        if($AD[3]/( $AD[0]+$AD[3])>=0.5){
          push(@{$seq{$head_fields[$variant_idx]}},$alt_al[2]);
          print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[4],"alt","\n";
        } else {
          push(@{$seq{$head_fields[$variant_idx]}},$fields[3]);
          print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[3],"ref","\n";
        }       
      }
      if ($gtype_fields[0]=~/1\/2/){

        if($AD[2]/( $AD[1]+$AD[2])>=0.5){
          push(@{$seq{$head_fields[$variant_idx]}},$alt_al[1]);
          print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[4],"alt2","\n";
        } else {
          push(@{$seq{$head_fields[$variant_idx]}},$alt_al[0]);
          print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[3],"alt1","\n";
        }       
      }
      if ($gtype_fields[0]=~/1\/3/){

        if($AD[3]/( $AD[1]+$AD[3])>=0.5){
          push(@{$seq{$head_fields[$variant_idx]}},$alt_al[2]);
          print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[4],"alt3","\n";
        } else {
          push(@{$seq{$head_fields[$variant_idx]}},$alt_al[0]);
          print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[3],"alt1","\n";
        }       
      }
      if ($gtype_fields[0]=~/2\/3/){

        if($AD[3]/( $AD[2]+$AD[3])>=0.5){
          push(@{$seq{$head_fields[$variant_idx]}},$alt_al[2]);
          print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[4],"alt3","\n";
        } else {
          push(@{$seq{$head_fields[$variant_idx]}},$alt_al[1]);
          print $head_fields[$variant_idx],"\t",$gtype_fields[0],"\t",$fields[3],"alt2","\n";
        }
      }       
    }
  } 
}
my %rez;
my @seq_obj;
for my $var (sort keys %seq){
  $rez{$var}=join("",@{$seq{$var}});
  my $seq=new_sequence($rez{$var},$var);
  push (@seq_obj,$seq);
}
write_sequence(">$outname", 'fasta', @seq_obj);
print Dumper(\%rez);

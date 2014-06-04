#!/usr/bin/perl;
use strict;
use warnings;
use File::Basename;

print $ARGV[0],"\n";
my @blats=glob "crd2/*";
map {$_=~s/crd2\///} @blats;
#map {print $_,"\n"} @blats;
#exit;
for my $blat_name (@blats){	#jāmaina blat_name
	my %start_coords;
	my %end_coords;
	#print $blat_name,"\n";
	for my $end (qw/starts ends/){
	  print STDERR  "processing $ARGV[0]/blat/$ARGV[0].$blat_name"."_clean.$end.blat\n";
	  open FH, "<", "$ARGV[0]/blat/$ARGV[0].$blat_name"."_clean.$end.blat"
	    or die "Cannot open file\n";
	  my $last_read_name="dummy";
	  my %last_start_coords;
	  my %last_end_coords;
	  while (<FH>){
		
		#jāizdzēš atkārtojumi
	    chomp($_);
		my @fields=split("\t",$_);
	    if ($_=~/^\d/) {
		  
		  if ($last_read_name eq $fields[9]){
	#		print $last_read_name, $fields[9],"\n";
	#		%start_coords=%last_start_coords;
	#		%end_coords=%last_end_coords;
		    next; 
		  }
	#	  %last_start_coords=%start_coords;
	#	  %last_end_coords=%end_coords;
	      
	      if ($fields[8] eq "+" && $end eq "starts") {
			#saskaita IS6110 5' un 3' koordinātes un norāda tālāko nolasījuma pozīciju, kas iziet no IS6110 gala
	        $start_coords{$fields[16]}{"count"}+=1;
			$start_coords{$fields[16]}{"other"}=$fields[15] if $start_coords{$fields[16]}{"other"}>$fields[15] or !$start_coords{$fields[16]}{"other"};
	      } elsif ($fields[8] eq "-" && $end eq "starts"){
	        $start_coords{$fields[15]}{"count"}+=1;
			$start_coords{$fields[15]}{"other"}=$fields[16] if $start_coords{$fields[15]}{"other"}<$fields[16] or !$start_coords{$fields[15]}{"other"};
	      } elsif ($fields[8] eq "+" && $end eq "ends"){       
	        $end_coords{$fields[15]}{"count"}+=1;
			$end_coords{$fields[15]}{"other"}=$fields[16] if $end_coords{$fields[15]}{"other"}<$fields[16] or !$end_coords{$fields[15]}{"other"};
	      } elsif ($fields[8] eq "-" && $end eq "ends"){
	        $end_coords{$fields[16]}{"count"}+=1;
			$end_coords{$fields[16]}{"other"}=$fields[15] if $end_coords{$fields[16]}{"other"}>$fields[15] or !$end_coords{$fields[16]}{"other"};
	      } else {
	        warn "Not a standart case\n"; 
	      }
		  $last_read_name=$fields[9];
	    }
		
	  }
	  close FH;
	}
	open my $fh, ">>", "crd2/$blat_name/$ARGV[0].$blat_name.crd";
	#print "starts\tends\taccession\n";f
	#map {print "$_\n" if $_ > 1} sort {$a <=> $b} keys %start_coords;
	#print "ends\n";
	#map {print "$_\n" if $_ > 1} sort {$a <=> $b} keys %end_coords;
    #noņem potenciālās pielīdzināšanas kļūdas
	clean_hash(\%start_coords);
	clean_hash(\%end_coords);
	
	#jāmēģina atrast coordināšu pāri ņemot vērā genoma pārkārtojumus
	
	#pārējās koordinātes parāda tāpat
#	for my $end_key (sort {$a <=> $b} keys %end_coords){
#	  for my $start_key (sort {$a <=> $b } keys %start_coords){  
#	    my $value = abs ($end_key - $start_key);

#		if ($start_coords{$start_key}{"count"}>=5 && $end_coords{$end_key}{"count"}>=5) {
##		  if ($start_coords{$start_key+1}{"count"}+$end_coords{$end_key+1}{"count"}>$start_coords{$start_key}{"count"}+$end_coords{$end_key}{"count"}){			
#		    
##		  }
#	  	  if ($value<10){	#jaizmēģina bez value, ņemot vērā genoma pārkārtojumus
#	 	     print $fh "$start_coords{$start_key}{other}\t$start_key\t$end_key\t$end_coords{$end_key}{other}\t$value\t$ARGV[0]\t$start_coords{$start_key}{count}\t$end_coords{$end_key}{count}\t$blat_name\n";
#	 	     delete $start_coords{$start_key};
#	 	     delete $end_coords{$end_key};    
#	 	   }
#	 	}
#	  }
#	}
	#parāda to kas palicis
	map {print $fh "#start_unpaired\t$_\t$start_coords{$_}{other}\t$ARGV[0]\t$start_coords{$_}{count}\t$blat_name\n" if $start_coords{$_}{"count"}>=5} sort {$a <=> $b} keys %start_coords;	#parāda nemapētos nolasījumus, ja tiem ir vairāk par noteiktu slieksni readu skaits
	map {print $fh "#end_unpaired\t$_\t$end_coords{$_}{other}\t$ARGV[0]\t$end_coords{$_}{count}\t$blat_name\n" if $end_coords{$_}{"count"}>=5} sort {$a <=> $b} keys %end_coords;
	close $fh;
}
sub clean_hash {
  #Nonjem tas kordinates no key-value para kuras ir mazāk par nākamo koordināti un kuras ir tuvāk kādai koordinātei par 10 (Katram paraugam atsevišķi jādarbojas)
  my $hashref = shift;
  my @keys = sort {$a <=> $b} keys %{$hashref};
  my @bad_idx;
  my $thresh=5;
  #if (abs($keys[0] - $keys[1]) < 10 && $hashref->{$keys[0]} < 10){
  if (abs($keys[0] - $keys[1]) < $thresh && $hashref->{$keys[0]}{"count"} < $hashref->{$keys[1]}{"count"}){
    push(@bad_idx, 0);
  }
  for my $idx (1..$#keys-1){
    #print $hashref->{$keys[$idx]},"\n";
    #if ((abs($keys[$idx-1] - $keys[$idx]) < 10 || abs($keys[$idx+1] - $keys[$idx]) < 10) && $hashref->{$keys[$idx]} < 10){
	if ( (abs($keys[$idx-1] - $keys[$idx]) < $thresh && $hashref->{$keys[$idx]}{"count"} < $hashref->{$keys[$idx-1]}{"count"})|| ((abs($keys[$idx+1] - $keys[$idx]) < $thresh) && $hashref->{$keys[$idx]}{"count"} < $hashref->{$keys[$idx+1]}{"count"})){
      push(@bad_idx, $idx);
    }
  }
  if (abs($keys[$#keys-1] - $keys[$#keys]) < $thresh && $hashref->{$keys[$#keys]}{"count"} <$hashref->{$keys[$#keys-1]}{"count"} ){
    push(@bad_idx, $#keys);
  }

  for my $bdix (@bad_idx){
    print STDERR "Removing bad coordinate $keys[$bdix]\n";
    delete $hashref->{$keys[$bdix]};
  }
}

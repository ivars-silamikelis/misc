#!/usr/bin/perl;
use strict;
use threads;
use warnings;
use Bio::Perl;
use Getopt::Std;
use Bio::DB::Sam;
use Statistics::R;
use Switch;
my %opts;
getopt('irocQMFs', \%opts);
my %sams;

#$sams{"2"} = Bio::DB::Sam->new(-bam => $opts{"i"}, -fasta => $opts{"r"});
#$sams{"3"} = Bio::DB::Sam->new(-bam => $opts{"i"}, -fasta => $opts{"r"});
#$sams{"4"} = Bio::DB::Sam->new(-bam => $opts{"i"}, -fasta => $opts{"r"});

#$ref_length=
my $header=<<EOF;
##fileformat=VCFv4.1
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##reference=file:///$opts{"r"}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$opts{"s"}
EOF
print $header;
my @threads = initThreads();
my %thr;
for my $thread (@threads){
	my $start;
	my $end;
	switch ($thread){
		case 1 {$start=1;$end=10000}
		case 2 {$start=10001;$end=20000}
		case 3 {$start=20001;$end=30001}
		case 4 {$start=30001;$end=40001}
	}
	$sams{$thread} = Bio::DB::Sam->new(-bam => $opts{"i"}, -fasta => $opts{"r"});
	$thr{$thread}=threads->create(\&vcf,$thread,$sams{$thread},$start,$end);
	#=
}

for (1..4){
	$thr{"$_"}->join();
}	

sub vcf {
my $i = shift;
my $sam=shift;
my $start=shift;
my $end=shift;

print $start,"\t",$end,"\n";
#my $R = Statistics::R->new();

my @targets    = $sam->seq_ids;
for my $target (@targets){
	
	my $snp_caller = sub {

		my ($seqid,$pos,$p) = @_;
		#print STDERR $pos,"\n";
		my $refbase = $sam->segment($seqid,$pos,$pos)->dna;
        my $total=0;
		my $different=0;
		my $qbase;
		my $mq_avrg=0;
		my $qscore_total=0;
		my $fwd_count=0;
		my $rev_count=0;
		my $fisher_pvalue=1;
		my $ref_total=0;
		for my $pileup (@$p) {
		
		
			my $b = $pileup->alignment;
			my $mq = $b->qual;
			next if $pileup->indel or $pileup->is_refskip;      # don't deal with these
			$qbase  = substr($b->qseq,$pileup->qpos,1);
    	    next if $qbase =~ /[nN]/;

			my $qscore = $b->qscore->[$pileup->qpos];
			next if $qscore < 20;
		
			next if $mq < 20;
			$qscore_total+=$qscore;
      	  $total++;
     	   if ($refbase ne $qbase){
				$different++;
			if ($b->strand == 1){
				$fwd_count+=1;
			}
			if ($b->strand == -1){
				$rev_count+=1;
			}		
			} else {
				$ref_total++;		
			}
		}
		#print $fwd_count,"\t",$rev_count,"\n";
		if ($total >= 10 && $different/$total >= 0.6) {
			# fisher exact test for strand bias
		#	$R->set("fwd",$fwd_count);
		#	$R->set("rev",$rev_count);
		#	$R->run(q'result <- fisher.test(matrix(c(fwd,50,rev,50),ncol=2))');
		#	$R->run(q'pval <- result[["p.value"]]');
			#		q`result <- fisher.test(mat)`);
		#	$fisher_pvalue=$R->get("pval");
		    #next if $fisher_pvalue<0.01;
			my $avg_qscore=int($qscore_total/$total);
           
		   print "$seqid\t$pos\t.\t$refbase\t$qbase\t$qscore_total\t.\tFisherPvalueSB=$fisher_pvalue\tGT:AD:DP\t1:$different,".$ref_total.":$total\n";

			
        }
    };
	
#	$sam->pileup("$target:$start-$end",$snp_caller);

}


#$R->stop();

}




sub initThreads{
	my @initThreads;
	for (my $i=1; $i<=4;$i++){
		push(@initThreads,$i);
	}
	return @initThreads;
}

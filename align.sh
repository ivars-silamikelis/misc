#!/bin/bash;
i=$1
  name=${i/.fastq/}
  echo "aligning $name"
  tmap map2 -f ~/Data/Mtb/references/H37Rv/H37Rv.fasta \
  -r $i \
  -i fastq \
  -o 2 \
  -s $name.bam

  samtools sort $name.bam $name.sorted
  samtools view -H  $name.sorted.bam|sed "s/SM:NOSM/SM:$name/"|samtools reheader - $name.sorted.bam >$name.sorted.reheaded.bam
  samtools index $name.sorted.reheaded.bam
 rm  $name.sorted.bam
 rm  $name.bam

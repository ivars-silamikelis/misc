#!/bin/bash;
picard="/home/ivars/programs/picard-tools-1.51/picard-tools-1.51"
  i=$1
  name=${i/.fastq/}
  echo "aligning $name"
  tmap map3 -f $2 \
  -r $i \
  -i fastq \
  -o 2 \
  -s $name.bam

#  samtools sort $name.bam $name.sorted
#  samtools view -H  $name.sorted.bam|sed "s/SM:NOSM/SM:$name/"|samtools reheader - $name.sorted.bam >$name.sorted.reheaded.bam
#  samtools index $name.sorted.reheaded.bam
#  rm  $name.bam
#  rm  $name.sorted.bam
  samtools sort $name.bam $name.sorted
  echo "Marking for duplicates $name"
  java -jar $picard/MarkDuplicates.jar INPUT=$name.sorted.bam OUTPUT=$name.marked.bam METRICS_FILE=$name.metrics REMOVE_DUPLICATES=true
  samtools sort $name.marked.bam $name.marked.sorted
  samtools view -H  $name.marked.sorted.bam|sed "s/SM:NOSM/SM:$name/"|samtools reheader - $name.marked.sorted.bam >$name.marked.sorted.reheaded.bam
  samtools index $name.marked.sorted.reheaded.bam
  rm  $name.bam
  rm  $name.sorted.bam
  rm $name.marked.sorted.bam


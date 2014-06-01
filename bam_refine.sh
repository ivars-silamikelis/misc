#!/bin/bash
bam_files=("$@")

for i in "${bam_files[@]}";do
  echo "sorting $i.bam"
  samtools sort $i.bam $i.sorted
  echo "indexing $i.sorted.bam"
  samtools index $i.sorted.bam
  echo "removing $i.bam"
  rm $i.bam
done

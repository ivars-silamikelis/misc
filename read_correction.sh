picard="/home/ivars/programs/picard-tools-1.113/picard-tools-1.113"
name=${1/.bam/}
echo "Marking for duplicates $name"
java -jar $picard/MarkDuplicates.jar INPUT=$1 OUTPUT=$name.marked.bam METRICS_FILE=$name.metrics
#samtools index $name.marked.bam
i=$name.marked.bam
name=${i/.marked.bam/}

echo "processing $i"
java -Xmx2g -jar $gatk \
  -T RealignerTargetCreator \
  -R $ref \
  -I $i \
  -o $name.intervals


java -Xmx2g -jar $gatk \
   -T IndelRealigner \
   -R $ref \
   -I $i \
   -targetIntervals $name.intervals \
   -o $name.marked.realigned.bam
   samtools index $name.marked.realigned.bam
i=$name.marked.realigned.bam

#pair endiem vajag fixet readus

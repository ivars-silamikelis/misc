gatk=/home/ivars/programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar
vcfutils=/usr/local/bin/vcfutils.pl
varscan=/home/ivars/programs/VarScan/VarScan.v2.3.6.jar
ref=/home/ivars/Data/Mtb/references/H37Rv/H37Rv.fasta
shopt -s nullglob
bams=(*realigned.bam)
for i in ${bams[@]};do 
	gatk_ins+="-I $i "
	sam_ins+=" "$i
done
echo ${gatk_ins[@]}
threads=4

	  samtools mpileup -Iuf $ref ${sam_ins[@]} | bcftools view -bvcg - > all.raw.bcf  
	  bcftools view all.raw.bcf | perl $vcfutils varFilter -D100 > all.samtools.vcf
	  
	echo ${sam_ins[@]} > all.name
	
	  samtools mpileup -B -f $ref ${sam_ins[@]} | java -jar $varscan mpileup2snp --min-coverage 10 --min-var-freq 0.8 --output-vcf 1 --vcf-sample-list all.name >all.varscan.vcf
	  rm all.name
	  
	  java -jar $gatk \
	     -T UnifiedGenotyper \
	     -ploidy 1 \
	     -glm SNP \
	     -R $ref \
	  	-nt $threads \
	     ${gatk_ins[@]} \
	     -stand_emit_conf 10 \
	     -stand_call_conf 30 \
	     -o all.gatk_ug.vcf

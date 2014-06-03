for i in *.samtools.vcf;do
	name=${i/.samtools.vcf}
	echo $name >>progress.txt
	python process_vcf.py $name
done

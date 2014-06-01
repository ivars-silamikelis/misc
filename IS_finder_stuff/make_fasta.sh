for i in $(ls -d */);do 
	name=${i///}	
	echo $name
	cat $name/*.crd >$name/all.txt
	perl crd2fasta_revcom_extended.pl $name $name/all.txt >dump
	
done

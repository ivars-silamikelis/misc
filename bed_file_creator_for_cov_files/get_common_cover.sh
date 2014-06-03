#!/bin/bash
echo -e "AL123456\t1\t5000000" >all.txt
for i in *.cov;do
	bedtools intersect -a $i -b all.txt >all.tmp.txt
	mv all.tmp.txt all.txt
done
bedtools sort -i all.txt >sorted.all.txt
bedtools merge -i sorted.all.txt  > all.txt
sed -i "1i accession\tstart\tend" all.txt

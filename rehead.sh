#!/bin/bash

cd SN2-1_amp
samtools view -H SN2-1_amp.trimmed.Q20.P50.H37Rv.sorted.bam|sed 's/SM:NOSM/SM:SN2-1_amp/'|samtools reheader - SN2-1_amp.trimmed.Q20.P50.H37Rv.sorted.bam >SN2-1_amp_trimmed.Q20.P50.H37Rv.sorted.reheaded.bam
samtools index SN2-1_amp_trimmed.Q20.P50.H37Rv.sorted.reheaded.bam


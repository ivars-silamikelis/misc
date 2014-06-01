#!/bin/bash
name=${1/_in.iontor.fastq/}
echo argv: $name
mira --project=$name --job=denovo,genome,accurate,iontor -notraceinfo


#!/bin/bash

url=https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
filename=$(basename ${url})

if [ ! -e ${filename} ]; then
	wget --quiet -O - ${url} | gzip -d | bgzip > ${filename}
	samtools faidx ${filename}
fi

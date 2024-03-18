#!/bin/bash

url=http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
filename=$(basename ${url})

if [ ! -e ${filename} ]; then
	wget --quiet -O - ${url} | gzip -d | bgzip > ${filename}
	samtools faidx ${filename}
fi

#!/bin/bash

set -e

VG_DIR=$(dirname $(dirname $(dirname $(dirname ${BASH_SOURCE[0]}))))

GFF3_TO_GENEPRED=$(which gff3ToGenePred)
if [ -z ${GFF3_TO_GENEPRED} ]; then
  echo "Downloading gff3ToGenePred command line tool"
  wget hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred
  chmod a+x gff3ToGenePred
  GFF3_TO_GENEPRED=./gff3ToGenePred
fi


# Ensembl
echo "Downloading ENSEMBL"
mkdir -p ensembl
cd ensembl

if false; then
  for release in 102 103 104; do
    gff=Homo_sapiens.GRCh38.${release}.gff3.gz;
    GENEPRED="$(basename ${gff} .gff.gz).genePred"
    echo "GenePred = ${GENEPRED}";

    if [[ ! -e ${gff} ]]; then
      wget ftp://ftp.ensembl.org/pub/release-${release}/gff3/homo_sapiens/${gff}
    fi

    if [[ ! -e ${GENEPRED} ]]; then
      ${GFF3_TO_GENEPRED} -processAllGeneChildren ${gff} ${GENEPRED}
    fi

    echo "Inserting gene annotation"

    python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=Ensembl \
      --gff ${gff} \
      --genePred ${GENEPRED}

  done
fi

cd ..

# RefSeq
echo "Downloading RefSeq"

mkdir -p refseq
cd refseq

for release in 109.20210226 109.20210514; do
  gff=GCF_000001405.39_GRCh38.p13_genomic.${release}.gff.gz
  GENEPRED="$(basename ${gff} .gff.gz).genePred"

  if [[ ! -e ${gff} ]]; then
    echo "Downloading '${gff}'"
    # FTP is corrupt, trying http

    wget http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/${release}/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz
    mv GCF_000001405.39_GRCh38.p13_genomic.gff.gz ${gff}
  fi

  if [[ ! -e ${GENEPRED} ]]; then
      ${GFF3_TO_GENEPRED} -processAllGeneChildren -maxParseErrors=-1 -geneNameAttr=Name -rnaNameAttr=transcript_id ${gff} ${GENEPRED}
  fi

  python3.8 ${VG_DIR}/manage.py import_gene_annotation --genome-build=GRCh38 --replace --annotation-consortium=RefSeq \
    --gff ${gff} \
    --genePred ${GENEPRED}

done

cd ..
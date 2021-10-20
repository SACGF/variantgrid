#!/bin/bash

# v81 (points to 75) and earlier at GTFs that don't have transcript versions - just skip them

#82 is first GFF3 for GRCh37
#83 has no data
#84 is 82 again
#86 is 85 again
pyreference_args=()
for release in 82 85 87; do
  filename=Homo_sapiens.GRCh37.${release}.gff3.gz
  url=ftp://ftp.ensembl.org/pub/grch37/release-${release}/gff3/homo_sapiens/${filename}
  pyreference_file=${filename}.json.gz
  if [[ ! -e ${filename} ]]; then
    wget ${url}
  fi
  if [[ ! -e ${pyreference_file} ]]; then
    pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
  fi
  pyreference_args+=(--pyreference-json ${pyreference_file})
done

VG_DIR=$(dirname $(dirname $(dirname $(dirname ${BASH_SOURCE[0]}))))
merged_file="vg_gene_annotation.ensembl.grch37.json.gz"
if [[ ! -e ${merged_file} ]]; then
  python3 ${VG_DIR}/manage.py import_gene_annotation --annotation-consortium=Ensembl --genome-build=GRCh37 \
    --save-merged-file=${merged_file} ${pyreference_args[@]}
fi
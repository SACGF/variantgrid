# Need this in your path

GFF3_TO_GENEPRED=$(which gff3ToGenePred)
if [ -z ${GFF3_TO_GENEPRED} ]; then
  wget hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred
  chmod a+x gff3ToGenePred
  GFF3_TO_GENEPRED=./gff3ToGenePred
fi

for gtf in *.gtf.gz
do
    GENEPRED="$(basename ${gtf} .gtf.gz).genePred"
    ${GFF3_TO_GENEPRED} -processAllGeneChildren -geneNameAttr=Name -rnaNameAttr=transcript_id ${gtf} ${GENEPRED}
done

for gff3 in Homo_sapiens.*.gff3.gz
do
    GENEPRED="$(basename ${gff3} .gff3.gz).genePred"
    ${GFF3_TO_GENEPRED} -processAllGeneChildren ${gff3} ${GENEPRED}
done

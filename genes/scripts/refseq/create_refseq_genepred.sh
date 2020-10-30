# Need this in your path

GFF3_TO_GENEPRED=$(which gff3ToGenePred)
if [ -z ${GFF3_TO_GENEPRED} ]; then
  wget hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred
  chmod a+x gff3ToGenePred
  GFF3_TO_GENEPRED=./gff3ToGenePred
fi


for gff3 in *.gff3.gz
do
    GENEPRED="$(basename ${gff3} .gff3.gz).genePred"
    ${GFF3_TO_GENEPRED} -processAllGeneChildren -geneNameAttr=Name -rnaNameAttr=transcript_id ${gff3} ${GENEPRED}
done

for gff in *.gff.gz
do
    GENEPRED="$(basename ${gff} .gff.gz).genePred"
    ${GFF3_TO_GENEPRED} -processAllGeneChildren -geneNameAttr=Name -rnaNameAttr=transcript_id ${gff} ${GENEPRED}
done

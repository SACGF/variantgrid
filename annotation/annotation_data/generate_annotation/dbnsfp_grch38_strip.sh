#!/bin/bash

set -e

# Download 4.8 from https://sites.google.com/site/jpopgen/dbNSFP

# https://m.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#dbnsfp

# zcat dbNSFP4.8a_variant.chr1.gz | head -n1 > header.txt
# mkdir -p /tmp/dbnsfp38
# zgrep -h -v ^#chr dbNSFP4.8a_variant.chr* | sort -T /tmp/dbnsfp38 -k1,1 -k2,2n - | cat header.txt - | bgzip -c > dbNSFP4.8a_grch38.gz
# tabix -s 1 -b 2 -e 2 dbNSFP4.8a_grch38.gz


# All of this python is just to get the columns used in cut and tabix args

# Get dbNSFP fields used by VariantGrid - run python3 manage.py shell
# In [12]: ",".join(ColumnVEPField.get_source_fields(vep_plugin='d'))                                                                                                                                     

# Get column names from dbNSFP data file
# df = pd.read_csv("header.txt", sep='\t', index_col=None, nrows=0)
# vep_fields = 'GERP++_RS,Interpro_domain,CADD_raw_rankscore,REVEL_rankscore,BayesDel_noAF_rankscore,ClinPred_rankscore,VEST4_rankscore,MetaLR_rankscore,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence,AlphaMissense_rankscore,AlphaMissense_pred'
# new_vep_fields = "CADD_raw,REVEL_score,BayesDel_noAF_score,ClinPred_score,ClinPred_pred,VEST4_score,MetaLR_score,MetaLR_pred,AlphaMissense_score,MutPred_score,MutPred_rankscore,MutPred_protID,MutPred_AAchange,MutPred_Top5features"
# columns = ['#chr', 'pos(1-based)', 'ref', 'alt', 'aaref', 'aaalt', 'Ensembl_transcriptid'] + vep_fields.split(",") + new_vep_fields.split(",")
# cols = []
# for i in columns:
#    cols.append(list(df.columns).index(i) + 1)
# print(",".join([str(c) for c in sorted(cols)]))
# columns are: '1,2,3,4,5,6,15,68,69,73,74,75,83,84,85,86,87,88,89,105,106,108,109,110,138,139,140,146,147,148,149,150,151,152,189,447'

TMP_DIR=/tmp/dbnsfp38
CUT_COLUMNS="1,2,3,4,5,6,15,68,69,73,74,75,83,84,85,86,87,88,89,105,106,108,109,110,138,139,140,146,147,148,149,150,151,152,189,447"
SEQ_COL=1  # chr
POS_COL=2  # pos(1-based)

version=4.8a
out_file=dbNSFP${version}_grch38.stripped

mkdir -p ${TMP_DIR}

# Sort chromosomes individually as that's much more efficient
cat header.txt | cut -f ${CUT_COLUMNS} > ${OUT_FILE}
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
    zgrep -h -v ^#chr dbNSFP${version}_variant.chr${chrom}.gz | cut -f ${CUT_COLUMNS} | sort -T ${TMP_DIR} -k${SEQ_COL},${SEQ_COL} -k${POS_COL},${POS_COL}n - >> ${out_file}
done

bgzip ${out_file}
tabix -s ${SEQ_COL} -b ${POS_COL} -e 2 ${out_file}.gz # cols are: 1=chr, 2=pos

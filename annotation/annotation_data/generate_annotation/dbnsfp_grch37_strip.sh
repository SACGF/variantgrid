#!/bin/bash

set -e

# Download 4.5 from https://sites.google.com/site/jpopgen/dbNSFP

# https://m.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#dbnsfp
# zcat dbNSFP4.5a_variant.chr1.gz | head -n1 > header.txt
# mkdir /tmp/dbsnp37
# zgrep -h -v ^#chr dbNSFP4.5a_variant.chr* | awk '$8 != "." ' | sort -T /tmp/dbsnp37 -k8,8 -k9,9n - | cat header.txt - | bgzip -c > dbNSFP4.5a_grch37.gz
# tabix -s 8 -b 9 -e 9 dbNSFP4.5a.grch37.gz


# All of this python is just to get the columns used in cut and tabix args at bottom of this file

# Get dbNSFP fields used by VariantGrid - run python3 manage.py shell
# In [12]: ",".join(ColumnVEPField.get_source_fields(vep_plugin='d'))                                                                                                                                     

# Get column names from dbNSFP data file
# import pandas as pd
# df = pd.read_csv("header.txt", sep='\t', index_col=None, nrows=0)
# vep_fields = 'GERP++_RS,Interpro_domain,CADD_raw_rankscore,REVEL_rankscore,BayesDel_noAF_rankscore,ClinPred_rankscore,VEST4_rankscore,MetaLR_rankscore,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence,AlphaMissense_rankscore,AlphaMissense_pred'
# columns = ['ref', 'alt', 'aaref', 'aaalt', 'hg19_chr', 'hg19_pos(1-based)', 'Ensembl_transcriptid'] + vep_fields.split(",")
# cols = []
# for i in columns:
#    cols.append(list(df.columns).index(i) + 1)
# print(",".join([str(c) for c in sorted(cols)]))
# columns are: '3,4,5,6,8,9,15,69,74,84,106,109,139,140,142,143,144,145,146,148,185,705'

CUT_COLUMNS="3,4,5,6,8,9,15,69,74,84,106,109,139,140,142,143,144,145,146,148,185,705"
SEQ_COL=5 # hg19_chr
POS_COL=6 # hg19_pos(1-based)
OUT_FILE=dbNSFP4.5a.grch37.stripped
TMP_DIR=/tmp # /hpcfs/groups/phoenix-hpc-sacgf/scratch/dbnsfp_GRCh37
mkdir -p ${TMP_DIR}

# Sort chromosomes individually as that's much more efficient
cat header.txt | cut -f ${CUT_COLUMNS} > ${OUT_FILE}
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
    zgrep -h -v ^#chr dbNSFP4.5a_variant.chr${chrom}.gz | awk '$8 != "." ' | cut -f ${CUT_COLUMNS} | sort -T ${TMP_DIR} -k${SEQ_COL},${SEQ_COL} -k${POS_COL},${POS_COL}n - >> ${OUT_FILE}
done

bgzip ${OUT_FILE}
tabix -s ${SEQ_COL} -b ${POS_COL} -e ${POS_COL} ${OUT_FILE}.gz


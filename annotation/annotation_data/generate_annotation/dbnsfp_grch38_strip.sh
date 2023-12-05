#!/bin/bash

set -e

# Download 4.5 from https://sites.google.com/site/jpopgen/dbNSFP

# https://m.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#dbnsfp

# zcat dbNSFP4.5a_variant.chr1.gz | head -n1 > header.txt
# mkdir /tmp/dbsnp38
# zgrep -h -v ^#chr dbNSFP4.5a_variant.chr* | sort -T /tmp/dbsnp38 -k1,1 -k2,2n - | cat header.txt - | bgzip -c > dbNSFP4.5a_grch38.gz
# tabix -s 1 -b 2 -e 2 dbNSFP4.5a_grch38.gz


# All of this python is just to get the columns used in cut and tabix args

# Get dbNSFP fields used by VariantGrid - run python3 manage.py shell
# In [12]: ",".join(ColumnVEPField.get_source_fields(vep_plugin='d'))                                                                                                                                     

# Get column names from dbNSFP data file
# df = pd.read_csv("header.txt", sep='\t', index_col=None, nrows=0)
# vep_fields = 'GERP++_RS,Interpro_domain,CADD_raw_rankscore,REVEL_rankscore,BayesDel_noAF_rankscore,ClinPred_rankscore,VEST4_rankscore,MetaLR_rankscore,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence,AlphaMissense_rankscore,AlphaMissense_pred'
# columns = ['#chr', 'pos(1-based)', 'ref', 'alt', 'aaref', 'aaalt', 'Ensembl_transcriptid'] + vep_fields.split(",")
# cols = []
# for i in columns:
#    cols.append(list(df.columns).index(i) + 1)
# print(",".join([str(c) for c in sorted(cols)]))
# columns are: '1,2,3,4,5,6,15,69,74,84,106,109,139,140,142,143,144,145,146,148,185,705'

CUT_COLUMNS="1,2,3,4,5,6,15,69,74,84,106,109,139,140,142,143,144,145,146,148,185,705"
OUT_FILE=dbNSFP4.5a.grch38.stripped.gz
TMP_DIR=/tmp # /hpcfs/groups/phoenix-hpc-sacgf/scratch/dbnsfp4.5_GRCh38
mkdir -p ${TMP_DIR}

# Sort chromosomes individually as that's much more efficient
cat header.txt | cut -f ${CUT_COLUMNS} | bgzip > ${OUT_FILE}
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
    zgrep -h -v ^#chr dbNSFP4.5a_variant.chr${chrom}.gz | cut -f ${CUT_COLUMNS} | sort -T ${TMP_DIR} -k1,1 -k2,2n - | bgzip >> ${OUT_FILE}
done

tabix -s 1 -b 2 -e 2 ${OUT_FILE} # cols are: 1=chr, 2=pos

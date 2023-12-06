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

# Note: We can't do this per-contig then join them, as some variants switch contigs between builds
CUT_COLUMNS="3,4,5,6,8,9,15,69,74,84,106,109,139,140,142,143,144,145,146,148,185,705"
SEQ_COL=3  # hg19_chr was col 5 (but 3rd after cut)
POS_COL=4  # hg19_pos(1-based) was 6 (but 4th after cut)

version=4.5a
out_vcf=dbNSFP${version}_grch37.gz

# cd /hpcfs/groups/phoenix-hpc-sacgf/reference/annotation/dbnsfp/dbnsfp4.5

zcat dbNSFP${version}_variant.chr1.gz | head -n1 > h
zgrep -h -v ^#chr dbNSFP${version}_variant.chr* | awk '$8 != "." ' | sort -T ${TMP_DIR} -k8,8 -k9,9n - | cat h - | bgzip -c > ${out_vcf}
zcat ${out_vcf} | cut -f  ${CUT_COLUMNS}  > dbNSFP${version}_grch37.stripped
bgzip dbNSFP${version}_grch37.stripped

tabix -s ${SEQ_COL} -b ${POS_COL} -e ${POS_COL} dbNSFP${version}_grch37.stripped.gz

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



IN_FILE=dbNSFP4.5a.grch37.gz
OUT_FILE=dbNSFP4.5a.grch37.stripped.gz

# Header needs to start with #
(echo -n "#" ; zcat ${IN_FILE} | cut -f 3,4,5,6,8,9,15,69,74,84,106,109,139,140,142,143,144,145,146,148,185,705 ) | bgzip > ${OUT_FILE}
tabix -s 5 -b 6 -e 6 ${OUT_FILE} # cols are: 1=ref, 2=alt, 3=chr, 4=pos


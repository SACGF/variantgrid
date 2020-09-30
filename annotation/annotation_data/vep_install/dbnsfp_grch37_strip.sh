#!/bin/bash

set -e

# All of this python is just to get the columns used in cut and tabix args

# Get dbNSFP fields used by VariantGrid - run python3 manage.py shell
# In [12]: ",".join(ColumnVEPField.get_source_fields(vep_plugin='d'))                                                                                                                                     
# Out[12]: 'MutationAssessor_pred,Polyphen2_HVAR_pred,REVEL_score,MutationTaster_pred,GERP++_RS,FATHMM_pred,Interpro_domain,CADD_raw,CADD_phred'

# Get column names from dbNSFP data file
# df = pd.read_csv("./dbNSFP4.0a.grch37.gz", sep='\t', index_col=None, nrows=0)
# vep_fields = 'MutationAssessor_pred,Polyphen2_HVAR_pred,REVEL_score,MutationTaster_pred,GERP++_RS,FATHMM_pred,Interpro_domain,CADD_raw,CADD_phred'
# columns = ['ref', 'alt',  'hg19_chr', 'hg19_pos(1-based)', 'Ensembl_transcriptid'] + vep_fields.split(",")
# cols = []
# for i in columns:
#    cols.append(list(df.columns).index(i) + 1)
# ",".join([str(c) for c in sorted(cols)])
# Out[7]: '3,4,8,9,15,48,55,60,63,79,102,104,137,373'


IN_FILE=dbNSFP4.0a.grch37.gz
OUT_FILE=dbNSFP4.0a.grch37.stripped.gz

# Header needs to start with #
(echo -n "#" ; zcat ${IN_FILE} | cut -f 3,4,8,9,15,48,55,60,63,79,102,104,137,373 ) | bgzip > ${OUT_FILE}
tabix -s 3 -b 4 -e 4 ${OUT_FILE} # cols are: 1=ref, 2=alt, 3=chr, 4=pos


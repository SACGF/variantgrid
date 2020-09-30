#!/bin/bash

set -e

# All of this python is just to get the columns used in cut and tabix args

# Get dbNSFP fields used by VariantGrid - run python3 manage.py shell
# In [12]: ",".join(ColumnVEPField.get_source_fields(vep_plugin='d'))                                                                                                                                     
# Out[12]: 'MutationAssessor_pred,Polyphen2_HVAR_pred,REVEL_score,MutationTaster_pred,GERP++_RS,FATHMM_pred,Interpro_domain,CADD_raw,CADD_phred'

# Get column names from dbNSFP data file
# df = pd.read_csv("./dbNSFP4.0a.grch38.gz", sep='\t', index_col=None, nrows=0)
# vep_fields = 'MutationAssessor_pred,Polyphen2_HVAR_pred,REVEL_score,MutationTaster_pred,GERP++_RS,FATHMM_pred,Interpro_domain,CADD_raw,CADD_phred'
# columns = ['#chr', 'pos(1-based)', 'ref', 'alt', 'Ensembl_transcriptid'] + vep_fields.split(",")
# cols = []
# for i in columns:
#    cols.append(list(df.columns).index(i) + 1)
# ",".join([str(c) for c in sorted(cols)])
# Out[18]: '1,2,3,4,15,48,55,60,63,79,102,104,137,373'

IN_FILE=dbNSFP4.0a.grch38.gz
OUT_FILE=dbNSFP4.0a.grch38.stripped.gz

# Header needs to start with #
(echo -n "#" ; zcat ${IN_FILE} | cut -f 1,2,3,4,15,48,55,60,63,79,102,104,137,373 ) | bgzip > ${OUT_FILE}
tabix -s 1 -b 2 -e 2 ${OUT_FILE} # cols are: 1=chr, 2=pos


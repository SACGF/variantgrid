#!/bin/bash

set -e

# All of this python is just to get the columns used in cut and tabix args

# Get dbNSFP fields used by VariantGrid - run python3 manage.py shell
# In [12]: ",".join(ColumnVEPField.get_source_fields(vep_plugin='d'))                                                                                                                                     

# Get column names from dbNSFP data file
# df = pd.read_csv("./dbNSFP4.3a_grch38.gz", sep='\t', index_col=None, nrows=0)
# vep_fields = 'GERP++_RS,Interpro_domain,CADD_raw_rankscore,REVEL_rankscore,BayesDel_noAF_rankscore,ClinPred_rankscore,VEST4_rankscore,MetaLR_rankscore'
# columns = ['#chr', 'pos(1-based)', 'ref', 'alt', 'aaref', 'aaalt', 'Ensembl_transcriptid'] + vep_fields.split(",")
# cols = []
# for i in columns:
#    cols.append(list(df.columns).index(i) + 1)
# ",".join([str(c) for c in sorted(cols)])
# columns are: '1,2,3,4,5,6,15,69,74,84,104,107,119,156,640'

# Download 4.3 from https://sites.google.com/site/jpopgen/dbNSFP

# https://m.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#dbnsfp

# zcat dbNSFP4.3a_variant.chr1.gz | head -n1 > h
# zgrep -h -v ^#chr dbNSFP4.3a_variant.chr* | sort -T /path/to/tmp_folder -k1,1 -k2,2n - | cat h - | bgzip -c > dbNSFP4.3a_grch38.gz
# tabix -s 1 -b 2 -e 2 dbNSFP4.3a_grch38.gz


IN_FILE=dbNSFP4.3a_grch38.gz
OUT_FILE=dbNSFP4.3a_grch38.stripped.gz

# Header needs to start with #
(echo -n "#" ; zcat ${IN_FILE} | cut -f 1,2,3,4,5,6,15,69,74,84,104,107,119,156,640 ) | bgzip > ${OUT_FILE}
tabix -s 1 -b 2 -e 2 ${OUT_FILE} # cols are: 1=chr, 2=pos


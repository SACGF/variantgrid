#!/bin/bash

set -e

# All of this python is just to get the columns used in cut and tabix args

# Get dbNSFP fields used by VariantGrid - run python3 manage.py shell
# In [12]: ",".join(ColumnVEPField.get_source_fields(vep_plugin='d'))                                                                                                                                     

# Get column names from dbNSFP data file
# df = pd.read_csv("./dbNSFP4.3a_grch37.gz", sep='\t', index_col=None, nrows=0)
# vep_fields = 'GERP++_RS,Interpro_domain,CADD_raw_rankscore,REVEL_rankscore,BayesDel_noAF_rankscore,ClinPred_rankscore,VEST4_rankscore,MetaLR_rankscore'
# columns = ['ref', 'alt', 'aaref', 'aaalt', 'hg19_chr', 'hg19_pos(1-based)', 'Ensembl_transcriptid'] + vep_fields.split(",")
# cols = []
# for i in columns:
#    cols.append(list(df.columns).index(i) + 1)
# ",".join([str(c) for c in sorted(cols)])
# columns are: '3,4,5,6,8,9,15,69,74,84,104,107,119,156,640'


# Download 4.3 from https://sites.google.com/site/jpopgen/dbNSFP

# https://m.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#dbnsfp
# zcat dbNSFP4.3a_variant.chr1.gz | head -n1 > h
# zgrep -h -v ^#chr dbNSFP4.3a_variant.chr* | awk '$8 != "." ' | sort -T /path/to/tmp_folder -k8,8 -k9,9n - | cat h - | bgzip -c > dbNSFP4.3a_grch37.gz
# tabix -s 8 -b 9 -e 9 dbNSFP4.3a_grch37.gz

IN_FILE=dbNSFP4.3a.grch37.gz
OUT_FILE=dbNSFP4.3a.grch37.stripped.gz

# Header needs to start with #
(echo -n "#" ; zcat ${IN_FILE} | cut -f 3,4,5,6,8,9,15,69,74,84,104,107,119,156,640 ) | bgzip > ${OUT_FILE}
tabix -s 5 -b 6 -e 6 ${OUT_FILE} # cols are: 1=ref, 2=alt, 3=chr, 4=pos


update expression_cuffdiffrecord cdr
set gene_id = eurs.ensembl_gene_id
from snpdb_ensemblucscrefseq eurs
inner join snpdb_ucscalias ua on (ua.ucsc_id = eurs.ucsc_id)
where upper(cdr.reference_id) = upper(ua.alias)
and cdr.cuff_diff_file_id = %(cuff_diff_file_id)s;

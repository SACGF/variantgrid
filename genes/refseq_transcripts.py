# Copied from pyhgvs.models.hgvs_name.get_refseq_type
# @see https://github.com/counsyl/hgvs/blob/master/pyhgvs/__init__.py

REFSEQ_PREFIXES = [
    ('AC_', 'genomic',
     'Complete genomic molecule, usually alternate assembly'),
    ('NC_', 'genomic',
     'Complete genomic molecule, usually reference assembly'),
    ('NG_', 'genomic', 'Incomplete genomic region'),
    ('NT_', 'genomic', 'Contig or scaffold, clone-based or WGS'),
    ('NW_', 'genomic', 'Contig or scaffold, primarily WGS'),
    ('NS_', 'genomic', 'Environmental sequence'),
    ('NZ_', 'genomic', 'Unfinished WGS'),
    ('NM_', 'mRNA', ''),
    ('NR_', 'RNA', ''),
    ('XM_', 'mRNA', 'Predicted model'),
    ('XR_', 'RNA', 'Predicted model'),
    ('AP_', 'Protein', 'Annotated on AC_ alternate assembly'),
    ('NP_', 'Protein', 'Associated with an NM_ or NC_ accession'),
    ('YP_', 'Protein', ''),
    ('XP_', 'Protein', 'Predicted model, associated with an XM_ accession'),
    ('ZP_', 'Protein', 'Predicted model, annotated on NZ_ genomic records'),
]


REFSEQ_PREFIX_LOOKUP = dict(
    (prefix, (kind, description))
    for prefix, kind, description in REFSEQ_PREFIXES
)


def get_refseq_type(name):
    """
    Return the RefSeq type for a refseq name.
    """
    prefix = name[:3]
    return REFSEQ_PREFIX_LOOKUP.get(prefix, (None, ''))[0]


def transcript_is_lrg(transcript_accession: str):
    return transcript_accession.startswith("LRG_")

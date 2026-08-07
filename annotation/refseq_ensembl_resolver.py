"""
RefSeq <-> Ensembl transcript resolution for VEP annotation insertion.

Some VEP plugin outputs (currently dbNSFP per-transcript scores) come as
&-separated arrays parallel to an Ensembl transcript-id array. To pick the
value matching the picked VEP transcript on a RefSeq pipeline we need to
translate NM_xxx -> ENSTxxx (unversioned, since dbNSFP's array is unversioned).

This module hides that translation behind a `RefSeqEnsemblTranscriptResolver`
protocol so the inserter doesn't bake in the lookup strategy. The current
implementation backs onto `DBNSFPGeneAnnotation` (gene-level mapping, GRCh37+38).
A MANE-backed implementation could be added later for GRCh38, or a
TranscriptVersion-join implementation for full per-transcript mapping.

The resolver name is persisted on `VariantAnnotationVersion.transcript_resolver`
so we can tell, after the fact, which strategy was used to populate the
per-transcript dbNSFP scores for a given annotation run.
"""

from typing import Optional, Protocol


class RefSeqEnsemblTranscriptResolver(Protocol):
    """ Translate between RefSeq and Ensembl transcript identifiers (unversioned).
        Implementations should be cheap to call repeatedly — cache as needed. """

    name: str  # stable identifier persisted on VariantAnnotationVersion.transcript_resolver

    def refseq_to_ensembl(self, gene_symbol: str, refseq: str) -> Optional[str]:
        """ Returns an unversioned Ensembl transcript id (e.g. 'ENST00000262340')
            for the given gene_symbol + RefSeq accession (e.g. 'NM_005845'),
            or None if no mapping is available. """
        ...

    def ensembl_to_refseq(self, gene_symbol: str, ensembl: str) -> Optional[str]:
        """ Inverse of refseq_to_ensembl. """
        ...


class DBNSFPGeneResolver:
    """ Backed by `annotation.models.DBNSFPGeneAnnotation`, which carries one
        (refseq_transcript, ensembl_transcript) pair per gene_symbol — i.e. dbNSFP's
        chosen canonical pair. Build-agnostic (GRCh37 + GRCh38).

        Trade-off: gene-level only, so all RefSeq accessions on a gene resolve to
        the same ENST. For the current consumers (dbNSFP per-transcript scores)
        this is acceptable because score variation across transcripts of the same
        gene is usually small, and falling back to "no match" leaves the inserter
        free to use a representative-value collapse instead.
    """

    name = "dbnsfp_gene"

    def __init__(self):
        # Imported lazily to avoid Django app-loading order issues at module import time.
        from annotation.models.models import DBNSFPGeneAnnotation, DBNSFPGeneAnnotationVersion
        self._refseq_to_ensembl: dict[tuple[str, str], str] = {}
        self._ensembl_to_refseq: dict[tuple[str, str], str] = {}
        version = DBNSFPGeneAnnotationVersion.latest()
        if version is None:
            return
        qs = DBNSFPGeneAnnotation.objects.filter(version=version).values_list(
            "gene_symbol_id", "refseq_transcript_id", "ensembl_transcript_id",
        )
        for symbol, refseq, ensembl in qs:
            if refseq and ensembl:
                self._refseq_to_ensembl[(symbol, refseq)] = ensembl
                self._ensembl_to_refseq[(symbol, ensembl)] = refseq

    def refseq_to_ensembl(self, gene_symbol: str, refseq: str) -> Optional[str]:
        return self._refseq_to_ensembl.get((gene_symbol, refseq))

    def ensembl_to_refseq(self, gene_symbol: str, ensembl: str) -> Optional[str]:
        return self._ensembl_to_refseq.get((gene_symbol, ensembl))

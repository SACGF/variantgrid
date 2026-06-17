from typing import Optional

from cdot.hgvs.dataproviders import LocalDataProvider, FastaSeqFetcher, ChainedSeqFetcher
from django.conf import settings
from django.db.models import Q
from hgvs.exceptions import HGVSDataNotAvailableError

from genes.models import TranscriptVersion, TranscriptVersionSequenceInfo, NoTranscript, MANE
from genes.models_enums import MANEStatus
from genes.transcripts_utils import get_refseq_type
from snpdb.models import Contig

# cdot tag spellings expected by cdot's DEFAULT_TAG_PRIORITY list
_MANE_STATUS_TO_CDOT_TAG = {
    MANEStatus.MANE_SELECT: "MANE_Select",
    MANEStatus.MANE_PLUS_CLINICAL: "MANE_Plus_Clinical",
}


def is_contig(ac):
    if refseq_type := get_refseq_type(ac):
        return refseq_type == "genomic"
    return False


class DBTranscriptSeqFetcher:
    """ We store some transcripts in the database, use them if we have them   """
    source = "VG DB Transcript sequences (transcript fasta/API)"

    def __init__(self, retrieve_transcripts: Optional[bool] = None):
        if retrieve_transcripts is None:
            retrieve_transcripts = settings.HGVS_RETRIEVE_TRANSCRIPT_SEQUENCE
        self.retrieve_transcripts = retrieve_transcripts

    def fetch_seq(self, ac, start_i=None, end_i=None):
        if is_contig(ac):
            raise HGVSDataNotAvailableError("This SeqFetcher implementation doesn't handle contigs")

        tvsi = None
        try:
            tvsi = TranscriptVersionSequenceInfo.get(ac, retrieve=self.retrieve_transcripts)
        except NoTranscript:
            pass

        if tvsi is None:
            raise HGVSDataNotAvailableError(f"No database sequence for {ac}")
        transcript_seq = tvsi.sequence
        if start_i is None:
            start_i = 0
        if end_i is None:
            end_i = len(transcript_seq)
        return transcript_seq[start_i:end_i]


class SingleBuildFastaSeqFetcher(FastaSeqFetcher):
    """ We want this to fail with ContigNotInBuildError if asked for a contig not in the
        the genome build provided """
    def __init__(self, genome_build):
        super().__init__(genome_build.genome_fasta.filename)
        self.genome_build = genome_build

    def fetch_seq(self, ac, start_i=None, end_i=None):
        # If FastaSeqFetcher is called with contigs from another build it will try and lookup transcripts
        # Just fail ASAP
        if is_contig(ac):
            if ac not in self.contig_fastas:
                raise Contig.ContigNotInBuildError(f"Contig '{ac}' not in genome build '{self.genome_build.name}'")

        return super().fetch_seq(ac, start_i=start_i, end_i=end_i)


class DjangoTranscriptDataProvider(LocalDataProvider):

    def __init__(self, genome_build, retrieve_transcripts: Optional[bool] = None):
        """ retrieve_transcripts: attempt to retrieve transcript if not present (slower) """
        self.db_transcript_seqfetcher = DBTranscriptSeqFetcher(retrieve_transcripts)
        self.fasta_seqfetcher = SingleBuildFastaSeqFetcher(genome_build)
        seqfetcher = ChainedSeqFetcher(self.db_transcript_seqfetcher, self.fasta_seqfetcher)

        super().__init__(assemblies=[genome_build.name], seqfetcher=seqfetcher)
        self.genome_build = genome_build

    def _get_transcript_ids_for_gene(self, gene):
        tv_qs = TranscriptVersion.objects.filter(genome_build=self.genome_build,
                                                 gene_version__gene_symbol__symbol=gene)
        return list(tv_qs.values_list("transcript", flat=True))

    def _get_contig_interval_tree(self, alt_ac):
        raise NotImplementedError()

    def _get_transcript(self, tx_ac):
        tv = TranscriptVersion.get_transcript_version(self.genome_build, tx_ac)
        return tv.data

    def _get_tags_by_tx_ac(self, tx_acs: list[str], genome_build: str) -> dict[str, list[str]]:
        """ Batch tag lookup for gene-symbol resolution.

            Start from cdot's JSON tags (the per-transcript default), then for any
            transcript that has none, supplement from VG's MANE table in a single
            query (older VG imports may have no JSON `tag` data). MANE matching is
            version-insensitive: VG stores transcript versions separately, so the
            incoming accessions are typically versionless (e.g. "NM_000059") and we
            match on transcript id, ignoring any ".N" version. """
        tags_by_ac = super()._get_tags_by_tx_ac(tx_acs, genome_build)

        missing = [ac for ac in tx_acs if not tags_by_ac.get(ac)]
        if missing:
            # Map versionless transcript id -> the accession the caller asked for
            versionless_to_ac = {ac.split(".", 1)[0]: ac for ac in missing}
            mane_qs = MANE.objects.filter(
                Q(refseq_transcript_version__transcript_id__in=versionless_to_ac) |
                Q(ensembl_transcript_version__transcript_id__in=versionless_to_ac)
            ).values_list("refseq_transcript_version__transcript_id",
                          "ensembl_transcript_version__transcript_id", "status")
            for refseq_tid, ensembl_tid, status in mane_qs:
                if tag := _MANE_STATUS_TO_CDOT_TAG.get(MANEStatus(status)):
                    for tid in (refseq_tid, ensembl_tid):
                        if ac := versionless_to_ac.get(tid):
                            if not tags_by_ac.get(ac):
                                tags_by_ac[ac] = [tag]
        return tags_by_ac

    def _get_gene(self, gene):
        # TODO: Do we need this? We could always store this from cdot in the DB
        raise NotImplementedError()

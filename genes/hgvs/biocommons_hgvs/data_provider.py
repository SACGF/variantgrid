from cdot.hgvs.dataproviders import LocalDataProvider, FastaSeqFetcher, ChainedSeqFetcher
from hgvs.exceptions import HGVSDataNotAvailableError

from genes.models import TranscriptVersion, TranscriptVersionSequenceInfo, NoTranscript
from genes.transcripts_utils import get_refseq_type
from snpdb.models import Contig

def is_contig(ac):
    if refseq_type := get_refseq_type(ac):
        return refseq_type == "genomic"
    return False


class DBTranscriptSeqFetcher:
    """ We store some transcripts in the database, use them if we have them   """
    source = "VG DB Transcript sequences (transcript fasta/API)"

    def fetch_seq(self, ac, start_i=None, end_i=None):
        if is_contig(ac):
            raise HGVSDataNotAvailableError("This SeqFetcher implementation doesn't handle contigs")

        tvsi = TranscriptVersionSequenceInfo.get(ac, retrieve=False)
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

    def __init__(self, genome_build):
        db_transcript_seqfetcher = DBTranscriptSeqFetcher()
        fasta_seqfetcher = SingleBuildFastaSeqFetcher(genome_build)
        seqfetcher = ChainedSeqFetcher(db_transcript_seqfetcher, fasta_seqfetcher)

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

    def _get_gene(self, gene):
        # TODO: Do we need this? We could always store this from cdot in the DB
        raise NotImplementedError()

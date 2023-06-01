from cdot.hgvs.dataproviders import LocalDataProvider, FastaSeqFetcher

from genes.models import TranscriptVersion
from genes.refseq_transcripts import get_refseq_type
from snpdb.models import Contig


class SingleBuildFastaSeqFetcher(FastaSeqFetcher):
    """ We want this to fail with ContigNotInBuildError if asked for a contig not in the
        the genome build provided """
    def __init__(self, genome_build):
        super().__init__(genome_build.genome_fasta.filename)
        self.genome_build = genome_build

    def fetch_seq(self, ac, start_i=None, end_i=None):
        # If FastaSeqFetcher is called with contigs from another build it will try and lookup transcripts
        # Just fail ASAP
        if refseq_type := get_refseq_type(ac):
            if refseq_type == "genomic":
                if ac not in self.contig_fastas:
                    raise Contig.ContigNotInBuildError(f"Contig '{ac}' not in genome build '{self.genome_build.name}'")
        return super().fetch_seq(ac, start_i=start_i, end_i=end_i)


class DjangoTranscriptDataProvider(LocalDataProvider):

    def __init__(self, genome_build):
        seqfetcher = SingleBuildFastaSeqFetcher(genome_build)
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

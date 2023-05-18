from cdot.hgvs.dataproviders import LocalDataProvider, FastaSeqFetcher

from genes.models import TranscriptVersion


class DjangoTranscriptDataProvider(LocalDataProvider):

    def __init__(self, genome_build):
        seqfetcher = FastaSeqFetcher(genome_build.genome_fasta.filename)
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

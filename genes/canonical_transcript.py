import abc

from genes.models import GeneSymbol, Gene
from genes.models_enums import MANEStatus, AnnotationConsortium
from snpdb.models import GenomeBuild


class CanonicalTranscriptManager(abc.ABC):
    def __init__(self, genome_build: GenomeBuild, annotation_consortium: AnnotationConsortium):
        self.genome_build = genome_build
        self.annotation_consortium = annotation_consortium

    @abc.abstractmethod
    def get_for_gene_symbol(self, gene_symbol: GeneSymbol):
        pass

    @abc.abstractmethod
    def get_for_gene(self, gene: Gene):
        pass


class MANECanonicalTranscriptManager(CanonicalTranscriptManager):
    def get_for_gene_symbol(self, gene_symbol: GeneSymbol):
        pass

    def get_for_gene(self, gene: Gene):
        pass


class RefSeqSelectCanonicalTranscriptManager(CanonicalTranscriptManager):
    def get_for_gene_symbol(self, gene_symbol: GeneSymbol):
        pass

    def get_for_gene(self, gene: Gene):
        pass


class VEPCanonicalTranscriptManager(CanonicalTranscriptManager):
    """ This is for when transcripts don't define a
        and should only ever be used as a fallback """
    def get_for_gene_symbol(self, gene_symbol: GeneSymbol):
        pass

    def get_for_gene(self, gene: Gene):
        vav = None
        gene.get_vep_canonical_transcript(variant_annotation_version=vav)


def canonical_transcript_manager_factory(genome_build: GenomeBuild, annotation_consortium: AnnotationConsortium):
    if genome_build.name == "GRCh37":
        # MANE
        pass
    elif genome_build.name == "GRCh38":
        # MANE
        pass


# TODO:
# * Have a way to do it via enrichment kit (using CanonicalTranscriptCollection)
# * Build something that takes a list of transcripts then gets the most canonical one (ie refseq select etc)

def get_canonical_transcript_for_gene_symbol(genome_build: GenomeBuild, annotation_consortium, gene_symbol: GeneSymbol):
    if genome_build.name != "GRCh38":
        raise ValueError("Only know how to handle GRCh38")

    transcript_version = None
    if mane := gene_symbol.mane_set.filter(status=MANEStatus.MANE_SELECT).first():
        transcript_version = mane.get_transcript_version(annotation_consortium)
    return transcript_version

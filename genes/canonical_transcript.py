from genes.models import GeneSymbol
from genes.models_enums import MANEStatus, AnnotationConsortium
from snpdb.models import GenomeBuild

# Perhaps make abstract classes that do it for us??


def get_canonical_transcript_for_gene_symbol(genome_build: GenomeBuild, annotation_consortium, gene_symbol: GeneSymbol):
    if genome_build.name != "GRCh38":
        raise ValueError("Only know how to handle GRCh38")

    transcript_version = None
    if mane := gene_symbol.mane_set.filter(status=MANEStatus.MANE_SELECT).first():
        transcript_version = mane.get_transcript_version(annotation_consortium)
    return transcript_version

from genes.models import GeneSymbol, Transcript, Gene, GeneAnnotationImport, GeneVersion, TranscriptVersion
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild


def create_fake_transcript_version(genome_build: GenomeBuild) -> TranscriptVersion:
    runx1_symbol = GeneSymbol.objects.get_or_create(symbol='RUNX1')[0]
    transcript, _ = Transcript.objects.get_or_create(identifier="ENST00000300305",
                                                     annotation_consortium=AnnotationConsortium.ENSEMBL)
    gene, _ = Gene.objects.get_or_create(identifier="ENSG00000159216",
                                         annotation_consortium=AnnotationConsortium.ENSEMBL)
    import_source, _ = GeneAnnotationImport.objects.get_or_create(filename="fake", genome_build=genome_build,
                                                                  annotation_consortium=AnnotationConsortium.ENSEMBL)
    gene_version, _ = GeneVersion.objects.get_or_create(gene=gene, gene_symbol=runx1_symbol, version=1,
                                                        genome_build=genome_build, import_source=import_source)
    transcript_version, _ = TranscriptVersion.objects.get_or_create(transcript=transcript, gene_version=gene_version,
                                                                    version=1, genome_build=genome_build,
                                                                    import_source=import_source)
    return transcript_version

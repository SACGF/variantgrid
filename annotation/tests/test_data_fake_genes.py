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

    data = {'id': 'ENST00000300305.7', 'end': 35049344, 'chrom': '21',
            'exons': [[34787800, 34792610], [34799300, 34799462], [34834409, 34834601], [34859473, 34859578],
                      [34880556, 34880713], [34886842, 34887096], [34892924, 34892963], [35048841, 35049344]],
            'start': 34787801, 'strand': '-', 'cds_end': 35048899, 'cds_start': 34792134}

    transcript_version, _ = TranscriptVersion.objects.get_or_create(transcript=transcript, gene_version=gene_version,
                                                                    version=7, genome_build=genome_build,
                                                                    import_source=import_source, data=data)
    return transcript_version

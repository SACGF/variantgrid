from genes.models import GeneSymbol, Transcript, Gene, GeneAnnotationImport, GeneVersion, TranscriptVersion, \
    GeneAnnotationRelease, ReleaseGeneSymbol, ReleaseGeneSymbolGene
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild


def create_fake_transcript_version(genome_build: GenomeBuild,
                                   release: GeneAnnotationRelease = None) -> TranscriptVersion:
    runx1_symbol = GeneSymbol.objects.get_or_create(symbol='RUNX1')[0]
    transcript, _ = Transcript.objects.get_or_create(identifier="ENST00000300305",
                                                     annotation_consortium=AnnotationConsortium.ENSEMBL)
    gene, _ = Gene.objects.get_or_create(identifier="ENSG00000159216",
                                         annotation_consortium=AnnotationConsortium.ENSEMBL)
    import_source, _ = GeneAnnotationImport.objects.get_or_create(url="fake", genome_build=genome_build,
                                                                  annotation_consortium=AnnotationConsortium.ENSEMBL)
    gene_version, _ = GeneVersion.objects.get_or_create(gene=gene, gene_symbol=runx1_symbol, version=1,
                                                        genome_build=genome_build, import_source=import_source)

    data = {
        "id": "ENST00000300305.7",
        "hgnc": "10471",
        "biotype": [],
        "gene_name": "RUNX1",
        "stop_codon": 1888,
        "start_codon": 445,
    }

    # This is from GRCh38 (doesn't exist for 37) but we'll use it for both
    build_data = {
        "url": "ftp://ftp.ensembl.org/pub/release-106/gff3/homo_sapiens/Homo_sapiens.GRCh38.106.gff3.gz",
        "exons": [[34787800, 34792610, 7, 1413, 6222, None], [34799300, 34799462, 6, 1251, 1412, None],
                  [34834409, 34834601, 5, 1059, 1250, None], [34859473, 34859578, 4, 954, 1058, None],
                  [34880556, 34880713, 3, 797, 953, None], [34886842, 34887096, 2, 543, 796, None],
                  [34892924, 34892963, 1, 504, 542, None], [35048841, 35049344, 0, 1, 503, None]],
        "contig": "21",
        "strand": "-",
        "cds_end": 35048899,
        "cds_start": 34792134,
    }

    data["genome_builds"] = {genome_build.name: build_data}
    contig = genome_build.chrom_contig_mappings[build_data["contig"]]
    tv_defaults = {
        "gene_version": gene_version,
        "contig": contig,
        "import_source": import_source,
        "data": data,
    }
    transcript_version, _ = TranscriptVersion.objects.update_or_create(transcript=transcript,
                                                                       version=7,
                                                                       genome_build=genome_build,
                                                                       defaults=tv_defaults)

    if release:
        release_symbol = ReleaseGeneSymbol.objects.create(release=release, gene_symbol=runx1_symbol)
        ReleaseGeneSymbolGene.objects.create(release_gene_symbol=release_symbol, gene=gene)

    return transcript_version

from genes.models import GeneSymbol, Transcript, Gene, GeneAnnotationImport, GeneVersion, TranscriptVersion, \
    GeneAnnotationRelease, ReleaseGeneSymbol, ReleaseGeneSymbolGene
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild


def _create_fake_gene_version(genome_build: GenomeBuild, gene_id, gene_symbol_str, annotation_consortium):
    gene_symbol = GeneSymbol.objects.get_or_create(symbol=gene_symbol_str)[0]
    gene, _ = Gene.objects.get_or_create(identifier=gene_id,
                                         annotation_consortium=annotation_consortium)
    import_source, _ = GeneAnnotationImport.objects.get_or_create(url="fake", genome_build=genome_build,
                                                                  annotation_consortium=AnnotationConsortium.ENSEMBL)
    gene_version, _ = GeneVersion.objects.get_or_create(gene=gene, gene_symbol=gene_symbol, version=1,
                                                        genome_build=genome_build, import_source=import_source)
    return gene_version


def _insert_transcript_data(genome_build, data: dict, gene_version: GeneVersion,
                            release: GeneAnnotationRelease = None):
    build_data = data["genome_builds"][genome_build.name]
    contig = genome_build.chrom_contig_mappings[build_data["contig"]]
    tv_defaults = {
        "gene_version": gene_version,
        "contig": contig,
        "import_source": gene_version.import_source,
        "data": data,
    }

    transcript_id, version = TranscriptVersion.get_transcript_id_and_version(data["id"])
    transcript, _ = Transcript.objects.get_or_create(identifier=transcript_id,
                                                     annotation_consortium=gene_version.gene.annotation_consortium)
    transcript_version, _ = TranscriptVersion.objects.update_or_create(transcript=transcript,
                                                                       version=version,
                                                                       genome_build=genome_build,
                                                                       defaults=tv_defaults)

    if release:
        release_symbol = ReleaseGeneSymbol.objects.create(release=release, gene_symbol=gene_version.gene_symbol)
        ReleaseGeneSymbolGene.objects.create(release_gene_symbol=release_symbol, gene=gene_version.gene)

    return transcript_version


def create_fake_transcript_version(genome_build: GenomeBuild,
                                   release: GeneAnnotationRelease = None) -> TranscriptVersion:

    gene_version = _create_fake_gene_version(genome_build, "ENSG00000159216", "RUNX1", AnnotationConsortium.ENSEMBL)

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
    return _insert_transcript_data(genome_build, data, gene_version, release)


def create_gata2_transcript_version(genome_build) -> TranscriptVersion:

    gene_version = _create_fake_gene_version(genome_build, "2624", "GATA2", AnnotationConsortium.REFSEQ)
    nm_001145661_2 = {"id": "NM_001145661.2",
                      "cdot": "0.2.17",
                      "hgnc": "4171",
                      "biotype": ["protein_coding"],
                      "gene_name": "GATA2",
                      "stop_codon": 1878,
                      "start_codon": 435,
                      "genome_builds": {
                          "GRCh37": {
                              "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
                              "exons": [[128198269, 128200161, 6, 1579, 3470, None],
                                        [128200661, 128200787, 5, 1453, 1578, None],
                                        [128202702, 128202848, 4, 1307, 1452, None],
                                        [128204569, 128205211, 3, 665, 1306, None],
                                        [128205645, 128205919, 2, 391, 664, None],
                                        [128206553, 128206764, 1, 180, 390, None],
                                        [128207194, 128207373, 0, 1, 179, None]],
                              "contig": "NC_000003.11",
                              "strand": "-",
                              "cds_end": 128205874,
                              "cds_start": 128199861},
                          "GRCh38": {
                              "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.40-RS_2023_03/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                              "exons": [[128479426, 128481318, 6, 1579, 3470, None],
                                        [128481818, 128481944, 5, 1453, 1578, None],
                                        [128483859, 128484005, 4, 1307, 1452, None],
                                        [128485726, 128486368, 3, 665, 1306, None],
                                        [128486802, 128487076, 2, 391, 664, None],
                                        [128487710, 128487921, 1, 180, 390, None],
                                        [128488351, 128488530, 0, 1, 179, None]],
                              "contig": "NC_000003.12",
                              "strand": "-",
                              "cds_end": 128487031,
                              "cds_start": 128481018}}}

    return _insert_transcript_data(genome_build, nm_001145661_2, gene_version)

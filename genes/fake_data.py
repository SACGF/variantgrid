"""
Fake genes and transcripts: RUNX1 (Ensembl, minus strand), GATA2 (RefSeq, minus) and PTEN (RefSeq, plus), each with
real exon data for GRCh37 and GRCh38 - the transcripts tests build on, and the 'create_fake_data genes' step, which
makes them for a build that has no real gene annotation to place fake variants in.

The 'create_fake_data gene_level' step makes the gene-level events the search page gives as examples (fusion, copy
number, splice), resolved from their written strings as a classification import would.
"""
from genes.gene_copy_number import create_gene_copy_number_events_for_variants
from genes.gene_fusions import create_gene_fusions_for_variants
from genes.gene_level_strings import resolve_gene_level_string
from genes.models import (
    Gene,
    GeneAnnotationImport,
    GeneAnnotationRelease,
    GeneSymbol,
    GeneVersion,
    ReleaseGeneSymbol,
    ReleaseGeneSymbolGene,
    ReleaseGeneVersion,
    Transcript,
    TranscriptVersion,
)
from genes.models_enums import AnnotationConsortium
from library.fake_data import FakeData, FakeDataContext, register
from library.utils import sha256sum_str
from snpdb.gene_level_variants import GENE_LEVEL_SVLEN
from snpdb.models import Contig, GenomeBuild, Locus, Sequence, Variant, VariantCoordinate

# The examples on the search page's gene-level receivers (snpdb/signals/variant_search.py)
GENE_LEVEL_EXAMPLES = ["BCR::ABL1", "CD74--ROS1", "EGFR amplification", "PTEN loss", "AR V7", "MET exon 14 skipping"]


def create_fake_gene_version(genome_build: GenomeBuild, gene_id, gene_symbol_str, annotation_consortium):
    gene_symbol = GeneSymbol.objects.get_or_create(symbol=gene_symbol_str)[0]
    gene, _ = Gene.objects.get_or_create(identifier=gene_id,
                                         annotation_consortium=annotation_consortium)
    import_source, _ = GeneAnnotationImport.objects.get_or_create(url="fake", genome_build=genome_build,
                                                                  annotation_consortium=AnnotationConsortium.ENSEMBL)
    gene_version, _ = GeneVersion.objects.get_or_create(gene=gene, gene_symbol=gene_symbol, version=1,
                                                        genome_build=genome_build, import_source=import_source)
    return gene_version


def insert_transcript_data(genome_build, data: dict, gene_version: GeneVersion,
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
        # What import_cdot_gene_annotation_release writes - which gene versions the release holds
        ReleaseGeneVersion.objects.get_or_create(release=release, gene_version=gene_version)

    return transcript_version


def create_fake_transcript_version(genome_build: GenomeBuild,
                                   release: GeneAnnotationRelease = None) -> TranscriptVersion:

    gene_version = create_fake_gene_version(genome_build, "ENSG00000159216", "RUNX1", AnnotationConsortium.ENSEMBL)

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
    return insert_transcript_data(genome_build, data, gene_version, release)


def create_gata2_transcript_version(genome_build) -> TranscriptVersion:

    gene_version = create_fake_gene_version(genome_build, "2624", "GATA2", AnnotationConsortium.REFSEQ)
    nm_001145661_2 = {"id": "NM_001145661.2",
                      "hgnc": "4171",
                      "biotype": ["protein_coding"],
                      "gene_name": "GATA2",
                      "protein": "NP_001139133.1",
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

    return insert_transcript_data(genome_build, nm_001145661_2, gene_version)


def create_col18a1_transcript_version(genome_build) -> TranscriptVersion:
    """ Plus strand, with a 9 base alignment gap in exon 33 (the RefSeq mRNA lacks 9 genomic bases),
        so c. positions after it are 9 lower than gap-unaware numbering gives """

    gene_version = create_fake_gene_version(genome_build, "80781", "COL18A1", AnnotationConsortium.REFSEQ)
    nm_130444_2 = {"id": "NM_130444.2",
                   "hgnc": "2195",
                   "biotype": ["mRNA"],
                   "gene_name": "COL18A1",
                   "stop_codon": 5277,
                   "start_codon": 21,
                   "genome_builds": {
                       "GRCh37": {
                           "url": "postgresql://uta.biocommons.org/uta_20241220",
                           "exons": [[46875423, 46876795, 0, 1, 1372, None],
                                     [46888155, 46888700, 1, 1373, 1917, None],
                                     [46893808, 46893895, 2, 1918, 2004, None],
                                     [46895389, 46895449, 3, 2005, 2064, None],
                                     [46896264, 46896394, 4, 2065, 2194, None],
                                     [46897324, 46897401, 5, 2195, 2271, None],
                                     [46897663, 46897879, 6, 2272, 2487, None],
                                     [46898240, 46898267, 7, 2488, 2514, None],
                                     [46899815, 46899878, 8, 2515, 2577, None],
                                     [46899983, 46900070, 9, 2578, 2664, None],
                                     [46900380, 46900434, 10, 2665, 2718, None],
                                     [46900613, 46900772, 11, 2719, 2877, None],
                                     [46901876, 46901939, 12, 2878, 2940, None],
                                     [46902708, 46902735, 13, 2941, 2967, None],
                                     [46906774, 46906906, 14, 2968, 3099, None],
                                     [46907360, 46907423, 15, 3100, 3162, None],
                                     [46908331, 46908358, 16, 3163, 3189, None],
                                     [46909399, 46909435, 17, 3190, 3225, None],
                                     [46910188, 46910260, 18, 3226, 3297, None],
                                     [46910749, 46910785, 19, 3298, 3333, None],
                                     [46911138, 46911228, 20, 3334, 3423, None],
                                     [46912448, 46912478, 21, 3424, 3453, None],
                                     [46912600, 46912627, 22, 3454, 3480, None],
                                     [46913076, 46913139, 23, 3481, 3543, None],
                                     [46913414, 46913489, 24, 3544, 3618, None],
                                     [46914458, 46914485, 25, 3619, 3645, None],
                                     [46914775, 46914829, 26, 3646, 3699, None],
                                     [46915271, 46915346, 27, 3700, 3774, None],
                                     [46916413, 46916482, 28, 3775, 3843, None],
                                     [46916963, 46917006, 29, 3844, 3886, None],
                                     [46917512, 46917575, 30, 3887, 3949, None],
                                     [46923924, 46923968, 31, 3950, 3993, None],
                                     [46924329, 46924470, 32, 3994, 4125, "M96 D9 M36"],
                                     [46925047, 46925192, 33, 4126, 4270, None],
                                     [46925271, 46925345, 34, 4271, 4344, None],
                                     [46925751, 46925880, 35, 4345, 4473, None],
                                     [46927474, 46927507, 36, 4474, 4506, None],
                                     [46929269, 46929515, 37, 4507, 4752, None],
                                     [46929977, 46930175, 38, 4753, 4950, None],
                                     [46931024, 46931140, 39, 4951, 5066, None],
                                     [46932101, 46933634, 40, 5067, 6599, None]],
                           "contig": "NC_000021.8",
                           "strand": "+",
                           "cds_end": 46932312,
                           "cds_start": 46875444},
                       "GRCh38": {
                           "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/109/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
                           "exons": [[45455509, 45456881, 0, 1, 1372, None],
                                     [45468241, 45468786, 1, 1373, 1917, None],
                                     [45473894, 45473981, 2, 1918, 2004, None],
                                     [45475475, 45475535, 3, 2005, 2064, None],
                                     [45476350, 45476480, 4, 2065, 2194, None],
                                     [45477410, 45477487, 5, 2195, 2271, None],
                                     [45477749, 45477965, 6, 2272, 2487, None],
                                     [45478326, 45478353, 7, 2488, 2514, None],
                                     [45479901, 45479964, 8, 2515, 2577, None],
                                     [45480069, 45480156, 9, 2578, 2664, None],
                                     [45480466, 45480520, 10, 2665, 2718, None],
                                     [45480699, 45480858, 11, 2719, 2877, None],
                                     [45481962, 45482025, 12, 2878, 2940, None],
                                     [45482794, 45482821, 13, 2941, 2967, None],
                                     [45486860, 45486992, 14, 2968, 3099, None],
                                     [45487446, 45487509, 15, 3100, 3162, None],
                                     [45488417, 45488444, 16, 3163, 3189, None],
                                     [45489485, 45489521, 17, 3190, 3225, None],
                                     [45490274, 45490346, 18, 3226, 3297, None],
                                     [45490835, 45490871, 19, 3298, 3333, None],
                                     [45491224, 45491314, 20, 3334, 3423, None],
                                     [45492534, 45492564, 21, 3424, 3453, None],
                                     [45492686, 45492713, 22, 3454, 3480, None],
                                     [45493162, 45493225, 23, 3481, 3543, None],
                                     [45493500, 45493575, 24, 3544, 3618, None],
                                     [45494544, 45494571, 25, 3619, 3645, None],
                                     [45494861, 45494915, 26, 3646, 3699, None],
                                     [45495357, 45495432, 27, 3700, 3774, None],
                                     [45496499, 45496568, 28, 3775, 3843, None],
                                     [45497049, 45497092, 29, 3844, 3886, None],
                                     [45497598, 45497661, 30, 3887, 3949, None],
                                     [45504010, 45504054, 31, 3950, 3993, None],
                                     [45504415, 45504556, 32, 3994, 4125, "M110 D9 M22"],
                                     [45505133, 45505278, 33, 4126, 4270, None],
                                     [45505357, 45505431, 34, 4271, 4344, None],
                                     [45505837, 45505966, 35, 4345, 4473, None],
                                     [45507560, 45507593, 36, 4474, 4506, None],
                                     [45509355, 45509601, 37, 4507, 4752, None],
                                     [45510063, 45510261, 38, 4753, 4950, None],
                                     [45511110, 45511226, 39, 4951, 5066, None],
                                     [45512187, 45513720, 40, 5067, 6599, None]],
                           "contig": "NC_000021.9",
                           "strand": "+",
                           "cds_end": 45512398,
                           "cds_start": 45455530}}}

    return insert_transcript_data(genome_build, nm_130444_2, gene_version)


def create_gata2_as1_transcript_version(genome_build) -> TranscriptVersion:
    """ Non-coding (lnc_RNA) transcript, GRCh37 only """

    gene_version = create_fake_gene_version(genome_build, "101927167", "GATA2-AS1", AnnotationConsortium.REFSEQ)
    nr_125398_1 = {"id": "NR_125398.1",
                   "hgnc": "51108",
                   "biotype": ["lnc_RNA", "ncRNA"],
                   "gene_name": "GATA2-AS1",
                   "genome_builds": {
                       "GRCh37": {
                           "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.25-RS_2024_09/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
                           "exons": [[128208045, 128208513, 0, 1, 468, None],
                                     [128208621, 128208912, 1, 469, 759, None],
                                     [128220451, 128222050, 2, 760, 2358, None]],
                           "contig": "NC_000003.11",
                           "strand": "+"}}}

    return insert_transcript_data(genome_build, nr_125398_1, gene_version)


def create_pten_transcript_version(genome_build) -> TranscriptVersion:
    """ Plus strand - GATA2 above is minus, so the pair covers both orientations """

    gene_version = create_fake_gene_version(genome_build, "5728", "PTEN", AnnotationConsortium.REFSEQ)
    nm_000314_8 = {"id": "NM_000314.8",
                   "hgnc": "9588",
                   "biotype": ["protein_coding"],
                   "gene_name": "PTEN",
                   "stop_codon": 2057,
                   "start_codon": 845,
                   "genome_builds": {
                       "GRCh37": {
                           "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.25-RS_2024_09/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
                           "exons": [[89623381, 89624305, 0, 1, 924, None],
                                     [89653781, 89653866, 1, 925, 1009, None],
                                     [89685269, 89685314, 2, 1010, 1054, None],
                                     [89690802, 89690846, 3, 1055, 1098, None],
                                     [89692769, 89693008, 4, 1099, 1337, None],
                                     [89711874, 89712016, 5, 1338, 1479, None],
                                     [89717609, 89717776, 6, 1480, 1646, None],
                                     [89720650, 89720875, 7, 1647, 1871, None],
                                     [89725043, 89731687, 8, 1872, 8515, None]],
                           "contig": "NC_000010.10",
                           "strand": "+",
                           "cds_end": 89725229,
                           "cds_start": 89624226},
                       "GRCh38": {
                           "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                           "exons": [[87863624, 87864548, 0, 1, 924, None],
                                     [87894024, 87894109, 1, 925, 1009, None],
                                     [87925512, 87925557, 2, 1010, 1054, None],
                                     [87931045, 87931089, 3, 1055, 1098, None],
                                     [87933012, 87933251, 4, 1099, 1337, None],
                                     [87952117, 87952259, 5, 1338, 1479, None],
                                     [87957852, 87958019, 6, 1480, 1646, None],
                                     [87960893, 87961118, 7, 1647, 1871, None],
                                     [87965286, 87971930, 8, 1872, 8515, None]],
                           "contig": "NC_000010.11",
                           "strand": "+",
                           "cds_end": 87965472,
                           "cds_start": 87864469}}}

    return insert_transcript_data(genome_build, nm_000314_8, gene_version)


# When the build has real gene annotation, fake variants go in real genes that tagging and classification
# are busiest in (the steps that use them each prefer their own subset)
WELL_KNOWN_GENES = [
    "TP53", "KRAS", "NRAS", "BRAF", "EGFR", "PIK3CA", "KIT", "JAK2", "DNMT3A", "TET2", "ASXL1", "RUNX1", "IDH1",
    "IDH2", "FLT3", "NPM1", "SF3B1", "CALR", "PTEN", "SMAD4", "BRCA1", "BRCA2", "ATM", "MSH2", "MLH1", "MSH6",
    "PMS2", "PALB2", "CHEK2", "APC", "CFTR", "TTN", "MYBPC3", "RYR1", "DMD", "FBN1", "NF1", "LDLR", "SCN1A",
    "COL4A5", "STK11", "CDH1", "RB1", "VHL", "RET", "SDHB", "MUTYH", "BMPR1A",
]


def has_real_gene_annotation(genome_build: GenomeBuild) -> bool:
    return GeneAnnotationRelease.objects.filter(genome_build=genome_build) \
        .exclude(gene_annotation_import__url="fake").exists()


@register
class FakeGenes(FakeData):
    name = "genes"
    help = ("RUNX1, GATA2 and PTEN transcripts, when the build has no real gene annotation - "
            "otherwise later steps use real, well known genes")

    def create(self, context: FakeDataContext, **options):
        genome_build = context.genome_build
        if has_real_gene_annotation(genome_build):
            context.genes = WELL_KNOWN_GENES
            context.stdout.write(f"{genome_build} has real gene annotation - using {len(context.genes)} real genes")
            return

        if genome_build.name not in ("GRCh37", "GRCh38"):
            raise ValueError(f"The fake transcripts only have exons for GRCh37 and GRCh38, not {genome_build}")
        transcript_versions = [create_fake_transcript_version(genome_build),
                               create_gata2_transcript_version(genome_build),
                               create_pten_transcript_version(genome_build)]
        context.genes = [tv.gene_version.gene_symbol_id for tv in transcript_versions]
        context.stdout.write(f"Fake transcripts for {', '.join(context.genes)}")

    def delete(self, context: FakeDataContext, **options):
        context.stdout.write("Leaving the fake transcripts - tests and the fake annotation version build on them")


def get_sequence(seq: str) -> Sequence:
    """ Upper-cased, as every path that inserts a Sequence in production does """
    seq = seq.upper()
    sequence, _ = Sequence.objects.get_or_create(seq=seq, defaults={"seq_sha256_hash": sha256sum_str(seq)})
    return sequence


def create_gene_level_variant(variant_coordinate: VariantCoordinate) -> Variant:
    """ The Variant the insert pipeline would have created for a resolved gene-level identity. Production has one
        way in - the VCF insert pipeline (@see snpdb.gene_level_variants) - which is far more machinery than one
        event needs; this makes the same rows from the same resolved identity """
    locus, _ = Locus.objects.get_or_create(contig=Contig.get_gene_level(),
                                           position=variant_coordinate.position,
                                           ref=get_sequence(variant_coordinate.ref))
    variant, _ = Variant.objects.get_or_create(locus=locus, alt=get_sequence(variant_coordinate.alt),
                                               svlen=GENE_LEVEL_SVLEN,
                                               defaults={"end": variant_coordinate.position})
    return variant


def _resolve_gene_level_examples(context: FakeDataContext) -> dict[str, VariantCoordinate]:
    """ An example whose genes this database doesn't have (only the fake genes, say) is skipped """
    variant_coordinates = {}
    for example in GENE_LEVEL_EXAMPLES:
        resolution = resolve_gene_level_string(example, context.genome_build)
        if resolution:
            variant_coordinates[example] = resolution.resolved.variant_coordinate
        else:
            context.stdout.write(f"Skipping '{example}': {resolution.reason}")
    return variant_coordinates


@register
class FakeGeneLevelVariants(FakeData):
    name = "gene_level"
    help = ("The gene fusions, copy number and splice events the search page gives as examples, "
            "resolved from their strings as a classification import would")
    requires = ("genes",)

    def create(self, context: FakeDataContext, **options):
        variant_ids = []
        for example, variant_coordinate in _resolve_gene_level_examples(context).items():
            variant = create_gene_level_variant(variant_coordinate)
            variant_ids.append(variant.pk)
            context.stdout.write(f"{example} -> {variant}")
        variant_qs = Variant.objects.filter(pk__in=variant_ids)
        create_gene_fusions_for_variants(variant_qs)
        create_gene_copy_number_events_for_variants(variant_qs)

    def delete(self, context: FakeDataContext, **options):
        for example, variant_coordinate in _resolve_gene_level_examples(context).items():
            deleted, _ = Variant.objects.filter(locus__contig=Contig.get_gene_level(),
                                                locus__position=variant_coordinate.position,
                                                alt__seq=variant_coordinate.alt.upper()).delete()
            context.stdout.write(f"{example}: deleted {deleted} rows")

"""
Helpers shared by the 'manage.py create_fake_data' subcommands, which live in the app whose data they
create (@see snpdb.management.commands.create_fake_data).
"""
from annotation.models import VariantAnnotation, VariantAnnotationVersion
from genes.models import TranscriptVersion
from snpdb.models import GenomeBuild


def get_variant_ids_by_gene(genome_build: GenomeBuild, genes: list[str],
                            without_alleles: bool = False) -> dict[str, list[int]]:
    """ without_alleles restricts to variants no allele has claimed yet, so fake data can make its own
        alleles and take them away again without touching anything real """
    gene_symbol_field = "transcript_version__gene_version__gene_symbol_id"
    transcript_versions_qs = TranscriptVersion.objects.filter(genome_build=genome_build,
                                                              gene_version__gene_symbol__in=genes)
    variant_annotation_qs = VariantAnnotation.objects.filter(version=VariantAnnotationVersion.latest(genome_build),
                                                             transcript_version__in=transcript_versions_qs)
    if without_alleles:
        variant_annotation_qs = variant_annotation_qs.filter(variant__variantallele__isnull=True)
    variant_ids_by_gene = {}
    for gene_symbol, variant_id in variant_annotation_qs.values_list(gene_symbol_field, "variant_id"):
        variant_ids_by_gene.setdefault(gene_symbol, []).append(variant_id)
    return variant_ids_by_gene


def zipf_weight(rank: int) -> float:
    """ Long tail - a few genes/records get most of the action, the rest get a little """
    return 1 / (rank + 1) ** 0.8

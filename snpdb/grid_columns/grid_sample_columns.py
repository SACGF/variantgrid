"""
Helper functions for grid to join to cohort genotype

This is in snpdb not analysis as it needs to be called by snpdb.tasks.cohort_genotype_tasks
and I don't want cross-dependencies between those apps
"""
from collections.abc import Iterable

from django.db.models import Value
from django.db.models.functions import Coalesce

from snpdb.models import VCF, Cohort, CohortGenotype

# The sample's whole call is drawn in the one zygosity cell - glyph, frequency, depths, then GQ/PL
# as quality marks and the sample filters. The rest ride along hidden and stay in the CSV.
# @see VariantGridFormat.sampleZygosity
SAMPLE_COMPOSITE_COLUMNS = ['samples_allele_depth', 'samples_allele_frequency', 'samples_read_depth',
                            'samples_genotype_quality', 'samples_phred_likelihood', 'samples_filters']
# What the zygosity header's sort menu offers - whichever of them the sample's VCF actually has
SAMPLE_SORT_KEY_LABELS = {
    'samples_zygosity': 'Zygosity',
    'samples_allele_frequency': 'Allele frequency',
    'samples_allele_depth': 'Allele depth',
    'samples_read_depth': 'Read depth',
    'samples_genotype_quality': 'Genotype quality',
    'samples_phred_likelihood': 'Phred likelihood',
    'samples_filters': 'Filters',
}


def get_available_format_columns(cohorts):
    available_format_columns = {
        # We want to always show some fields
        "samples_zygosity": True,
        "samples_allele_depth": True,
        "samples_allele_frequency": True,
        "samples_read_depth": True,
        "samples_genotype_quality": False,
        "samples_phred_likelihood": False,
        "samples_filters": False,
    }

    vcf_qs = VCF.objects.filter(sample__cohortsample__cohort__in=cohorts).distinct()
    for gq, pl, ft in vcf_qs.values_list("genotype_quality_field", "phred_likelihood_field", "sample_filters_field"):
        if gq:
            available_format_columns["samples_genotype_quality"] = True
        if pl:
            available_format_columns["samples_phred_likelihood"] = True
        if ft:
            available_format_columns["samples_filters"] = True
    return available_format_columns


def get_variantgrid_zygosity_annotation_kwargs(cohorts: Iterable[Cohort], common_variants: bool,
                                               annotation_gnomad_version=None):
    available_format_columns = get_available_format_columns(cohorts)
    annotation_kwargs = {}

    cgc_fields = {f.name: f for f in CohortGenotype._meta.fields}
    for cohort in cohorts:
        # How did this ever work with multiple cohorts - did it overwrite??
        # TODO: After we've done this - try and optimise to only doing rare if we can
        cgc = cohort.cohort_genotype_collection
        annotation_kwargs.update(cgc.get_annotation_kwargs(common_variants=common_variants,
                                                           annotation_gnomad_version=annotation_gnomad_version))

        for column, (is_array, empty_value) in CohortGenotype.COLUMN_IS_ARRAY_EMPTY_VALUE.items():
            if not available_format_columns[column]:
                continue  # No data, so don't show

            if is_array:
                empty_data = [empty_value] * cohort.sample_count
            else:
                empty_data = empty_value * cohort.sample_count

            output_field = cgc_fields[column]
            packed_column = cgc.get_packed_column_alias(column)
            alias = cgc.cohortgenotype_alias
            annotation_kwargs[packed_column] = Coalesce(f"{alias}__{column}", Value(empty_data),
                                                        output_field=output_field)

    return annotation_kwargs

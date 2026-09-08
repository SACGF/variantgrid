"""
Helper functions for grid to join to cohort genotype

This is in snpdb not analysis as it needs to be called by snpdb.tasks.cohort_genotype_tasks
and I don't want cross-dependencies between those apps
"""
from collections.abc import Iterable
from typing import Optional

from django.db.models import TextField, Value
from django.db.models.fields.json import KeyTextTransform
from django.db.models.functions import Coalesce

from snpdb.models import VCF, Cohort, CohortGenotype, CohortGenotypeCollection, Sample

# The caller's copy number / copy ratio. Unlike its neighbours it has no packed CohortGenotype column -
# which key means it is per VCF (@see VCF.copy_number_field), so it is read out of the stored JSON
COPY_NUMBER_COLUMN = 'samples_copy_number'
# The sample's whole call is drawn in the one zygosity cell - glyph, frequency, depths, then GQ/PL
# as quality marks, the copy number and the sample filters. The rest ride along hidden and stay in
# the CSV. @see VariantGridFormat.sampleZygosity
SAMPLE_COMPOSITE_COLUMNS = ['samples_allele_depth', 'samples_allele_frequency', 'samples_read_depth',
                            'samples_genotype_quality', 'samples_phred_likelihood', 'samples_filters',
                            COPY_NUMBER_COLUMN]
# What the zygosity header's sort menu offers - whichever of them the sample's VCF actually has
SAMPLE_SORT_KEY_LABELS = {
    'samples_zygosity': 'Zygosity',
    'samples_allele_frequency': 'Allele frequency',
    'samples_allele_depth': 'Allele depth',
    'samples_read_depth': 'Read depth',
    'samples_genotype_quality': 'Genotype quality',
    'samples_phred_likelihood': 'Phred likelihood',
    'samples_filters': 'Filters',
    COPY_NUMBER_COLUMN: 'Copy number',
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
        COPY_NUMBER_COLUMN: False,
    }

    vcf_qs = VCF.objects.filter(sample__cohortsample__cohort__in=cohorts).distinct()
    for gq, pl, ft, cn in vcf_qs.values_list("genotype_quality_field", "phred_likelihood_field",
                                             "sample_filters_field", "copy_number_field"):
        if gq:
            available_format_columns["samples_genotype_quality"] = True
        if pl:
            available_format_columns["samples_phred_likelihood"] = True
        if ft:
            available_format_columns["samples_filters"] = True
        if cn:
            available_format_columns[COPY_NUMBER_COLUMN] = True
    return available_format_columns


def get_copy_number_alias(cgc: CohortGenotypeCollection, sample_id: int) -> str:
    """ One alias per sample, not per cohort - there is no packed column to index into """
    return f"{cgc.cohortgenotype_alias}_{COPY_NUMBER_COLUMN}_{sample_id}"


def get_copy_number_annotation(cgc: CohortGenotypeCollection, sample: Sample) -> Optional[Coalesce]:
    """ The caller's copy number or copy ratio, read out of the stored CohortGenotype JSON: this
        sample's dict in the per-sample FORMAT list then the field's one-element array, falling back
        to INFO for the single-sample VCFs that put it there. None when the VCF declares no such
        field. @see VCF.copy_number_field """
    field = sample.vcf.copy_number_field
    if not field:
        return None
    alias = cgc.cohortgenotype_alias
    sample_index = cgc.get_array_index_for_sample_id(sample.pk)
    return Coalesce(KeyTextTransform("0", f"{alias}__format__{sample_index}__{field}"),
                    KeyTextTransform(field, f"{alias}__info"),
                    output_field=TextField())


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

        if available_format_columns[COPY_NUMBER_COLUMN]:
            for sample in cohort.get_samples():
                if (copy_number := get_copy_number_annotation(cgc, sample)) is not None:
                    annotation_kwargs[get_copy_number_alias(cgc, sample.pk)] = copy_number

    return annotation_kwargs

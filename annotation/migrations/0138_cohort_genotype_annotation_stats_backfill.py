"""
    Backfill CohortGenotypeVariantAnnotationStats /
    CohortGenotypeGeneAnnotationStats / CohortGenotypeClinVarAnnotationStats
    from the legacy Sample*AnnotationStats(+PassingFilter) tables. Per-sample
    rows only; aggregate rows (sample IS NULL) are populated lazily.
"""
from django.db import migrations


VARIANT_ANN_FIELDS = (
    "variant_dbsnp_count", "snp_dbsnp_count", "insertions_dbsnp_count", "deletions_dbsnp_count",
    "ref_high_or_moderate_count", "het_high_or_moderate_count",
    "hom_high_or_moderate_count", "unk_high_or_moderate_count",
)
GENE_ANN_FIELDS = (
    "gene_count",
    "ref_omim_phenotype_count", "het_omim_phenotype_count",
    "hom_omim_phenotype_count", "unk_omim_phenotype_count",
)
CLINVAR_ANN_FIELDS = (
    "clinvar_count",
    "ref_clinvar_pathogenic_count", "het_clinvar_pathogenic_count",
    "hom_clinvar_pathogenic_count", "unk_clinvar_pathogenic_count",
)


def _backfill(apps, schema_editor):
    SVAS = apps.get_model("annotation", "SampleVariantAnnotationStats")
    SVASPF = apps.get_model("annotation", "SampleVariantAnnotationStatsPassingFilter")
    SGAS = apps.get_model("annotation", "SampleGeneAnnotationStats")
    SGASPF = apps.get_model("annotation", "SampleGeneAnnotationStatsPassingFilter")
    SCVAS = apps.get_model("annotation", "SampleClinVarAnnotationStats")
    SCVASPF = apps.get_model("annotation", "SampleClinVarAnnotationStatsPassingFilter")
    NewVar = apps.get_model("annotation", "CohortGenotypeVariantAnnotationStats")
    NewGene = apps.get_model("annotation", "CohortGenotypeGeneAnnotationStats")
    NewCV = apps.get_model("annotation", "CohortGenotypeClinVarAnnotationStats")
    CohortGenotypeCollection = apps.get_model("snpdb", "CohortGenotypeCollection")

    def _cgc_for_sample(ss):
        cohort = ss.sample.vcf.cohort
        try:
            return CohortGenotypeCollection.objects.get(
                cohort=cohort, cohort_version=cohort.version, collection_type="U")
        except CohortGenotypeCollection.DoesNotExist:
            return None

    def _copy(old_qs, NewModel, version_field, version_attr, fields, passing_filter):
        rows = []
        for ss in old_qs.select_related(
            "sample__vcf__cohort", "code_version", version_attr,
        ).iterator():
            cgc = _cgc_for_sample(ss)
            if cgc is None:
                continue
            kwargs = {
                "cohort_genotype_collection": cgc,
                "sample": ss.sample,
                "filter_key": None,
                "passing_filter": passing_filter,
                "code_version": ss.code_version,
                version_field: getattr(ss, version_attr),
            }
            row = NewModel(**kwargs)
            for f in fields:
                setattr(row, f, getattr(ss, f, 0))
            rows.append(row)
            if len(rows) >= 2000:
                NewModel.objects.bulk_create(rows, batch_size=2000)
                rows = []
        if rows:
            NewModel.objects.bulk_create(rows, batch_size=2000)

    _copy(SVAS.objects.all(), NewVar, "variant_annotation_version", "variant_annotation_version",
          VARIANT_ANN_FIELDS, passing_filter=False)
    _copy(SVASPF.objects.all(), NewVar, "variant_annotation_version", "variant_annotation_version",
          VARIANT_ANN_FIELDS, passing_filter=True)
    _copy(SGAS.objects.all(), NewGene, "gene_annotation_version", "gene_annotation_version",
          GENE_ANN_FIELDS, passing_filter=False)
    _copy(SGASPF.objects.all(), NewGene, "gene_annotation_version", "gene_annotation_version",
          GENE_ANN_FIELDS, passing_filter=True)
    _copy(SCVAS.objects.all(), NewCV, "clinvar_version", "clinvar_version",
          CLINVAR_ANN_FIELDS, passing_filter=False)
    _copy(SCVASPF.objects.all(), NewCV, "clinvar_version", "clinvar_version",
          CLINVAR_ANN_FIELDS, passing_filter=True)


def _noop(apps, schema_editor):
    pass


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0137_cohortgenotypeclinvarannotationstats_and_more"),
        ("snpdb", "0179_cohort_genotype_stats_backfill"),
    ]

    operations = [
        migrations.RunPython(_backfill, _noop),
    ]

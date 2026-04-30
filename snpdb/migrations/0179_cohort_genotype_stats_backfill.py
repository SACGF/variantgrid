"""
    Backfill CohortGenotypeStats from existing SampleStats / SampleStatsPassingFilter
    rows. Per-sample rows only — aggregate rows (sample IS NULL) are populated
    lazily by calculate_cohort_stats when a CohortNode/TrioNode/PedigreeNode
    asks for them.

    Each old row is copied to a new CohortGenotypeStats row keyed by
        (cohort_genotype_collection=sample.vcf.cohort.cohort_genotype_collection,
         sample=sample,
         filter_key=NULL,
         passing_filter=False/True).
"""
from django.db import migrations


GENOTYPE_FIELDS = (
    "variant_count", "snp_count", "insertions_count", "deletions_count",
    "ref_count", "het_count", "hom_count", "unk_count",
    "x_hom_count", "x_het_count", "x_unk_count",
)


def _backfill(apps, schema_editor):
    SampleStats = apps.get_model("snpdb", "SampleStats")
    SampleStatsPassingFilter = apps.get_model("snpdb", "SampleStatsPassingFilter")
    CohortGenotypeStats = apps.get_model("snpdb", "CohortGenotypeStats")
    CohortGenotypeCollection = apps.get_model("snpdb", "CohortGenotypeCollection")

    def _cgc_for_sample(ss):
        # The Cohort.cohort_genotype_collection cached_property is a real
        # method on the runtime model and isn't available via apps.get_model,
        # so we re-implement the lookup using the underlying fields.
        cohort = ss.sample.vcf.cohort
        try:
            return CohortGenotypeCollection.objects.get(
                cohort=cohort, cohort_version=cohort.version, collection_type="U")
        except CohortGenotypeCollection.DoesNotExist:
            return None

    def _copy(old_qs, passing_filter):
        rows = []
        for ss in old_qs.select_related("sample__vcf__cohort", "code_version").iterator():
            cgc = _cgc_for_sample(ss)
            if cgc is None:
                continue  # Sample's cohort has no current CGC — skip; recompute will recreate.
            row = CohortGenotypeStats(
                cohort_genotype_collection=cgc,
                sample=ss.sample,
                filter_key=None,
                passing_filter=passing_filter,
                code_version=ss.code_version,
                import_status=ss.import_status,
            )
            for f in GENOTYPE_FIELDS:
                setattr(row, f, getattr(ss, f))
            rows.append(row)
            if len(rows) >= 2000:
                CohortGenotypeStats.objects.bulk_create(rows, batch_size=2000)
                rows = []
        if rows:
            CohortGenotypeStats.objects.bulk_create(rows, batch_size=2000)

    _copy(SampleStats.objects.all(), passing_filter=False)
    _copy(SampleStatsPassingFilter.objects.all(), passing_filter=True)


def _noop(apps, schema_editor):
    # The schema-create migration's reverse drops the table entirely.
    # No backfill rollback needed.
    pass


class Migration(migrations.Migration):

    dependencies = [
        ("snpdb", "0178_cohortgenotypestats"),
    ]

    operations = [
        migrations.RunPython(_backfill, _noop),
    ]

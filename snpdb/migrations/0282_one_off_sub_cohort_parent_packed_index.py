from collections import defaultdict

from django.db import migrations
from django.db.models import F, Value


def _sub_cohort_parent_packed_index(apps, _schema_editor):
    """ Sub cohorts read their parent's CohortGenotype arrays, but were renumbered 0..n-1 like a custom
        cohort (#2155). Put back the parent's index for each sample """
    Cohort = apps.get_model("snpdb", "Cohort")
    CohortSample = apps.get_model("snpdb", "CohortSample")

    parent_cohort_id_by_sub_cohort_id = dict(Cohort.objects.filter(parent_cohort__isnull=False)
                                             .values_list("pk", "parent_cohort_id"))
    if not parent_cohort_id_by_sub_cohort_id:
        return

    packed_index_by_cohort_and_sample = {}
    cohort_ids = set(parent_cohort_id_by_sub_cohort_id) | set(parent_cohort_id_by_sub_cohort_id.values())
    for cohort_id, sample_id, packed_index in CohortSample.objects.filter(cohort_id__in=cohort_ids).values_list(
            "cohort_id", "sample_id", "cohort_genotype_packed_field_index"):
        packed_index_by_cohort_and_sample[(cohort_id, sample_id)] = packed_index

    cohort_samples_by_sub_cohort_id = defaultdict(list)
    for cs in CohortSample.objects.filter(cohort_id__in=parent_cohort_id_by_sub_cohort_id):
        cohort_samples_by_sub_cohort_id[cs.cohort_id].append(cs)

    for sub_cohort_id, cohort_samples in cohort_samples_by_sub_cohort_id.items():
        parent_cohort_id = parent_cohort_id_by_sub_cohort_id[sub_cohort_id]
        parent_indexes = [packed_index_by_cohort_and_sample.get((parent_cohort_id, cs.sample_id))
                          for cs in cohort_samples]
        if None in parent_indexes:
            print(f"Sub cohort {sub_cohort_id} has samples outside parent cohort {parent_cohort_id} - skipping")
            continue
        if all(cs.cohort_genotype_packed_field_index == i for cs, i in zip(cohort_samples, parent_indexes)):
            continue

        print(f"Sub cohort {sub_cohort_id} packed indexes set from parent cohort {parent_cohort_id}")
        # packed fields are unique per cohort so need to move them sufficiently high to be able to rearrange
        highest_packed_index = max(*parent_indexes, *(cs.cohort_genotype_packed_field_index for cs in cohort_samples))
        out_of_range = F("cohort_genotype_packed_field_index") + Value(highest_packed_index + 1)
        CohortSample.objects.filter(cohort_id=sub_cohort_id).update(cohort_genotype_packed_field_index=out_of_range)
        for cs, parent_index in zip(cohort_samples, parent_indexes):
            cs.cohort_genotype_packed_field_index = parent_index
        CohortSample.objects.bulk_update(cohort_samples, ["cohort_genotype_packed_field_index"], batch_size=1000)


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0281_genomic_intervals_collection_error_message'),
    ]

    operations = [
        migrations.RunPython(_sub_cohort_parent_packed_index, reverse_code=migrations.RunPython.noop),
    ]

from collections import Counter

import celery

from patients.models_enums import Zygosity
from snpdb.models import VCF, SampleLocusCount


@celery.shared_task
def do_sample_locus_count_for_vcf_id(vcf_id):
    vcf = VCF.objects.get(pk=vcf_id)
    cohort = vcf.cohort
    cgc = cohort.cohort_genotype_collection
    qs = vcf.get_variant_qs().order_by("locus_id")

    total_sample_locus_count = [Counter()] * cohort.sample_count
    sample_locus_count = [0] * cohort.sample_count
    last_locus_id = 0
    values_qs = qs.values_list("locus_id", f"{cgc.cohortgenotype_alias}__samples_zygosity")
    for locus_id, samples_zygosity in values_qs.iterator():
        if last_locus_id != locus_id:
            for i, count in enumerate(sample_locus_count):
                total_sample_locus_count[i][count] += 1
            sample_locus_count = [0] * cohort.sample_count
        last_locus_id = locus_id

        for i, sz in enumerate(samples_zygosity):
            if sz in Zygosity.VARIANT:
                sample_locus_count[i] += 1

    # Last locus
    for i, count in enumerate(sample_locus_count):
        total_sample_locus_count[i][count] += 1

    for cs in cohort.cohortsample_set.all():
        locus_counts = total_sample_locus_count[cs.cohort_genotype_packed_field_index]
        keys = list(locus_counts)
        keys.remove(0)  # Don't care about missing
        for locus_count in keys:
            count = locus_counts[locus_count]
            SampleLocusCount.objects.update_or_create(sample=cs.sample, locus_count=locus_count,
                                                      defaults={"count": count})

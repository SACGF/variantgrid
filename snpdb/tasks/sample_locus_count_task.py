from collections import Counter

import celery

from patients.models_enums import Zygosity
from snpdb.models import VCF, SampleLocusCount


@celery.shared_task
def do_sample_locus_count_for_vcf_id(vcf_id):
    vcf = VCF.objects.get(pk=vcf_id)
    cgc = vcf.cohort.cohort_genotype_collection
    samples = list(vcf.cohort.get_samples())
    num_samples = len(samples)
    qs = vcf.get_variant_qs().order_by("locus_id")

    total_sample_locus_count = [Counter() for _ in samples]
    sample_locus_count = [0] * num_samples
    last_locus_id = 0
    for locus_id, samples_zygosity in qs.values_list("locus_id", f"{cgc.cohortgenotype_alias}__samples_zygosity"):
        if last_locus_id != locus_id:
            for i, count in enumerate(sample_locus_count):
                total_sample_locus_count[i][count] += 1
            sample_locus_count = [0] * num_samples
        last_locus_id = locus_id

        for i, sz in enumerate(samples_zygosity):
            if sz in Zygosity.VARIANT:
                sample_locus_count[i] += 1

    # Last locus
    for i, count in enumerate(sample_locus_count):
        total_sample_locus_count[i][count] += 1

    for i, sample in enumerate(samples):
        locus_counts = total_sample_locus_count[i]
        keys = list(locus_counts)
        keys.remove(0)  # Don't care about missing
        for locus_count in keys:
            count = locus_counts[locus_count]
            SampleLocusCount.objects.update_or_create(sample=sample, locus_count=locus_count,
                                                      defaults={"count": count})

from django.contrib.auth.models import User
from django.utils import timezone

from snpdb.models import CohortGenotypeCollection, Trio, CohortSample, ImportStatus, Sample, VCF, GenomeBuild, Cohort, \
    assign_permission_to_user_and_groups


def create_fake_trio(user: User, genome_build: GenomeBuild) -> Trio:
    vcf = VCF.objects.create(name="test_urls_vcf", genotype_samples=1, genome_build=genome_build,
                             import_status=ImportStatus.SUCCESS,
                             user=user, date=timezone.now())
    sample = Sample.objects.create(name="sample1", vcf=vcf, import_status=ImportStatus.SUCCESS)
    assign_permission_to_user_and_groups(user, vcf)
    assign_permission_to_user_and_groups(user, sample)

    mother_sample = Sample.objects.create(name="mother", vcf=vcf)
    father_sample = Sample.objects.create(name="father", vcf=vcf)
    cohort = Cohort.objects.create(name="test_urls_cohort", vcf=vcf, genome_build=genome_build,
                                   import_status=ImportStatus.SUCCESS)

    proband_cs = CohortSample.objects.create(cohort=cohort, sample=sample,
                                             cohort_genotype_packed_field_index=1, sort_order=1)
    mother_cs = CohortSample.objects.create(cohort=cohort, sample=mother_sample,
                                            cohort_genotype_packed_field_index=2, sort_order=2)
    father_cs = CohortSample.objects.create(cohort=cohort, sample=father_sample,
                                            cohort_genotype_packed_field_index=3, sort_order=3)
    trio = Trio.objects.create(name="test_urls_trio",
                               cohort=cohort,
                               mother=mother_cs,
                               mother_affected=True,
                               father=father_cs,
                               father_affected=False,
                               proband=proband_cs)
    assign_permission_to_user_and_groups(user, cohort)

    # Cohort version has been bumped every time a cohort sample has been added
    CohortGenotypeCollection.objects.create(cohort=cohort,
                                            cohort_version=cohort.version,
                                            num_samples=cohort.cohortsample_set.count())

    return trio

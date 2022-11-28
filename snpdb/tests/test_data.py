from django.contrib.auth.models import User
from django.utils import timezone

from library.guardian_utils import assign_permission_to_user_and_groups
from pedigree.models import Pedigree, PedFileFamily, PedFile
from snpdb.models import CohortGenotypeCollection, Trio, CohortSample, ImportStatus, Sample, VCF, GenomeBuild, Cohort, \
    VCFFilter

def create_fake_cohort(user: User, genome_build: GenomeBuild) -> Cohort:
    vcf = VCF.objects.create(name="test_urls_vcf", genotype_samples=1, genome_build=genome_build,
                             import_status=ImportStatus.SUCCESS,
                             user=user, date=timezone.now())
    VCFFilter.objects.create(vcf=vcf, filter_code="X", filter_id='YOUSHALLNOTPASS', description="fdas")
    sample = Sample.objects.create(name="proband", vcf=vcf, import_status=ImportStatus.SUCCESS)
    assign_permission_to_user_and_groups(user, vcf)
    assign_permission_to_user_and_groups(user, sample)

    mother_sample = Sample.objects.create(name="mother", vcf=vcf)
    father_sample = Sample.objects.create(name="father", vcf=vcf)
    cohort = Cohort.objects.create(name="test_urls_cohort", user=user, vcf=vcf, genome_build=genome_build,
                                   import_status=ImportStatus.SUCCESS)

    CohortSample.objects.create(cohort=cohort, sample=sample,
                                cohort_genotype_packed_field_index=0, sort_order=0)
    CohortSample.objects.create(cohort=cohort, sample=mother_sample,
                                cohort_genotype_packed_field_index=1, sort_order=1)
    CohortSample.objects.create(cohort=cohort, sample=father_sample,
                                cohort_genotype_packed_field_index=2, sort_order=2)

    assign_permission_to_user_and_groups(user, cohort)

    # Cohort version has been bumped every time a cohort sample has been added
    CohortGenotypeCollection.objects.create(cohort=cohort,
                                            cohort_version=cohort.version,
                                            num_samples=cohort.cohortsample_set.count())
    return cohort


def create_fake_trio(user: User, genome_build: GenomeBuild) -> Trio:
    cohort = create_fake_cohort(user, genome_build)
    proband_cs = cohort.cohortsample_set.get(sample__name='proband')
    mother_cs = cohort.cohortsample_set.get(sample__name='mother')
    father_cs = cohort.cohortsample_set.get(sample__name='father')

    trio = Trio.objects.create(name="test_urls_trio",
                               user=user,
                               cohort=cohort,
                               mother=mother_cs,
                               mother_affected=True,
                               father=father_cs,
                               father_affected=False,
                               proband=proband_cs)

    return trio


def create_fake_pedigree(user: User, genome_build: GenomeBuild) -> Pedigree:
    cohort = create_fake_cohort(user, genome_build)

    ped_file = PedFile.objects.get_or_create(name="fakepf", user=user,
                                             import_status=ImportStatus.SUCCESS)[0]
    assign_permission_to_user_and_groups(user, ped_file)
    ped_file_family = PedFileFamily.objects.get_or_create(name="fake family", ped_file=ped_file)[0]
    pedigree = Pedigree.objects.get_or_create(user=user, name="fake pedigree",
                                              cohort=cohort, ped_file_family=ped_file_family)[0]
    return pedigree

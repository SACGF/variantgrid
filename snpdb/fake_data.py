"""
Fake samples, cohorts and people: the VCF / cohort / trio / quad / duo / pedigree builders tests use, and the
'create_fake_data' steps "people" (the fake organisation, labs and users every later step shares) and "trio" (a trio
with genotypes over the fake variants). fake_alleles_for_variants / delete_unused_fake_alleles are the alleles the
tag and classification steps share.
"""
import random
import secrets
from dataclasses import dataclass

from django.contrib.auth.models import Group, User
from django.db import transaction
from django.db.models import ProtectedError
from django.utils import timezone

from library.fake_data import FakeData, FakeDataContext, register
from library.guardian_utils import all_users_group, assign_permission_to_user_and_groups
from pedigree.models import PedFile, PedFileFamily, Pedigree
from snpdb.models import (
    VCF,
    Allele,
    AlleleConversionTool,
    AlleleOrigin,
    Cohort,
    CohortGenotypeCollection,
    CohortSample,
    Duo,
    DuoRelationship,
    GenomeBuild,
    ImportStatus,
    Organization,
    Quad,
    Sample,
    Trio,
    Variant,
    VariantAllele,
    VCFFilter,
)
from snpdb.models.models import Lab
from snpdb.models.models_cohort import CohortGenotype


def create_fake_cohort(user: User, genome_build: GenomeBuild, name: str = "test_urls") -> Cohort:
    """ 3 samples (proband, mother, father) on one VCF, with an empty CohortGenotypeCollection """
    vcf = VCF.objects.create(name=f"{name}_vcf", genotype_samples=1, genotype_field="GT", allele_depth_field="AD", genome_build=genome_build,
                             import_status=ImportStatus.SUCCESS,
                             user=user, date=timezone.now())
    VCFFilter.objects.create(vcf=vcf, filter_code="X", filter_id='YOUSHALLNOTPASS', description="fdas")
    sample = Sample.objects.create(name="proband", vcf=vcf, import_status=ImportStatus.SUCCESS)
    assign_permission_to_user_and_groups(user, vcf)
    assign_permission_to_user_and_groups(user, sample)

    mother_sample = Sample.objects.create(name="mother", vcf=vcf)
    father_sample = Sample.objects.create(name="father", vcf=vcf)
    cohort = Cohort.objects.create(name=f"{name}_cohort", user=user, vcf=vcf, genome_build=genome_build,
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


def create_fake_trio(user: User, genome_build: GenomeBuild, name: str = "test_urls") -> Trio:
    """ Mother affected, father not, over its own create_fake_cohort """
    cohort = create_fake_cohort(user, genome_build, name=name)
    proband_cs = cohort.cohortsample_set.get(sample__name='proband')
    mother_cs = cohort.cohortsample_set.get(sample__name='mother')
    father_cs = cohort.cohortsample_set.get(sample__name='father')

    trio = Trio.objects.create(name=f"{name}_trio",
                               user=user,
                               cohort=cohort,
                               mother=mother_cs,
                               mother_affected=True,
                               father=father_cs,
                               father_affected=False,
                               proband=proband_cs)

    return trio


def create_fake_quad(user: User, genome_build: GenomeBuild, sibling_affected: bool = False,
                     name: str = "test_quad") -> Quad:
    """4-sample Cohort (proband, mother, father, sibling) + a Quad."""
    vcf = VCF.objects.create(
        name=f"{name}_vcf", genotype_samples=1, genotype_field="GT", allele_depth_field="AD", genome_build=genome_build,
        import_status=ImportStatus.SUCCESS, user=user, date=timezone.now()
    )
    proband_sample = Sample.objects.create(name="proband", vcf=vcf, import_status=ImportStatus.SUCCESS)
    mother_sample = Sample.objects.create(name="mother", vcf=vcf)
    father_sample = Sample.objects.create(name="father", vcf=vcf)
    sibling_sample = Sample.objects.create(name="sibling", vcf=vcf)

    assign_permission_to_user_and_groups(user, vcf)
    assign_permission_to_user_and_groups(user, proband_sample)

    cohort = Cohort.objects.create(
        name=f"{name}_cohort", user=user, vcf=vcf,
        genome_build=genome_build, import_status=ImportStatus.SUCCESS
    )
    for i, sample in enumerate([proband_sample, mother_sample, father_sample, sibling_sample]):
        CohortSample.objects.create(
            cohort=cohort, sample=sample,
            cohort_genotype_packed_field_index=i, sort_order=i
        )
    assign_permission_to_user_and_groups(user, cohort)

    CohortGenotypeCollection.objects.create(
        cohort=cohort, cohort_version=cohort.version,
        num_samples=cohort.cohortsample_set.count()
    )

    proband_cs = cohort.cohortsample_set.get(sample__name='proband')
    mother_cs = cohort.cohortsample_set.get(sample__name='mother')
    father_cs = cohort.cohortsample_set.get(sample__name='father')
    sibling_cs = cohort.cohortsample_set.get(sample__name='sibling')

    return Quad.objects.create(
        name=name,
        user=user,
        cohort=cohort,
        mother=mother_cs, mother_affected=False,
        father=father_cs, father_affected=False,
        proband=proband_cs,
        sibling=sibling_cs, sibling_affected=sibling_affected,
    )


def create_fake_duo(user: User, genome_build: GenomeBuild,
                    relationship: str = DuoRelationship.MOTHER,
                    relative_affected: bool = False, name: str = "test_duo") -> Duo:
    """2-sample Cohort (proband, relative) + a Duo - the relative is named after the relationship."""
    vcf = VCF.objects.create(
        name=f"{name}_vcf", genotype_samples=1, genotype_field="GT", allele_depth_field="AD", genome_build=genome_build,
        import_status=ImportStatus.SUCCESS, user=user, date=timezone.now()
    )
    relative_name = DuoRelationship(relationship).label.lower()
    proband_sample = Sample.objects.create(name="proband", vcf=vcf, import_status=ImportStatus.SUCCESS)
    relative_sample = Sample.objects.create(name=relative_name, vcf=vcf)

    assign_permission_to_user_and_groups(user, vcf)
    assign_permission_to_user_and_groups(user, proband_sample)

    cohort = Cohort.objects.create(
        name=f"{name}_cohort", user=user, vcf=vcf,
        genome_build=genome_build, import_status=ImportStatus.SUCCESS
    )
    for i, sample in enumerate([proband_sample, relative_sample]):
        CohortSample.objects.create(
            cohort=cohort, sample=sample,
            cohort_genotype_packed_field_index=i, sort_order=i
        )
    assign_permission_to_user_and_groups(user, cohort)

    CohortGenotypeCollection.objects.create(
        cohort=cohort, cohort_version=cohort.version,
        num_samples=cohort.cohortsample_set.count()
    )

    return Duo.objects.create(
        name=name,
        user=user,
        cohort=cohort,
        proband=cohort.cohortsample_set.get(sample__name='proband'),
        relative=cohort.cohortsample_set.get(sample__name=relative_name),
        relationship=relationship,
        relative_affected=relative_affected,
    )


def create_fake_pedigree(user: User, genome_build: GenomeBuild, cohort: Cohort = None,
                         name: str = "fake pedigree") -> Pedigree:
    """ Over cohort, or its own create_fake_cohort """
    if cohort is None:
        cohort = create_fake_cohort(user, genome_build)

    ped_file = PedFile.objects.get_or_create(name="fakepf", user=user,
                                             import_status=ImportStatus.SUCCESS)[0]
    assign_permission_to_user_and_groups(user, ped_file)
    ped_file_family = PedFileFamily.objects.get_or_create(name="fake family", ped_file=ped_file)[0]
    pedigree = Pedigree.objects.get_or_create(user=user, name=name,
                                              cohort=cohort, ped_file_family=ped_file_family)[0]
    return pedigree


DEFAULT_GENOTYPE_VALUES = object()  # so a test can ask for a NULL array, eg a VCF with no AF field


def make_cohort_genotype(cgc: CohortGenotypeCollection, variant, samples_zygosity: str,
                         allele_depth=DEFAULT_GENOTYPE_VALUES, allele_frequency=DEFAULT_GENOTYPE_VALUES):
    """ One CohortGenotype row. samples_zygosity has a letter per sample in packed order - E=HET, R=HOM_REF,
        O=HOM_ALT, U=UNKNOWN, .=MISSING - so a trio "ERR" is a het de novo.
        allele_depth/allele_frequency are per sample, in packed order - the mosaic modes read them """
    n = len(samples_zygosity)
    CohortGenotype.objects.create(
        collection=cgc,
        variant=variant,
        ref_count=samples_zygosity.count('R'),
        het_count=samples_zygosity.count('E'),
        hom_count=samples_zygosity.count('O'),
        samples_zygosity=samples_zygosity,
        samples_allele_depth=[20] * n if allele_depth is DEFAULT_GENOTYPE_VALUES else allele_depth,
        samples_allele_frequency=[100] * n if allele_frequency is DEFAULT_GENOTYPE_VALUES else allele_frequency,
        samples_read_depth=[30] * n,
        samples_genotype_quality=[30] * n,
        samples_phred_likelihood=[0] * n,
    )


FAKE_ORGANIZATION_GROUP = "fake_health"
GERMLINE = "germline"
SOMATIC = "somatic"
BATCH_SIZE = 5000


@dataclass(frozen=True)
class FakeLab:
    name: str
    group_name: str
    focus: str


FAKE_LABS = [
    FakeLab("Fake Germline Unit", f"{FAKE_ORGANIZATION_GROUP}/fake_germline_unit", GERMLINE),
    FakeLab("Fake Somatic Unit", f"{FAKE_ORGANIZATION_GROUP}/fake_somatic_unit", SOMATIC),
]

# Username and the focus of the lab(s) they're in - None is both
FAKE_USERS = {
    "fake_freya_germline": GERMLINE,
    "fake_gus_germline": GERMLINE,
    "fake_hana_germline": GERMLINE,
    "fake_alice_somatic": SOMATIC,
    "fake_bruno_somatic": SOMATIC,
    "fake_chen_somatic": SOMATIC,
    "fake_rotating_registrar": None,
    "fake_import_bot": None,
}


@register
class FakePeople(FakeData):
    name = "people"
    help = ("The 'Fake Health Network' organisation, a germline and a somatic lab, and the fake_* users "
            "(password printed when they are created) every later step uses")

    def create(self, context: FakeDataContext, **options):
        organization, _ = Organization.objects.get_or_create(
            group_name=FAKE_ORGANIZATION_GROUP, defaults={"name": "Fake Health Network", "short_name": "FakeHealth"})
        labs = {}
        for fake_lab in FAKE_LABS:
            labs[fake_lab.focus], _ = Lab.objects.get_or_create(
                group_name=fake_lab.group_name,
                defaults={"name": fake_lab.name, "city": "Faketown", "organization": organization})

        password = secrets.token_urlsafe(12)
        created_usernames = []
        for username, focus in FAKE_USERS.items():
            user, created = User.objects.get_or_create(username=username,
                                                       defaults={"first_name": "Fake", "last_name": username[5:]})
            if created:
                user.set_password(password)
                user.save()
                lab_groups = [lab.group for lab_focus, lab in labs.items() if focus in (None, lab_focus)]
                user.groups.add(all_users_group(), *lab_groups)
                created_usernames.append(username)
            context.users.append(user)
        context.labs = [labs[GERMLINE], labs[SOMATIC]]

        if created_usernames:
            context.stdout.write(f"Created users {', '.join(created_usernames)} with password: {password}")
        else:
            context.stdout.write("Fake users already exist")

    def delete(self, context: FakeDataContext, **options):
        User.objects.filter(username__in=FAKE_USERS).delete()
        Lab.objects.filter(group_name__in=[fake_lab.group_name for fake_lab in FAKE_LABS]).delete()
        Organization.objects.filter(group_name=FAKE_ORGANIZATION_GROUP).delete()
        Group.objects.filter(name__startswith=FAKE_ORGANIZATION_GROUP).delete()
        context.stdout.write("Deleted the fake users, labs and organisation")


FAKE_TRIO_NAME = "fake"

# samples_zygosity in cohort order - proband, mother (affected), father - and how often each turns up
TRIO_GENOTYPE_WEIGHTS = {
    "EER": 30,  # inherited from mother
    "ERE": 25,  # inherited from father
    "REE": 10,  # parents carry it, proband doesn't
    "EEE": 10,
    "OOO": 7,
    "OEE": 8,  # recessive
    "ERR": 5,  # de novo
    "ORR": 2,
    "E.R": 3,  # a missing call
}


@register
class FakeTrio(FakeData):
    name = "trio"
    help = "A VCF with proband / mother / father samples, cohort, genotypes over the fake variants, a trio and pedigree"
    requires = ("people", "variants")

    @classmethod
    def add_arguments(cls, parser):
        parser.add_argument("--max-genotypes", type=int, default=1000,
                            help="Genotype at most this many of the variants")

    def create(self, context: FakeDataContext, **options):
        owner = context.users[0]
        if trio := Trio.objects.filter(name=f"{FAKE_TRIO_NAME}_trio", user=owner).first():
            context.trio = trio
            context.stdout.write(f"{trio} already exists")
            return

        rng = random.Random(context.seed)
        trio = create_fake_trio(owner, context.genome_build, name=FAKE_TRIO_NAME)
        create_fake_pedigree(owner, context.genome_build, cohort=trio.cohort, name=f"{FAKE_TRIO_NAME} pedigree")

        variant_ids = sorted(variant_id for variant_ids in context.variant_ids_by_gene.values()
                             for variant_id in variant_ids)
        if len(variant_ids) > options["max_genotypes"]:
            variant_ids = sorted(rng.sample(variant_ids, options["max_genotypes"]))
        genotypes = list(TRIO_GENOTYPE_WEIGHTS)
        weights = list(TRIO_GENOTYPE_WEIGHTS.values())
        cgc = trio.cohort.cohort_genotype_collection
        for variant_id in variant_ids:
            make_cohort_genotype(cgc, Variant(pk=variant_id), rng.choices(genotypes, weights=weights)[0])
        context.trio = trio
        context.stdout.write(f"Created {trio} with {len(variant_ids)} genotypes")

    def delete(self, context: FakeDataContext, **options):
        trios = list(Trio.objects.filter(name=f"{FAKE_TRIO_NAME}_trio", user__username__in=FAKE_USERS)
                     .select_related("cohort__vcf"))
        for trio in trios:
            cohort = trio.cohort
            Pedigree.objects.filter(cohort=cohort).delete()
            Trio.objects.filter(pk=trio.pk).delete()
            Cohort.objects.filter(pk=cohort.pk).delete()
            VCF.objects.filter(pk=cohort.vcf_id).delete()
        PedFile.objects.filter(name="fakepf", user__username__in=FAKE_USERS).delete()
        context.stdout.write(f"Deleted {len(trios)} fake trios")


def fake_alleles_for_variants(genome_build: GenomeBuild, variant_ids: list[int]) -> dict[int, int]:
    """ allele_id per variant_id, making an allele where a variant has none. The fake steps share them, so tags
        and classifications on a variant land on the one allele, as real ones do """
    allele_ids = dict(VariantAllele.objects.filter(genome_build=genome_build, variant_id__in=variant_ids)
                      .values_list("variant_id", "allele_id"))
    without_allele = [variant_id for variant_id in dict.fromkeys(variant_ids) if variant_id not in allele_ids]
    alleles = Allele.objects.bulk_create([Allele() for _ in without_allele], batch_size=BATCH_SIZE)
    VariantAllele.objects.bulk_create([
        VariantAllele(variant_id=variant_id, genome_build=genome_build, allele=allele,
                      origin=AlleleOrigin.IMPORTED_TO_DATABASE, allele_linking_tool=AlleleConversionTool.SAME_CONTIG)
        for variant_id, allele in zip(without_allele, alleles)], batch_size=BATCH_SIZE)
    allele_ids.update((variant_id, allele.pk) for variant_id, allele in zip(without_allele, alleles))
    return allele_ids


def delete_unused_fake_alleles(allele_ids) -> int:
    """ Deletes the alleles no other fake data still uses - VariantTag and Classification protect their alleles,
        so those wait for the step that made them to be deleted. Returns how many were deleted """
    alleles_qs = Allele.objects.filter(pk__in=allele_ids, clingen_allele__isnull=True)
    try:
        with transaction.atomic():
            _, deleted = alleles_qs.delete()
    except ProtectedError as protected_error:
        in_use = {getattr(obj, "allele_id", None) for obj in protected_error.protected_objects}
        _, deleted = alleles_qs.exclude(pk__in=in_use).delete()
    return deleted.get(Allele._meta.label, 0)

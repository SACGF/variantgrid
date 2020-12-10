from django.contrib.auth.models import User
from django.db import models
from django.db.models.deletion import CASCADE
from django.urls.base import reverse

from library.django_utils import SortByPKMixin
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsAutoInitialSaveMixin, \
    GuardianPermissionsMixin
from patients.models_enums import Sex
from snpdb.models import ImportStatus, Cohort, CohortSample, Sample


class PedFile(GuardianPermissionsMixin, models.Model):
    """ See http://www.helsinki.fi/~tsjuntun/autogscan/pedigreefile.html """
    name = models.TextField()
    user = models.ForeignKey(User, on_delete=CASCADE)
    import_status = models.CharField(max_length=1, choices=ImportStatus.CHOICES, default=ImportStatus.CREATED)

    def __str__(self):
        name = f"({self.pk}) {self.name}"
        if self.import_status != ImportStatus.SUCCESS:
            name += f" (status={self.import_status})"
        return name

    def get_absolute_url(self):
        return reverse('view_ped_file', kwargs={"ped_file_id": self.pk})


class PedFileFamily(models.Model):
    """ Collects records together by family/pedigree (1st column)
        there can be multiple pedigrees per PED file """
    ped_file = models.ForeignKey(PedFile, on_delete=CASCADE)
    name = models.TextField(blank=False)

    def get_records_count(self):
        return self.pedfilerecord_set.count()

    @property
    def errors(self):
        records = list(self.pedfilerecord_set.all())
        return validate(records)

    def is_valid(self):
        return not self.errors

    @staticmethod
    def filter_for_user(user):
        pf_qs = PedFile.filter_for_user(user)
        pf_qs = pf_qs.filter(import_status=ImportStatus.SUCCESS)
        return PedFileFamily.objects.filter(ped_file__in=pf_qs)

    def __str__(self):
        return f"Family {self.name} (ped file: {self.ped_file.name})"


class PedFileRecord(models.Model):
    family = models.ForeignKey(PedFileFamily, on_delete=CASCADE)
    sample = models.TextField(blank=False, null=False)
    father = models.ForeignKey("self", related_name='father_rel', null=True, on_delete=CASCADE)
    mother = models.ForeignKey("self", related_name='mother_rel', null=True, on_delete=CASCADE)
    sex = models.CharField(max_length=1, choices=Sex.choices, null=True)
    affection = models.BooleanField(null=True)

    @property
    def name(self):
        return self.sample

    @property
    def affected(self):
        if self.affection:
            description = "affected"
        else:
            description = "unaffected"
        return description

    def __str__(self):
        description = f"{self.sample} ({self.get_sex_display()}/{self.affected})"
        if self.father:
            description += f" father: {self.father.sample}"
        if self.mother:
            description += f" mother: {self.mother.sample}"
        return description


def validate(records):
    """ Should work on pedfile and pedigree """
    errors_list = []
    # Make sure there is at least 1 sample
    if len(records) < 1:
        errors_list.apped("There must be at least 1 record")

    any_affected = False
    for r in records:
        # Validate sexes
        errors = []
        if r.father and r.father.sex != Sex.MALE:
            errors.append(f"Father ({r.father.name}) is not male")
        if r.mother and r.mother.sex != Sex.FEMALE:
            errors.append(f"Mother ({r.mother.name}) is not female")
        if r.affection:
            any_affected = True

        if errors:
            msg = f"{r.name} : {','.join(errors)}"
            errors_list.append(msg)
    if not any_affected:
        errors_list.append("Must have at least 1 affected individual.")
    return errors_list


class Pedigree(GuardianPermissionsAutoInitialSaveMixin, SortByPKMixin, models.Model):
    user = models.ForeignKey(User, on_delete=CASCADE)
    name = models.TextField(blank=False)
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    ped_file_family = models.ForeignKey(PedFileFamily, on_delete=CASCADE)

    def get_samples(self, affected=None):
        """ returns a Sample queryset """

        pfr_qs = self.ped_file_family.pedfilerecord_set.all()
        if affected is not None:
            pfr_qs = pfr_qs.filter(affection=affected)
        samples = pfr_qs.values_list("cohortsamplepedfilerecord__cohort_sample__sample", flat=True)
        return Sample.objects.filter(pk__in=samples).order_by("pk")

    @property
    def genome_build(self):
        return self.cohort.genome_build

    def get_absolute_url(self):
        return reverse('view_pedigree', kwargs={"pedigree_id": self.pk})

    @classmethod
    def get_listing_url(cls):
        return reverse('pedigrees')

    def __str__(self):
        return self.name


class CohortSamplePedFileRecord(models.Model):
    """ There should be one of these for every ped_file_record to make a valid pedigree """
    pedigree = models.ForeignKey(Pedigree, on_delete=CASCADE)
    cohort_sample = models.ForeignKey(CohortSample, on_delete=CASCADE)
    ped_file_record = models.ForeignKey(PedFileRecord, on_delete=CASCADE)


class PedigreeInheritance:
    AUTOSOMAL_RECESSIVE = 'R'
    AUTOSOMAL_DOMINANT = 'D'
    CHOICES = (
        (AUTOSOMAL_RECESSIVE, 'Auto. Recessive'),
        (AUTOSOMAL_DOMINANT, 'Auto. Dominant'),
    )


def create_automatch_pedigree(user, ped_file_family, cohort):
    name = f"Auto Pedigree for {ped_file_family}/{cohort}"
    pedigree = Pedigree.objects.create(user=user,
                                       name=name,
                                       cohort=cohort,
                                       ped_file_family=ped_file_family)

    cohort_samples_by_sample_name = {cs.sample.name: cs for cs in cohort.cohortsample_set.all()}
    for pfr in ped_file_family.pedfilerecord_set.all():
        cohort_sample = cohort_samples_by_sample_name.get(pfr.sample)
        if cohort_sample:
            CohortSamplePedFileRecord.objects.create(pedigree=pedigree,
                                                     cohort_sample=cohort_sample,
                                                     ped_file_record=pfr)

    return pedigree

from typing import Union

from django.contrib.auth.models import Group, User
from django.core.exceptions import PermissionDenied
from django.db import models
from django.db.models import Q, QuerySet
from django.db.models.deletion import CASCADE, PROTECT, SET_NULL
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel

from genes.models import GeneList, GeneSymbol
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.enums import ModificationOperation
from library.preview_request import PreviewModelMixin
from pathtests.models_enums import (
    CaseState,
    CaseWorkflowStatus,
    InvestigationType,
    PathologyTestGeneModificationOutcome,
)
from patients.models import (
    TEST_PATIENT_KWARGS,
    Clinician,
    ExternallyManagedModel,
    Patient,
    get_lead_scientist_users_for_user,
)
from patients.models_enums import PopulationGroup
from seqauto.models import EnrichmentKit, Experiment, SequencingRun
from snpdb.models import Sample, Wiki
from snpdb.models.models_enums import ImportStatus


class PathologyTest(TimeStampedModel):
    name = models.TextField(primary_key=True)
    curator = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    deleted = models.BooleanField(default=False)
    empty_test = models.BooleanField(default=False)  # For custom tests

    def can_write(self, user_or_group: Union[User, Group]) -> bool:
        return self.is_curator(user_or_group)

    def is_curator(self, user_or_group: Union[User, Group]):
        return self.curator == user_or_group

    def check_is_curator(self, user_or_group: Union[User, Group]):
        if not self.is_curator(user_or_group):
            raise PermissionDenied(f"You are not the curator of pathology test {self}")

    def get_active_test_version(self):
        active_test_version = None
        try:
            active_test_version = self.activepathologytestversion.pathology_test_version
        except Exception:
            pass
        return active_test_version

    def get_latest_confirmed_version(self):
        return self.pathologytestversion_set.filter(confirmed_date__isnull=False).order_by("-version").first()

    def delete_test(self):
        """ Soft delete - the API stops serving it until it is restored """
        self.deleted = True
        self.save()
        ActivePathologyTestVersion.objects.filter(pathology_test=self).delete()

    def restore_test(self):
        """ Inverse of delete_test: the latest confirmed version becomes active again """
        self.deleted = False
        self.save()
        if latest_confirmed := self.get_latest_confirmed_version():
            latest_confirmed.set_as_active_test()

    def get_absolute_url(self):
        return reverse("view_pathology_test", kwargs={"name": self.name})

    def __str__(self):
        return self.name


class PathologyTestWiki(Wiki):
    pathology_test = models.OneToOneField(PathologyTest, on_delete=CASCADE)

    def _get_restricted_object(self):
        return self.pathology_test


class PathologyTestSynonyms(models.Model):
    pathology_test = models.ForeignKey(PathologyTest, on_delete=CASCADE)
    synonym_name = models.TextField()


class PathologyTestVersion(TimeStampedModel):
    pathology_test = models.ForeignKey(PathologyTest, on_delete=CASCADE)
    version = models.IntegerField(default=1)
    confirmed_date = models.DateTimeField(null=True)
    gene_list = models.ForeignKey(GeneList, on_delete=PROTECT)
    enrichment_kit = models.ForeignKey(EnrichmentKit, null=True, blank=False, on_delete=PROTECT)

    class Meta:
        unique_together = ("pathology_test", "version")

    @property
    def can_modify(self):
        return self.confirmed_date is None

    @property
    def can_confirm(self):
        return self.can_modify and self.enrichment_kit

    @property
    def is_active_test(self):
        try:
            return self.activepathologytestversion is not None
        except ActivePathologyTestVersion.DoesNotExist:
            return False

    def set_as_active_test(self):
        manager = ActivePathologyTestVersion.objects
        manager.update_or_create(pathology_test=self.pathology_test,
                                 defaults={"pathology_test_version": self})

    def is_curator(self, user_or_group: Union[User, Group]):
        return self.pathology_test.is_curator(user_or_group)

    def check_is_curator(self, user_or_group: Union[User, Group]):
        self.pathology_test.check_is_curator(user_or_group)

    def next_version(self):
        """ Clones w/new version """

        copy = self
        gene_list = copy.gene_list.clone()
        copy.pk = None
        copy.version += 1
        copy.confirmed_date = None
        gene_list.locked = False
        gene_list.save()
        copy.gene_list = gene_list
        copy.save()
        return copy

    def save(self, *args, **kwargs):
        super().save(*args, **kwargs)

        if self.gene_list.import_status != ImportStatus.SUCCESS:
            import_status = self.gene_list.get_import_status_display()
            msg = f"{self} assigned gene_list {self.gene_list} ({self.gene_list.pk}) with invalid status of {import_status} (not success)"
            raise ValueError(msg)

        self.pathology_test.save()  # Updated last modified

        # Lock gene list if confirmed
        if self.confirmed_date and self.gene_list:
            self.gene_list.locked = True
            self.gene_list.save()

    def get_absolute_url(self):
        return reverse("view_pathology_test_version", kwargs={"pk": self.pk})

    def __str__(self):
        return f"{self.pathology_test} (v{self.version})"


class ActivePathologyTestVersion(models.Model):
    pathology_test = models.OneToOneField(PathologyTest, on_delete=CASCADE)
    pathology_test_version = models.OneToOneField(PathologyTestVersion, on_delete=CASCADE)


class PathologyTestGeneModificationRequest(TimeStampedModel):
    pathology_test_version = models.ForeignKey(PathologyTestVersion, on_delete=CASCADE)
    outcome = models.CharField(max_length=1, choices=PathologyTestGeneModificationOutcome.choices, default=PathologyTestGeneModificationOutcome.PENDING)
    operation = models.CharField(max_length=1, choices=ModificationOperation.CHOICES)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    user = models.ForeignKey(User, on_delete=CASCADE)
    comments = models.TextField(blank=True)

    def __str__(self):
        return f"{self.pathology_test_version} {self.get_operation_display()} {self.gene_symbol}: {self.get_outcome_display()}"


class Case(GuardianPermissionsMixin, PreviewModelMixin, ExternallyManagedModel):
    name = models.TextField(null=True, blank=True)
    lead_scientist = models.ForeignKey(User, null=True, blank=True, on_delete=SET_NULL)
    result_required_date = models.DateTimeField(null=True, blank=True)
    patient = models.ForeignKey(Patient, on_delete=CASCADE)
    report_date = models.DateTimeField(null=True, blank=True)
    details = models.TextField(blank=True)
    status = models.CharField(max_length=1, choices=CaseState.choices, default=CaseState.OPEN)
    workflow_status = models.CharField(max_length=2, choices=CaseWorkflowStatus.choices, default=CaseWorkflowStatus.NA)
    investigation_type = models.CharField(max_length=1, choices=InvestigationType.choices, default=InvestigationType.SINGLE_SAMPLE)

    @classmethod
    def get_permission_class(cls):
        return Patient

    def get_permission_object(self):
        # A case holds clinical details about its patient, so is as confidential as the patient
        return self.patient

    @classmethod
    def _filter_from_permission_object_qs(cls, queryset):
        return cls.objects.filter(patient__in=queryset)

    def can_write(self, user) -> bool:
        return ExternallyManagedModel.can_write(self, user) and GuardianPermissionsMixin.can_write(self, user)

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-file"

    def get_absolute_url(self):
        return reverse("view_case", kwargs={"pk": self.pk})

    @property
    def is_open(self):
        return not self.is_closed

    @property
    def is_closed(self):
        return self.status in CaseState.CLOSED_STATES

    def __str__(self):
        if self.name:
            return self.name
        if self.external_pk:
            return str(self.external_pk)
        return str(self.pk)


class CaseClinician(models.Model):
    case = models.ForeignKey(Case, on_delete=CASCADE)
    clinician = models.ForeignKey(Clinician, on_delete=CASCADE)
    specified_on_clinical_grounds = models.BooleanField(default=False)

    def __str__(self):
        return f"{self.clinician} for {self.case}"


def get_cases_qs(user):
    cases_qs = Case.filter_for_user(user)
    try:
        test_patient = Patient.objects.get(**TEST_PATIENT_KWARGS)
        cases_qs = cases_qs.exclude(patient=test_patient)
    except Exception:
        pass
    return cases_qs


def cases_for_user(user):
    """ The cases the "My cases" page shows: led by the user or a lead scientist they follow """
    users_list = get_lead_scientist_users_for_user(user)
    return get_cases_qs(user).filter(lead_scientist__in=users_list)


class PathologyTestOrder(PreviewModelMixin, ExternallyManagedModel):
    """ external_pk = usually comes from LIMS
        custom_gene_list = override of pathology_test_version
    """
    case = models.ForeignKey(Case, null=True, on_delete=CASCADE)
    pathology_test_version = models.ForeignKey(PathologyTestVersion, null=True, on_delete=CASCADE)
    custom_gene_list = models.ForeignKey(GeneList, null=True, on_delete=CASCADE)
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)
    started_library = models.DateTimeField(null=True)
    finished_library = models.DateTimeField(null=True)
    started_sequencing = models.DateTimeField(null=True)
    finished_sequencing = models.DateTimeField(null=True)
    order_completed = models.DateTimeField(null=True)
    experiment = models.ForeignKey(Experiment, null=True, on_delete=SET_NULL)
    sequencing_run = models.ForeignKey(SequencingRun, null=True, on_delete=SET_NULL)

    @classmethod
    def filter_for_user(cls, user) -> QuerySet['PathologyTestOrder']:
        if user.is_superuser:
            return cls.objects.all()
        return cls.objects.filter(Q(case__in=Case.filter_for_user(user)) | Q(case__isnull=True, user=user))

    def can_view(self, user) -> bool:
        if self.case:
            return self.case.can_view(user)
        return user.is_superuser or self.user == user

    def check_can_view(self, user):
        if not self.can_view(user):
            raise PermissionDenied(f"You do not have READ permission to view {self}")

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-clipboard-list"

    def get_absolute_url(self):
        return reverse("view_pathology_test_order", kwargs={"pk": self.pk})

    def __str__(self):
        if self.external_pk:
            return str(self.external_pk)
        return str(self.pk)


class PathologyTestOrderSample(models.Model):
    pathology_test_order = models.ForeignKey(PathologyTestOrder, on_delete=CASCADE)
    sample = models.OneToOneField(Sample, on_delete=CASCADE)


class PathologyTestOrderPopulation(models.Model):
    pathology_test_order = models.ForeignKey(PathologyTestOrder, on_delete=CASCADE)
    population = models.CharField(max_length=3, choices=PopulationGroup.choices)


def get_external_order_system_last_checked():
    # TODO: Make this more general - e.g. register somewhere in settings?
    try:
        from sapath.models.sapath_helix import HelixNGSOrdersImport
        return HelixNGSOrdersImport.get_last_checked()
    except Exception:  # App not registered?
        return None

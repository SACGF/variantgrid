from django.contrib.auth.models import User
from django.db import models
from django.db.models.deletion import CASCADE, SET_NULL, PROTECT
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel

from annotation.models.models_mim_hpo import HPOSynonym, MIMMorbidAlias
from genes.models import GeneList, Gene
from library.enums import ModificationOperation
from pathtests.models_enums import PathologyTestGeneModificationOutcome, \
    CaseState, InvestigationType, CaseWorkflowStatus
from patients.models import Patient, Clinician, ExternallyManagedModel, \
    TEST_PATIENT_KWARGS
from patients.models_enums import PopulationGroup
from seqauto.models import Experiment, SequencingRun, EnrichmentKit
from snpdb.models import Wiki, Sample
from snpdb.models.models_enums import ImportStatus


class PathologyTest(TimeStampedModel):
    name = models.TextField(primary_key=True)
    curator = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    deleted = models.BooleanField(default=False)
    empty_test = models.BooleanField(default=False)  # For custom tests

    def can_write(self, user):
        return self.is_curator(user)

    def is_curator(self, user):
        return self.curator == user

    def get_active_test_version(self):
        active_test_version = None
        try:
            active_test_version = self.activepathologytestversion.pathology_test_version
        except:
            pass
        return active_test_version

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

    def is_curator(self, user):
        return self.pathology_test.is_curator(user)

    def next_version(self):
        """ Clones w/new version """

        copy = self
        copy.pk = None
        copy.version += 1
        copy.confirmed_date = None
        gene_list = copy.gene_list.clone()
        gene_list.locked = False
        gene_list.save()
        copy.gene_list = gene_list
        copy.save()
        return copy

    def save(self, **kwargs):
        super().save(**kwargs)

        if self.gene_list.import_status != ImportStatus.SUCCESS:
            import_status = self.gene_list.get_import_status_display()
            msg = f"{self} assigned gene_list {self.gene_list} ({self.gene_list.pk}) with invalid status of {import_status} (not success)"
            raise ValueError(msg)

        self.pathology_test.save()  # Updated last modified

        # Lock gene list if confirmed
        if self.confirmed_date and self.gene_list:
            self.gene_list.locked = True
            self.gene_list.save()

    def replace_gene_list(self, gene_list, clone=True):
        """ Set new list, set name/category from old one """
        if clone:
            gene_list = gene_list.clone()

        existing_gene_list = self.gene_list

        # overwrite pk of new list with data from existing gene list
        existing_gene_list.pk = gene_list.pk
        existing_gene_list.save()
        self.gene_list = existing_gene_list
        self.save()

    def get_absolute_url(self):
        return reverse("view_pathology_test_version", kwargs={"pk": self.pk})

    def __str__(self):
        return f"{self.pathology_test} (v{self.version})"


class ActivePathologyTestVersion(models.Model):
    pathology_test = models.OneToOneField(PathologyTest, on_delete=CASCADE)
    pathology_test_version = models.OneToOneField(PathologyTestVersion, on_delete=CASCADE)


class PathologyTestGeneModificationRequest(TimeStampedModel):
    pathology_test_version = models.ForeignKey(PathologyTestVersion, on_delete=CASCADE)
    outcome = models.CharField(max_length=1, choices=PathologyTestGeneModificationOutcome.CHOICES, default=PathologyTestGeneModificationOutcome.PENDING)
    operation = models.CharField(max_length=1, choices=ModificationOperation.CHOICES)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    user = models.ForeignKey(User, on_delete=CASCADE)
    comments = models.TextField(blank=True)

    def __str__(self):
        operation = self.get_operation_display()
        name = " ".join((self.pathology_test_version, operation, self.gene.gene_symbol))
        outcome = self.get_outcome_display()
        return f"{name}: {outcome}"


class RelatedGeneLists(models.ForeignKey):
    """ A way to link eg competitor tests to yours """
    pathology_test = models.ForeignKey(PathologyTest, on_delete=CASCADE)
    gene_list = models.ForeignKey(GeneList, on_delete=CASCADE)
    comments = models.TextField(blank=True)


class Case(ExternallyManagedModel):
    name = models.TextField(null=True, blank=True)
    lead_scientist = models.ForeignKey(User, null=True, blank=True, on_delete=SET_NULL)
    result_required_date = models.DateTimeField(null=True, blank=True)
    patient = models.ForeignKey(Patient, on_delete=CASCADE)
    report_date = models.DateTimeField(null=True, blank=True)
    details = models.TextField(blank=True)
    status = models.CharField(max_length=1, choices=CaseState.CHOICES, default=CaseState.OPEN)
    workflow_status = models.CharField(max_length=2, choices=CaseWorkflowStatus.CHOICES, default=CaseWorkflowStatus.NA)
    investigation_type = models.CharField(max_length=1, choices=InvestigationType.CHOICES, default=InvestigationType.SINGLE_SAMPLE)

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


def get_cases_qs():
    cases_qs = Case.objects.all()
    try:
        test_patient = Patient.objects.get(**TEST_PATIENT_KWARGS)
        cases_qs = cases_qs.exclude(patient=test_patient)
    except:
        pass
    return cases_qs


def cases_for_user(user):
    # TODO: this is disabled - need to re-implement show_cases()
    return Case.objects.none()
    #users_list = get_lead_scientist_users_for_user(user)
    #cases_qs = get_cases_qs().filter(lead_scientist__in=users_list)


class PathologyTestOrder(ExternallyManagedModel):
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


class PathologyTestOrderHPO(models.Model):
    pathology_test_order = models.ForeignKey(PathologyTestOrder, on_delete=CASCADE)
    hpo_synonym = models.ForeignKey(HPOSynonym, on_delete=CASCADE)


class PathologyTestOrderOMIM(models.Model):
    pathology_test_order = models.ForeignKey(PathologyTestOrder, on_delete=CASCADE)
    mim_morbid_alias = models.ForeignKey(MIMMorbidAlias, on_delete=CASCADE)


def get_external_order_system_last_checked():
    # TODO: Make this more general - eg register somewhere in settings?
    try:
        from sapath.models.sapath_helix import HelixNGSOrdersImport
        return HelixNGSOrdersImport.get_last_checked()
    except:  # App not registered?
        return None

"""

Find/surface potential variants / classifications for human review

Note on naming:
* Review, Flag are already used in VG code base
* "Finding" in clinical contexts means real observed/reported result

"""
from celery.canvas import Signature
from django.contrib.auth.models import User
from django.db import models
from django.db.models import CASCADE
from django.db.models.deletion import SET_NULL
from django.urls import reverse
from model_utils.models import TimeStampedModel

from analysis.models import Analysis
from annotation.models import ClinVar, AnnotationVersion
from classification.models import Classification
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsAutoInitialSaveMixin
from library.utils.django_utils import get_cached_project_git_hash
from patients.models_enums import Zygosity
from snpdb.models import Variant, Sample, ProcessingStatus


class CandidateSearchType(models.TextChoices):
    REANALYSIS_NEW_ANNOTATION = 'R', 'Reanalysis - new annotation'
    CROSS_SAMPLE_CLASSIFICATION = 'S', 'Cross Sample Classification'
    CLASSIFICATION_EVIDENCE_UPDATE = 'E', 'Classification Evidence Update'

    def get_methods(self) -> str:
        """ This is the description of the latest methods, which will be saved to the DB at execution time
            So you can update it over time and it'll keep historical methods """
        if self is CandidateSearchType.REANALYSIS_NEW_ANNOTATION:
            methods = """Find analyses for P/LP ClinVar variants in latest version that are not in analysis version."""
        elif self is CandidateSearchType.CROSS_SAMPLE_CLASSIFICATION:
            methods = """Find samples that have a variant classified in another sample, but not themselves"""
        elif self is CandidateSearchType.CLASSIFICATION_EVIDENCE_UPDATE:
            methods = """Find Classifications that have new annotation information (gnomAD and pathogenicity predictions)"""
        else:
            raise ValueError(f"Method not defined for {self}")
        return methods

    def is_analysis_type(self):
        return self in (CandidateSearchType.REANALYSIS_NEW_ANNOTATION, )

    def is_classification_type(self):
        return self in (CandidateSearchType.CROSS_SAMPLE_CLASSIFICATION, CandidateSearchType.CLASSIFICATION_EVIDENCE_UPDATE)


class CandidateStatus(models.TextChoices):
    OPEN = 'O', 'OPEN'
    RESOLVED = 'R', 'RESOLVED'
    HIDDEN = 'D', 'HIDDEN'
    HIGHLIGHTED = 'L', 'HIGHLIGHTED'


class CandidateSearchVersion(TimeStampedModel):
    search_type = models.CharField(choices=CandidateSearchType.choices, max_length=1)
    code_version = models.IntegerField()
    celery_task_name = models.TextField()
    methods = models.TextField()

    class Meta:
        unique_together = ('search_type', 'code_version')

    def __str__(self) -> str:
        description = self.get_search_type_display()
        if self.code_version > 1:
            description += f"/v{self.code_version}"
        return description

    def get_type_base_page_url(self) -> str | None:
        url = None
        search_type = CandidateSearchType(self.search_type)
        if search_type.is_analysis_type():
            url = reverse("reanalysis_candidate_search")
        elif search_type.is_classification_type():
            url = reverse("classification_candidate_search")
        return url



class CandidateSearchRun(GuardianPermissionsAutoInitialSaveMixin, TimeStampedModel):
    search_version = models.ForeignKey(CandidateSearchVersion, on_delete=CASCADE)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    celery_task = models.CharField(max_length=36, null=True)
    status = models.CharField(max_length=1, choices=ProcessingStatus.choices, default=ProcessingStatus.CREATED)
    error_exception = models.TextField(null=True, blank=True)
    config_snapshot = models.JSONField(default=dict)
    git_hash = models.TextField()

    TASKS_AND_VERSIONS = {
        CandidateSearchType.REANALYSIS_NEW_ANNOTATION: ("analysis.tasks.reanalysis_tasks.ReAnalysisNewAnnotationTask", 1),
        CandidateSearchType.CROSS_SAMPLE_CLASSIFICATION: ("classification.tasks.classification_candidate_search_tasks.CrossSampleClassificationCandidateSearchTask", 1),
        CandidateSearchType.CLASSIFICATION_EVIDENCE_UPDATE: ("classification.tasks.classification_candidate_search_tasks.ClassificationEvidenceUpdateCandidateSearchTask", 1),
    }

    CANDIDATE_GRID_COLUMNS = {
        CandidateSearchType.REANALYSIS_NEW_ANNOTATION: ["sample__name", "analysis", "annotation_version", "zygosity"],
        CandidateSearchType.CROSS_SAMPLE_CLASSIFICATION: ["classification", "sample__name", "zygosity"],
        CandidateSearchType.CLASSIFICATION_EVIDENCE_UPDATE: ["classification", "annotation_version"],
    }

    def __str__(self):
        return f"{self.search_version}: {self.pk}"

    def is_running(self) -> bool:
        return self.status in ProcessingStatus.RUNNING_STATES

    def get_absolute_url(self) -> str:
        return reverse("view_candidate_search_run", kwargs={"pk": self.pk})

    @staticmethod
    def create_and_launch_job(user, search_type, config_snapshot: dict) -> 'CandidateSearchRun':
        celery_task_name, code_version = CandidateSearchRun.TASKS_AND_VERSIONS[search_type]
        methods = search_type.get_methods()
        search_version, _ = CandidateSearchVersion.objects.update_or_create(search_type=search_type,
                                                                            code_version=code_version,
                                                                            defaults={
                                                                                "methods": methods,
                                                                                "celery_task_name": celery_task_name
                                                                            })
        csr = CandidateSearchRun.objects.create(
            search_version=search_version,
            user=user,
            config_snapshot=config_snapshot,
            git_hash=get_cached_project_git_hash()
        )

        task = Signature(celery_task_name, args=(csr.pk, ))
        result = task.apply_async()
        CandidateSearchRun.objects.filter(pk=csr.pk).update(celery_task=result.id)
        return csr

    def get_zygosities_from_config(self) -> list[Zygosity]:
        """ Helper method if you have stored SampleCandidatesSearchForm in config_snapshot"""
        ZYG_NAMES = {
            "hom_ref": Zygosity.HOM_REF,
            "het": Zygosity.HET,
            "hom_alt": Zygosity.HOM_ALT,
        }

        zygosities = []
        for zyg, code in ZYG_NAMES.items():
            if self.config_snapshot.get(zyg):
                zygosities.append(code)

        return zygosities


class Candidate(TimeStampedModel):
    search_run = models.ForeignKey(CandidateSearchRun, on_delete=CASCADE)
    status = models.CharField(choices=CandidateStatus.choices, max_length=1, default=CandidateStatus.OPEN)
    notes = models.TextField(null=True, blank=True)
    evidence = models.JSONField(default=dict)
    reviewer = models.ForeignKey(User, null=True, blank=True, on_delete=SET_NULL)
    reviewer_comment = models.TextField(null=True, blank=True)
    variant = models.ForeignKey(Variant, null=True, blank=True, on_delete=CASCADE)
    classification = models.ForeignKey(Classification, null=True, blank=True, on_delete=CASCADE)
    analysis = models.ForeignKey(Analysis, null=True, blank=True, on_delete=CASCADE)
    annotation_version = models.ForeignKey(AnnotationVersion, null=True, blank=True, on_delete=CASCADE)
    # We can't point to ClinVar directly - so if you need to retrieve it
    sample = models.ForeignKey(Sample, null=True, on_delete=CASCADE)
    zygosity = models.CharField(choices=Zygosity.CHOICES, null=True, blank=True, max_length=1)

    class Meta:
        ordering = ('-created',)

    def get_clinvar(self) -> ClinVar|None:
        clinvar = None
        if self.variant and self.annotation_version and self.annotation_version.clinvar_version:
            clinvar = ClinVar.objects.filter(version=self.annotation_version.clinvar_version, variant=self.variant).first()
        return clinvar

    @staticmethod
    def get_permission_check(pk, user, write=False):
        """ Checks user has permission on CandidateSearchRun """
        c = Candidate.objects.get(pk=pk)
        c.search_run.check_permission(user, write)
        return c

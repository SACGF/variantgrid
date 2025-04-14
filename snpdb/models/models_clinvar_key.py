#####
# CLINVAR KEY
# All of this really should be moved to its own app, only put in snpdb so Labs could reference ClinVarKey
# which puts too much functionality into snpdb for optional functionality wrapped up in classifications
#####

from functools import cached_property
from typing import Optional
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied, ValidationError
from django.db.models import TextChoices, QuerySet, CASCADE
from django.db import models
from django_extensions.db.models import TimeStampedModel
from snpdb.models import Lab
import re

from library.utils import JsonObjType


class ClinVarAffectedStatus(TextChoices):
    # "yes", "no", "unknown", "not provided", "not applicable"
    yes = "yes"
    no = "no"
    unknown = "unknown"
    not_provided = "not provided"
    not_applicable = "not applicable"


class ClinVarCitationsModes(TextChoices):
    all = "all"
    interpretation_summary_only = "interpret"


class ClinVarCitationSource(models.TextChoices):
    PUBMED = 'PUBMED', 'PubMed'
    NCBI_BOOKSHELF = 'NCBI_BOOKSHELF', 'BookShelf'
    PUBMED_CENTRAL = 'PMC', 'pmc'

    @staticmethod
    def from_legacy_code(code: str) -> Optional['ClinVarCitationSource']:
        """
        Turns any kind of possible citation prefix into a CitationSource
        :param code: String representing PubMed, PubMedCentral or NCBI's Bookshelf
        :return: A CitationSource if valid, None otherwise
        """
        match code:
            case "P": return ClinVarCitationSource.PUBMED
            case "PMID": return ClinVarCitationSource.PUBMED
            case "PUBMED": return ClinVarCitationSource.PUBMED

            case "C": return ClinVarCitationSource.PUBMED_CENTRAL
            case "PMCID": return ClinVarCitationSource.PUBMED_CENTRAL
            case "PUBMEDCENTRAL": return ClinVarCitationSource.PUBMED_CENTRAL
            case "PMC": return ClinVarCitationSource.PUBMED_CENTRAL

            case "N": return ClinVarCitationSource.NCBI_BOOKSHELF
            case "NBK": return ClinVarCitationSource.NCBI_BOOKSHELF
            case "NCBIBOOKSHELF": return ClinVarCitationSource.NCBI_BOOKSHELF
            case "BOOKSHELF": return ClinVarCitationSource.NCBI_BOOKSHELF
            case "BOOKSHELF ID": return ClinVarCitationSource.NCBI_BOOKSHELF


class ClinVarKey(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar key"

    id = models.TextField(primary_key=True)
    api_key = models.TextField(null=True, blank=True)
    org_id = models.TextField(null=False, blank=True, default='')  # maybe this should have been the id?
    name = models.TextField(blank=True, default='')
    last_full_run = models.DateTimeField(null=True)

    default_affected_status = models.TextField(choices=ClinVarAffectedStatus.choices, null=True, blank=True)
    inject_acmg_description = models.BooleanField(blank=True, default=False)
    include_interpretation_summary = models.BooleanField(blank=True, default=True)
    citations_mode = models.TextField(choices=ClinVarCitationsModes.choices, default=ClinVarCitationsModes.all)

    @property
    def label(self) -> str:
        return self.name if self.name else self.id

    def __str__(self):
        return f"ClinVarKey ({self.label})"

    def __lt__(self, other):
        return self.id < other.id

    @staticmethod
    def clinvar_keys_for_user(user: User) -> QuerySet:
        """
        Ideally this would be on ClinVarKey but can't be due to ordering
        """
        if user.is_superuser:
            return ClinVarKey.objects.all().order_by('pk')

        labs = Lab.valid_labs_qs(user).filter(clinvar_key__isnull=False)
        if labs:
            return ClinVarKey.objects.filter(pk__in=labs.values_list('clinvar_key', flat=True)).order_by('pk')
        else:
            return ClinVarKey.objects.none()

    def check_user_can_access(self, user: User) -> None:
        """
        :throw PermissionDenied if user is not associated to ClinVarKey
        """
        if not user.is_superuser:
            allowed_clinvar_keys = ClinVarKey.clinvar_keys_for_user(user)
            if not allowed_clinvar_keys.filter(pk=self.pk).exists():
                raise PermissionDenied("User does not belong to a lab that uses the submission key")


class ClinVarKeyExcludePatternMode(TextChoices):
    EXCLUDE_IF_MATCH = 'X'
    EXCLUDE_IF_NO_MATCH = 'N'


class ClinVarKeyExcludePattern(TimeStampedModel):
    clinvar_key = models.ForeignKey(ClinVarKey, on_delete=CASCADE)
    evidence_key = models.TextField()  # Actually EvidenceKey but can't link to that from snpdb
    pattern = models.TextField()
    name = models.TextField(blank=True)
    case_insensitive = models.BooleanField(default=True, blank=True)
    mode = models.TextField(choices=ClinVarKeyExcludePatternMode.choices, default=ClinVarKeyExcludePatternMode.EXCLUDE_IF_MATCH)

    @cached_property
    def regex(self):
        flags = 0
        if self.case_insensitive:
            flags = re.RegexFlag.IGNORECASE
        return re.compile(self.pattern, flags=flags)

    def should_exclude(self, text: str) -> bool:
        matches = bool(self.regex.search(text))
        if self.mode == ClinVarKeyExcludePatternMode.EXCLUDE_IF_MATCH:
            return matches
        else:
            return not matches

    def clean(self):
        try:
            re.compile(self.pattern)
        except:
            raise ValidationError({'pattern': ValidationError(f'{self.pattern} is not a valid regular expression')})

        from classification.models import EvidenceKeyMap
        if EvidenceKeyMap.cached_key(self.evidence_key).is_dummy:
            raise ValidationError({'evidence_key': ValidationError(f'{self.evidence_key} is not a valid EvidenceKey')})

    def __str__(self):
        from classification.models import EvidenceKeyMap
        return EvidenceKeyMap.cached_key(self.evidence_key).pretty_label + " : " + (self.name or self.pattern)


class ClinVarExportTypeBucket(TextChoices):
    GERMLINE = "G", "Germline"
    ONCOGENIC = "O", "Oncogenic"
    CLINICAL_IMPACT = "C", "Clinical Impact"


class _Patterned(models.Model):
    pattern = models.TextField()

    @cached_property
    def _re_pattern(self):
        return re.compile(self.pattern, re.IGNORECASE)

    def matches(self, text: str) -> bool:
        if self.pattern:
            return bool(self._re_pattern.match(text))
        return False

    class Meta:
        abstract = True


class ClinVarExportAssertionMethod(_Patterned, TimeStampedModel):
    id = models.TextField(primary_key=True)
    label = models.TextField()
    reference_db = models.TextField(max_length=20, choices=ClinVarCitationSource.choices, null=True, blank=True)
    reference_id = models.TextField(null=True, blank=True)
    reference_url = models.TextField(null=True, blank=True)
    export_type = models.TextField(max_length=1, choices=ClinVarExportTypeBucket.choices)

    def clean(self):
        if not (self.reference_db and self.reference_id) and not self.reference_url:
            raise ValidationError("Must specify either reference_db & reference_id OR reference_url")

    @property
    def as_clinvar_json(self) -> JsonObjType:
        if url := self.reference_url:
            return {"url": url}
        else:
            return {"db": ClinVarCitationSource(self.reference_db).value, "id": self.reference_id}


class ClinVarExportAssertionMethodMapping(_Patterned, TimeStampedModel):
    clinvar_key = models.ForeignKey(ClinVarKey, on_delete=models.CASCADE)
    order = models.IntegerField(help_text="Lower order is more important")
    pattern = models.TextField()
    assertion_method = models.ForeignKey(ClinVarExportAssertionMethod, on_delete=models.CASCADE)

    class Meta:
        unique_together = (
                ("clinvar_key", "order"),
                ("clinvar_key", "pattern")
        )
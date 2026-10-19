"""
ClinVarKeyExcludePattern: a per-ClinVarKey regex over one evidence key that keeps matching (or non-matching)
classifications out of that key's ClinVar submissions. Applied by
classification/models/clinvar_export_exclude_utils.py:ClinVarExcludePatternUtil, edited as an inline on the
ClinVarKey admin (classification/admin/clinvar_export_admin.py).
"""
import re
from functools import cached_property
from re import RegexFlag

from django.core.exceptions import ValidationError
from django.db import models
from django.db.models import TextChoices
from django.db.models.deletion import CASCADE
from django_extensions.db.models import TimeStampedModel

from classification.models.clinvar_key_models import ClinVarKey
from classification.models.evidence_key import EvidenceKeyMap


class ClinVarKeyExcludePatternMode(TextChoices):
    EXCLUDE_IF_MATCH = 'X'
    EXCLUDE_IF_NO_MATCH = 'N'


class ClinVarKeyExcludePattern(TimeStampedModel):
    clinvar_key = models.ForeignKey(ClinVarKey, on_delete=CASCADE)
    evidence_key = models.TextField()  # EvidenceKey.key - not a FK as existing rows aren't guaranteed to be valid keys
    pattern = models.TextField()
    name = models.TextField(blank=True)
    case_insensitive = models.BooleanField(default=True, blank=True)
    mode = models.TextField(choices=ClinVarKeyExcludePatternMode.choices, default=ClinVarKeyExcludePatternMode.EXCLUDE_IF_MATCH)

    class Meta:
        db_table = "snpdb_clinvarkeyexcludepattern"

    @cached_property
    def regex(self):
        flags = 0
        if self.case_insensitive:
            flags = RegexFlag.IGNORECASE
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
        except Exception as exc:
            raise ValidationError({'pattern': ValidationError(f'{self.pattern} is not a valid regular expression')}) from exc

        if EvidenceKeyMap.cached_key(self.evidence_key).is_dummy:
            raise ValidationError({'evidence_key': ValidationError(f'{self.evidence_key} is not a valid EvidenceKey')})

    def __str__(self):
        return EvidenceKeyMap.cached_key(self.evidence_key).pretty_label + " : " + (self.name or self.pattern)

"""
ClinVarKey: one ClinVar submitter account (API key, org id and how records are converted for it). Labs submit under
the key on their classification/models/classification_settings_models.py:LabClassificationSettings, and
ClinVarAllele / ClinVarExportBatch (classification/models/clinvar_export_models.py) hang off it.
"""
import re
from re import RegexFlag
from typing import Optional

from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db import models
from django.db.models import QuerySet, TextChoices
from django_extensions.db.models import TimeStampedModel

from library.utils import JsonObjType
from snpdb.models.models import Lab


class ClinVarAssertionMethods(TextChoices):
    # "yes", "no", "unknown", "not provided", "not applicable"
    yes = "yes"
    no = "no"
    unknown = "unknown"
    not_provided = "not provided"
    not_applicable = "not applicable"


class ClinVarCitationsModes(TextChoices):
    all = "all"
    interpretation_summary_only = "interpret"


class ClinVarKey(TimeStampedModel):
    class Meta:
        verbose_name = "ClinVar key"

    id = models.TextField(primary_key=True)
    api_key = models.TextField(null=True, blank=True)
    org_id = models.TextField(null=False, blank=True, default='')  # maybe this should have been the id?
    name = models.TextField(blank=True, default='')
    last_full_run = models.DateTimeField(null=True)

    default_affected_status = models.TextField(choices=ClinVarAssertionMethods.choices, null=True, blank=True)
    inject_acmg_description = models.BooleanField(blank=True, default=False)
    include_interpretation_summary = models.BooleanField(blank=True, default=True)
    assertion_method_lookup = models.JSONField(null=False, default=dict)
    citations_mode = models.TextField(choices=ClinVarCitationsModes.choices, default=ClinVarCitationsModes.all)

    def assertion_criteria_vg_to_code(self, vg_value: str) -> Optional[JsonObjType]:
        """
        :param vg_value: Value as stored in the evidence key
        :return: The code we map to, None indicates we don't have a mapping
        """
        def map_value(raw_value: str):
            if not raw_value:
                raw_value = ""

            if self.assertion_method_lookup:
                if lookups := self.assertion_method_lookup.get("lookups"):
                    for lookup in lookups:
                        if match_text := lookup.get("match"):
                            if re.compile(match_text, RegexFlag.IGNORECASE).match(raw_value):
                                return lookup.get("citation")
                        else:
                            raise ValueError("Assertion Method Lookup, lacking a 'match' field")
                else:
                    raise ValueError("Assertion Method Lookup, needs a key 'lookups', with dict entries of 'match' and 'citation'")
            elif raw_value == "acmg":
                # if no method lookups have been setup, only accept "acmg"
                return raw_value
            else:
                raise ValueError(f"Unable to map {raw_value=}")

        mapped_value = map_value(vg_value)
        if mapped_value == "acmg":
            return {
                "db": "PubMed",
                "id": "PMID:25741868"
            }
        return mapped_value

    @property
    def label(self) -> str:
        return self.name if self.name else self.id

    def __str__(self):
        return f"ClinVarKey ({self.label})"

    def __lt__(self, other):
        return self.id < other.id

    @property
    def labs(self) -> QuerySet[Lab]:
        """ The labs that submit under this key """
        return Lab.objects.filter(classification_settings__clinvar_key=self).order_by('name')

    @staticmethod
    def clinvar_keys_for_user(user: User) -> QuerySet:
        if user.is_superuser:
            return ClinVarKey.objects.all().order_by('pk')
        labs = Lab.valid_labs_qs(user)
        return ClinVarKey.objects.filter(lab_classification_settings__lab__in=labs).distinct().order_by('pk')

    def check_user_can_access(self, user: User) -> None:
        """
        :throw PermissionDenied if user is not associated to ClinVarKey
        """
        if not user.is_superuser:
            allowed_clinvar_keys = ClinVarKey.clinvar_keys_for_user(user)
            if not allowed_clinvar_keys.filter(pk=self.pk).exists():
                raise PermissionDenied("User does not belong to a lab that uses the submission key")

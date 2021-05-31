from dataclasses import dataclass
from enum import Enum

from django.conf import settings
from django.contrib.auth.models import User
from django.db import models, transaction
from django.db.models.deletion import CASCADE
from django.db.models.query import QuerySet
import django.dispatch
from django.dispatch.dispatcher import receiver
from django.utils.timezone import now
from django_extensions.db.models import TimeStampedModel
from typing import List, Optional, Iterable

from lazy import lazy

from flags.models.models import FlagsMixin, FlagCollection, FlagTypeContext, \
    flag_collection_extra_info_signal, FlagInfos
from library.log_utils import report_message
from snpdb.models import Variant
from snpdb.models.models_variant import Allele
from classification.enums import ShareLevel, SpecialEKeys
from classification.enums.clinical_context_enums import ClinicalContextStatus
from classification.models.flag_types import classification_flag_types
from classification.models.classification import Classification, \
    ClassificationModification

clinical_context_signal = django.dispatch.Signal(providing_args=["clinical_context", "status", "is_significance_change", "cause"])

CS_TO_NUMBER = {
    'B': 1,
    'LB': 1,
    'VUS': 2,
    'VUS_A': 2,
    'VUS_B': 2,
    'VUS_C': 2,
    'LP': 3,
    'P': 3,
    'R': 4
}


class DiscordanceLevel(str, Enum):
    """
    Values are assigned based on how big the differences are considered when it comes to discordance
    e.g. difference between Benign and Likely Benign isn't as important as the difference between Likely Benign and VUS etc
    """

    NO_ENTRIES = 'no_entries'
    CONCORDANT_AGREEMENT = 'concordant_agreement'  # complete agreement
    CONCORDANT_CONFIDENCE = 'concordant_confidence'  # Benign vs Likely Benign
    DISCORDANT = 'discordant'

    @property
    def label(self) -> str:
        if self == DiscordanceLevel.CONCORDANT_AGREEMENT:
            return "Concordant (Agreement)"
        if self == DiscordanceLevel.CONCORDANT_CONFIDENCE:
            return "Concordant (Confidence)"
        if self == DiscordanceLevel.NO_ENTRIES:
            return "No Shared Submissions"
        if self == DiscordanceLevel.DISCORDANT:
            return "Discordant"
        return "Unknown"

    @property
    def bs_status(self) -> str:
        if self == DiscordanceLevel.CONCORDANT_AGREEMENT:
            return "success"
        if self == DiscordanceLevel.CONCORDANT_CONFIDENCE:
            return "warning"
        if self == DiscordanceLevel.NO_ENTRIES:
            return "secondary"
        return "danger"


@dataclass
class DiscordanceStatus:
    level: DiscordanceLevel
    lab_count: int

    @staticmethod
    def calculate(modifications: Iterable[ClassificationModification]) -> 'DiscordanceStatus':
        cs_scores = set()
        cs_values = set()
        labs = set()
        level: Optional[DiscordanceLevel]
        for vcm in modifications:
            labs.add(vcm.classification.lab)
            clin_sig = vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            if vcm.share_level_enum.is_discordant_level and vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE) is not None:
                strength = CS_TO_NUMBER.get(clin_sig)
                if strength:
                    cs_scores.add(strength)
                    cs_values.add(clin_sig)
        if len(cs_scores) > 1:
            level = DiscordanceLevel.DISCORDANT
        elif len(cs_values) > 1:
            level = DiscordanceLevel.CONCORDANT_CONFIDENCE
        elif len(labs) == 0:
            level = DiscordanceLevel.NO_ENTRIES
        else:
            level = DiscordanceLevel.CONCORDANT_AGREEMENT
        return DiscordanceStatus(level=level, lab_count=len(labs))


class ClinicalContext(FlagsMixin, TimeStampedModel):
    default_name = 'default'

    # variant = models.ForeignKey(Variant, on_delete=CASCADE)
    allele = models.ForeignKey(Allele, on_delete=CASCADE)
    name = models.TextField(null=True, blank=True)
    status = models.TextField(choices=ClinicalContextStatus.CHOICES, null=True, blank=True)
    last_evaluation = models.JSONField(null=True, blank=True)

    class Meta:
        unique_together = ('allele', 'name')

    def flag_type_context(self) -> FlagTypeContext:
        return FlagTypeContext.objects.get(pk='clinical_context')

    def calculate_status(self) -> str:
        """
        deprecated, try to use DiscordanceLevels instead
        """
        status = self.calculate_discordance_status()
        level = status.level
        if level == DiscordanceLevel.DISCORDANT:
            return ClinicalContextStatus.DISCORDANT
        else:
            return ClinicalContextStatus.CONCORDANT

    def calculate_discordance_status(self) -> DiscordanceStatus:
        return DiscordanceStatus.calculate(self.classification_modifications)

    @lazy
    def relevant_classification_count(self) -> int:
        return len([vcm for vcm in self.classification_modifications if CS_TO_NUMBER.get(vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))])

    def is_discordant(self):
        return self.status == ClinicalContextStatus.DISCORDANT

    @lazy
    def latest_report(self) -> Optional['DiscordanceReport']:
        from classification.models import DiscordanceReport
        return DiscordanceReport.latest_report(self)

    @transaction.atomic
    def recalc_and_save(self, cause: str):
        """
        Updates this ClinicalContext with the new status and applies flags where appropriate.
        """
        old_status = self.status
        new_status = self.calculate_status()

        is_significance_change = old_status != new_status

        self.last_evaluation = {
            "date": now().timestamp(),
            "trigger": cause,
            "considers_vcms": [vcm.id_str for vcm in self.classification_modifications],
            "old_status": old_status,
            "new_status": new_status
        }

        self.status = new_status

        self.save()

        clinical_context_signal.send(sender=ClinicalContext, clinical_context=self, status=new_status, is_significance_change=is_significance_change, cause=cause)
        # might be worth checking to see if the clinical context has changed against post signal

        if settings.DISCORDANCE_ENABLED:
            if new_status == ClinicalContextStatus.DISCORDANT:
                self.flag_collection_safe.get_or_create_open_flag_of_type(
                    flag_type=classification_flag_types.clinical_context_discordance,
                    reopen=True
                )

            else:
                self.flag_collection_safe.close_open_flags_of_type(
                    flag_type=classification_flag_types.clinical_context_discordance
                )

            if is_significance_change and (old_status or new_status == ClinicalContextStatus.DISCORDANT):
                report_message(f'Allele ID {self.allele_id} clinical grouping {old_status} -> {new_status}', extra_data={
                    'allele_id': self.allele_id
                })

        else:
            self.flag_collection_safe.close_open_flags_of_type(
                flag_type=classification_flag_types.clinical_context_discordance,
                comment='Discordance functionality has been disabled'
            )

    @property
    def is_default(self) -> bool:
        return self.name == ClinicalContext.default_name

    @staticmethod
    def default_group_for(variant: Variant) -> Optional['ClinicalContext']:
        if variant.allele:
            dg, _ = ClinicalContext.objects.get_or_create(allele=variant.allele, name=ClinicalContext.default_name)
            return dg
        print("Warning, classification variant doesn't have allele")
        return None

    @staticmethod
    def for_variant(variant: Variant) -> QuerySet:
        ClinicalContext.default_group_for(variant)  # ensure we always have the default group
        return ClinicalContext.objects.filter(allele=variant.allele).order_by('created')

    @property
    def classifications_qs(self) -> QuerySet:
        return Classification.objects.filter(clinical_context=self, share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS, withdrawn=False)

    @property
    def classifications_associated_qs(self) -> QuerySet:
        return Classification.objects.filter(clinical_context=self)

    @property
    def classification_modifications(self) -> List[ClassificationModification]:
        vcms = [vcm for vcm in (vc.last_published_version for vc in self.classifications_qs) if vcm is not None]
        return vcms

    def to_json(self) -> dict:
        return {
            'id': self.id,
            'name': self.name,
            'is_default': self.is_default
        }

    def __str__(self) -> str:
        if self.is_default:
            return 'Default for allele'
        return self.name


@receiver(flag_collection_extra_info_signal, sender=FlagCollection)
def get_extra_info(flag_infos: FlagInfos, user: User, **kwargs):  # pylint: disable=unused-argument
    ccs = ClinicalContext.objects.filter(flag_collection__in=flag_infos.ids).select_related('allele')
    for cc in ccs:
        flag_infos.set_extra_info(cc.flag_collection_id, {
            'label': f'{cc.allele} - "{cc.name}"',
            'cc_id': cc.id
        }, source_object=cc)

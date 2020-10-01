from django.conf import settings
from django.contrib.auth.models import User
from django.db import models, transaction
from django.db.models.deletion import CASCADE
from django.db.models.query import QuerySet
import django.dispatch
from django.dispatch.dispatcher import receiver
from django.utils.timezone import now
from django_extensions.db.models import TimeStampedModel
from typing import List, Optional

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
        min_score = None
        max_score = None
        has_vcs = False

        for vcm in self.classification_modifications:
            has_vcs = True
            cs = vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            # database jsonb requires string indexes only
            score = CS_TO_NUMBER.get(cs)
            if score is None:
                pass
            else:
                if min_score is None or score < min_score:
                    min_score = score
                if max_score is None or score > max_score:
                    max_score = score

        if not has_vcs:
            return ClinicalContextStatus.CONCORDANT
        if min_score is None or max_score is None or min_score == max_score:
            return ClinicalContextStatus.CONCORDANT
        return ClinicalContextStatus.DISCORDANT

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
def get_extra_info(flag_infos: FlagInfos, user: User, **kwargs):
    ccs = ClinicalContext.objects.filter(flag_collection__in=flag_infos.ids).select_related('allele')
    for cc in ccs:
        flag_infos.set_extra_info(cc.flag_collection_id, {
            'label': f'{cc.allele} - "{cc.name}"',
            'cc_id': cc.id
        }, source_object=cc)

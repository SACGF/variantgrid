from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from typing import Dict, List, Optional

from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.dispatch import receiver

from classification.models import ClassificationModification, EvidenceMixin, classification_flag_types, \
    classification_post_publish_signal, Classification
from flags.models import FlagStatus, Flag
from library.guardian_utils import admin_bot
from library.utils import DebugTimer
from snpdb.models import ClinVarKey, ClinVarKeyExcludePattern


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender,
              classification: Classification,
              previously_published: ClassificationModification,
              newly_published: ClassificationModification,
              user: User,
              debug_timer: DebugTimer,
              **kwargs):  # pylint: disable=unused-argument
    """
    On publishing of a classification, see if if it matches the criteria of the clinvar key
    and update it's ignore clinvar properties accordingly
    (will have no effect if ClinVarKey is not configured)
    """

    if clinvar_key := classification.lab.clinvar_key:
        ClinVarExcludePatternUtil(clinvar_key).ignore_pattern_for(newly_published).apply()

    debug_timer.tick("Reviewed ClinVar exclusion")


@dataclass(frozen=True)
class ExcludeStatus:
    """
    Quick class that describes what a record's current exclusion status is and what it is calculated it should be
    """
    old_exclude: bool
    new_exclude: bool
    manual_edit: bool

    def __str__(self):
        old_str = "Exclude" if self.old_exclude else "Consider"
        new_str = "Exclude" if self.new_exclude else "Consider"
        result: str
        if old_str == new_str:
            result = f"{old_str} (No Change)"
        else:
            result = f"{old_str} -> {new_str}"
        if self.manual_edit:
            result += " - but flag manually edited so not changing"
        return result


class ExcludeRecord:
    """
    Performs the calculation to see if a record should be excluded (and can update the flag status if asked)
    """

    def __init__(self, record: ClassificationModification, matches_ignores: List[ClinVarKeyExcludePattern]):
        self.record = record
        self.matches_ignores = matches_ignores

    @property
    def has_ignore_matches(self) -> bool:
        return bool(self.matches_ignores)

    @cached_property
    def existing_flag(self) -> Flag:
        return self.record.classification.flag_collection_safe.get_flag_of_type(classification_flag_types.classification_not_public, open_only=False)

    @cached_property
    def is_currently_ignored(self) -> bool:
        if flag := self.existing_flag:
            return flag.resolution.status == FlagStatus.OPEN
        return False

    @cached_property
    def existing_flag_manual(self) -> bool:
        """
        Is there an exclude from ClinVar flag, where the last action was taken by a user
        (If so, let's not muck around with it)
        """
        if flag := self.existing_flag:
            if last_comment := flag.flagcomment_set.order_by('-created').first():
                return last_comment.user != admin_bot()
        return False

    @property
    def status(self):
        return ExcludeStatus(
            old_exclude=self.is_currently_ignored,
            new_exclude=self.has_ignore_matches,
            manual_edit=self.existing_flag_manual
        )

    def apply(self):
        """
        Updates the exclude from ClinVar flag to the new status
        Unless existing_flag_manual in which case don't touch it
        """
        if self.existing_flag_manual or self.is_currently_ignored == self.has_ignore_matches:
            return

        elif not self.is_currently_ignored and self.has_ignore_matches:
            comment = "Not submitting to ClinVar as the record matched the following pattern(s):\n"
            comment += "\n".join(str(ip) for ip in self.matches_ignores)

            self.record.classification.flag_collection_safe.ensure_resolution(
                flag_type=classification_flag_types.classification_not_public,
                resolution="np_dont_share",
                comment=comment
            )

        elif self.is_currently_ignored and not self.has_ignore_matches:
            self.record.classification.flag_collection_safe.ensure_resolution(
                flag_type=classification_flag_types.classification_not_public,
                resolution="np_share",
                comment="Record no longer matches any ClinVar ignore patterns"
            )

        else:
            raise ValueError(f"Unable to handle status somehow {self.status}")


class ClinVarExcludePatternUtil:
    """
    Utility for updating one or all records (for a single ClinVarKey) exclude from ClinVarFlags
    """

    def __init__(self, clinvar_key: ClinVarKey):
        self.clinvar_key = clinvar_key
        self.exclude_patterns: List[ClinVarKeyExcludePattern] = list(clinvar_key.clinvarkeyexcludepattern_set.all())

    def matching_exclude_patterns(self, record: EvidenceMixin) -> Optional[List[ClinVarKeyExcludePattern]]:
        """
        Return list of exclude patterns that the record has failed
        """
        if not self.exclude_patterns:
            return None

        matching_patterns: List[ClinVarKeyExcludePattern] = []
        for pattern in self.exclude_patterns:
            value = record.get(pattern.evidence_key)
            if pattern.should_exclude(str(value)):
                matching_patterns.append(pattern)
        return matching_patterns

    def _all_classifications_for_key(self) -> QuerySet[ClassificationModification]:
        labs = self.clinvar_key.lab_set.all()
        return ClassificationModification.objects.filter(
            classification__lab__in=labs,
            is_last_published=True
        )

    def ignore_pattern_for(self, record: ClassificationModification) -> ExcludeRecord:
        return ExcludeRecord(record, self.matching_exclude_patterns(record))

    def run_all(self, apply: bool) -> Dict[ExcludeStatus, List[int]]:
        """
        Evaluate all the classifications for the ClinVarKey
        @param apply: If true update the flags, if false just report on what changes would take place
        """
        counter = defaultdict(list)
        for cm in self._all_classifications_for_key():
            ir = self.ignore_pattern_for(cm)
            counter[ir.status].append(cm.classification_id)
            if apply:
                ir.apply()

        # show most recent IDs first as probably most interesting thing for testing
        for ids in counter.values():
            ids.sort(reverse=True)

        return counter

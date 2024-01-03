from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from typing import Optional, Iterable, List, Union, Set, Dict, Tuple
from django.utils.safestring import SafeString
from classification.enums import CriteriaEvaluation
from library.utils import first


_UNDERSTOOD_NAMESPACES = {"acmg", "horak"}


@dataclass(frozen=True)
class CriteriaStrength:
    ekey: 'EvidenceKey'
    strength: Optional[str] = None

    @property
    def acmg_point(self) -> Optional[int]:
        if self.strength.endswith("X") and self.is_expected_direction:
            return CriteriaEvaluation.POINTS.get(self.ekey.default_crit_evaluation)
        return CriteriaEvaluation.POINTS.get(self.strength)

    @property
    def strength_direction(self) -> str:
        if not self.strength:
            return "N"
        else:
            return self.strength[0]

    @property
    def is_met(self):
        return CriteriaEvaluation.is_met(self.strength)

    @property
    def is_unspecified(self) -> bool:
        return self.strength and self.strength.endswith("X")

    def __hash__(self):
        return hash(self.ekey) + hash(self.strength)

    def __eq__(self, other):
        return self.ekey == other.ekey and self.strength == other.strength

    @property
    def is_default_strength(self):
        return self.strength == self.ekey.default_crit_evaluation

    @property
    def is_expected_direction(self) -> bool:
        if self.strength:
            return self.ekey.key[0].upper() == self.strength_direction
        else:
            return True

    @staticmethod
    def strength_suffix_for(strength: str, short: bool = False):
        if strength is None:
            return "NM"
        elif strength == "X":
            if short:
                return "X"
            else:
                return "unspecified"
        elif strength.endswith("X"):
            if short:
                return strength[0] + "_X"
            else:
                return strength[0] + "_unspecified"
        return strength

    @property
    def short(self):
        return format(self, 'short')

    def __format__(self, format_spec):
        # Make sure criteria are in camel case so removing spaces still leaves it readable
        pretty_label = self.ekey.pretty_label.replace(" ", "")
        suffix = self.strength
        if suffix:
            if not bool(self.ekey.namespace) or self.ekey.namespace in _UNDERSTOOD_NAMESPACES:
                if self.is_default_strength:
                    return pretty_label

            if self.ekey.crit_uses_points:
                score = self.acmg_point
                if score is not None:
                    suffix = f"{score}"
            elif self.is_expected_direction:
                suffix = suffix[1:]

        return f"{pretty_label}_{CriteriaStrength.strength_suffix_for(suffix, short=format_spec=='short')}"

    def __repr__(self):
        return format(self)

    @property
    def _sort_key(self) -> Tuple[int, str]:
        index = 0
        try:
            index = CriteriaEvaluation.ALL_STRENGTHS.index(self.strength)
        except:
            pass
        return index, self.ekey.pretty_label

    def __lt__(self, other) -> bool:
        return self._sort_key < other._sort_key


@dataclass(frozen=True)
class AcmgPointScore:
    points: Optional[int]
    has_criteria: bool
    is_acmg_standard: bool
    has_unspecified: bool = False

    def __bool__(self):
        return self.has_criteria

    def __lt__(self, other):
        return self.points < other.points

    def __str__(self):
        if not self.has_criteria:
            return "-"
        return f"{self.points}{'*' if not self.is_acmg_standard else ''}"

    @property
    def sort_string(self) -> str:
        adjusted_score = (self.points or 0) + 100
        return f"{adjusted_score:03}"

    @staticmethod
    def most_extreme_point_score(scores: Iterable['AcmgPointScore']) -> 'AcmgPointScore':
        """
        Return the ACMG Criteria with the largest absolute points (e.g. -4 is more extreme than 3)
        If any record with criteria has non-standard ACMG criteria, mark the result as non-standard
        :param scores: A list of scores to inspect for the most extreme score
        :return: The most extreme score, or NO_CRITERIA if none of the scores has_criteria
        """
        biggest_score: Optional[AcmgPointScore] = None
        lowest_score: Optional[AcmgPointScore] = None
        unspecified: Optional[AcmgPointScore] = None
        has_non_standard = False
        for score in scores:
            if score.has_unspecified:
                unspecified = score
            elif score.has_criteria:
                if not biggest_score or biggest_score < score:
                    biggest_score = score
                if not lowest_score or score < lowest_score:
                    lowest_score = score
                has_non_standard = has_non_standard or not score.is_acmg_standard

        if not biggest_score:
            if unspecified:
                return unspecified
            return AcmgPointScore.NO_CRITERIA

        most_extreme = biggest_score if abs(biggest_score.points) >= abs(lowest_score.points) else lowest_score
        return most_extreme


AcmgPointScore.NO_CRITERIA = AcmgPointScore(points=None, has_criteria=False, is_acmg_standard=True)


class CriteriaStrengths:

    def __init__(self, strengths: Iterable[CriteriaStrength], is_acmg_standard: bool = True):
        self.strength_map = {strength.ekey.key.lower(): strength for strength in strengths}
        self.is_acmg_standard = is_acmg_standard

    @property
    def strengths(self) -> Iterable[CriteriaStrength]:
        return self.strength_map.values()

    @property
    def strength_list_met(self) -> List[CriteriaStrength]:
        return [cs for cs in self.strengths if CriteriaEvaluation.is_met(cs.strength)]

    def summary_string(self, acmg_only: bool = True):
        def report_me(cs: CriteriaStrength):
            nonlocal acmg_only
            return CriteriaEvaluation.is_met(cs.strength) and (not acmg_only or not cs.ekey.namespace)

        return ", ".join(str(cs) for cs in sorted(self.strength_list_met))

    def summary_string_short(self) -> str:
        return ", ".join(f"{cs:short}" for cs in sorted(strength for strength in self.strength_list_met if not strength.ekey.namespace))

    @property
    def has_criteria(self):
        return any(s.is_met for s in self.strengths)

    def __getitem__(self, item) -> Union[None, CriteriaStrength, List[CriteriaStrength]]:

        # this bit is to not confuse Django templates by having it call getitem for regular properties
        if hasattr(self, item):
            return getattr(self, item)

        if isinstance(item, str):
            from classification.models import EvidenceKeyMap
            parts = [item]
            if '__' in item:
                parts = item.split('__')
            prefixed = []
            for part in parts:
                part = part.replace("_", ":")
                prefixed.append(part)

            return [self.strength_map.get(part.lower()) or CriteriaStrength(EvidenceKeyMap.cached_key(part), None) for part in prefixed]

    def __contains__(self, item) -> bool:
        if isinstance(item, str):
            return item.lower() in self.strength_map.keys()

    @property
    def has_non_standard_strengths(self):
        return any(CriteriaEvaluation.POINTS.get(s.strength) is None for s in self.strengths if s.is_met)

    @cached_property
    def acmg_point_score(self) -> AcmgPointScore:
        has_unspecified = False
        if not self.has_criteria:
            return AcmgPointScore.NO_CRITERIA
        points = 0
        for cs in self.strength_list_met:
            if cs.is_unspecified:
                has_unspecified = True
            if bit_score := cs.acmg_point:
                points += bit_score
        return AcmgPointScore(points=points, has_criteria=True, is_acmg_standard=self.is_acmg_standard, has_unspecified=has_unspecified)


@dataclass(frozen=True)
class CriteriaCompare:
    any_not_met: bool
    strength_counts: Dict[CriteriaStrength, int]

    @property
    def has_value(self):
        return bool(self.strength_counts)

    @property
    def direction(self) -> str:
        directions: Set[str] = set()
        for strength in self.strength_counts.keys():
            directions.add(strength.strength_direction)
        if len(directions) != 1:
            # multiple or no directions, use "N"
            return "N"
        else:
            return first(directions)

    @property
    def has_difference(self) -> bool:
        return (len(self.strength_counts.keys()) + (1 if self.any_not_met else 0)) > 1

    @property
    def has_difference_unspecified(self) -> bool:
        if not self.any_not_met:
            if len(self.strength_counts.keys()) == 2 and any(strength.is_unspecified for strength in self.strength_counts.keys()):
                return True
            # two unspecified strengths are still has_difference_unspecified as they may be different strengths (whole point of unspecified)
            if len(self.strength_counts.keys()) == 1 and first(self.strength_counts.keys()).is_unspecified and first(self.strength_counts.values()) > 1:
                return True
        return False

    @property
    def html(self):
        if not self.has_value:
            return ""

        hidden_text: str = ""
        icon_name: str
        style: str = ""
        title: str = ""

        icon_name = "fa-solid fa-circle"
        title = "This allele has the same strength across all records for this criteria."

        if self.has_difference_unspecified:
            icon_name = "fa-solid fa-circle-question"
            title = "This allele has an unspecified strength for the given criteria. Unable to determine if this represents a difference or not."
        elif self.has_difference:
            icon_name = "fa-solid fa-circle-half-stroke"
            title = "This allele has different strengths for the given criteria."

        if self.direction == "N":
            title = "This allele has both benign and pathogenic strengths for the given criteria."
            icon_name = "fa-solid fa-circle-half-stroke"
            style += "color:#d8d;"
            hidden_text = "!!"
        else:
            if self.direction == "P":
                style += "color:#d88;"
                hidden_text = "P"
            elif self.direction == "B":
                style += "color:#88f;"
                hidden_text = "B"

            if self.has_difference_unspecified:
                hidden_text += "?"
            elif self.has_difference:
                hidden_text += "!"

        return SafeString(f'<i class="{icon_name} hover-detail" style="{style}" title="{title}"></i><span style="opacity:0;width:0">{hidden_text}</span>')


class CriteriaSummarizer:
    def __init__(self, strengths: List[CriteriaStrengths]):
        self.strengths = strengths

    def __getitem__(self, item: str):
        strength_counts: Dict[CriteriaStrength, int] = defaultdict(int)
        any_not_met = False

        for strengths in self.strengths:
            # if we only want to compare when we have criteria
            values = strengths[item]
            if not isinstance(values, list):
                values = [values]
            if len(values) > 1 and (met_values := [value for value in values if value.is_met]):
                values = met_values
            for strength in values:
                if not strength.is_met:
                    any_not_met = True
                else:
                    strength_counts[strength] += 1

        return CriteriaCompare(
            any_not_met=any_not_met,
            strength_counts=strength_counts
        )

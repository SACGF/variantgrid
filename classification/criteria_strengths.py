from dataclasses import dataclass
from typing import Optional, Iterable, List, Any, Union

from lazy import lazy

from classification.enums import CriteriaEvaluation


@dataclass(frozen=True)
class CriteriaStrength:
    ekey: 'EvidenceKey'
    strength: Optional[str] = None

    @property
    def acmg_point(self) -> Optional[int]:
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
            if not self.ekey.namespace:
                if self.is_expected_direction and self.is_default_strength:
                    return pretty_label

            if self.is_expected_direction:
                suffix = suffix[1:]

        return f"{pretty_label}_{CriteriaStrength.strength_suffix_for(suffix, short=format_spec=='short')}"

    def __repr__(self):
        return format(self)

    @property
    def _sort_key(self) -> int:
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

    def __bool__(self):
        return self.has_criteria

    def __lt__(self, other):
        return self.points < other.points

    def with_non_standard(self) -> 'AcmgPointScore':
        return AcmgPointScore(points=self.points, has_criteria=self.has_criteria, is_acmg_standard=False)

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
        If any record with criteria has non standard ACMG criteria, mark the result as non-standard
        :param scores: A list of scores to inspect for the most extreme score
        :return: The most extreme score, or NO_CRITERIA if none of the scores has_criteria
        """
        biggest_score: Optional[AcmgPointScore] = None
        lowest_score: Optional[AcmgPointScore] = None
        has_non_standard = False
        for score in scores:
            if score.has_criteria:
                if not biggest_score or biggest_score < score:
                    biggest_score = score
                if not lowest_score or score < lowest_score:
                    lowest_score = score
                has_non_standard = has_non_standard or not score.is_acmg_standard

        if not biggest_score:
            return AcmgPointScore.NO_CRITERIA

        most_extreme = biggest_score if abs(biggest_score.points) >= abs(lowest_score.points) else lowest_score
        if has_non_standard:
            most_extreme = most_extreme.with_non_standard()
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

    @property
    def has_criteria(self):
        return any(s.is_met for s in self.strengths)

    def __getitem__(self, item) -> Union[None, CriteriaStrength, List[CriteriaStrength]]:
        if hasattr(self, item):
            return getattr(self, item)
        if isinstance(item, str):
            from classification.models import EvidenceKeyMap
            if '_' in item:
                parts = item.split('_')
                return [self.strength_map.get(part.lower()) or CriteriaStrength(EvidenceKeyMap.cached_key(part), None) for part in parts]

            return self.strength_map.get(item.lower()) or CriteriaStrength(EvidenceKeyMap.cached_key(item), None)

    def __contains__(self, item) -> bool:
        if isinstance(item, str):
            return item.lower() in self.strength_map.keys()

    @property
    def has_non_standard_strengths(self):
        return any(CriteriaEvaluation.POINTS.get(s.strength) is None for s in self.strengths if s.is_met)

    @lazy
    def acmg_point_score(self) -> AcmgPointScore:
        if not self.has_criteria:
            return AcmgPointScore.NO_CRITERIA
        points = 0
        for cs in self.strength_list_met:
            if bit_score := cs.acmg_point:
                points += bit_score
        return AcmgPointScore(points=points, has_criteria=True, is_acmg_standard=self.is_acmg_standard)

from dataclasses import dataclass
from typing import Optional, Iterable, List, Any

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
    def strength_suffix_for(strength: str):
        if strength is None:
            return "NM"
        elif strength == "X":
            return "unspecified"
        elif strength.endswith("X"):
            return strength[0] + "_unspecified"
        return strength

    @property
    def strength_suffix(self):
        suffix = self.strength
        if self.is_expected_direction and suffix:
            suffix = suffix[1:]

        return self.strength_suffix_for(suffix)

    def __str__(self) -> str:
        # Make sure criteria are in camel case so removing spaces still leaves it readable
        pretty_label = self.ekey.pretty_label.replace(" ", "")
        suffix = self.strength
        if suffix:
            matches_direction = False
            if not self.ekey.namespace:
                if self.is_expected_direction:
                    return pretty_label

            if self.is_expected_direction:
                suffix = suffix[1:]

        return f"{pretty_label}_{CriteriaStrength.strength_suffix_for(suffix)}"

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


class CriteriaStrengths:

    def __init__(self, strengths: Iterable[CriteriaStrength], source: Any):
        self.strength_map = {strength.ekey.key.lower(): strength for strength in strengths}
        self.source = source

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

    def __getitem__(self, item) -> Optional[CriteriaStrength]:
        if isinstance(item, str):
            from classification.models import EvidenceKeyMap
            return self.strength_map.get(item.lower()) or CriteriaStrength(EvidenceKeyMap.cached_key(item), None)

    def __contains__(self, item) -> bool:
        if isinstance(item, str):
            return item.lower() in self.strength_map.keys()

    @property
    def has_non_standard_strengths(self):
        return any(CriteriaEvaluation.POINTS.get(s.strength) is None for s in self.strengths if s.is_met)

    @property
    def acmg_point_score(self) -> int:
        score = 0
        for cs in self.strength_list_met:
            if bit_score := cs.acmg_point:
                score += bit_score
        return score

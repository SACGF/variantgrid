"""Convert the legacy values of a select / multi-select evidence key to its option keys.

Every write normalises a select value (Classification.process_entry), but records imported before that
can still hold an option's label ("Autosomal recessive") or a comma separated string of labels. This
rewrites the current evidence of internal-lab records as the sorted list of option keys - note, explain
and immutability kept - and republishes a record whose latest version was already published. History is
left as it was submitted; the new modification by admin_bot is the audit trail.

Dry run by default: prints the stored value -> normalised tally, per lab outcome counts and the values
that match no option. --apply writes. Re-running is a no-op for records already clean.
"""
import json
import logging
from collections import Counter, defaultdict
from enum import StrEnum
from typing import Any

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from classification.enums import EvidenceKeyValueType, SpecialEKeys, SubmissionSource
from classification.models.classification import Classification
from classification.models.condition_text_matching import ConditionTextMatch
from classification.models.evidence_key import EvidenceKey, EvidenceKeyMap
from library.django_utils.model_utils import ArrayLength
from library.guardian_utils import admin_bot

OPTION_VALUE_TYPES = (EvidenceKeyValueType.SELECT, EvidenceKeyValueType.MULTISELECT)


class Outcome(StrEnum):
    CLEAN = "clean"
    PUBLISHED = "normalised and published"
    UNSUBMITTED_EDITS = "normalised, not published (unsubmitted edits)"
    WITHDRAWN = "normalised, not published (withdrawn)"
    LEFT = "left (unmatched)"
    FAILED = "failed"


def normalise_option_value(e_key: EvidenceKey, value: Any) -> tuple[Any, list[str]]:
    """ (the value process_entry would store, the parts that stop it being stored) - the second is empty
        when every part matches an option, or the key allows custom values """
    normalised: set[str] = set()
    unmatched: list[str] = []
    options = e_key.virtual_options or []
    for part in Classification.option_value_parts(e_key, value):
        if option := Classification.match_option(options, str(part)):
            normalised.add(option["key"])
        elif e_key.allow_custom_values:
            normalised.add(str(part))
        else:
            unmatched.append(str(part))

    normalised_value: Any = sorted(normalised)
    if e_key.value_type != EvidenceKeyValueType.MULTISELECT:
        if len(normalised_value) > 1:
            unmatched.append("more than one value for a single select")
        elif normalised_value:
            normalised_value = normalised_value[0]
    return normalised_value, unmatched


class Command(BaseCommand):
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument("key", help="A select or multi-select evidence key, e.g. mode_of_inheritance")
        parser.add_argument("--apply", action="store_true", help="Write the normalised values (default: report only)")

    def handle(self, *args, **options):
        key = options["key"]
        apply = options["apply"]
        e_key = EvidenceKeyMap.cached_key(key)
        if e_key.is_dummy or e_key.value_type not in OPTION_VALUE_TYPES:
            raise CommandError(f"'{key}' is not a select or multi-select evidence key")

        user = admin_bot()
        value_tally: Counter[tuple[str, str]] = Counter()
        lab_outcomes: dict[str, Counter[Outcome]] = defaultdict(Counter)
        left: list[tuple[int, str, list[str]]] = []

        candidate_pks = Classification.objects.filter(lab__external=False, **{f"evidence__{key}__value__isnull": False}) \
                                              .order_by("pk").values_list("pk", flat=True)
        for pk in list(candidate_pks):
            classification = Classification.objects.select_related("lab").get(pk=pk)
            try:
                with transaction.atomic():
                    outcome = self._normalise(classification, key, user, apply, value_tally, left)
            except Exception:
                logging.exception("Classification %d: normalising %s failed", pk, key)
                outcome = Outcome.FAILED
            lab_outcomes[str(classification.lab)][outcome] += 1

        self._report(key, apply, value_tally, lab_outcomes, left)

    @staticmethod
    def _normalise(classification: Classification, key: str, user, apply: bool,
                   value_tally: Counter, left: list) -> Outcome:
        blob = classification.evidence[key]
        stored = blob["value"]
        normalised, unmatched = normalise_option_value(classification.evidence_keys.get(key), stored)
        if unmatched:
            left.append((classification.pk, json.dumps(stored), unmatched))
            return Outcome.LEFT
        if normalised == stored:
            return Outcome.CLEAN

        value_tally[(json.dumps(stored), json.dumps(normalised))] += 1
        if classification.withdrawn:
            outcome = Outcome.WITHDRAWN
        elif classification.has_outstanding_changes():
            outcome = Outcome.UNSUBMITTED_EDITS
        else:
            outcome = Outcome.PUBLISHED

        if apply:
            # VARIANT_GRID can change an API-immutable value; note / explain go in so the patch keeps them
            classification.patch_value({key: {"value": normalised, "note": blob.get("note"),
                                              "explain": blob.get("explain")}},
                                       user=user, source=SubmissionSource.VARIANT_GRID, save=True)
            if outcome == Outcome.PUBLISHED:
                classification.publish_latest(user=user)
        return outcome

    def _report(self, key: str, apply: bool, value_tally: Counter, lab_outcomes: dict[str, Counter],
                left: list[tuple[int, str, list[str]]]):
        write = self.stdout.write
        write(f"{key} - {'applied' if apply else 'dry run, nothing written'}")
        write("\nStored value -> normalised (count)")
        for (stored, normalised), count in value_tally.most_common():
            write(f"  {stored} -> {normalised}\t{count}")

        write("\nPer lab")
        totals: Counter[Outcome] = Counter()
        for lab_name in sorted(lab_outcomes):
            outcomes = lab_outcomes[lab_name]
            totals.update(outcomes)
            write(f"  {lab_name}: " + ", ".join(f"{o.value} {outcomes[o]}" for o in Outcome if outcomes[o]))
        write("  Total: " + ", ".join(f"{o.value} {totals[o]}" for o in Outcome))

        if left:
            write("\nLeft untouched - parts that match no option")
            for pk, stored, unmatched in left:
                write(f"  {pk}\t{stored}\t{unmatched}")

        if key == SpecialEKeys.MODE_OF_INHERITANCE:
            write(f"\nMOI-level ConditionTextMatch rows with a resolution on a non-normalised MOI: "
                  f"{self._dirty_moi_condition_matches()}")

    @staticmethod
    def _dirty_moi_condition_matches() -> int:
        """ Publishing moves a record under the MOI node for its normalised value, so a resolution made on the
            node for the legacy value does not follow it (ConditionTextMatch.sync_condition_text_classification) """
        e_key = EvidenceKeyMap.cached_key(SpecialEKeys.MODE_OF_INHERITANCE)
        moi_matches = ConditionTextMatch.objects.filter(classification__isnull=True, mode_of_inheritance__isnull=False) \
                                                .alias(xrefs_length=ArrayLength("condition_xrefs")) \
                                                .filter(xrefs_length__gt=0)
        return sum(1 for moi in moi_matches.values_list("mode_of_inheritance", flat=True)
                   if normalise_option_value(e_key, moi)[0] != moi)

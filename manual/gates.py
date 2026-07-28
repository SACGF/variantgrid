"""
Gates are named prerequisites a manual `manage` task can declare via
`ManualOperation(..., requires=["gate-name"])`. The migrator only auto-runs a task once all
of its gates are satisfied, so a step can't run before its prerequisite (e.g. an external
data upgrade) is done.

Gate *definitions* (the check code) live here, keyed by gate name. The *association* between a
command and its gates lives on the task in the DB (see ManualMigrationTask.requires), declared
at the ManualOperation call site.

A requires entry is resolved as:
  - "after:<task_id>"  -> AfterTaskGate (depends on another manual task completing)
  - a name in GATES    -> that registered gate
  - anything else      -> UnknownGate (always blocks, so a typo is noticed rather than skipped)
"""
from dataclasses import dataclass
from typing import Callable

from genes.models import TranscriptVersion
from manual.models import ManualGateSatisfied, ManualMigrationOutstanding, ManualMigrationTask
from ontology.models import OntologyVersion

AFTER_PREFIX = "after:"


@dataclass
class GateStatus:
    name: str
    satisfied: bool
    kind: str  # "auto" | "manual" | "after" | "unknown"
    description: str


class Gate:
    kind = "auto"

    def __init__(self, description: str):
        self.description = description

    def satisfied(self, name: str) -> bool:
        raise NotImplementedError()

    def status(self, name: str) -> GateStatus:
        return GateStatus(name=name, satisfied=self.satisfied(name), kind=self.kind,
                          description=self.description)


class AutoGate(Gate):
    """ Satisfied when a predicate over live models returns True. """
    kind = "auto"

    def __init__(self, check: Callable[[], bool], description: str):
        super().__init__(description)
        self.check = check

    def satisfied(self, name: str) -> bool:
        return bool(self.check())


class ManualGate(Gate):
    """ Satisfied once an operator confirms it via `manage.py manual_gate --satisfy <name>`. """
    kind = "manual"

    def satisfied(self, name: str) -> bool:
        return ManualGateSatisfied.is_satisfied(name)


class AfterTaskGate(Gate):
    """ Satisfied when another manual task is no longer outstanding (done, or never requested). """
    kind = "after"

    def __init__(self, task_id: str):
        super().__init__(f"After manual task: {ManualMigrationTask.describe_manual(task_id)}")
        self.task_id = task_id

    def satisfied(self, name: str) -> bool:
        task = ManualMigrationTask.objects.filter(pk=self.task_id).first()
        if task is None:
            return True  # never requested -> nothing to wait for
        return ManualMigrationOutstanding.outstanding_task(task) is None


class UnknownGate(Gate):
    kind = "unknown"

    def __init__(self):
        super().__init__("Unknown gate (typo in requires?) - blocking to be safe")

    def satisfied(self, name: str) -> bool:
        return False


# Gate definitions, keyed by name. AfterTaskGate is resolved dynamically from "after:<task_id>".
GATES: dict[str, Gate] = {
    "ontology-imported": AutoGate(
        lambda: OntologyVersion.latest(validate=False) is not None,
        "Ontology (MONDO/HPO/OMIM/HGNC) has been imported",
    ),
    "cdot-current": AutoGate(
        TranscriptVersion.data_is_current_cdot_format,
        "Transcripts are in current cdot format (run import_cdot_latest)",
    ),
    "variant-annotation-current": ManualGate(
        "Variant annotation upgraded to the latest version",
    ),
}


def resolve_gate(name: str) -> Gate:
    if name.startswith(AFTER_PREFIX):
        return AfterTaskGate(task_id=name[len(AFTER_PREFIX):])
    return GATES.get(name, UnknownGate())


def gate_status(name: str) -> GateStatus:
    return resolve_gate(name).status(name)


def blocked_by(requires: list[str]) -> list[str]:
    """ Return the subset of required gate names that are not currently satisfied. """
    return [name for name in (requires or []) if not resolve_gate(name).satisfied(name)]

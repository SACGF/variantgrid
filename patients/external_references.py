"""
How a client names a Patient / Specimen / Extraction, and how VG turns that into a row.

An externally managed object is found by its ExternalPK, which is unique on
(code, external_type, external_manager) - so a code needs its type to mean anything, and reaching one
takes both. An object that is not externally managed falls back to the local reference beside it,
which is what a bare string means. The caller supplies whichever it has, and each narrows the query
where present, so no deployment is forced into a canonical form it does not run.

Resolution has three outcomes rather than two. Nothing found is Pending, not an error: an extraction
legitimately arrives after the VCF that names it (@see patients.tasks.extraction_matching_tasks).
More than one found is Needs attention - reference_id is unique per parent rather than globally, so a
bare reference can match several rows, and that is ambiguity rather than a precedence rule.
"""
from dataclasses import dataclass
from typing import Optional

from django.contrib.auth.models import User
from django.db.models import Model, Q, QuerySet

from patients.models_enums import MatchStatus

REFERENCE_KEYS = frozenset({"code", "external_type", "external_manager", "reference_id"})
# Not a reference key - it says how the reference was arrived at, and only VG ever sets it
DERIVED_KEY = "derived"


class ExternalReferenceError(ValueError):
    """ A reference we can reject on shape alone, before any lookup """


@dataclass(frozen=True)
class ExternalReference:
    reference_id: Optional[str] = None      # the local reference, whatever the model calls it
    code: Optional[str] = None              # ExternalPK.code - meaningless without external_type
    external_type: Optional[str] = None     # which identifier scheme the code belongs to
    external_manager: Optional[str] = None  # ExternalModelManager.name
    # Inferred from a sample name rather than asserted by a client (@see the sample-name fallback in
    # upload.vcf.vcf_import). Takes no part in the query - it only says how much to trust the row
    derived: bool = False

    @classmethod
    def from_data(cls, data) -> 'ExternalReference':
        """ Accepts the bare-string and object forms - shape only, no database access """
        if isinstance(data, str):
            data = {"reference_id": data}
        if not isinstance(data, dict):
            raise ExternalReferenceError(f"Reference must be a string or an object, got {type(data).__name__}")
        if unknown := set(data) - REFERENCE_KEYS - {DERIVED_KEY}:
            raise ExternalReferenceError(f"Unknown reference key(s): {', '.join(sorted(unknown))}. "
                                         f"Accepted keys: {', '.join(sorted(REFERENCE_KEYS))}")
        values = {k: str(v) for k, v in data.items() if k in REFERENCE_KEYS and v is not None}
        reference = cls(derived=bool(data.get(DERIVED_KEY, False)), **values)
        if not (reference.reference_id or reference.code):
            raise ExternalReferenceError("Reference must supply 'reference_id' or 'code'")
        if bool(reference.code) != bool(reference.external_type):
            # ExternalPK is unique on (code, external_type, external_manager) - one without the other
            # names nothing, and is worth saying now rather than parking a row nobody can resolve
            raise ExternalReferenceError("'code' and 'external_type' must be supplied together")
        return reference

    @property
    def has_external_pk(self) -> bool:
        return bool(self.code)

    def as_json(self) -> dict:
        """ What gets parked on a row - round-trips back through from_data """
        data = {k: v for k, v in self.__dict__.items() if k in REFERENCE_KEYS and v is not None}
        if self.derived:
            data[DERIVED_KEY] = True
        return data

    def q(self, model: type[Model]) -> Q:
        q = Q()
        if self.code:
            q &= Q(external_pk__code=self.code, external_pk__external_type=self.external_type)
        if self.external_manager:
            q &= Q(external_pk__external_manager=self.external_manager)
        if self.reference_id:
            q &= Q(**{model.LOCAL_REFERENCE_FIELD: self.reference_id})
        return q

    def halves(self) -> list['ExternalReference']:
        """ The external and local halves on their own, for telling 'neither matched' apart from
            'each matched something different' """
        parts = []
        if self.has_external_pk:
            parts.append(ExternalReference(code=self.code, external_type=self.external_type,
                                           external_manager=self.external_manager))
        if self.reference_id:
            parts.append(ExternalReference(reference_id=self.reference_id))
        return parts

    def __str__(self):
        return ", ".join(f"{k}={v}" for k, v in sorted(self.as_json().items()))


@dataclass
class ResolvedReference:
    reference: ExternalReference
    status: str                     # MatchStatus
    obj: Optional[Model] = None
    error: Optional[str] = None

    @property
    def matched(self) -> bool:
        return self.status == MatchStatus.MATCHED


def resolve_reference(model: type[Model], reference: ExternalReference,
                      user: Optional[User] = None) -> ResolvedReference:
    """ user scopes the search through filter_for_user - a client can only attach objects it can see,
        and one it cannot reads as Pending rather than Needs attention """
    qs: QuerySet = model.filter_for_user(user) if user else model.objects.all()
    matches = list(qs.filter(reference.q(model))[:2])
    if len(matches) == 1:
        return ResolvedReference(reference, MatchStatus.MATCHED, obj=matches[0])
    if len(matches) > 1:
        return ResolvedReference(reference, MatchStatus.NEEDS_ATTENTION,
                                 error=f"{reference} matches more than one {model.__name__}")

    halves = reference.halves()
    if len(halves) > 1:
        # An ExternalPK and a local reference that each resolve, but to different rows. The combined
        # query just returns nothing, which reads as Pending - but this will never settle on its own
        found = {h: qs.filter(h.q(model)).first() for h in halves}
        distinct = {obj.pk for obj in found.values() if obj}
        if len(distinct) > 1:
            described = "; ".join(f"{h} -> {obj}" for h, obj in found.items() if obj)
            return ResolvedReference(reference, MatchStatus.NEEDS_ATTENTION,
                                     error=f"{reference} names more than one {model.__name__}: {described}")
    return ResolvedReference(reference, MatchStatus.PENDING,
                             error=f"No {model.__name__} found for {reference}")

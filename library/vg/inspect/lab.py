"""
`vg inspect lab <pk|group_name|name>`: a Lab - organization, group and members (heads marked),
classification counts by share level and withdrawn state, recent activity, ClinVar key, upload
location and the settings that shape its classifications.
"""
from typing import Any

from django.db.models import Count

from classification.models.classification import Classification
from classification.models.classification_grouping import ClassificationGrouping
from library.vg.inspect import capped, ref
from snpdb.models import Lab, LabHead


def load(key: str) -> Lab:
    key = key.strip()
    lab = None
    if key.isdigit():
        lab = Lab.objects.filter(pk=int(key)).first()
    lab = lab or Lab.objects.filter(group_name=key).first() or Lab.objects.filter(name__iexact=key).first()
    if lab is None:
        raise LookupError(f"No Lab with pk, group_name or name {key!r}")
    return lab


def inspect(key: str, depth: int) -> dict[str, Any]:
    lab = load(key)
    heads = set(LabHead.objects.filter(lab=lab).values_list("user_id", flat=True))
    members = list(lab.active_users.order_by("username"))
    classifications = Classification.objects.filter(lab=lab)
    data: dict[str, Any] = {
        "id": lab.pk,
        "name": lab.name,
        "group_name": lab.group_name,
        "organization": {"id": lab.organization_id, "name": lab.organization.name, "group_name": lab.organization.group_name},
        "external": lab.external,
        "research": lab.research,
        "location": ", ".join(p for p in (lab.city, lab.state.name if lab.state_id else None, lab.country.name if lab.country_id else None) if p) or None,
        "clinvar_key": lab.clinvar_key_id,
        "upload_location": lab.upload_location,
        "upload_automatic": lab.upload_automatic,
        "consolidates_variant_classifications": lab.consolidates_variant_classifications,
        "classification_config": bool(lab.classification_config),
        "members": capped(members, lambda u: u.username + (" (head)" if u.pk in heads else ""), cap=25),
        "classifications": {
            "total": classifications.count(),
            "withdrawn": classifications.filter(withdrawn=True).count(),
            "by_share_level": dict(classifications.values_list("share_level").annotate(n=Count("pk")).order_by("share_level")),
            "groupings": ClassificationGrouping.objects.filter(lab=lab).count(),
            "dirty_groupings": ClassificationGrouping.objects.filter(lab=lab, dirty=True).count(),
        },
        "url": lab.get_absolute_url(),
    }
    if depth >= 2:
        data["recent_classifications"] = capped(classifications.order_by("-modified"),
                                                lambda c: {**ref("classification", c, c.lab_record_id or ""), "modified": c.modified.date().isoformat()})
    return data

"""
`vg inspect user <username|pk>`: a User - flags, groups, the labs and organizations they belong
to, their settings (default build, columns) and what they own: VCFs, analyses, classifications,
gene lists, plus last login and recent page views.
"""
from typing import Any

from django.contrib.auth.models import User

from analysis.models.models_analysis import Analysis
from classification.models.classification import Classification
from eventlog.models import ViewEvent
from genes.models.models_gene_list import GeneList
from library.vg.inspect import capped, ref
from library.vg.inspect.common import lab_summary
from snpdb.models import VCF, Lab, UserSettings


def load(key: str) -> User:
    user = User.objects.filter(username=key).first() or (User.objects.filter(pk=int(key)).first() if key.isdigit() else None)
    if user is None:
        raise LookupError(f"No User {key!r}")
    return user


def inspect(key: str, depth: int) -> dict[str, Any]:
    user = load(key)
    labs = list(Lab.valid_labs_qs(user).select_related("organization").order_by("group_name"))
    user_settings = UserSettings.get_for_user(user)
    data: dict[str, Any] = {
        "id": user.pk,
        "username": user.username,
        "email": user.email,
        "name": user.get_full_name() or None,
        "active": user.is_active,
        "superuser": user.is_superuser,
        "staff": user.is_staff,
        "joined": user.date_joined.date().isoformat(),
        "last_login": user.last_login.strftime("%Y-%m-%d %H:%M") if user.last_login else None,
        "groups": list(user.groups.order_by("name").values_list("name", flat=True)),
        "labs": [{**lab_summary(lab), "organization": lab.organization.name} for lab in labs],
        "settings": {"default_build": user_settings.default_genome_build.name if user_settings.default_genome_build else None,
                     "lab": str(user_settings.get_lab()) if labs else None},
        "owns": {"vcfs": VCF.objects.filter(user=user).count(), "analyses": Analysis.objects.filter(user=user).count(),
                 "classifications": Classification.objects.filter(user=user).count(), "gene_lists": GeneList.objects.filter(user=user).count()},
    }
    if depth >= 2:
        data["recent_vcfs"] = capped(VCF.objects.filter(user=user).order_by("-pk"), lambda v: {**ref("vcf", v, v.name or ""), "date": v.date.date().isoformat()})
        data["recent_analyses"] = capped(Analysis.objects.filter(user=user).order_by("-pk"), lambda a: ref("analysis", a, a.name))
        data["recent_views"] = capped(ViewEvent.objects.filter(user=user).order_by("-pk"),
                                      lambda e: f"{e.created.strftime('%Y-%m-%d %H:%M')} {e.view_name} {e.args or ''}"[:160], cap=5)
    return data

import re
from datetime import datetime
from typing import Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.db import transaction
from django.http.request import HttpRequest
from django.shortcuts import render, redirect, get_object_or_404
from django.urls import reverse
from requests.models import Response

from classification.models import ClassificationRef, ClinicalContextRecalcTrigger
from classification.models.allele_overlap import AlleleOverlap, OverlapCounts, OverlapSet
from classification.models.classification import Classification
from classification.models.clinical_context_models import ClinicalContext
from classification.models.flag_types import classification_flag_types
from library.django_utils import require_superuser
from snpdb.models import Allele, Lab
from snpdb.models.models_variant import Variant


def view_overlaps(request: HttpRequest, lab_id: Optional[int] = None) -> Response:
    user = request.user
    user: User = request.user

    all_labs = list(Lab.valid_labs_qs(request.user, admin_check=True))
    if len(all_labs) == 1 and not lab_id:
        return redirect(reverse('view_overlaps', kwargs={'lab_id': all_labs[0].pk}))

    context = {
        "selected_lab": lab_id or 0,
        "labs": all_labs
    }

    return render(request, "classification/overlaps.html", context)


def view_overlaps_detail(request: HttpRequest, lab_id: Optional[int] = None) -> Response:
    user = request.user
    allele_and_vcs = AlleleOverlap.overlaps_for_user(user, lab_id = lab_id)
    overlap_counts = OverlapCounts(allele_and_vcs)

    context = {
        "overlap_sets": OverlapSet.as_sets(allele_and_vcs),
        "overlap_counts": overlap_counts
    }

    return render(request, "classification/overlaps_detail.html", context)


# POST only

@transaction.atomic()
def post_clinical_context(request: HttpRequest) -> Response:

    came_from_variant_id = request.POST.get('variant')
    came_from_allele_id = request.POST.get('allele')
    came_from_obj = None
    if came_from_variant_id:
        came_from_obj = Variant.objects.get(pk=came_from_variant_id)
    else:
        came_from_obj = Allele.objects.get(pk=came_from_allele_id)

    ccre = re.compile('cc_([0-9]+)')
    affected_ccs = set()
    for key, cc_name in request.POST.items():
        match = ccre.match(key)
        if match:
            cc_name = cc_name.strip()
            vcid = int(match[1])
            vc = Classification.objects.get(pk=vcid)
            old_cc = vc.clinical_context

            if old_cc is None or old_cc.name != cc_name:
                if old_cc:
                    affected_ccs.add(old_cc)
                updated_cc, _ = ClinicalContext.objects.get_or_create(allele=vc.variant.allele, name=cc_name)
                affected_ccs.add(updated_cc)
                vc.clinical_context = updated_cc
                vc.save()

                old_cc_name = None
                if old_cc:
                    old_cc_name = old_cc.name

                vc.flag_collection_safe.add_flag(
                    flag_type=classification_flag_types.classification_clinical_context_changed,
                    user=request.user,
                    comment=f'Clinical grouping changed from "{old_cc_name}" to "{cc_name}"',
                    permission_check=False
                )

    for cc in affected_ccs:
        cc.recalc_and_save(cause='Record(s) moved between clinical groupings', cause_code=ClinicalContextRecalcTrigger.CLINICAL_GROUPING_SET)

    #pretend to do something
    return redirect(came_from_obj)


@require_superuser
def view_clinical_context(request: HttpRequest, pk: int) -> Response:
    cc = get_object_or_404(ClinicalContext, pk=pk)
    user = request.user
    context = {"cc": cc, "timezone": settings.TIME_ZONE}
    last_evaluation = cc.last_evaluation
    context["has_details"] = True
    if last_evaluation:
        context["date"] = datetime.fromtimestamp(last_evaluation.get('date'))
        for copy_me in ["trigger", "old_status", "new_status"]:
            context[copy_me] = last_evaluation.get(copy_me)
        considers_vcms = last_evaluation.get('considers_vcms')
        if considers_vcms is not None:
            try:
                records = [ClassificationRef.init_from_str(user=user, id_str=vcm_id).modification for vcm_id in considers_vcms]
                context["records"] = records
            except:
                context["records_error"] = f"Could not resolve ids to classifications {', '.join(considers_vcms)}"
    else:
        context["date"] = cc.modified
        context["trigger"] = "Not recalculated since upgrade"
        context["old_status"] = "?"
        context["new_status"] = cc.status
        context["records_error"] = "Unknown"

    return render(request, 'classification/clinical_context.html', context)

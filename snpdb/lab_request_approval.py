"""
Approving a LabRequest (snpdb/models/models_lab_request.py): approve_lab_request turns the requester's private lab
into the requested lab in place, so everything they uploaded or classified there comes with it.
Lives outside snpdb.models because it fixes classification permissions, and snpdb's models sit below classification's.
"""
from typing import Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.db import transaction

from classification.models.classification import Classification
from classification.models.classification_grouping import ClassificationGrouping
from classification.models.overlaps_model import OverlapContributionNextStep
from classification.services.overlaps_services import OverlapServices
from snpdb.enums import LabRequestStatus
from snpdb.models.models import Lab, LabHead, Organization
from snpdb.models.models_lab_request import (
    LabRequest,
    get_user_private_lab,
    group_name_slug,
    lab_group_name,
)


@transaction.atomic
def approve_lab_request(lab_request: LabRequest, reviewed_by: User, organization: Optional[Organization],
                        organization_name: str, lab_name: str) -> Lab:
    """ Moves the requester's private lab into organization (created from organization_name when None) and renames
        it lab_name, making the requester its head.

        The guardian Group row is renamed rather than replaced: object permissions and membership hang off its pk.
        Lab.group looks the group up by name, so it is fetched before Lab.group_name changes """
    user = lab_request.user
    lab = get_user_private_lab(user)
    if lab is None:
        raise ValueError(f"{user} has no private lab to turn into '{lab_name}'")

    if organization is None:
        organization = Organization.objects.create(name=organization_name, group_name=group_name_slug(organization_name))

    group = lab.group
    group.name = lab_group_name(organization.group_name, lab_name)
    group.save()

    lab.name = lab_name
    lab.group_name = group.name
    lab.organization = organization
    lab.url = lab_request.website
    lab.contact_email = lab_request.email
    lab.save()

    organization.add_member(user)
    LabHead.objects.get_or_create(lab=lab, user=user)
    _fix_classifications(lab)

    lab_request.set_reviewed(LabRequestStatus.APPROVED, reviewed_by)
    lab_request.lab = lab
    lab_request.save()

    transaction.on_commit(lambda: lab_request.notify_requester(
        f"{settings.SITE_NAME} lab request approved",
        f"Your lab [{lab}]({lab.get_absolute_url()}) is set up, and you are its lab head. Everything you had in "
        f"your private lab is now in it. Add colleagues from the lab's Members tab."))
    return lab


def _fix_classifications(lab: Lab):
    """ Organisation-level shares go to Lab.group_institution, which comes from the group name, so they move to the
        new organisation. Overlap snapshots hold lab group names: recalculating them without notifications stops the
        rename looking like a lab joining each overlap """
    for classification in Classification.objects.filter(lab=lab):
        classification.fix_permissions(fix_modifications=True)

    ClassificationGrouping.objects.filter(lab=lab).update(dirty=True)

    overlap_ids = set(OverlapContributionNextStep.objects.filter(
        contribution__classification_grouping__lab=lab).values_list("overlap_id", flat=True))
    with OverlapServices.discordance_notifications_suppressed():
        OverlapServices.recalc_overlaps_batch(overlap_ids)

from django.core.management import BaseCommand
from django.db.models import Q

from classification.enums import SubmissionSource
from classification.models import Classification, SpecialEKeys
from library.guardian_utils import admin_bot
from library.utils import Value


class Command(BaseCommand):
    @staticmethod
    def _update_allele_origin(user, qs, old_value, allele_origin_value, allele_origin_confirmation_value):
        note = f"Auto set from '{old_value}' (legacy)"

        patch = {
            SpecialEKeys.ALLELE_ORIGIN: {"value": allele_origin_value, "note": note},
            'allele_origin_confirmation': {"value": allele_origin_confirmation_value, "note": note},
        }

        classification_ids = []
        for classification in qs:
            classification_ids.append(classification.pk)

            # keep track of if there are outstanding changes before we patch the value
            # as we don't want to publish changes a user hasn't published yet
            has_outstanding_changes = classification.has_outstanding_changes()
            classification.patch_value(
                patch=patch,
                user=user,
                source=SubmissionSource.FORM, # set this to form so as to not make it immutable
                # leave_existing_values=True, # safety just in case you set this on a record that already has allele origin
                save=True
            )
            if not has_outstanding_changes:
                # by not providing share level, it will default to publishing at the existing share level
                classification.publish_latest(user=user)

        if classification_ids:
            print(f"Updated allele_origin from {old_value} => {allele_origin_value}")
            print(", ".join(classification_ids))


    def handle(self, *args, **options):
        raise ValueError("Don't run this - haven't confirmed it's the right thing to do yet. "
                         "See https://github.com/SACGF/variantgrid_private/issues/2926")

        user = admin_bot()
        # Only alter ones that are from this system (not external)
        qs_internal = Classification.objects.filter(lab__external=False)

        Q_OLD_VALUE_NEW_VALUE = [
            (Q(evidence__allele_origin__isnull=True), "blank", "germline"),
            (Q(evidence__allele_origin__value="likely_germline"), "likely_germline", "germline"),
            # It was decided to handle these manually - https://github.com/SACGF/variantgrid_private/issues/2926
            #(Q(evidence__allele_origin__value="unknown"), "unknown", "germline"),
            #(Q(evidence__allele_origin__value="other"), "other", "germline"),
            #(Q(evidence__allele_origin__value="tested_inconclusive"), "tested_inconclusive", "germline"),
            (Q(evidence__allele_origin__value="likely_somatic"), "likely_somatic", "somatic"),
        ]

        for q, old_value, allele_origin_value in Q_OLD_VALUE_NEW_VALUE:
            qs = qs_internal.filter(q)
            self._update_allele_origin(user, qs,
                                       old_value=old_value,
                                       allele_origin_value=allele_origin_value,
                                       allele_origin_confirmation_value="unconfirmed")


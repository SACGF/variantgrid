from typing import Optional

from django.core.management import BaseCommand

from classification.enums import SubmissionSource
from classification.models import ClassificationModification, EvidenceKeyMap
from library.guardian_utils import admin_bot
from snpdb.models import Lab


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--key', type='str', required=True)
        parser.add_argument('--lab', type=str, required=False)
        parser.add_argument('--apply', action="store_true")

    def handle(self, *args, **options):
        lab_id: str = options["lab"]
        lab: Optional[Lab] = None
        if lab_id:
            if lab_id.isnumeric():
                lab = Lab.objects.get(id=lab_id)
            else:
                lab = Lab.objects.get(group_name=lab_id)
            print(f"Lab = {str(lab)}")

        qs = ClassificationModification.objects.filter(is_last_published=True, is_last_edited=True)
        if lab:
            qs = qs.filter(classification__lab=lab)

        valid_evidence_keys = set(k.key for k in EvidenceKeyMap.cached().all_keys)

        apply = options["apply"]

        mod: ClassificationModification
        for mod in qs.iterator():
            if bonus_keys := set(mod.published_evidence.keys()).difference(valid_evidence_keys):
                print(f"{mod.classification.pk} has illegal evidence keys = {bonus_keys}")

                if apply:
                    patch = {k: None for k in bonus_keys}
                    classification = mod.classification
                    classification.patch_value(patch=patch,
                                               source=SubmissionSource.API,
                                               save=True,
                                               user=admin_bot())
                    classification.publish_latest(user=admin_bot())
                    print("Updated")

from django.core.management import BaseCommand

from classification.models import Classification
from library.guardian_utils import admin_bot


def remove_evidence_key(key: str):
    # could either try to wipe it from history, but taking the safer route of just wiping out the value
    key_in_evidence = {f"evidence__{key}__isnull": False}
    affected_qs = Classification.objects.filter(**key_in_evidence)

    print(f"Number of keys with evidence key {key} = {affected_qs.count()}")
    fixed_records = 0

    user = admin_bot()
    for c in affected_qs:
        c.patch_value(
            patch={
                key: None
            },
            user=user,
            save=True
        )
        c.revalidate(user=user)
        c.publish_latest(user=user)
        fixed_records += 1
        if fixed_records % 10 == 0:
            print(f"Fixed {fixed_records}")

    print(f"Complete")


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--key', type=str, default=False)

    def handle(self, *args, **options):
        key = options["key"]
        remove_evidence_key(key)

from django.core.management import BaseCommand
from django.db.models import QuerySet
from classification.classification_import import reattempt_variant_matching
from classification.models import Classification
from library.guardian_utils import admin_bot
import re


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--file', type=str, required=True)

    def handle(self, *args, **options):
        filename = options["file"]
        qs: QuerySet[Classification] = Classification.objects.none()

        with open(filename) as file:
            all_ids = set(int(cid) for cid in re.compile(r"\d+").findall(file.read()))
            print(f"Found {len(all_ids)} IDs in file")
            qs = Classification.objects.filter(pk__in=all_ids)
            if qs.count() != len(all_ids):
                for c in qs:
                    all_ids.remove(c.pk)

                print(f"Number of classification IDs not present in the database = {len(all_ids)}")
                print(f"Example missing IDs = {list(all_ids)[0:5]}")
                return

        if qs:
            user = admin_bot()
            reattempt_variant_matching(user, qs)
            print("Re-matching has been queued")

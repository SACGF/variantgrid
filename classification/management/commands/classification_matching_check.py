from collections import Counter

from django.core.management import BaseCommand
from django.db.models import Max

from classification.models import Classification


class Command(BaseCommand):
    """
        Classifications are matched to variants, and adding/changing transcript data (cdot), or matching algorithms
        can alter which transcript/variant is used.

        This can be re-triggered, but that may alter the classification's variant. This tool allows checking without
        altering the classification
    """

    def _get_last_modified(self):
        qs = Classification.objects.all()
        data = qs.aggregate(Max("modified"))
        return data["modified__max"]

    def handle(self, *args, **options):
        start_last_modified = self._get_last_modified()

        diff = Counter()
        for c in Classification.objects.all():
            existing_tuple = tuple()
            if v := c.variant:
                existing_tuple = v.as_tuple()
            variant_tuple = c.get_variant_coordinates_from_evidence(update_flags=False)
            if existing_tuple == variant_tuple:
                diff["no change"] += 1
            elif not existing_tuple:
                diff["gained match"] += 1
            elif not variant_tuple:
                diff["lost match"] += 1
                print(f"{c.pk} old: {existing_tuple} -> {variant_tuple}")
            else:
                diff["match changed"] += 1
                print(f"{c.pk} old: {existing_tuple} -> {variant_tuple}")

        print(diff)
        end_last_modified = self._get_last_modified()
        if start_last_modified != end_last_modified:
            print(f"Beware - looks like someone edited a classification!!! {start_last_modified=} vs {end_last_modified=}")


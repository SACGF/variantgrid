import bisect
from functools import total_ordering

from django.core.management import BaseCommand

from classification.enums import SpecialEKeys
from classification.models import Classification


@total_ordering
class ConversionSize:

    def __init__(self, vc: Classification):
        self.ref_length = vc.update_cached_c_hgvs()
        self.vc_id = vc.id
        self.chgvs = vc.get(SpecialEKeys.C_HGVS)
        vc.save()

    def __lt__(self, other):
        return self.ref_length < other.ref_length


class Command(BaseCommand):

    def handle(self, *args, **options):
        conversions = []
        update_count = 0

        for vc in Classification.objects.all():
            conversion = ConversionSize(vc)

            bisect.insort(conversions, conversion)
            if len(conversions) > 10:
                conversions.pop(0)

            update_count += 1
            if update_count % 100 == 0:
                print(f"Completed {update_count}")
        print(f"Bulk Update of Cached c.hgvs - completed")
        print(f"Biggest ref lengths are:")
        for conversion in conversions[::-1]:
            print(f"{conversion.ref_length} from vc.id {conversion.vc_id} {conversion.chgvs}")

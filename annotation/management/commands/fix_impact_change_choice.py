from django.core.management.base import BaseCommand
from django.db.models import Case, When, Value

from analysis.models import DamageNode
from annotation.models import VariantAnnotationVersion, VariantAnnotation, VariantTranscriptAnnotation


class Command(BaseCommand):
    """ Do this as a management command not migration so its not in a transaction (which got too big) """

    def handle(self, *args, **options):
        DamageNode.objects.update(impact_min=self._get_case("impact_min"))

        impact_case = self._get_case("impact")
        for vav in VariantAnnotationVersion.objects.all():
            print(f"Updating {vav} VariantAnnotation... (may take a while)")
            VariantAnnotation.objects.filter(version=vav).update(impact=impact_case)
            print(f"Updating {vav} VariantTranscriptAnnotation... (may take a really long while)")
            VariantTranscriptAnnotation.objects.filter(version=vav).update(impact=impact_case)

    @staticmethod
    def _get_case(field_name) -> Case:
        IMPACT_OLD_NEW = {
            "O": "1",
            "L": "2",
            "M": "3",
            "H": "4",
        }
        return Case(*[When(**{field_name: k}, then=Value(v)) for k, v in IMPACT_OLD_NEW.items()])

from django.core.management.base import BaseCommand
from django.db.models import Case, When, Value

from analysis.models import DamageNode
from annotation.models import VariantAnnotationVersion, VariantAnnotation, VariantTranscriptAnnotation


class Command(BaseCommand):
    """ Do this as a management command not migration so its not in a transaction (which got too big) """
    def add_arguments(self, parser):
        parser.add_argument('--small_updates', action='store_true', help="Do update in small chunks (to reduce transaction size)")

    def handle(self, *args, **options):
        small_updates = options["small_updates"]
        IMPACT_OLD_NEW = {
            "O": "1",
            "L": "2",
            "M": "3",
            "H": "4",
        }

        DamageNode.objects.update(impact_min=self._get_case(IMPACT_OLD_NEW, "impact_min"))

        if small_updates:
            for vav in VariantAnnotationVersion.objects.all():
                for old, new in IMPACT_OLD_NEW.items():
                    print(f"Updating {vav} VariantAnnotation {old} -> {new}")
                    VariantAnnotation.objects.filter(version=vav, impact=old).update(impact=new)
                for old, new in IMPACT_OLD_NEW.items():
                    print(f"Updating {vav} VariantTranscriptAnnotation... {old} -> {new}")
                    VariantTranscriptAnnotation.objects.filter(version=vav, impact=old).update(impact=new)
        else:
            # Do as 1 big case statement
            impact_case = self._get_case(IMPACT_OLD_NEW, "impact")
            for vav in VariantAnnotationVersion.objects.all():
                print(f"Updating {vav} VariantAnnotation... (may take a while)")
                VariantAnnotation.objects.filter(version=vav).update(impact=impact_case)
                print(f"Updating {vav} VariantTranscriptAnnotation... (may take a really long while)")
                VariantTranscriptAnnotation.objects.filter(version=vav).update(impact=impact_case)

    @staticmethod
    def _get_case(old_new, field_name) -> Case:
        return Case(*[When(**{field_name: k}, then=Value(v)) for k, v in old_new.items()])

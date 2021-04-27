from django.core.management import BaseCommand

from classification.management.commands.evidence_key_to_unit import EvidenceKeyToUnit
from classification.models import Classification


class Command(BaseCommand):

    def handle(self, *args, **options):
        print("Reviewing every classification, this may take some time")
        EvidenceKeyToUnit(key_names=["allele_frequency"]).migrate(qs=Classification.objects.all())
        print("Done")
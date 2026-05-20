from django.core.management.base import BaseCommand

from annotation.models.models_phenotype_match import TextPhenotype, PatientTextPhenotype
from annotation.phenotype_matching import bulk_patient_phenotype_matching


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--clear', action='store_true')
        parser.add_argument('--cores', type=int, default=1,
                            help='Number of parallel workers for sentence NLP matching (default 1)')

    def handle(self, *args, **options):
        if options.get("clear"):
            TextPhenotype.objects.all().delete()
            PatientTextPhenotype.objects.all().delete()

        bulk_patient_phenotype_matching(cores=options["cores"])

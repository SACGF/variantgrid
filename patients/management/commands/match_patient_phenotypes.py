import logging
from collections import Counter

from django.core.management.base import BaseCommand
from django.db.models import Count

from annotation.models.models_phenotype_match import TextPhenotype, PatientTextPhenotype, TextPhenotypeMatch
from annotation.phenotype_matching import bulk_patient_phenotype_matching


def _get_ontology_text_match_counts() -> dict:
    ontology_counts = Counter()
    tpm_qs = TextPhenotypeMatch.objects.values("ontology_term__ontology_service")
    tpm_qs = tpm_qs.annotate(count=Count("pk")).values_list("ontology_term__ontology_service", "count")
    for ontology_service, count in tpm_qs:
        ontology_counts[ontology_service] += count
    return ontology_counts


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--clear', action='store_true')

    def handle(self, *args, **options):
        before_counts = _get_ontology_text_match_counts()

        if options.get("clear"):
            TextPhenotype.objects.all().delete()
            TextPhenotypeMatch.objects.all().delete()  # This removes the cache
            PatientTextPhenotype.objects.all().delete()

        bulk_patient_phenotype_matching()
        after_counts = _get_ontology_text_match_counts()
        if before_counts:
            if before_counts == after_counts:
                logging.info("No change")
                for ontology_service, count in after_counts.items():
                    print(f"{ontology_service}: {count:,}")

            else:
                all_services = set(before_counts.keys()) | set(after_counts.keys())
                for ontology_service in sorted(all_services):
                    before = before_counts.get(ontology_service, 0)
                    after = after_counts.get(ontology_service, 0)
                    print(f"{ontology_service}: {before:,} -> {after:,}")

from collections import Counter

from django.core.management.base import BaseCommand
from django.db.models import Count

from patients.models.models_patient import Patient
from patients.models.models_phenotype import TextPhenotype, TextPhenotypeMatch
from patients.phenotype_matching import bulk_patient_phenotype_matching, requeue_sentences


def _get_ontology_text_match_counts() -> dict:
    ontology_counts = Counter()
    tpm_qs = TextPhenotypeMatch.objects.values("ontology_term__ontology_service")
    tpm_qs = tpm_qs.annotate(count=Count("pk")).values_list("ontology_term__ontology_service", "count")
    for ontology_service, count in tpm_qs:
        ontology_counts[ontology_service] += count
    return ontology_counts


class Command(BaseCommand):
    category = "maintenance"

    def add_arguments(self, parser):
        group = parser.add_mutually_exclusive_group()
        group.add_argument('--stale', action='store_true',
                           help='Rematch sentences matched with an older matcher or ontology version')
        group.add_argument('--clear', action='store_true',
                           help='Rematch every sentence (patient/cohort links and approvals are kept)')
        parser.add_argument('--cores', type=int, default=1,
                            help='Number of parallel workers for sentence NLP matching (default 1)')

    def handle(self, *args, **options):
        before_counts = _get_ontology_text_match_counts()

        num_requeued = 0
        if options["clear"]:
            num_requeued = requeue_sentences(TextPhenotype.objects.filter(match_version__isnull=False))
        elif options["stale"]:
            num_requeued = requeue_sentences(TextPhenotype.stale_qs())

        bulk_patient_phenotype_matching(Patient.with_phenotype_text(), cores=options["cores"])

        print(f"Sentences requeued: {num_requeued:,}")
        # This is a very blunt count (ie individual stuff may have changed)
        after_counts = _get_ontology_text_match_counts()
        all_services = set(before_counts.keys()) | set(after_counts.keys())
        for ontology_service in sorted(all_services):
            before = before_counts.get(ontology_service, 0)
            after = after_counts.get(ontology_service, 0)
            print(f"{ontology_service}: {before:,} -> {after:,}")

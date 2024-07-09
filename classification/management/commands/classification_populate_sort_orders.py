from django.core.management import BaseCommand

from classification.models import ClassificationModification, calculate_somatic_clinical_significance_order


class Command(BaseCommand):

    def handle(self, *args, **options):
        for cm in ClassificationModification.objects.filter(is_last_published=True).iterator():
            cm.somatic_clinical_significance_sort = calculate_somatic_clinical_significance_order(cm)
            cm.save(update_fields=['somatic_clinical_significance_sort'])

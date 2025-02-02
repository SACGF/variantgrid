from django.core.management import BaseCommand

from classification.models import ClassificationModification


class Command(BaseCommand):

    def handle(self, *args, **options):
        for cm in ClassificationModification.objects.filter(is_last_published=True).iterator():
            cm: ClassificationModification
            if clin_sig_value := cm.somatic_clinical_significance_value:
                cm.somatic_clinical_significance_sort = clin_sig_value.sort_value
            else:
                cm.somatic_clinical_significance_sort = None
            cm.save(update_fields=['somatic_clinical_significance_sort'])

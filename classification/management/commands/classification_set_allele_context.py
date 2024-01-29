from django.core.management import BaseCommand
from classification.models import Classification, ClinicalContext
from classification.models.clinical_context_utils import update_clinical_contexts
from library.utils import iter_fixed_chunks


class Command(BaseCommand):

    def handle(self, *args, **options):
        print("About to update all classifications to put them in the right allele origin clinical context")
        for count, chunk in enumerate(iter_fixed_chunks(Classification.objects.iterator(chunk_size=100), chunk_size=100)):
            classify: Classification
            for classify in chunk:
                calculated_allele_origin = classify.calc_allele_origin_bucket()
                if calculated_allele_origin != classify.allele_origin_bucket:
                    print(f"Setting {classify.cr_lab_id} to {calculated_allele_origin}")
                    classify.allele_origin_bucket = calculated_allele_origin
                    classify.save(update_fields=['allele_origin_bucket'])

            update_clinical_contexts(chunk)
            print(f"Processed {count*100} Classifications")
        print("Done")

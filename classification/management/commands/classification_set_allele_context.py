from django.conf import settings
from django.core.management import BaseCommand
from django.db.models import QuerySet

from classification.enums import SpecialEKeys, AlleleOriginBucket
from classification.models import Classification, ClinicalContext, EvidenceKeyMap, ClassificationModification
from classification.models.clinical_context_utils import update_clinical_contexts, update_clinical_context
from library.utils import iter_fixed_chunks


class Command(BaseCommand):

    @staticmethod
    def update_bucket(qs: QuerySet[Classification], allele_origin_bucket: AlleleOriginBucket):
        for count, chunk in enumerate(iter_fixed_chunks(Classification.objects.iterator(chunk_size=100), chunk_size=100)):
            for classification in chunk:
                classification.allele_origin_bucket = allele_origin_bucket
                classification.save(update_fields=["allele_origin_bucket"])
            update_clinical_contexts(chunk)

    def handle(self, *args, **options):
        print("About to update all classifications to put them in the correct allele origin clinical context")

        option_key_to_bucket = {}
        for allele_origin_option in EvidenceKeyMap.cached_key(SpecialEKeys.ALLELE_ORIGIN).virtual_options:
            allele_origin_option_key = allele_origin_option.get("key")
            option_key_to_bucket[allele_origin_option_key] = AlleleOriginBucket.bucket_for_allele_origin(allele_origin_option_key)

        cm_qs = ClassificationModification.objects.filter(is_last_published=True)
        for option_key, bucket in option_key_to_bucket.items():
            misaligned_qs = cm_qs.filter(published_evidence__allele_origin__value=option_key).exclude(classification__allele_origin_bucket=bucket)
            if misaligned_qs:
                print(f"*** For allele origin \"{option_key}\" updating {misaligned_qs.count()} records to be in bucket \"{bucket}\"")
                Command.update_bucket(
                    Classification.objects.filter(pk__in=misaligned_qs.values_list("classification_id", flat=True)),
                    allele_origin_bucket=bucket
                )
            elif count := cm_qs.filter(published_evidence__allele_origin__value=option_key).count():
                print(f"There are { count } records of \"{option_key}\", all with the correct bucket \"{bucket}\"")
            else:
                print(f"There are no records with an allele origin of \"{option_key}\"")

        no_allele_origin = cm_qs.filter(classification__allele_origin_bucket__isnull=True)
        if no_allele_origin_c := Classification.objects.filter(pk__in=no_allele_origin.values_list("classification_id", flat=True)):
            print(f"Still {no_allele_origin_c.count()} records with no allele origin bucket, setting them to default of \"{settings.ALLELE_ORIGIN_NOT_PROVIDED_BUCKET}\"")
            Command.update_bucket(
                no_allele_origin_c,
                allele_origin_bucket=settings.ALLELE_ORIGIN_NOT_PROVIDED_BUCKET
            )
        else:
            print("All records now have an allele origin bucket")

        # for count, chunk in enumerate(iter_fixed_chunks(Classification.objects.iterator(chunk_size=100), chunk_size=100)):
        #     classify: Classification
        #     for classify in chunk:
        #         calculated_allele_origin = classify.calc_allele_origin_bucket()
        #         if calculated_allele_origin != classify.allele_origin_bucket:
        #             print(f"Setting {classify.cr_lab_id} to {calculated_allele_origin}")
        #             classify.allele_origin_bucket = calculated_allele_origin
        #             classify.save(update_fields=['allele_origin_bucket'])
        #
        #     update_clinical_contexts(chunk)
        #     print(f"Processed {count*100} Classifications")
        print("Done")

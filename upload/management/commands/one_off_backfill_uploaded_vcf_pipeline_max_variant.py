import time

from django.core.management.base import BaseCommand
from django.db.models import Max

from annotation.annotation_version_querysets import pipeline_type_variant_q
from annotation.models.models_enums import VariantAnnotationPipelineType
from snpdb.archive import DataArchivedError
from snpdb.models import Variant, CohortGenotypeCollection
from upload.models import UploadedVCF, UploadedVCFPipelineMaxVariant


class Command(BaseCommand):
    """ Backfill per-pipeline-type max-variant rows for VCFs imported before UploadedVCFPipelineMaxVariant
        existed, so the next annotation upgrade evaluates completeness correctly (issue #1656).

        A VCF with no rows is treated as "fully annotated" (nothing to wait for), so existing VCFs must
        get their rows before that check matters. Only VCFs that don't already have rows are processed. """

    def handle(self, *args, **options):
        # Only VCFs still lacking rows - newly imported VCFs already create them via the import path
        uploaded_vcf_qs = UploadedVCF.objects.filter(vcf__isnull=False, pipeline_max_variants__isnull=True)
        total = uploaded_vcf_qs.count()
        self.stdout.write(f"Backfilling pipeline max-variant rows for {total} UploadedVCF(s)")

        processed = 0
        tally_interval_secs = 30
        last_tally = time.monotonic()
        for uploaded_vcf in uploaded_vcf_qs.iterator():
            processed += 1
            try:
                cgc = uploaded_vcf.vcf.cohort.cohort_genotype_collection
            except (CohortGenotypeCollection.DoesNotExist, DataArchivedError, AttributeError):
                continue  # No genotype data (or archived) - nothing to enumerate

            collection_ids = [cgc.pk]
            if cgc.common_collection_id:
                collection_ids.append(cgc.common_collection_id)

            # Only non-reference variants are annotated
            variants_qs = Variant.objects.filter(cohortgenotype__collection_id__in=collection_ids) \
                .filter(Variant.get_no_reference_q())

            max_by_pipeline_type = {}
            sv_q = pipeline_type_variant_q(VariantAnnotationPipelineType.STRUCTURAL_VARIANT)
            if variants_qs.filter(sv_q).exists():
                # Mixed - pay the per-type scan
                for pipeline_type in VariantAnnotationPipelineType:
                    q = pipeline_type_variant_q(pipeline_type)
                    max_id = variants_qs.filter(q).aggregate(m=Max("pk"))["m"]
                    if max_id is not None:
                        max_by_pipeline_type[pipeline_type.value] = max_id
            else:
                # Common case: short-only - overall max is the STANDARD max, and there's no SV row
                max_id = variants_qs.aggregate(m=Max("pk"))["m"]
                if max_id is not None:
                    max_by_pipeline_type[VariantAnnotationPipelineType.STANDARD.value] = max_id

            UploadedVCFPipelineMaxVariant.objects.bulk_create([
                UploadedVCFPipelineMaxVariant(uploaded_vcf=uploaded_vcf, pipeline_type=pipeline_type,
                                              max_variant_id=max_id)
                for pipeline_type, max_id in max_by_pipeline_type.items()
            ])

            now = time.monotonic()
            if now - last_tally >= tally_interval_secs:
                self.stdout.write(f"... {processed}/{total}")
                last_tally = now

        self.stdout.write(f"Done ({processed}/{total})")

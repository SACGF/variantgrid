import logging

from django.conf import settings

from snpdb.models import AlleleOrigin, VariantAllele, Allele, SequenceRole, Variant, AlleleLiftover, ProcessingStatus
from upload.models import ModifiedImportedVariant
from upload.vcf.bulk_minimal_vcf_processor import BulkMinimalVCFProcessor


class BulkAlleleLinkingVCFProcessor(BulkMinimalVCFProcessor):
    """ Reads VCF with allele_id as ID column - then links to existing Alleles """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.allele_ids = []
        self.liftover = self.upload_pipeline.uploaded_file.uploadedliftover.liftover
        self.allele_liftovers_by_allele_id = {}
        for al in AlleleLiftover.objects.filter(liftover=self.liftover):
            self.allele_liftovers_by_allele_id[al.allele_id] = al

    @property
    def genome_build(self):
        return self.upload_step.genome_build

    def process_entry(self, variant):
        super().process_entry(variant)
        self.allele_ids.append(variant.ID)

    def batch_handle_variant_ids(self, variant_ids):
        super().batch_handle_variant_ids(variant_ids)

        variants_with_allele_qs = VariantAllele.objects.filter(variant__in=variant_ids,
                                                               genome_build=self.genome_build)
        variants_with_existing_allele = dict(variants_with_allele_qs.values_list("variant", "allele"))

        miv_qs = ModifiedImportedVariant.objects.filter(variant__in=variant_ids)

        invalid_contig_variant_ids = set()
        if settings.LIFTOVER_TO_CHROMOSOMES_ONLY:
            invalid_contigs = self.genome_build.contigs.exclude(role=SequenceRole.ASSEMBLED_MOLECULE)
            variant_qs = Variant.objects.filter(locus__contig__in=invalid_contigs, pk__in=variant_ids)
            invalid_contig_variant_ids = dict(variant_qs.values_list("pk", "locus__contig__name"))

        normalized_variants = set(miv_qs.values_list("variant_id", flat=True))
        variant_alleles = []
        updated_allele_liftovers = []
        for variant_id, allele_id in zip(variant_ids, self.allele_ids):
            variant_id = int(variant_id)
            allele_id = int(allele_id)

            al = self.allele_liftovers_by_allele_id[allele_id]
            updated_allele_liftovers.append(al)

            if invalid_contig := invalid_contig_variant_ids.get(variant_id):
                error_message = f"settings.LIFTOVER_TO_CHROMOSOMES_ONLY=True disabled liftover to non-chrom contig: {invalid_contig}"
                al.status = ProcessingStatus.ERROR
                al.error = {"message": error_message}
                continue

            existing_allele_id = variants_with_existing_allele.get(variant_id)
            if existing_allele_id:
                al.status = ProcessingStatus.SKIPPED
                error_message = f"{variant_id=} already linked to AlleleID={existing_allele_id}"
                if allele_id != existing_allele_id:
                    self.merge_alleles(allele_id, existing_allele_id)
                    error_message += f" (merged {allele_id=} + {existing_allele_id=})"
                al.error = {"message": error_message}
                continue  # Only one VariantAllele allowed, bulk_create below would have not inserted anything

            if variant_id in normalized_variants:
                origin = AlleleOrigin.LIFTOVER_NORMALIZED
            else:
                origin = AlleleOrigin.LIFTOVER

            va = VariantAllele(variant_id=variant_id,
                               allele_id=allele_id,
                               genome_build=self.genome_build,
                               origin=origin,
                               allele_linking_tool=self.liftover.conversion_tool)
            variant_alleles.append(va)
            al.status = ProcessingStatus.SUCCESS

        if updated_allele_liftovers:
            AlleleLiftover.objects.bulk_update(updated_allele_liftovers, fields=["status", "error"])

        if variant_alleles:
            logging.info("Inserting %d variant_alleles", len(variant_alleles))
            VariantAllele.objects.bulk_create(variant_alleles, ignore_conflicts=True)

        self.allele_ids = []

    def merge_alleles(self, allele_1_id, allele_2_id):
        # We always merge into the lowest PK so that things are consistent if merge(a, b) and merge(b, a) both run
        # in a race condition
        allele_1_id, allele_2_id = sorted((allele_1_id, allele_2_id))
        allele_1 = Allele.objects.get(pk=allele_1_id)
        allele_2 = Allele.objects.get(pk=allele_2_id)
        allele_1.merge(self.liftover.conversion_tool, allele_2)


class FailedLiftoverVCFProcessor(BulkMinimalVCFProcessor):
    """ Reads VCF with allele_id as ID column - """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.allele_ids = []
        self.liftover = self.upload_pipeline.uploaded_file.uploadedliftover.liftover

    @property
    def genome_build(self):
        return self.upload_step.genome_build

    def process_entry(self, variant):
        # In May 2024 I raised an issue about adding rejection details
        # maybe this has been implemented, and we can store it, https://github.com/freeseek/score/issues/10
        self.allele_ids.append(variant.ID)

    def finish(self):
        error = {"message": "BCFTools +liftover rejected variant"}
        qs = AlleleLiftover.objects.filter(liftover=self.liftover, allele__in=self.allele_ids)
        self.rows_processed = qs.update(status=ProcessingStatus.ERROR, error=error)

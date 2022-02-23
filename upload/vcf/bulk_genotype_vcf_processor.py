import csv
import logging
import os
from typing import Optional, Set

import cyvcf2
import numpy as np
from django.conf import settings

from library.django_utils import thread_safe_unique_together_get_or_create
from library.django_utils.django_file_utils import get_import_processing_filename
from library.git import Git
from library.postgres_utils import postgres_arrays
from library.utils import double_quote
from library.vcf_utils import VCFConstant
from patients.models_enums import Zygosity
from snpdb.common_variants import get_classified_high_frequency_variants_qs
from snpdb.models import CohortGenotype
from snpdb.models.models_enums import ProcessingStatus
from upload.models import UploadPipeline, PipelineFailedJobTerminateEarlyException, \
    VCFImporter, UploadStep, UploadStepTaskType, VCFPipelineStage
from upload.tasks.vcf.import_sql_copy_task import ImportCohortGenotypeSQLCopyTask
from upload.vcf.abstract_bulk_vcf_processor import AbstractBulkVCFProcessor
from upload.vcf.sql_copy_files import write_sql_copy_csv, COHORT_GENOTYPE_HEADER


class BulkGenotypeVCFProcessor(AbstractBulkVCFProcessor):
    # v1. Refactored to use Normalise/separate known/unknown variants
    # v2. Split read header (create vcf) from parsing VCF (so can be in parallel)
    # v3. PL/AD getters plus big refactoring
    # v5. Somatic, not inserting ref variants
    # v6. Used to mark CohortGenotypeCollection as fixed by fix_cohort_genotype_task
    # v7. Insert HOM_REF properly (March 23 2020)
    # v8. VCF filters are in CohortGenotype now
    # v9. No ObservedVariant (sample level data)
    # v10. Shift fields from UploadedVCF to VCF. Fix various bugs
    # v11. ?
    # v12. Ensure missing data in FreeBayes is -1 not -2147483648 (CyVCF2 returns this from format)
    # v13. Support CLCAD2 from CLC workbench
    # v14. AF in unit not percent. Use AF from VCF if provided. Support for sample level filters (FT)
    # v15. Any genotype as '.' -> unknown zygosity
    # v16. Don't insert a reference variant for each unknown ALT (only if provided) - version never deployed
    # v17. Use refactored RedisVariantPKLookup
    # v18. Handle Mixed diploid/haploid calls eg male chrX GT=1 - see https://github.com/brentp/cyvcf2/issues/227
    # v19. Go back to using CyVCF code as of version >= 0.30.14
    # v20. Split into common/uncommon genotype collections
    VCF_IMPORTER_VERSION = 20  # Change this if you make a major change to the code.
    # Need to distinguish between no entry and 0, can't use None w/postgres command line inserts
    DEFAULT_AD_FIELD = 'AD'  # What CyVCF2 uses
    # GL = Genotype Likelihood - used by freeBayes v1.2.0: log10-scaled likelihoods of the data given the called
    # genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy
    GENOTYPE_LIKELIHOOD = 'GL'
    # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
    ALT_CYVCF_GT_ZYGOSITIES = [Zygosity.HOM_REF, Zygosity.HET, Zygosity.UNKNOWN_ZYGOSITY, Zygosity.HOM_ALT]
    DOESNT_MATTER = 0  # For unknown, anything should work
    # PL field will contain three numbers, corresponding to the three possible genotypes (0/0, 0/1, and 1/1)
    # Need to work out how to map CyVCF2GTType to PL index.
    # This ends up being the same for both ref/alt
    CYVCF_HAPLOID_PL_INDEX = [0, DOESNT_MATTER, DOESNT_MATTER, 1]
    CYVCF_DIPLOID_PL_INDEX = [0, 1, DOESNT_MATTER, 2]
    CYVCF_PL_INDEX_FOR_PLOIDY = [None, CYVCF_HAPLOID_PL_INDEX, CYVCF_DIPLOID_PL_INDEX]
    EMPTY_PL_ARRAY = "{-1,-1,-1}"  # Shortcut for PL = '.'
    # Used to adjust cohort_gt_vaf_index, numbers come from fields at top of process_cohort_genotypes
    COHORT_GT_NUM_ADDED_FIELDS = 3

    @staticmethod
    def get_vcf_importer_version():
        """ Get identifier for this importer """

        vcf_importer, _ = thread_safe_unique_together_get_or_create(VCFImporter,
                                                                    name="PythonKnownVariantsImporter",
                                                                    version=BulkGenotypeVCFProcessor.VCF_IMPORTER_VERSION,
                                                                    vcf_parser="cyvcf2",
                                                                    vcf_parser_version=cyvcf2.__version__,
                                                                    code_git_hash=Git(settings.BASE_DIR).hash)
        return vcf_importer

    def __init__(self, upload_step, cohort_genotype_collection,
                 uploaded_vcf, preprocess_vcf_import_info,
                 batch_size=settings.SQL_BATCH_INSERT_SIZE):
        super().__init__(upload_step, preprocess_vcf_import_info, batch_size=batch_size)

        self.checked_pipeline_for_failures_this_row = False
        self.vcf = uploaded_vcf.vcf
        self.cohort_genotype_collection = cohort_genotype_collection
        self.cohort_genotype_file_id = 0

        num_samples = uploaded_vcf.vcf.genotype_samples
        # We need to sum up the ADs of decomposed (single allele per row) VCF
        # So have to store variants til finished with locus
        self.last_locus_tuple = None
        self.locus_ad_sum = np.zeros(num_samples, dtype='i')
        self.locus_variant_hashes = []
        self.locus_filters = []
        self.locus_cohort_genotypes = []
        self.locus_gnomad_af = []
        self.locus_allele_depths = []
        self.locus_modified_imported_variant_hashes = []
        self.locus_modified_imported_variants = []

        # Variant lists are kept in sync with super().variant_hashes
        self.variant_filters = []
        self.cohort_genotypes = []
        self.variant_gnomad_af = []

        # These are variant IDs that won't be kept in "common" CGC regardless of frequency
        self.uncommon_variant_ids = self._get_uncommon_variant_ids()

        if num_samples:
            # Only need this if we have genotypes (otherwise will be NoGenotypeProcessor)
            self.get_ref_alt_allele_depth = self.get_ref_alt_allele_depth_function(self.vcf)
        self.vcf_filter_map = uploaded_vcf.vcf.get_filter_dict()
        self.last_read_depth_str = None
        self.child_processes = []
        self.cohort_gt_vaf_index = COHORT_GENOTYPE_HEADER.index("samples_allele_frequency")
        if self.cohort_gt_vaf_index < 0:
            raise ValueError(f"Could not find 'samples_allele_frequency' in {COHORT_GENOTYPE_HEADER}")
        self.cohort_gt_vaf_index -= self.COHORT_GT_NUM_ADDED_FIELDS

    def _get_uncommon_variant_ids(self) -> Set[int]:
        uncommon_variant_ids = set()
        if common_collection := self.cohort_genotype_collection.common_collection:
            uncommon_qs = get_classified_high_frequency_variants_qs(common_collection.common_filter)
            uncommon_variant_ids.update(uncommon_qs.values_list("pk", flat=True))
        return uncommon_variant_ids

    def finish(self):
        """ This is called at the very end so we can collect any remaining items to process """

        logging.debug("VCFSplitToTableCSVs - finish")

        self.finished_locus()
        self.batch_process_check(0)  # Insert anything that is there
        # Make sure child processes have finished
        for p in self.child_processes:
            p.join()

    def convert_filters(self, vcf_filter) -> Optional[str]:
        filters = None
        if vcf_filter and vcf_filter != '.':  # FORMAT filters can be '.'
            filters_iter = (self.vcf_filter_map[f] for f in vcf_filter.split(";"))
            filters = ''.join(sorted(filters_iter))
        return filters

    @staticmethod
    def get_ref_alt_allele_depth_default(variant):
        return variant.gt_ref_depths, variant.gt_alt_depths

    @staticmethod
    def get_ref_alt_allele_depth_function(vcf):
        if vcf.allele_depth_field == BulkGenotypeVCFProcessor.DEFAULT_AD_FIELD:
            return BulkGenotypeVCFProcessor.get_ref_alt_allele_depth_default  # use cyvcf2 defaults
        if vcf.allele_depth_field == VCFConstant.CLCAD2:
            # @see http://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/700/index.php?manual=Annotation_variant_formats.html
            # CLCAD2 tag follow the order of REF and ALT, with one value for the REF and for each ALT
            def _get_clcad2_ref_alt_depths(variant):
                # VCF is decomposed so will only have 0 or 1 alt
                clcad2 = variant.format(vcf.allele_depth_field).flatten()
                ref = clcad2[0]
                if len(clcad2) == 2:
                    alt = clcad2[1]
                else:
                    alt = CohortGenotype.MISSING_NUMBER_VALUE
                return np.array([ref]), np.array([alt])
            return _get_clcad2_ref_alt_depths
        if vcf.ref_depth_field and vcf.alt_depth_field:  # explicitly set
            def get_ref_alt_depths(variant):
                # CyVCF2 returns -1 for empty in gt_ref_depths / gt_alt_depths but this value in format()
                # So convert to make consistent. @see https://github.com/brentp/cyvcf2/issues/172
                _NUMPY_INT_NAN_VALUE = -2147483648

                # As VCF is decomposed these arrays will always be length of 1, so return in same shape as default above
                ref_depth = variant.format(vcf.ref_depth_field).flatten()
                ref_depth[ref_depth == _NUMPY_INT_NAN_VALUE] = CohortGenotype.MISSING_NUMBER_VALUE
                alt_depth = variant.format(vcf.alt_depth_field).flatten()
                alt_depth[alt_depth == _NUMPY_INT_NAN_VALUE] = CohortGenotype.MISSING_NUMBER_VALUE
                return ref_depth, alt_depth
            return get_ref_alt_depths
        raise ValueError(f"Don't know how to get ref and alt allele depth for {vcf}")

    @staticmethod
    def get_format_array_str(variant, field, as_type=None) -> Optional[str]:
        format_array_str = None
        if field:
            format_array = variant.format(field)
            if format_array is not None:
                format_array[np.isnan(format_array) | (format_array < 0)] = CohortGenotype.MISSING_NUMBER_VALUE
                if as_type:
                    format_array = format_array.astype(as_type)
                format_array_str = postgres_arrays(format_array.flat)
        return format_array_str

    def get_phred_likelihood_str(self, variant, gt_types) -> Optional[str]:
        phred_likelihood_str = None
        if self.vcf.phred_likelihood_field:
            phred_likelihood = []
            pl = variant.format(self.vcf.phred_likelihood_field)
            if pl is not None:
                if self.vcf.phred_likelihood_field == BulkGenotypeVCFProcessor.GENOTYPE_LIKELIHOOD:
                    pl *= -10  # GQ = log10(P), PL = -10 log10(P)

                # Handle PL = '.' (all missing if unknown genotype)
                missing = np.isnan(pl) | (pl < 0)
                if missing.all():
                    return self.EMPTY_PL_ARRAY

                pl[missing] = CohortGenotype.MISSING_NUMBER_VALUE  # Handle individual PL array entry missing
                for i, (gt, genotype) in enumerate(zip(gt_types, variant.genotypes)):
                    # CyVCF2 genotypes entry for sample is eg [0, 1, True] (last element = is phased)
                    # Handle case where PL is not variant.ploidy + 1 - eg where male chrX given as haploid, ie GT=1
                    # see https://github.com/brentp/cyvcf2/issues/227
                    ploidy = len(genotype) - 1
                    pl_index_lookup = BulkGenotypeVCFProcessor.CYVCF_PL_INDEX_FOR_PLOIDY[ploidy]
                    pl_index = pl_index_lookup[gt]
                    pl_value = int(pl[i][pl_index])
                    phred_likelihood.append(pl_value)

                phred_likelihood_str = postgres_arrays(phred_likelihood)
        return phred_likelihood_str

    def get_samples_filters_str(self, variant) -> Optional[str]:
        samples_filters_str = None
        if self.vcf.sample_filters_field:
            sample_filters = []
            s_filters = variant.format(self.vcf.sample_filters_field)
            if s_filters is not None:
                for sf in s_filters:
                    sample_filters.append(self.convert_filters(sf))
                if any([s is not None for s in sample_filters]):  # Skip if all were PASS
                    # https://www.postgresql.org/docs/9.1/arrays.html - Postgres arrays use double quotes
                    cleaned_sample_filters = []
                    for sf in sample_filters:
                        if sf is not None:
                            sf = double_quote(sf)
                        cleaned_sample_filters.append(sf)
                    samples_filters_str = postgres_arrays(cleaned_sample_filters)
        return samples_filters_str

    def finished_locus(self):
        """ sum(AD) for this locus and add data to arrays """

        self.variant_hashes.extend(self.locus_variant_hashes)
        self.variant_filters.extend(self.locus_filters)
        self.variant_gnomad_af.extend(self.locus_gnomad_af)
        self.modified_imported_variant_hashes.extend(self.locus_modified_imported_variant_hashes)
        self.modified_imported_variants.extend(self.locus_modified_imported_variants)

        if self.vcf.allele_frequency_field:
            # AF will already be there
            self.cohort_genotypes.extend(self.locus_cohort_genotypes)
        else:
            # Calculate ourselves across locus
            for cgt, ad in zip(self.locus_cohort_genotypes, self.locus_allele_depths):
                vaf = ad / self.locus_ad_sum
                vaf[np.isnan(vaf)] = CohortGenotype.MISSING_NUMBER_VALUE
                cgt[self.cohort_gt_vaf_index] = postgres_arrays(vaf)
                self.cohort_genotypes.append(cgt)

        self.locus_variant_hashes = []
        self.locus_filters = []
        self.locus_gnomad_af = []
        self.locus_cohort_genotypes = []
        self.locus_allele_depths = []
        self.locus_modified_imported_variant_hashes = []
        self.locus_modified_imported_variants = []
        self.locus_ad_sum.fill(0)

    def process_entry(self, variant: cyvcf2.Variant):
        # Pre-processed by vcf_filter_unknown_contigs so only recognised contigs present
        ref, alt = self.get_ref_alt(variant)
        locus_tuple = (variant.CHROM, variant.POS, ref)

        # These ADs come out with empty value as -1 - that's what we want to store
        ref_allele_depth, alt_allele_depth = self.get_ref_alt_allele_depth(variant)
        # We want empty values as 0 so that adding them is ok
        empty_as_zero_ref_allele_depth = ref_allele_depth.copy()
        empty_as_zero_ref_allele_depth[empty_as_zero_ref_allele_depth < 0] = 0
        empty_as_zero_alt_allele_depth = alt_allele_depth.copy()
        empty_as_zero_alt_allele_depth[empty_as_zero_alt_allele_depth < 0] = 0

        if self.vcf.allele_frequency_field:
            allele_frequency_str = self.get_format_array_str(variant, self.vcf.allele_frequency_field, as_type=float)
        else:
            # We'll calculate ourselves across locus
            allele_frequency_str = None

        read_depth_str = self.get_format_array_str(variant, self.vcf.read_depth_field, as_type=int)
        genotype_quality_str = self.get_format_array_str(variant, self.vcf.genotype_quality_field, as_type=int)
        phred_likelihood_str = self.get_phred_likelihood_str(variant, variant.gt_types)
        samples_filters_str = self.get_samples_filters_str(variant)

        if self.last_locus_tuple:
            if self.last_locus_tuple != locus_tuple:
                self.finished_locus()

        alt_hash = self.variant_pk_lookup.get_variant_coordinate_hash(variant.CHROM, variant.POS, ref, alt)
        alt_zygosity = [BulkGenotypeVCFProcessor.ALT_CYVCF_GT_ZYGOSITIES[i] for i in variant.gt_types]
        alt_allele_depth_str = postgres_arrays(alt_allele_depth)

        cohort_gt = [str(variant.num_hom_ref),
                     str(variant.num_het),
                     str(variant.num_hom_alt),
                     str(variant.num_unknown),
                     ''.join(alt_zygosity),
                     alt_allele_depth_str,
                     allele_frequency_str,
                     read_depth_str,
                     genotype_quality_str,
                     phred_likelihood_str,
                     samples_filters_str]

        self.locus_variant_hashes.append(alt_hash)
        self.locus_filters.append(self.convert_filters(variant.FILTER))
        self.locus_cohort_genotypes.append(cohort_gt)
        self.locus_gnomad_af.append(variant.INFO.get("AF"))
        self.locus_allele_depths.append(empty_as_zero_alt_allele_depth)

        if self.preprocess_vcf_import_info:
            self.add_modified_imported_variant(variant, alt_hash,
                                               miv_hash_list=self.locus_modified_imported_variant_hashes,
                                               miv_list=self.locus_modified_imported_variants)

        self.locus_ad_sum += empty_as_zero_ref_allele_depth + empty_as_zero_alt_allele_depth
        self.last_read_depth_str = read_depth_str
        self.last_locus_tuple = locus_tuple

        self.rows_processed += 1
        self.batch_process_check()

    def check_pipeline_for_failures(self):
        """ If a pipeline has already died, quit the Reading early to save doing useless work
            Check it before we write files / create new jobs, and once every settings.VCF_IMPORT_PROCESS_CHECK_PIPELINE_FAIL_ROWS
            This stores checked_pipeline_for_failures_this_row so it only makes one SQL query per row at most
        """
        if not self.checked_pipeline_for_failures_this_row:
            logging.info("Reloading Pipeline to check for failures")
            pipeline = UploadPipeline.objects.get(pk=self.upload_pipeline.pk)
            if pipeline.status == ProcessingStatus.ERROR:
                raise PipelineFailedJobTerminateEarlyException()
            self.checked_pipeline_for_failures_this_row = True

    def batch_process_check(self, minimum_insert_size=None):
        if minimum_insert_size is None:
            minimum_insert_size = self.batch_size

        self.checked_pipeline_for_failures_this_row = False

        if self.variant_hashes and len(self.variant_hashes) >= minimum_insert_size:
            variant_ids = self.variant_pk_lookup.get_variant_ids(self.variant_hashes)
            self.set_max_variant(self.variant_hashes, variant_ids)

            if self.modified_imported_variants:
                variant_ids_by_hash = dict(zip(self.variant_hashes, variant_ids))
                self.process_modified_imported_variants(variant_ids_by_hash)

            self.process_cohort_genotypes(variant_ids)

    def process_cohort_genotypes(self, variant_ids):
        cohort_genotypes_common = []
        cohort_genotypes_rare = []
        # If you add any columns here, need to adjust COHORT_GT_NUM_ADDED_FIELDS
        for variant_id, filters, cohort_gt, gnomad_af in zip(variant_ids, self.variant_filters,
                                                             self.cohort_genotypes, self.variant_gnomad_af):
            common = gnomad_af and variant_id not in self.uncommon_variant_ids
            if common:
                cgc_id = self.cohort_genotype_collection.common_collection_id
                cohort_genotypes = cohort_genotypes_common
            else:  # Rare
                cgc_id = self.cohort_genotype_collection.pk
                cohort_genotypes = cohort_genotypes_rare
            cohort_genotypes.append([cgc_id, variant_id, filters] + cohort_gt)

        cg_destinations = {
            "common": (cohort_genotypes_common, {"base_table_name": "common"}),
            "rare": (cohort_genotypes_rare, {}),
        }
        total_variants = len(variant_ids)
        for cg_dest, (cohort_genotypes, partition_kwargs) in cg_destinations.items():
            num_genotypes = len(cohort_genotypes)
            logging.info(f"{cg_dest}: {num_genotypes} ({100 * num_genotypes / total_variants}%)")
            if not cohort_genotypes:
                continue
            cg_basename = f"cohort_genotype_step_{self.upload_step.pk}_batch_{self.cohort_genotype_file_id}.csv"
            cohort_genotypes_filename = get_import_processing_filename(self.upload_pipeline.pk, cg_basename)
            write_sql_copy_csv(cohort_genotypes, cohort_genotypes_filename,
                               quoting=csv.QUOTE_NONE, escapechar='\\', quotechar='')
            table_name = self.cohort_genotype_collection.get_partition_table(**partition_kwargs)
            num_cohort_genotypes = len(cohort_genotypes)
            self.create_cohort_genotype_job(table_name, num_cohort_genotypes, cohort_genotypes_filename)
            self.cohort_genotype_file_id += 1

        self.variant_hashes = []
        self.cohort_genotypes = []
        self.variant_gnomad_af = []
        self.check_pipeline_for_failures()  # Need to do this every so often

    def create_cohort_genotype_job(self, table_name, items_to_process, input_filename):
        if not os.path.exists(input_filename):
            msg = f"create_cohort_genotype_job: input file: '{input_filename}' does not exist."
            raise ValueError(msg)

        name = "CohortGenotypeCollection SQL COPY"
        sort_order = self.upload_pipeline.get_max_step_sort_order() + 1
        sql_job = UploadStep.objects.create(upload_pipeline=self.upload_pipeline,
                                            name=name,
                                            sort_order=sort_order,
                                            task_type=UploadStepTaskType.SQL,
                                            pipeline_stage=VCFPipelineStage.DATA_INSERTION,
                                            input_filename=input_filename,
                                            items_to_process=items_to_process,
                                            import_variant_table=table_name)

        sql_job.launch_task(ImportCohortGenotypeSQLCopyTask)

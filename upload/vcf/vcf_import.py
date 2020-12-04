from typing import Any

import cyvcf2
from django.conf import settings
from django.db.models.query_utils import Q
from django.urls.base import reverse
from django.utils import timezone
from django_messages.models import Message
import logging
import os
import re
import sys
import traceback

from library.guardian_utils import assign_permission_to_user_and_groups
from library.vcf_utils import cyvcf2_header_types, cyvcf2_header_get, VCFConstant,\
    cyvcf2_get_contig_lengths_dict
from seqauto.models import SampleSheetCombinedVCFFile, VCFFile, \
    VCFFromSequencingRun, SampleFromSequencingSample, QCGeneList
from seqauto.signals import backend_vcf_import_start_signal
from snpdb.models import VCF, ImportStatus, Sample, VCFFilter, \
    Cohort, CohortSample, UserSettings, VCFSourceSettings
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_enums import ImportSource
from snpdb.tasks.cohort_genotype_tasks import create_cohort_genotype_collection
from upload.models import UploadedVCF, PipelineFailedJobTerminateEarlyException, \
    BackendVCF, UploadStep, ModifiedImportedVariants
from upload.models_enums import UploadStepTaskType, VCFPipelineStage
from upload.tasks.vcf.import_sql_copy_task import ImportModifiedImportedVariantSQLCopyTask
from upload.vcf.bulk_genotype_vcf_processor import BulkGenotypeVCFProcessor
from upload.vcf.bulk_no_genotype_vcf_processor import BulkNoGenotypeVCFProcessor


def get_format_field(vcf_formats, wanted_format_id):
    """ returns field if present, None if missing from format header """
    format_id = None
    if wanted_format_id in vcf_formats:
        format_id = wanted_format_id

    return format_id


def set_allele_depth_format_fields(vcf: VCF, vcf_formats, vcf_source, default_allele_field):
    # Use FreeBayes AO/RO fields due to AD field not being decomposed properly on multi-alts
    # @see https://github.com/SACGF/variantgrid/issues/2126
    if VCFConstant.FREEBAYES in vcf_source:
        if VCFConstant.ALT_DEPTH_FIELD in vcf_formats and VCFConstant.REF_DEPTH_FIELD in vcf_formats:
            vcf.allele_depth_field = None
            vcf.alt_depth_field = VCFConstant.ALT_DEPTH_FIELD
            vcf.ref_depth_field = VCFConstant.REF_DEPTH_FIELD
            return

    # Most callers pack depths for all alleles in AD field
    if VCFConstant.CLCAD2 in vcf_formats:
        ad = VCFConstant.CLCAD2
    else:  # Default (GATK)
        ad = get_format_field(vcf_formats, default_allele_field)

    if ad is None:
        msg = f"Couldn't determine allele depth format field, source: '{vcf_source}', formats: {vcf_formats}"
        raise ValueError(msg)

    vcf.allele_depth_field = ad


def get_phred_likelihood_field(vcf_formats, vcf_source, default_phred_likelihood_field):
    pl = get_format_field(vcf_formats, default_phred_likelihood_field)
    if pl is None:
        if VCFConstant.GENOTYPE_LIKELIHOOD in vcf_formats:
            pl = VCFConstant.GENOTYPE_LIKELIHOOD
    return pl


def create_vcf_filters(vcf, header_types):
    filter_code_id = VCFFilter.ASCII_MIN
    filters = header_types.get("FILTER", [])
    for filter_id, filter_dict in filters.items():
        if filter_id == "PASS":  # Special - don't store this as vcf.Reader will not return it
            continue
        filter_description = filter_dict["Description"]

        filter_code = chr(filter_code_id)
        VCFFilter.objects.create(vcf=vcf,
                                 filter_code=filter_code,
                                 filter_id=filter_id,
                                 description=filter_description)
        filter_code_id += 1
        if filter_code_id > VCFFilter.ASCII_MAX:
            num_filters = VCFFilter.ASCII_MAX - VCFFilter.ASCII_MIN
            logging.warning("Warning: Run out of characters to store filters! Only storing 1st %d.", num_filters)
            continue


def create_cohort_genotype_collection_from_vcf(vcf: VCF):
    """ Cohort will exist if reloading from existing VCF """
    assert vcf.genome_build, "VCF must have a genome build"

    vcf_name = f"{vcf}"
    vcf_name = vcf_name[:40]
    defaults = {"name": f"VCF: {vcf.pk} {vcf_name}",
                "import_status": ImportStatus.IMPORTING,
                "genome_build": vcf.genome_build}
    cohort, created = Cohort.objects.update_or_create(vcf=vcf, defaults=defaults)
    if created:
        assign_permission_to_user_and_groups(vcf.user, cohort)

        for i, sample in enumerate(vcf.sample_set.order_by("pk")):
            CohortSample.objects.create(cohort=cohort,
                                        sample=sample,
                                        cohort_genotype_packed_field_index=i,
                                        sort_order=i)

        create_cohort_genotype_collection(cohort)


def create_vcf_from_vcf(upload_step, vcf_filename) -> VCF:
    """ Reads header only - returns VCF object """

    upload_pipeline = upload_step.upload_pipeline
    uploaded_file = upload_pipeline.uploaded_file
    # There is no VCF file if this method is called, but it could be possible for UploadedVCF
    # to exist if there was an error before previous VCF was created, so re-use if there
    try:
        # Reload as VCF may have been deleted before calling this and still attached
        uploaded_vcf = UploadedVCF.objects.get(upload_pipeline=upload_pipeline)
        uploaded_vcf.uploaded_file = uploaded_file
        uploaded_vcf.vcf_importer = None
        uploaded_vcf.save()
    except UploadedVCF.DoesNotExist:
        uploaded_vcf = UploadedVCF.objects.create(uploaded_file=uploaded_file,
                                                  upload_pipeline=upload_pipeline)
        logging.info("import_vcf_file: CREATED uploaded_vcf: %d", uploaded_vcf.pk)

    vcf_reader = cyvcf2.VCF(vcf_filename)
    num_genotype_samples = len(vcf_reader.samples)
    vcf = create_vcf_from_uploaded_vcf(uploaded_vcf, num_genotype_samples)
    # Save header ASAP if case something goes wrong
    vcf.header = vcf_reader.raw_header

    try:
        vcf.genome_build = vcf_detect_genome_build_from_header(vcf_reader)
    except GenomeBuildDetectionException:
        pass
    vcf.save()

    header_types = cyvcf2_header_types(vcf_reader)
    create_vcf_filters(vcf, header_types)
    vcf_formats = set(header_types["FORMAT"])
    source = cyvcf2_header_get(vcf_reader, "source", "")
    vcf.source = source
    if num_genotype_samples:
        set_allele_depth_format_fields(vcf, vcf_formats, source, VCFConstant.DEFAULT_ALLELE_FIELD)
        vcf.read_depth_field = get_format_field(vcf_formats, VCFConstant.DEFAULT_READ_DEPTH_FIELD)
        vcf.genotype_quality_field = get_format_field(vcf_formats, VCFConstant.DEFAULT_GENOTYPE_QUALITY_FIELD)
        vcf.phred_likelihood_field = get_phred_likelihood_field(vcf_formats, source, VCFConstant.DEFAULT_PHRED_LIKILIHOOD_FIELD)
    vcf.save()
    uploaded_vcf.vcf = vcf
    uploaded_vcf.save()

    backend_vcf = create_backend_vcf_links(uploaded_vcf)
    if num_genotype_samples > 0:
        sample_names = vcf_reader.samples
    else:  # Need at least 1 sample per VCF
        default_name = uploaded_vcf.uploaded_file.name
        sample_names = [default_name]

    no_dna_control_sample_pattern = None
    if settings.VCF_IMPORT_NO_DNA_CONTROL_SAMPLE_REGEX:
        no_dna_control_sample_pattern = re.compile(settings.VCF_IMPORT_NO_DNA_CONTROL_SAMPLE_REGEX)

    logging.info("Creating samples")
    for sample_name in sample_names:
        no_dna_control = False
        if no_dna_control_sample_pattern and no_dna_control_sample_pattern.findall(sample_name):
            no_dna_control = True

        sample = Sample.objects.create(vcf=vcf,
                                       vcf_sample_name=sample_name,
                                       name=sample_name,
                                       no_dna_control=no_dna_control,
                                       import_status=ImportStatus.IMPORTING)

    if backend_vcf:
        link_samples_and_vcfs_to_sequencing(backend_vcf)
        backend_vcf_import_start_signal.send(sender=os.path.basename(__file__), backend_vcf=backend_vcf)

    handle_vcf_source(vcf)
    return vcf


def handle_vcf_source(vcf):
    if vcf.source:
        for vss in VCFSourceSettings.objects.all():
            if m := re.match(vss.source_regex, vcf.source):
                if vss.operation == VCFSourceSettings.SET_SAMPLE_VARIANTS_TYPE:
                    print(f"VCF: {vcf} has source {vcf.source} -> variants_type {vss.value}")
                    vcf.sample_set.all().update(variants_type=vss.value)


def genotype_vcf_processor_factory(upload_step, cohort_genotype_collection, uploaded_vcf, preprocess_vcf_import_info):
    if uploaded_vcf.vcf.has_genotype:
        klass = BulkGenotypeVCFProcessor
    else:
        klass = BulkNoGenotypeVCFProcessor
    print(f"{uploaded_vcf.vcf} has {uploaded_vcf.vcf.genotype_samples} genotype samples, klass: {klass}")
    return klass(upload_step, cohort_genotype_collection, uploaded_vcf, preprocess_vcf_import_info)


def import_vcf_file(upload_step):
    """ This can be run in parallel """

    upload_pipeline = upload_step.upload_pipeline
    uploaded_vcf = upload_pipeline.uploadedvcf

    vcf = uploaded_vcf.vcf
    cohort_genotype_collection = vcf.cohort.cohort_genotype_collection

    preprocess_vcf_import_info = get_preprocess_vcf_import_info(upload_pipeline)
    bulk_inserter = genotype_vcf_processor_factory(upload_step, cohort_genotype_collection, uploaded_vcf, preprocess_vcf_import_info)
    uploaded_vcf.vcf_importer = bulk_inserter.get_vcf_importer_version()
    uploaded_vcf.save()

    fast_vcf_reader = cyvcf2.VCF(upload_step.input_filename)

    record_id = 1
    v: Any = None
    try:
        for v in fast_vcf_reader:
            bulk_inserter.process_entry(v)
            record_id += 1
    except PipelineFailedJobTerminateEarlyException:
        raise
    except Exception as e:
        failed_vcf_line = str(v)

        exc_traceback = sys.exc_info()[2]
        tb = ''.join(traceback.format_tb(exc_traceback))
        params = (str(e), upload_step.input_filename, record_id, failed_vcf_line, tb)
        message = "Exception: '%s'\nVCF file '%s' record: %d\nLine: '%s'\nTraceback:\n%s" % params
        raise ValueError(message)

    bulk_inserter.finish()
    update_uploaded_vcf_max_variant(uploaded_vcf.pk, bulk_inserter.get_max_variant_id())
    return bulk_inserter.rows_processed


def get_preprocess_vcf_import_info(upload_pipeline):
    try:
        preprocess_vcf_sub_step = upload_pipeline.uploadstep_set.filter(name=UploadStep.PREPROCESS_VCF_NAME).get()
        # ModifiedImportedVariants object was created in Preprocess upload step
        return ModifiedImportedVariants.objects.get(upload_step=preprocess_vcf_sub_step)
    except:
        return None


def read_vcf_file_set_max_variant(upload_step, vcf_processor_klass):
    """ This can be run in parallel """

    upload_pipeline = upload_step.upload_pipeline
    uploaded_vcf = upload_pipeline.uploadedvcf

    preprocess_vcf_import_info = get_preprocess_vcf_import_info(upload_pipeline)
    bulk_inserter = vcf_processor_klass(upload_step, preprocess_vcf_import_info)
    uploaded_vcf.vcf_importer = bulk_inserter.get_vcf_importer_version()
    uploaded_vcf.save()

    fast_vcf_reader = cyvcf2.VCF(upload_step.input_filename)
    record_id = 1
    v: Any = None
    try:
        for v in fast_vcf_reader:
            bulk_inserter.process_entry(v)
            record_id += 1
    except PipelineFailedJobTerminateEarlyException:
        raise
    except Exception as e:
        failed_vcf_line = str(v)

        exc_traceback = sys.exc_info()[2]
        tb = ''.join(traceback.format_tb(exc_traceback))
        params = (str(e), upload_step.input_filename, record_id, failed_vcf_line, tb)
        message = "Exception: '%s'\nVCF file '%s' record: %d\nLine: '%s'\nTraceback:\n%s" % params
        raise ValueError(message)

    bulk_inserter.finish()
    update_uploaded_vcf_max_variant(uploaded_vcf.pk, bulk_inserter.get_max_variant_id())
    return bulk_inserter.rows_processed


def update_uploaded_vcf_max_variant(pk, max_inserted_variant_id):
    """ This can be run in parallel """

    if max_inserted_variant_id is not None:
        uploaded_vcf_qs = UploadedVCF.objects.filter(pk=pk)
        q_not_set_or_less_than = Q(max_variant__isnull=True) | Q(max_variant_id__lt=max_inserted_variant_id)
        uploaded_vcf_qs.filter(q_not_set_or_less_than).update(max_variant_id=max_inserted_variant_id)


def create_vcf_from_uploaded_vcf(uploaded_vcf, num_genotype_samples):
    vcf = VCF.objects.create(name=uploaded_vcf.uploaded_file.name,
                             date=timezone.now(),
                             user=uploaded_vcf.uploaded_file.user,
                             genotype_samples=num_genotype_samples,
                             import_status=ImportStatus.IMPORTING)
    logging.debug("Saved vcf %d from uploaded_vcf %d", vcf.pk, uploaded_vcf.pk)

    # Assign view permission to all of users groups
    # TODO: Make this a user option (auto share to groups?)
    user = uploaded_vcf.uploaded_file.user
    assign_permission_to_user_and_groups(user, vcf)

    return vcf


def create_backend_vcf_links(uploaded_vcf):
    """ returns backend_vcf if we can link to SeqAuto data (None if Web VCF) """

    backend_vcf = None
    uploaded_file = uploaded_vcf.uploaded_file
    if uploaded_file.import_source == ImportSource.SEQAUTO:
        path = uploaded_file.path
        if path:
            combo_vcf = None
            vcf_file = None
            try:
                combo_vcf = SampleSheetCombinedVCFFile.objects.get(path=path)
            except:
                try:
                    vcf_file = VCFFile.objects.get(path=path)
                except:
                    msg = f"Couldn't find Combo or single VCF for path '{path}'"
                    raise ValueError(msg)

            backend_vcf = BackendVCF.objects.create(uploaded_vcf=uploaded_vcf,
                                                    vcf_file=vcf_file,
                                                    combo_vcf=combo_vcf)
    return backend_vcf


def link_samples_and_vcfs_to_sequencing(backend_vcf, replace_existing=False):
    if backend_vcf:
        logging.info("link_samples_and_vcfs_to_sequencing backend_vcf: %s", backend_vcf.pk)
        sequencing_run = backend_vcf.sample_sheet.sequencing_run
        vcf = backend_vcf.vcf
        vcf.fake_data = sequencing_run.fake_data
        vcf.save()

        try:
            pfsr = VCFFromSequencingRun.objects.get(sequencing_run=sequencing_run)
            if replace_existing or pfsr.vcf.import_status in ImportStatus.DELETION_STATES:
                pfsr.vcf = vcf  # OK to replace
                pfsr.save()
            else:
                logging.warning("SR %s already linked to non-deleting vcf: %s (%s/%s)", sequencing_run, pfsr.vcf, pfsr.vcf.pk, pfsr.vcf.import_status)
        except VCFFromSequencingRun.DoesNotExist:
            VCFFromSequencingRun.objects.create(vcf=vcf,
                                                sequencing_run=sequencing_run)

        samples_by_sequencing_sample = backend_vcf.get_samples_by_sequencing_sample()

        for sequencing_sample, sample in samples_by_sequencing_sample.items():
            bam_file = sequencing_sample.get_single_bam()
            if bam_file:
                sample.bam_file_path = bam_file.path
                sample.save()

            try:
                q_existing = Q(sequencing_sample=sequencing_sample) | Q(sample=sample)
                sfss = SampleFromSequencingSample.objects.get(q_existing)
                if replace_existing or sfss.sample.import_status in ImportStatus.DELETION_STATES:
                    logging.info("Updating SampleFromSequencingSample")
                    sfss.sequencing_sample = sequencing_sample
                    sfss.sample = sample
                    sfss.save()
                else:
                    logging.warning("SS %s already linked to sample: %s (%s/%s)", sfss.sequencing_sample, sfss.sample, sfss.sample.pk, sfss.sample.import_status)
            except SampleFromSequencingSample.DoesNotExist:
                SampleFromSequencingSample.objects.create(sample=sample,
                                                          sequencing_sample=sequencing_sample)

            # Link any QCGeneLists
            for qcgl in QCGeneList.objects.filter(qc__bam_file__unaligned_reads__sequencing_sample=sequencing_sample,
                                                  custom_text_gene_list__gene_list__isnull=False,
                                                  sample_gene_list__isnull=True).distinct():
                qcgl.create_and_assign_sample_gene_list(sample)


def create_import_success_message(vcf):
    subject = f"VCF '{vcf.name}' import complete"
    url = reverse('view_vcf', kwargs={'vcf_id': vcf.pk})
    body = f"VCF {vcf.name} imported as <a href='{url}'>vcf #{vcf.pk}</a>"

    user = vcf.user
    user_settings = UserSettings.get_for_user(user)
    if user_settings.import_messages:
        Message.objects.create(subject=subject,
                               body=body,
                               sender=user,
                               recipient=user)


def create_modified_imported_variants_job(upload_pipeline, num_modified_imported_variants, input_filename):
    # Should this be it's own table?? Should only be 10s of thousands per file - and we don't need it fast?

    if not os.path.exists(input_filename):
        msg = f"create_modified_imported_variants_job: input file: '{input_filename}' does not exist."
        raise ValueError(msg)

    name = "ModifiedImportedVariant SQL COPY"
    sort_order = upload_pipeline.get_max_step_sort_order() + 1
    sql_job = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                        name=name,
                                        sort_order=sort_order,
                                        task_type=UploadStepTaskType.SQL,
                                        pipeline_stage=VCFPipelineStage.DATA_INSERTION,
                                        input_filename=input_filename,
                                        items_to_process=num_modified_imported_variants)

    sql_job.launch_task(ImportModifiedImportedVariantSQLCopyTask)


class GenomeBuildDetectionException(Exception):
    pass


class ContigMismatchException(GenomeBuildDetectionException):
    pass


def vcf_detect_genome_build(vcf_filename):
    vcf_reader = cyvcf2.VCF(vcf_filename)
    return vcf_detect_genome_build_from_header(vcf_reader)


def vcf_detect_genome_build_from_header(cyvcf2_reader: cyvcf2.VCF):
    load_failures = []

    contig_lengths = cyvcf2_get_contig_lengths_dict(cyvcf2_reader)
    if contig_lengths:
        try:
            return get_genome_build_from_contig_lengths(contig_lengths)
        except ContigMismatchException as cme:
            load_failures.append(str(cme))

    reference = cyvcf2_header_get(cyvcf2_reader, "reference")
    if reference:
        try:
            return GenomeBuild.get_name_or_alias(reference)
        except GenomeBuild.DoesNotExist:
            load_failures.append(f"System does not support GenomeBuild matching VCF header reference='{reference}' - must be one of {GenomeBuild.available_names_or_aliases()}")

    error_msg = "Could not load GenomeBuild using VCF header."

    if load_failures:
        load_failures_str = ", ".join(load_failures)
        msg = error_msg + f" Attempts: {load_failures_str}"
    else:
        msg = error_msg + " No 'reference' or 'contig' entries in header."
    raise GenomeBuildDetectionException(msg)


def get_genome_build_from_contig_lengths(contig_lengths: dict, min_contig_matches=1) -> GenomeBuild:
    """ Match based on contig length (main chromosomes only)
        This isn't fast, designed to be called once per file, not in a loop """

    build_contig_diffs = {}
    potential_genome_builds = []

    for genome_build in GenomeBuild.builds_with_annotation():
        diff = None
        num_matching_contigs = 0
        num_missing_contigs = 0
        for vcf_contig, vcf_contig_length in contig_lengths.items():
            try:
                contig = genome_build.chrom_contig_mappings[vcf_contig]
                if vcf_contig_length == contig.length:
                    num_matching_contigs += 1
                else:
                    diff = f"VCF contig '{vcf_contig}' length: {vcf_contig_length} not equal to build: {contig.length}"
                    break
            except KeyError:
                num_missing_contigs += 1  # OK to miss some strange contigs

        if diff:
            build_contig_diffs[genome_build.name] = diff
        else:
            if num_matching_contigs > min_contig_matches:
                potential_genome_builds.append(genome_build)
            else:
                matches = f"Matched: {num_matching_contigs} (min: {min_contig_matches}) Missing: {num_missing_contigs}"
                build_contig_diffs[genome_build.name] = matches

    num_matches = len(potential_genome_builds)
    if num_matches:
        if num_matches == 1:
            return potential_genome_builds[0]
        else:
            matches = ", ".join(potential_genome_builds)
            raise ContigMismatchException(f"Multiple genomes: {matches} matched supplied contigs {contig_lengths}")
    else:
        differences = ", ".join([f"{gb}: {diff}" for gb, diff in build_contig_diffs.items()])
        raise ContigMismatchException(f"No genomes matched supplied contigs. Differences were: {differences}")

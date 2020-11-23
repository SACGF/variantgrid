from django.conf import settings
from io import TextIOWrapper
import celery
import logging
import os
import shutil
import subprocess
import tempfile

from snpdb.variants_to_vcf import write_qs_to_vcf_file_sort_alphabetically
from library.log_utils import log_traceback, get_traceback
from snpdb.models import VariantCollection, Variant, VCFBedIntersection
from snpdb.models.models_enums import ProcessingStatus


def create_vcf_bed_intersection(vcf, enrichment_kit):
    name = f"{vcf.pk}_{enrichment_kit}"
    kwargs = {'name': name,
              'vcf': vcf,
              'genomic_intervals': enrichment_kit.genomic_intervals,
              'left_padding': settings.DEFAULT_ENRICHMENT_KIT_LEFT_PADDING,
              'right_padding': settings.DEFAULT_ENRICHMENT_KIT_RIGHT_PADDING}

    pbi, created = VCFBedIntersection.objects.get_or_create(**kwargs)
    if created:
        logging.warning("VCF %s created VCFBedIntersection %s for enrichment_kit %s", vcf, pbi.pk, enrichment_kit)
        vcf_bed_intersection_task.apply_async((pbi.pk,))
    else:
        logging.warning("VCF %s already has VCFBedIntersection for enrichment_kit %s", vcf, enrichment_kit)


def create_backend_vcf_bed_intersections(backend_vcf):
    """ Create VCFBedIntersection for enrichment_kits from sequencing samples, and launches async jobs """
    try:
        vcf = backend_vcf.vcf
        vcf.cohort.cohort_genotype_collection  # Error and skip if missing
        enrichment_kits = backend_vcf.sample_sheet.get_sample_enrichment_kits()
        print(f"enrichment_kits = {enrichment_kits}")

        for enrichment_kit in enrichment_kits:
            if enrichment_kit.genomic_intervals:
                create_vcf_bed_intersection(vcf, enrichment_kit)
            else:
                logging.warning("Enrichment Kit %s has no genomic intervals", enrichment_kit)
    except:
        log_traceback()


@celery.task
def vcf_bed_intersection_task(vcf_bed_intersection_id):
    vbi = VCFBedIntersection.objects.get(pk=vcf_bed_intersection_id)
    vbi.status = ProcessingStatus.PROCESSING
    vbi.save()
    temp_dir_to_delete = None

    try:
        cohort_genotype_collection = vbi.vcf.cohort.cohort_genotype_collection
        bed_filename = vbi.genomic_intervals.processed_file
        name = f"{vbi.vcf} intersect {bed_filename}"

        if vbi.left_padding or vbi.right_padding:
            temp_dir_to_delete = tempfile.mkdtemp(prefix="bed_padding")
            padding_text = f"l_{vbi.left_padding}_r_{vbi.right_padding}"
            name += " " + padding_text
            padded_name = ".".join([str(vcf_bed_intersection_id), padding_text, "bed"])
            padded_bed_filename = os.path.join(temp_dir_to_delete, padded_name)

            reference_fasta_index = vbi.vcf.genome_build.genome_fasta.index_filename
            with open(padded_bed_filename, "w") as padded_bed_file:
                args = ["slopBed", "-i", bed_filename, '-g', reference_fasta_index]
                if vbi.left_padding:
                    args.extend(['-l', f"{vbi.left_padding}"])
                if vbi.right_padding:
                    args.extend(['-r', f"{vbi.right_padding}"])

                logging.info("exec: '%s' > '%s'", str(args), padded_bed_filename)

                slop_bed_pipe = subprocess.Popen(args, stdout=padded_bed_file)
                slop_bed_pipe.communicate()

            bed_filename = padded_bed_filename

        variant_collection = VariantCollection(name=name)
        variant_collection.save()

        script_name = os.path.join(settings.BASE_DIR, 'scripts', 'intersect_bed_and_upload_variant_collection.sh')

        # Open a pipe to intersectBed, which is ALSO going to upload
        args = [script_name, bed_filename, str(variant_collection.pk)]
        intercept_bed_pipe = subprocess.Popen(args, stdin=subprocess.PIPE, stderr=subprocess.PIPE)

        try:
            vcf_queryset = Variant.objects.filter(cohortgenotype__collection=cohort_genotype_collection)
            text_pipe_in = TextIOWrapper(intercept_bed_pipe.stdin, encoding='utf-8', write_through=True)
            write_qs_to_vcf_file_sort_alphabetically(vcf_queryset, text_pipe_in)
            intercept_bed_pipe.communicate()
        except:
            poll = intercept_bed_pipe.poll()
            if poll is None:
                intercept_bed_pipe.communicate()

            stderr_output = intercept_bed_pipe.stderr.read()
            logging.error(stderr_output)
            vbi.status = ProcessingStatus.ERROR
            vbi.error_exception = stderr_output
            raise

        vbi.variant_collection = variant_collection
        vbi.status = ProcessingStatus.SUCCESS
    except:
        tb = get_traceback()
        vbi.error_exception = vbi.error_exception or '' + tb

    vbi.save()

    if temp_dir_to_delete:
        shutil.rmtree(temp_dir_to_delete)

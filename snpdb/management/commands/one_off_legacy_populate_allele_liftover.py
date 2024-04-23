import logging
from collections import defaultdict

from django.core.management.base import BaseCommand
from django.db.models import Q, Min, F
from django.db.models.functions import Length

from annotation.models import AnnotationRangeLock, ClinVar
from genes.hgvs import HGVSMatcher
from library.guardian_utils import admin_bot
from snpdb.models import Variant, Sequence, GenomeBuild, Locus, Allele, ClinGenAllele, Contig, VariantAllele, \
    AlleleLiftover, LiftoverRun, VariantAlleleCollectionSource, AlleleConversionTool, ProcessingStatus


class Command(BaseCommand):
    """
        We originally tried to avoid making a record for each Allele that was lifted over

        Now we want to make an AlleleLiftover for each Allele - that shows what was responsible for it

    """
    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true')

    def handle(self, *args, **options):
        # By taking the first VariantAllele we'll get the original ones (ie not lifted over)
        allele_non_liftover_build = defaultdict(set)
        for allele_id, genome_build in Allele.objects.all().annotate(first_va_id=Min("variantallele__pk")).filter(
                variantallele__pk=F("first_va_id")).values_list("pk", "variantallele__genome_build"):
            allele_non_liftover_build[allele_id].add(genome_build)

        genome_builds = [GenomeBuild.objects.get(name=name) for name in ["GRCh37", "GRCh38"]]
        records = []

        # Handle ClinGen API errors
        logging.info("Looking for ClinGen API errors...")
        clingen_api_errors = {}  # Store for later error lookup
        for va in VariantAllele.objects.filter(clingen_error__isnull=False):
            clingen_api_errors[va.allele] = va.clingen_error
            try:
                other_va = VariantAllele.objects.filter(allele=va.allele).exclude(pk=va.pk).get()
                al = AlleleLiftover(liftover=va.liftover, allele=va.allele,
                                    error=other_va.clingen_error, status=ProcessingStatus.ERROR)
                records.append(al)
            except VariantAllele.DoesNotExist:
                pass

        # We need to represent a "failed clingen" LiftoverRun - just make 1 and link all fails we don't actually run
        # going forward we'll
        failed_clingen_liftovers_by_build = {}
        for genome_build in genome_builds:
            # TODO: Should we only have 1 failed ClinGen allele liftover?
            # Should we have notes on liftover run?
            allele_source = VariantAlleleCollectionSource.objects.create(genome_build=genome_build)
            lr = LiftoverRun.objects.create(user=admin_bot(),
                                            allele_source=allele_source,
                                            genome_build=genome_build,
                                            conversion_tool=AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)
            failed_clingen_liftovers_by_build[genome_build] = lr


        # Get the ones that were never lifted over by ClinGen because they couldn't be (eg missing for that build)
        logging.info("Looking for those unable to be lifted over via ClinGen...")
        for allele in Allele.objects.all().select_related("clingen_allele"):
            for genome_build in genome_builds:
                if genome_build.name in allele_non_liftover_build[allele.pk]:
                    continue  # 1st varinat for allele, not lifted over
                failed_clingen_liftover = failed_clingen_liftovers_by_build[genome_build]

                if cga := allele.clingen_allele:
                    try:
                        _ = cga.get_g_hgvs(genome_build)  # Fails if missing
                    except (ClinGenAllele.ClinGenBuildNotInResponseError, Contig.ContigNotInBuildError) as e:
                        al = AlleleLiftover(liftover=failed_clingen_liftover, allele=allele,
                                            error={"message": str(e)}, status=ProcessingStatus.ERROR)
                        records.append(al)
                else:
                    error = clingen_api_errors.get(allele.pk)
                    if error is None:
                        error = {"message": "No ClinGenAllele"}
                    al = AlleleLiftover(liftover=failed_clingen_liftover, allele=allele,
                                        error=error, status=ProcessingStatus.ERROR)
                    records.append(al)

        known_missing_upload_pipelines = set(failed_clingen_liftovers_by_build.values())

        num_liftover_runs = LiftoverRun.objects.all().count()

        # Insert the 1st for each status (success/failure)
        # Or should we create a new FailedClinGen liftover run all the time? But we don't want to run it...
        allele_conversion_tool_initial_liftover = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        for i, lr in enumerate(LiftoverRun.objects.all().order_by("created")):
            if i % 10 == 0:
                logging.info("Assigning alleles to Liftover runs - %d/%d", i, num_liftover_runs)
            try:
                status = lr.uploadedliftover.uploaded_file.uploadpipeline.status
            except Exception as e:
                if lr not in known_missing_upload_pipelines:
                    logging.error("Couldn't get upload pipeline status for run: %s", lr)
                continue
            allele_qs = lr.get_allele_qs()

            # TODO: We could optimise here - find all the alleles that have already been handled by this
            # Run before/after to verify that it's the same...
            if False:
                handled_allele_ids = set()
                for allele_id, genome_build_set in allele_non_liftover_build.items():
                    if lr.genome_build_id in genome_build_set:
                        handled_allele_ids.add(allele_id)

                for allele, data in allele_conversion_tool_initial_liftover.items():
                    if data.get(lr.genome_build_id, {}).get(status, {}).get(lr.conversion_tool):
                        handled_allele_ids.add(allele.pk)

                allele_qs = allele_qs.exclude(pk__in=handled_allele_ids)

            for allele in allele_qs:
                if lr.genome_build_id not in allele_non_liftover_build[allele.pk]:
                    act = allele_conversion_tool_initial_liftover[allele][lr.genome_build_id][status]
                    if lr.conversion_tool not in act:
                        act[lr.conversion_tool] = lr
                        error = None
                        if status != ProcessingStatus.SUCCESS:
                            error = {"message": f"LiftoverRun {lr.pk} failed"}
                        al = AlleleLiftover(liftover=lr, allele=allele, error=error, status=status)
                        records.append(al)

        # The way things happened historically is:
        # * Liftover had an "allele source"
        # * We checked the alleles for whether it had already been lifted over - if not it was just skipped
        # * Then in liftover we called conversion_tool, variant_id_or_coordinate = allele.get_liftover_tuple(genome_build) to determine what to run
        # * If the liftover worked, we then populated things and set the 

        # How are we going to handle the AllClassifications - just do them all again and again?


        # What if there are some left that have been lifted over, but don't have anything with them??

        if records:
            logging.info("Creating %d records...", len(records))
            AlleleLiftover.objects.bulk_create(records, batch_size=1000)

import logging
from collections import defaultdict

from django.core.management.base import BaseCommand
from django.db.models import Min, F, Count
from library.guardian_utils import admin_bot
from snpdb.models import GenomeBuild, Allele, ClinGenAllele, Contig, VariantAllele, \
    AlleleLiftover, LiftoverRun, VariantAlleleCollectionSource, AlleleConversionTool, ProcessingStatus


class Command(BaseCommand):
    """
        We originally tried to avoid making a record for each Allele that was lifted over

        Now we want to make an AlleleLiftover for each Allele - that shows what was responsible for it

    """
    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true')

    def handle(self, *args, **options):
        use_optimisation_exclusion = False  # Optimisation attempt for lots of allele sources

        if num_non_errors := AlleleLiftover.objects.exclude(status=ProcessingStatus.ERROR).count():
            raise ValueError(f"Error: {num_non_errors} non-error AlleleLiftover records existing - was this run tiwce?")

        # By taking the first VariantAllele we'll get the original ones (ie not lifted over)
        allele_non_liftover_build = defaultdict(set)
        for allele_id, genome_build in Allele.objects.all().annotate(first_va_id=Min("variantallele__pk")).filter(
                variantallele__pk=F("first_va_id")).values_list("pk", "variantallele__genome_build"):
            allele_non_liftover_build[allele_id].add(genome_build)

        genome_builds = [GenomeBuild.objects.get(name=name) for name in ["GRCh37", "GRCh38"]]
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

        records = []

        # Handle ClinGen API errors
        logging.info("Looking for ClinGen API errors...")
        clingen_api_errors = {}  # Store for later error lookup
        for va in VariantAllele.objects.filter(clingen_error__isnull=False):
            clingen_api_errors[va.allele] = va.clingen_error

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
        # allele/genome_build/status/tool
        allele_conversion_tool_initial_liftover = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        existing_allele_id_liftovers = defaultdict(set)
        # Load existing AlleleLiftover - these are errors...
        for al in AlleleLiftover.objects.all():
            existing_allele_id_liftovers[al.allele_id].add(al.liftover)

        liftover_run_qs = LiftoverRun.objects.all().order_by("created")
        # liftover_run_qs = LiftoverRun.objects.none()

        for i, lr in enumerate(liftover_run_qs):
            if i % 10 == 0:
                logging.info("Assigning alleles to Liftover runs - %d/%d", i, num_liftover_runs)
            try:
                status = lr.uploadedliftover.uploaded_file.uploadpipeline.status
            except Exception as e:
                if lr not in known_missing_upload_pipelines:
                    logging.error("Couldn't get upload pipeline status for run: %s", lr)
                continue
            allele_qs = lr.get_allele_qs()

            # Performance really degrades when you have heaps of all allele sources
            # This optimisation is slower for most sources, but tries to reduce worst case when you have lots of
            # all allele sources (meaning you have to loop again and again)
            if use_optimisation_exclusion:
                handled_allele_ids = set()
                for allele_id, genome_build_set in allele_non_liftover_build.items():
                    if lr.genome_build in genome_build_set:
                        handled_allele_ids.add(allele_id)

                for allele_id, data in allele_conversion_tool_initial_liftover.items():
                    if data.get(lr.genome_build, {}).get(status, {}).get(lr.conversion_tool):
                        handled_allele_ids.add(allele_id)

                allele_qs = allele_qs.exclude(pk__in=handled_allele_ids)

            for allele_id in allele_qs.values_list("pk", flat=True):
                if lr in existing_allele_id_liftovers[allele_id]:
                    continue

                if lr.genome_build_id not in allele_non_liftover_build[allele_id]:
                    act = allele_conversion_tool_initial_liftover[allele_id][lr.genome_build][status]
                    if lr.conversion_tool not in act:
                        act[lr.conversion_tool] = lr
                        error = None
                        if status != ProcessingStatus.SUCCESS:
                            error = {"message": f"LiftoverRun {lr.pk} failed"}
                        al = AlleleLiftover(liftover=lr, allele_id=allele_id, error=error, status=status)
                        records.append(al)

        # The way things happened historically is:
        # * Liftover had an "allele source"
        # * We checked the alleles for whether it had already been lifted over - if not it was just skipped
        # * Then in liftover we called conversion_tool, variant_id_or_coordinate = allele.get_liftover_tuple(genome_build) to determine what to run
        # * If the liftover worked, we then populated things and set the 

        if records:
            logging.info("Creating %d records...", len(records))
            AlleleLiftover.objects.bulk_create(records, batch_size=1000)
            records = []

        # Leftover alleles

        logging.info("Assigning any remaining Alleles (linked via web)")

        # There are some Alleles that have VariantAllele for both - but no official liftover run
        # This is done via the web (variant/allele page) when variants from both builds are already there

        allele_wo_liftover = Allele.objects.filter(alleleliftover__isnull=True)
        allele_wo_liftover_both_builds = allele_wo_liftover.annotate(num_builds=Count("variantallele")).filter(num_builds__gt=1)

        allele_linked_first_last_build = defaultdict(lambda: defaultdict(list))
        for allele in allele_wo_liftover_both_builds:
            first_va, second_va, *_ = tuple(allele.variantallele_set.order_by("pk"))
            allele_linked_first_last_build[first_va.genome_build][second_va.genome_build].append(second_va)

        for source_build in genome_builds:
            for dest_build in genome_builds:
                if source_build == dest_build:
                    continue
                if source_to_dest := allele_linked_first_last_build[source_build][dest_build]:
                    logging.info("Leftovers - %s first, linked to %s: %d", source_build, dest_build, len(source_to_dest))
                    allele_source = VariantAlleleCollectionSource.objects.create(genome_build=source_build)
                    lr = LiftoverRun.objects.create(user=admin_bot(),
                                                    allele_source=allele_source,
                                                    source_genome_build=source_build,
                                                    genome_build=dest_build,
                                                    conversion_tool=AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)
                    for va in source_to_dest:
                        al = AlleleLiftover(liftover=lr, allele=va.allele, status=ProcessingStatus.SUCCESS)
                        records.append(al)

        if records:
            logging.info("Creating %d records...", len(records))
            AlleleLiftover.objects.bulk_create(records, batch_size=1000)

        if num_left := allele_wo_liftover.count():
            logging.info("%d Alleles left without any liftover (could be from variant page)", num_left)

        # Delete any LiftoverRuns
        unused_liftover_runs_qs = LiftoverRun.objects.all().annotate(num_alleles=Count("alleleliftover")).filter(num_alleles=0)
        unused_liftover_runs_qs.delete()
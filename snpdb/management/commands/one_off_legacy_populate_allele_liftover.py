from collections import defaultdict

from django.core.management.base import BaseCommand
from django.db.models import Q, Min, F
from django.db.models.functions import Length

from annotation.models import AnnotationRangeLock, ClinVar
from genes.hgvs import HGVSMatcher
from snpdb.models import Variant, Sequence, GenomeBuild, Locus, Allele, ClinGenAllele, Contig, VariantAllele, \
    AlleleLiftover, LiftoverRun


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

        # We need to represent a "failed clingen" LiftoverRun - just make 1 and link all fails we don't actually run
        # going forward we'll

        # Now find all alleles that would have failed clingen
        clingen_api_errors = {}
        for allele, error in VariantAllele.objects.filter(error__isnull=False).values_list("allele", "error"):
            clingen_api_errors[allele] = error

        genome_builds = [GenomeBuild.objects.get(name=name) for name in ["GRCh37", "GRCh38"]]
        clingen_failures = []

        for allele in Allele.objects.all():
            for genome_build in genome_builds:
                if genome_build.name in allele_non_liftover_build[allele.pk]:
                    continue  # 1st varinat for allele, not lifted over
                if cga := allele.clingen_allele:
                    try:
                        g_hgvs = cga.get_g_hgvs(genome_build)
                    except (ClinGenAllele.ClinGenBuildNotInResponseError, Contig.ContigNotInBuildError) as e:
                        clingen_failures.append((allele, genome_build.name, {"message": str(e)}))
                else:
                    error = clingen_api_errors.get(allele.pk)
                    if error is None:
                        error = "No ClinGenAllele"

                    clingen_failures.append((allele, genome_build.name, error))

        # Or should we create a new FailedClinGen liftover run all the time? But we don't want to run it...

        clingen_failures = []
        for va in VariantAllele.objects.filter(clingen_error__isnull=False):
            try:
                other_va = VariantAllele.objects.filter(allele=va.allele).exclude(pk=va.pk).get()

                error = {}  # TODO
                al = AlleleLiftover(liftover=va.liftover, allele=va.allele, error=error)
                clingen_failures.append(al)
            except VariantAllele.DoesNotExist:
                pass

        # going to use for different purpose
        handled_allele_builds = allele_non_liftover_build

        allele_liftovers = defaultdict(set)
        for lr in LiftoverRun.objects.all().order_by("created"):
            allele_qs = lr.get_allele_qs()
            for allele in allele_qs:
                if lr.genome_build_id not in handled_allele_builds[allele.pk]:
                    allele_liftovers[allele].add(lr)
                    handled_allele_builds[allele.pk].add(lr.genome_build_id)

        # The way things happened historically is:
        # * Liftover had an "allele source"
        # * We checked the alleles for whether it had already been lifted over - if not it was just skipped
        # * Then in liftover we called conversion_tool, variant_id_or_coordinate = allele.get_liftover_tuple(genome_build) to determine what to run
        # * If the liftover worked, we then populated things and set the 

        # How are we going to handle the AllClassifications - just do them all again and again?


        # What if there are some left that have been lifted over, but don't have anything with them??

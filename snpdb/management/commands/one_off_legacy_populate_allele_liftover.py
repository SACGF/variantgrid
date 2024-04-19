from django.core.management.base import BaseCommand
from django.db.models import Q
from django.db.models.functions import Length

from annotation.models import AnnotationRangeLock, ClinVar
from genes.hgvs import HGVSMatcher
from snpdb.models import Variant, Sequence, GenomeBuild, Locus


class Command(BaseCommand):
    """
        We originally tried to avoid making a record for each Allele that was lifted over

        Now we want to make an AlleleLiftover for each Allele - that shows what was responsible for it

    """
    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true')

    def handle(self, *args, **options):
        # Open questions:
        # How did we trigger NCBI liftover running again?


        # The way things happened historically is:
        # * Liftover had an "allele source"
        # * We checked the alleles for whether it had already been lifted over - if not it was just skipped
        # * If the liftover worked, we then populated things

        # How are we going to handle the AllClassifications - just do them all again and again?


        # What if there are some left that have been lifted over, but don't have anything with them??

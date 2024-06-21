import operator
from functools import reduce

from django.core.management import BaseCommand
from django.db.models import Q
from hgvs.exceptions import HGVSInvalidVariantError
from pyhgvs import InvalidHGVSName

from classification.classification_import import reattempt_variant_matching
from classification.models.classification_variant_info_models import ImportedAlleleInfo, \
    ImportedAlleleInfoStatus
from genes.hgvs import HGVSMatcher, VariantResolvingError
from library.guardian_utils import admin_bot
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    """
        We introduced a check, where "NM_000127.2(EXT1):c.1315_1320dupTTTTTTT"

        Now gives error:

        Calculated reference 'TCACGT' length (6) different from provided reference 'TTTTTTT' length (7)

        ImportedAlleleInfo caches hgvs resolution, so need to find these historical ones and re-validate them

    """

    def handle(self, *args, **options):
        genome_build_matchers = {}

        for genome_build in GenomeBuild.builds_with_annotation():
            genome_build_matchers[genome_build] = HGVSMatcher(genome_build)

        regex_span_plus_provided = r"\d+_\d+(del|dup)[GATC]+"
        filters = [
            Q(imported_c_hgvs__regex=regex_span_plus_provided),
            Q(imported_g_hgvs__regex=regex_span_plus_provided)
        ]
        q = ~Q(status=ImportedAlleleInfoStatus.FAILED) & reduce(operator.or_, filters)
        iai_qs = ImportedAlleleInfo.objects.all()
        iai_ids_with_hgvs_errors = set()
        for iai in iai_qs.filter(q).iterator():
            matcher = genome_build_matchers[iai.imported_genome_build]
            hgvs_name = iai.imported_c_hgvs or iai.imported_g_hgvs
            try:
                matcher.get_variant_coordinate(hgvs_name)
            except (HGVSInvalidVariantError, InvalidHGVSName) as e:
                print(f"{hgvs_name}: {e}")
                iai_ids_with_hgvs_errors.add(iai.pk)
            except:
                pass  # If things fail for any other reason, we can't help - so no point re-matching

        rematch_iai_qs = ImportedAlleleInfo.objects.filter(pk__in=iai_ids_with_hgvs_errors)
        print(f"Have {rematch_iai_qs.count()} ImportAlleleIDs to fix")
        user = admin_bot()
        reattempt_variant_matching(user, rematch_iai_qs)

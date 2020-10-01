from django.core.management.base import BaseCommand

from genes.hgvs import HGVSMatcher
from snpdb.models.models_genome import GenomeBuild
from classification.models.classification import Classification


class Command(BaseCommand):

    def handle(self, *args, **options):
        build_matchers = {b: HGVSMatcher(b) for b in GenomeBuild.builds_with_annotation()}

        vc: Classification
        for vc in Classification.objects.filter(variant__isnull=False):
            orig_allele = vc.variant.variantallele
            allele = orig_allele.allele
            transcript = vc.transcript
            if transcript:
                try:
                    matcher = build_matchers[orig_allele.genome_build]
                    original_hgvs = matcher.variant_to_c_hgvs_parts(vc.variant, transcript, throw_on_issue=True).full_c_hgvs
                    for va in allele.variantallele_set.exclude(pk=orig_allele.pk):
                        matcher = build_matchers[va.genome_build]
                        converted_hgvs = matcher.variant_to_c_hgvs_parts(va.variant, transcript, throw_on_issue=True).full_c_hgvs
                        if original_hgvs != converted_hgvs:
                            cols = [vc.pk, vc.variant.pk, orig_allele.genome_build.name, original_hgvs, converted_hgvs]
                            print("\t".join((str(c) for c in cols)))
                except Exception as e:
                    print(e)

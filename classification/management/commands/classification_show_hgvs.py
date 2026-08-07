from django.core.management.base import BaseCommand

from classification.models.classification import Classification
from genes.hgvs import HGVSMatcher


class Command(BaseCommand):

    def handle(self, *args, **options):
        vc: Classification
        for vc in Classification.objects.filter(variant__isnull=False):
            orig_allele = vc.variant.variantallele_set.filter(genome_build=vc.get_genome_build())
            allele = orig_allele.allele
            transcript = vc.transcript
            if transcript:
                try:
                    matcher = HGVSMatcher.instance(orig_allele.genome_build)
                    original_hgvs = matcher.variant_to_c_hgvs_parts(vc.variant, transcript, throw_on_issue=True).full_c_hgvs
                    for va in allele.variantallele_set.exclude(pk=orig_allele.pk):
                        matcher = HGVSMatcher.instance(va.genome_build)
                        converted_hgvs = matcher.variant_to_c_hgvs_parts(va.variant, transcript, throw_on_issue=True).full_c_hgvs
                        if original_hgvs != converted_hgvs:
                            cols = [vc.pk, vc.variant.pk, orig_allele.genome_build.name, original_hgvs, converted_hgvs]
                            print("\t".join((str(c) for c in cols)))
                except Exception as e:
                    print(e)

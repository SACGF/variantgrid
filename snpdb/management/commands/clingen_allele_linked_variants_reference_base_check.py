from datetime import date

from django.core.management import BaseCommand

from snpdb.models import Variant, Allele, ClinGenAllele, Contig


class Command(BaseCommand):
    """ Indel representation in g.HGVS doesn't have reference base - so we may have variants with
        different reference bases linked to an allele """

    def handle(self, *args, **options):
        variants_with_clingen = Variant.objects.filter(variantallele__allele__clingen_allele__isnull=False)
        indel_qs = variants_with_clingen.exclude(Variant.get_snp_q())
        indel_alleles = indel_qs.values_list("variantallele__allele_id", flat=True)

        alleles_qs = Allele.objects.filter(clingen_allele__isnull=False, pk__in=indel_alleles).distinct()
        print(f"{alleles_qs.count()} indel alleles to check...")
        for i, allele in enumerate(alleles_qs):
            if i % 100 == 0:
                print(f"processed {i} alleles")
            for va in allele.variantallele_set.all():
                existing_vc = va.variant.coordinate
                try:
                    clingen_vc = allele.clingen_allele.get_variant_coordinate(va.genome_build)
                    if existing_vc != clingen_vc:
                        print(f"{allele} has variant {repr(existing_vc)} not matching expected for build: {repr(clingen_vc)}")

                except (ClinGenAllele.ClinGenBuildNotInResponseError, Contig.ContigNotInBuildError):
                    pass
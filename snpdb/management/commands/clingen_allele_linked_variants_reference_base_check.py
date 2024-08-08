import logging

from django.core.management import BaseCommand

from snpdb.models import Variant, Allele, ClinGenAllele, Contig, AlleleLiftover


class Command(BaseCommand):
    """ Indel representation in g.HGVS doesn't have reference base - so we may have variants with
        different reference bases linked to an allele """
    def add_arguments(self, parser):
        parser.add_argument('--dry-run', help="Just report, don't unlink variants from allele if incorrect", action='store_true')

    def handle(self, *args, **options):
        dry_run = options["dry_run"]

        variants_with_clingen = Variant.objects.filter(variantallele__allele__clingen_allele__isnull=False)
        indel_qs = variants_with_clingen.exclude(Variant.get_snp_q())
        indel_alleles = indel_qs.values_list("variantallele__allele_id", flat=True)

        alleles_qs = Allele.objects.filter(clingen_allele__isnull=False, pk__in=indel_alleles).distinct()
        logging.info(f"{alleles_qs.count()} indel alleles to check...")
        if dry_run:
            logging.info("Dry-run only - will not unlink incorrect variants...")
        else:
            logging.info("Going to unlink incorrect variants...")

        num_unlinked = 0
        for i, allele in enumerate(alleles_qs):
            if i % 500 == 0:
                logging.info(f"processed {i} alleles")
            for va in allele.variantallele_set.all():
                existing_vc = va.variant.coordinate
                try:
                    clingen_vc = allele.clingen_allele.get_variant_coordinate(va.genome_build)
                    if existing_vc != clingen_vc:
                        logging.info(f"{allele} has variant {repr(existing_vc)} not matching expected for build {va.genome_build}: {repr(clingen_vc)}")
                        if not dry_run:
                            liftover_res = AlleleLiftover.objects.filter(variant_allele=allele,
                                                                         liftover__genome_build=va.genome_build).delete()
                            logging.info(f"Removing liftover record: %s", liftover_res)
                            va_res = va.delete()
                            logging.info(f"Unlinking variant: %s", va_res)

                except (ClinGenAllele.ClinGenBuildNotInResponseError, Contig.ContigNotInBuildError):
                    pass

        if num_unlinked:
            logging.info("Unlinked %d variants, you'll need to re-run liftover to re-build these", num_unlinked)

from datetime import date

from django.core.management import BaseCommand
from django.db.models import Q
from django.utils import timezone

from snpdb.clingen_allele import ClinGenAlleleRegistryAPI
from snpdb.models import ClinGenAllele, Variant


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--indels', action='store_true', help='Only retrieve indels')
        parser.add_argument('--before-date', type=date.fromisoformat, help='Date (ISO)')

    def handle(self, *args, **options):
        indels = options["indels"]
        before_date = options.get("before_date")

        filter_message = []
        cga_qs = ClinGenAllele.objects.all()
        if before_date:
            filter_message.append(f"Earlier than {before_date}")
            cga_qs = cga_qs.filter(modified__lte=before_date)

        if indels:
            filter_message.append(f"(Indels only)")
            variants_with_clingen = Variant.objects.filter(variantallele__allele__clingen_allele__in=cga_qs)
            indel_qs = variants_with_clingen.filter(Q(locus__ref__length__gt=1) | Q(alt__length__gt=1))
            values_qs = indel_qs.values_list("variantallele__allele__clingen_allele__id", flat=True)
            cga_qs = ClinGenAllele.objects.filter(pk__in=values_qs)

        clingen_allele_by_id = {str(ca): ca for ca in cga_qs}

        print(f"There are {len(clingen_allele_by_id)} CA IDs {', '.join(filter_message)}")

        clingen_api = ClinGenAlleleRegistryAPI()
        records = clingen_allele_by_id.keys()
        clingen_to_update = []
        now = timezone.now()
        for api_response in clingen_api.hgvs_put(records, file_type="id"):
            clingen_allele_id = ClinGenAllele.get_id_from_response(api_response)
            ca = clingen_allele_by_id[ClinGenAllele.format_clingen_allele(clingen_allele_id)]
            ca.api_response = api_response
            ca.modified = now
            clingen_to_update.append(ca)

        if clingen_to_update:
            print(f"Updating {len(clingen_to_update)} ClinGenAllele records...")
            ClinGenAllele.objects.bulk_update(clingen_to_update, ["api_response", "modified"], batch_size=2000)

from django.core.management import BaseCommand

from classification.models import ClinVarExport
from classification.models.clinvar_export_prepare import ClinvarAlleleExportPrepare
from snpdb.models import Allele


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--prepare',  action='store_true', default=False)
        parser.add_argument('--mondo', action='store_true', default=False)

    def handle(self, *args, **options):
        if options["prepare"]:
            for count, allele in enumerate(Allele.objects.all()):
                report = ClinvarAlleleExportPrepare(allele).update_export_records()
                if count % 100 == 0:
                    print(f"Processing Allele no {count}")
                # print(report)
            print(f"Completed {count}")

        if options["mondo"]:
            for clinvar_export in ClinVarExport.objects.filter(condition__display_text__istartswith="OMIM"):
                omim_condition = clinvar_export.condition_resolved
                mondo_condition = omim_condition.as_mondo_if_possible()
                if mondo_condition.mondo_term:
                    clinvar_export.condition_resolved = mondo_condition
                    clinvar_export.save()
                    print(f"ClinVarExport {clinvar_export.pk} Converted {omim_condition} -> {mondo_condition}")

from django.core.management import BaseCommand

from classification.models import ClinVarExport
from classification.models.clinvar_export_prepare import ClinvarExportPrepare
from snpdb.models import ClinVarKey


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--key', type=str, required=False)
        parser.add_argument('--all',  action='store_true', default=False)
        parser.add_argument('--mondo', action='store_true', default=False)

    def handle(self, *args, **options):
        if options["all"]:
            print("Preparing all ClinVarExports, this may take a while")
            report = ClinvarExportPrepare().update_export_records()
            print("\n".join(report))
        elif key_name := options["key"]:
            if clinvar_key := ClinVarKey.objects.filter(id=key_name).first():
                print(f"Updating ClinVarKey {key_name}")
                ClinvarExportPrepare.update_export_records_for_keys({clinvar_key})
            else:
                print(f"No ClinVar Key with id {key_name}")
                print("Values are:")
                for any_clinvar_key in ClinVarKey.objects.all():
                    print(any_clinvar_key.id)

        if options["mondo"]:
            for clinvar_export in ClinVarExport.objects.filter(condition__display_text__istartswith="OMIM"):
                omim_condition = clinvar_export.condition_resolved
                mondo_condition = omim_condition.as_mondo_if_possible()
                if mondo_condition.mondo_term:
                    clinvar_export.condition_resolved = mondo_condition
                    clinvar_export.save()
                    print(f"ClinVarExport {clinvar_export.pk} Converted {omim_condition} -> {mondo_condition}")

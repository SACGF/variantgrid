from django.db import migrations
from django.db.models import Q

from manual.operations.manual_operations import ManualOperation


def _has_clinvar_somatic_to_backfill(apps):
    """ ClinVar has stored SCI / ONC text since annotation 0102, but only imports from 0176 on derive
        somatic_tier / highest_oncogenicity from it - so loaded versions need the values calculated """
    GenomeBuild = apps.get_model("snpdb", "GenomeBuild")
    ClinVarVersion = apps.get_model("annotation", "ClinVarVersion")
    ClinVar = apps.get_model("annotation", "ClinVar")

    for genome_build in GenomeBuild.objects.filter(enabled=True):
        version_qs = ClinVarVersion.objects.filter(genome_build=genome_build).order_by('annotation_date')
        if clinvar_version := version_qs.last():
            derivable_qs = ClinVar.objects.filter(
                version=clinvar_version, somatic_tier__isnull=True, highest_oncogenicity=0
            ).filter(Q(somatic_clinical_significance__isnull=False) | Q(oncogenic_classification__isnull=False))
            if derivable_qs.exists():
                return True
    return False


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0176_clinvar_somatic_tier_and_highest_oncogenicity"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_backfill_clinvar_somatic"]),
                        note="Derive ClinVar somatic tier / oncogenicity from stored SCI and ONC text (#1791)",
                        test=_has_clinvar_somatic_to_backfill),
    ]

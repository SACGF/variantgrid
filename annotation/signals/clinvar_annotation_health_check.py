from datetime import timedelta

from django.dispatch import receiver

from annotation.models import ClinVarVersion
from library.health_check import health_check_signal, HealthCheckRequest, HealthCheckAge
from snpdb.models import GenomeBuild


@receiver(signal=health_check_signal)
def ontology_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    checks = list()
    for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
        latest_clinvar = ClinVarVersion.objects.filter(genome_build=genome_build).order_by(
                '-annotation_date').values_list('annotation_date', flat=True).first()
        checks.append(
            HealthCheckAge(
                name=f"ClinVar {genome_build}",
                now=health_request.now,
                warning_age=timedelta(days=60),
                last_performed=latest_clinvar
            )
        )
    return checks
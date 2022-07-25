from datetime import timedelta
from typing import Optional

from django.dispatch import receiver
from library.health_check import health_check_signal, HealthCheckRequest, HealthCheckAge
from ontology.models import OntologyImport


@receiver(signal=health_check_signal)
def ontology_health_check(sender, health_request: HealthCheckRequest, **kwargs):
    checks = list()
    warning_age = timedelta(days=60)
    for contexts in [
        ("mondo_file", "MONDO"),
        # OMIM can be loaded from 2 difference sources depending on OMIM license
        [("omim_file", "OMIM"), ("biomart_omim_aliases", "OMIM - via biomart")],
        ("hpo_file", "HPO"), ("gencc_file", "GenCC")]:
        if not isinstance(contexts, list):
            contexts = [contexts]

        check: Optional[HealthCheckAge] = None
        for context, label in contexts:
            if last_import := OntologyImport.objects.filter(
                    context=context, completed=True, processed_date__lte=health_request.now
            ).order_by('-processed_date').values_list('processed_date', flat=True).first():
                check = HealthCheckAge(
                    name=label,
                    now=health_request.now,
                    last_performed=last_import,
                    warning_age=warning_age
                )
                break
        if not check:
            check = HealthCheckAge(
                    name=contexts[0][1],
                    now=health_request.now,
                    last_performed=None,
                    warning_age=warning_age
                )
        checks.append(check)
    return checks

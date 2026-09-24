"""
deployment_check "Liftover": warns, per annotated build, about Alleles that have never been lifted over to it -
what adding a build (eg T2T) leaves behind until someone presses its button on the liftover page.
Entry point: check_alleles_never_lifted_over
"""
from django.urls import reverse

from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import Allele

# Counting stops here - a new build is missing every Allele, and "over N" says as much as the exact number
MAX_ALLELES_COUNTED = 10_000


def check_alleles_never_lifted_over() -> dict:
    liftover_checks = {}
    genome_builds = list(GenomeBuild.builds_with_annotation())
    if len(genome_builds) < 2:
        return liftover_checks

    liftover_url = reverse("liftover_runs")
    for genome_build in genome_builds:
        data = {"valid": True}  # Just a warning
        allele_qs = Allele.liftover_never_attempted_for_build(genome_build)
        if count := allele_qs[:MAX_ALLELES_COUNTED].count():
            if count == MAX_ALLELES_COUNTED:
                num_alleles = f"Over {MAX_ALLELES_COUNTED:,}"
            else:
                num_alleles = f"{count:,}"
            data["warning"] = f"{num_alleles} Alleles have no {genome_build} variant and have never been lifted " \
                              f"over to it. Go to the liftover page ({liftover_url}) and click " \
                              f"'Liftover variants' for {genome_build}"
        liftover_checks[f"Alleles lifted over to {genome_build}"] = data
    return liftover_checks

import logging

import celery
from django.contrib.auth.models import User

from library.log_utils import log_traceback
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import GenomeBuild, ImportSource, Allele


@celery.shared_task
def liftover_alleles(username, genome_build_name):
    user = User.objects.get(username=username)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    alleles = Allele.missing_variants_for_build(genome_build)

    for other_build in GenomeBuild.builds_with_annotation():
        if not GenomeBuild.is_equivalent(genome_build, other_build):
            try:
                logging.info("creating liftover pipelines from %s to %s", other_build, genome_build)
                create_liftover_pipelines(user, alleles, ImportSource.WEB, other_build, [genome_build])
                logging.info("/ finished creating pipelines")
            except Exception as e:
                log_traceback()

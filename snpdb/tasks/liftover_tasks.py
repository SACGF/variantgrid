import logging

import celery
from django.contrib.auth.models import User

from snpdb.liftover import create_liftover_pipelines
from snpdb.models import GenomeBuild, ImportSource, Allele


@celery.shared_task
def liftover_alleles(username, genome_build_name):
    user = User.objects.get(username=username)
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

    other_build = {
        "GRCh37": GenomeBuild.grch38(),
        "GRCh38": GenomeBuild.grch37(),
    }
    inserted_genome_build = other_build[genome_build_name]
    inserted_genome_build: GenomeBuild
    alleles = Allele.missing_variants_for_build(genome_build)
    logging.info("creating pipelines...")
    try:
        create_liftover_pipelines(user, alleles, ImportSource.WEB, inserted_genome_build, [genome_build])
        logging.info("/ finished creating pipelines")
    except Exception as e:
        logging.error(e)


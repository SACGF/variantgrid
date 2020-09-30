#!/usr/bin/env python3

"""
Monarch Disease Ontology

An ontology that harmonizes multiple disease resources.

http://www.obofoundry.org/ontology/mondo.html
"""

from django.core.management.base import BaseCommand, CommandError
import logging
import pronto
import re
import time

from annotation.models import MonarchDiseaseOntology

BULK_INSERT_SIZE = 10000
MATCH_PHENOTYPE_MIM = False


def insert_owl_mondo(ot, phenotype_mim_ids):
    start = time.time()

    logging.info("Deleting existing MonarchDiseaseOntology")
    MonarchDiseaseOntology.objects.all().delete()

    pattern = re.compile(r"http://identifiers.org/omim/(\d+)$")

    mondo_list = []

    for term in ot.terms():
        if term.obsolete:
            continue

        if term.id.startswith("MONDO"):
            try:
                mondo_id = MonarchDiseaseOntology.mondo_id_as_int(term.id)
            except ValueError:
                logging.info(f"Could not convert MONDO id '{term.id}' into int")
                continue

            mondo = MonarchDiseaseOntology(pk=mondo_id, name=term.name, definition=term.definition)
            mondo_list.append(mondo)

            if MATCH_PHENOTYPE_MIM:  # No phenotype mim yet
                for exact_match in term.other.get("exactMatch", []):
                    if "omim" in exact_match:
                        if m := pattern.match(exact_match):
                            phenotype_mim_id = m.group(1)
                            if phenotype_mim_id in phenotype_mim_ids:
                                mondo = MonarchDiseaseOntology(pk=mondo_id, phenotype_mim_id=phenotype_mim_id)
                                mondo_list.append(mondo)
                            else:
                                logging.warning("Mondo %s exact match phenotype_mim_id %s doesn't exist!", mondo_id, phenotype_mim_id)

    MonarchDiseaseOntology.objects.bulk_create(mondo_list)

    end = time.time()
    logging.info("Inserted %d Mondo records in %.2f seconds", len(mondo_list), end - start)


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--owl', required=True)

    def handle(self, *args, **options):
        owl_filename = options["owl"]
        if not owl_filename.endswith(".owl"):
            raise CommandError(f"Final arg: '{owl_filename}' must end in .owl")

        # TODO: Need to load phenotype mim ids
        # phenotype_mim_ids = list(PhenotypeMIM.objects.all().values_list("pk", flat=True))  # @UndefinedVariable

        phenotype_mim_ids = []

        ot = pronto.Ontology(owl_filename)
        insert_owl_mondo(ot, phenotype_mim_ids)

#!/usr/bin/env python3
from django.core.management.base import BaseCommand, CommandError
from django.db.models.expressions import F
import itertools
import logging
import pronto
import time

from annotation.management.commands.mim_import import load_diseases_to_genes_to_phenotypes
from annotation.models import HumanPhenotypeOntology, PhenotypeMIM, MIMMorbid
from annotation.models.models_enums import HPOSynonymScope
from annotation.models.models_mim_hpo import HPOSynonym, HPOEdge


DO_INSERT = True
BULK_INSERT_SIZE = 10000


def insert_human_phenotype_ontology(mim_ids, df):
    start = time.time()

    gb = df.groupby(["hpo_id", "hpo_name"])
    hpo_records = []
    phenotype_mim_records = []
    for (hpo_id, name), data in gb:
        hpo = HumanPhenotypeOntology(id=hpo_id, name=name)
        hpo_records.append(hpo)

        for disease_id in set(data["disease_id"]):
            if int(disease_id) not in mim_ids:
                continue
            phenotype_mim = PhenotypeMIM(hpo_id=hpo_id, mim_morbid_id=disease_id)
            phenotype_mim_records.append(phenotype_mim)

    logging.info("Inserting %d HumanPhenotypeOntology records", len(hpo_records))
    logging.info("Inserting %d PhenotypeMIM records", len(phenotype_mim_records))

    if DO_INSERT:
        if hpo_records:
            HumanPhenotypeOntology.objects.bulk_create(hpo_records, batch_size=BULK_INSERT_SIZE, ignore_conflicts=True)

        if phenotype_mim_records:
            PhenotypeMIM.objects.bulk_create(phenotype_mim_records, batch_size=BULK_INSERT_SIZE, ignore_conflicts=True)

    end = time.time()
    logging.info("Insert HPO took %.2f seconds", end - start)


def get_term_hpo_id(term):
    hpo_id = None
    if term.id.startswith(HumanPhenotypeOntology.PREFIX):
        hpo_id = int(term.id.replace(HumanPhenotypeOntology.PREFIX, ""))
    return hpo_id


def insert_owl_hpo(ot):
    start = time.time()
    existing_hpo = set()
    existing_hpo_no_definition = set()

    for hpo_id, definition in HumanPhenotypeOntology.objects.all().values_list("id", "definition"):
        hpo_id = int(hpo_id)
        existing_hpo.add(hpo_id)
        if not definition:
            existing_hpo_no_definition.add(hpo_id)

    hpo_records = []
    update_definitions = {}

    for term in ot.terms():
        hpo_id = get_term_hpo_id(term)
        if hpo_id:
            if hpo_id in existing_hpo:
                if term.definition and hpo_id in existing_hpo_no_definition:
                    update_definitions[hpo_id] = term.definition
            else:
                hpo = HumanPhenotypeOntology(id=hpo_id, name=term.name, definition=term.definition)
                hpo_records.append(hpo)

    logging.info("Inserting %d HumanPhenotypeOntology records", len(hpo_records))
    logging.info("Updating definitions for %d HumanPhenotypeOntology records", len(update_definitions))

    if DO_INSERT:
        if hpo_records:
            HumanPhenotypeOntology.objects.bulk_create(hpo_records)

        logging.info("Updating definitions via bulk update")
        update_hpo_definitions = []
        for (hpo_id, definition) in update_definitions.items():
            hpo = HumanPhenotypeOntology(pk=hpo_id, definition=definition)
            update_hpo_definitions.append(hpo)

        HumanPhenotypeOntology.objects.bulk_update(update_hpo_definitions,
                                                   fields=['definition'],
                                                   batch_size=BULK_INSERT_SIZE)

    end = time.time()
    logging.info("Insert HPO.owl took %.2f seconds", end - start)


def create_hpo_synonyms(ot):
    logging.info("Creating HPO Synonyms")
    start = time.time()

    scope_dict = {}
    for (k, v) in HPOSynonymScope.CHOICES:
        scope_dict[v.upper()] = k

    hpo_synonyms = []
    for term in ot.terms():
        hpo_id = get_term_hpo_id(term)
        if hpo_id:
            for synonym in term.synonyms:
                scope = scope_dict[synonym.scope]
                hpo_synonym = HPOSynonym(hpo_id=hpo_id, name=synonym.description, scope=scope)
                hpo_synonyms.append(hpo_synonym)

    logging.info("Inserting %d HumanPhenotypeOntology synonyms records", len(hpo_synonyms))
    if hpo_synonyms:
        HPOSynonym.objects.bulk_create(hpo_synonyms, batch_size=BULK_INSERT_SIZE, ignore_conflicts=True)
        hpo_synonyms = []

    logging.info("Ensuring there's a synonym for every HPO term")
    for hpo in HumanPhenotypeOntology.objects.exclude(hposynonym__name=F("name")):  # @UndefinedVariable
        hpo_synonym = HPOSynonym(hpo=hpo, name=hpo.name, scope=HPOSynonymScope.EXACT)
        hpo_synonyms.append(hpo_synonym)

    if hpo_synonyms:
        logging.info("Inserting %d Synonyms for exact HPOs", len(hpo_synonyms))
        HPOSynonym.objects.bulk_create(hpo_synonyms, batch_size=BULK_INSERT_SIZE, ignore_conflicts=True)

    end = time.time()
    logging.info("Create HPO synonyms took %.2f seconds", end - start)


def create_hpo_tree(ot):
    logging.info("Deleting existing edges")
    HPOEdge.objects.all().delete()

    logging.info("Creating Directed Acyclic Graph (this takes a while)")
    start = time.time()
    hpo_edge_list = []
    for term in ot.terms():
        hpo_id = get_term_hpo_id(term)
        if hpo_id:
            children = itertools.islice(term.subclasses(), 1, None)  # Skip 1st term (itself)
            for kid_term in children:
                kid_hpo_id = get_term_hpo_id(kid_term)
                if kid_hpo_id:
                    hpo_edge = HPOEdge(parent_id=hpo_id, child_id=kid_hpo_id)
                    hpo_edge_list.append(hpo_edge)

    HPOEdge.objects.bulk_create(hpo_edge_list, batch_size=BULK_INSERT_SIZE, ignore_conflicts=True)
    end = time.time()
    logging.info("Create HPO Graph took %.2f seconds", end - start)


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--disease_gene_txt', required=True)
        parser.add_argument('--owl', required=True)

    def handle(self, *args, **options):
        disease_gene_filename = options["disease_gene_txt"]
        owl_filename = options["owl"]
        if not owl_filename.endswith(".owl"):
            raise CommandError(f"Final arg: '{owl_filename}' must end in .owl")

        if not MIMMorbid.objects.exists():
            msg = "You have no OMIM records, please run 'mim_import' first"
            raise CommandError(msg)

        mim_ids = set(MIMMorbid.objects.all().values_list("pk", flat=True))
        df = load_diseases_to_genes_to_phenotypes(disease_gene_filename)
        insert_human_phenotype_ontology(mim_ids, df)

        ot = pronto.Ontology(owl_filename)
        insert_owl_hpo(ot)
        create_hpo_synonyms(ot)
        create_hpo_tree(ot)

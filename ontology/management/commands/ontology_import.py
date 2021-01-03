import json
import re
from collections import defaultdict
from typing import Optional, Dict, List, Iterable

from celery.worker.state import requests
from django.core.management import BaseCommand
from django.db import transaction

from genes.models import GeneSymbol, GeneSymbolAlias
from library.file_utils import file_md5sum
from ontology.models import OntologyImport, OntologySet, OntologyTerm, OntologyTermRelation, OntologyTermGeneRelation, \
    OntologyRelation

ID_EXTRACT_P = re.compile(r"^.*\/([A-Z]+)_([0-9]+)$")
HGNC_EXTRACT_P = re.compile(r"http://identifiers.org/hgnc/([0-9]+)")
OMIM_URL_P = re.compile(r"http://identifiers.org/omim/([0-9]+)")
GENE_SYMBOL_SEARCH = re.compile(r"([A-Z][A-Z,0-9]{2,})")

GENE_RELATIONS = {
    "http://purl.obolibrary.org/obo/RO_0004025": "disease causes dysfunction of",
    "http://purl.obolibrary.org/obo/RO_0004001": "has material basis in gain of function germline mutation in",
    "http://purl.obolibrary.org/obo/RO_0004021": "disease has basis in disruption of",
    "http://purl.obolibrary.org/obo/RO_0004020": "disease has basis in dysfunction of"
}

MATCH_TYPES = {
    "http://www.w3.org/2004/02/skos/core#exactMatch": OntologyRelation.EXACT,
    "http://www.w3.org/2004/02/skos/core#closeMatch": OntologyRelation.CLOSE,
    "http://www.w3.org/2004/02/skos/core#broadMatch": OntologyRelation.BROAD,
    "http://www.w3.org/2004/02/skos/core#narrowMatch": OntologyRelation.NARROW
}


class OntologyBuilder:
    """
    Not the most efficient way of getting data in (no bulk inserts are used)
    But simplest way to update existing values and create new ones
    """

    def __init__(self, filename: str, ontology_set: OntologySet, hash: str):
        self.ontology_import = OntologyImport.objects.create(ontology_set=OntologySet.MONDO, filename=filename, hash=hash)
        self.all_gene_symbols: Dict[str, GeneSymbol] = dict()
        for gene_symbol in GeneSymbol.objects.all():
            self.all_gene_symbols[gene_symbol.symbol] = gene_symbol

    def add_term(self, term_id: str, name: str = None, definition: str = None, extra: Dict = None, stub: bool = False) -> OntologyTerm:

        term_id = term_id.replace("_", ":")
        parts = term_id.split(":")
        ontology_set = parts[0]
        ontology_index = int(parts[1])

        if ontology_set != OntologySet.MONDO:
            if existing := OntologyTerm.objects.filter(id=term_id).exclude(from_import__ontology_set=OntologySet.MONDO).first():
                return existing # we've already imported from another file, don't import this stub

        if stub:
            if existing := OntologyTerm.objects.filter(id=term_id).first():
                return existing

        term, _ = OntologyTerm.objects.update_or_create(
            id=term_id,
            defaults={
                "ontology_set": ontology_set,
                "index": ontology_index,
                "name": name,
                "definition": definition,
                "extra": extra,
                "from_import": self.ontology_import
            }
        )
        return term

    def add_gene_relation(self, term_id: str, gene_symbol_str: str, relation: str) -> Optional[OntologyTermGeneRelation]:
        """
        Overrides any previous relationship - only one term -> gene relationship can exist regardless of the relation
        If the gene symbol doesn't exist, nothing will happen
        """
        if gene_symbol := self.all_gene_symbols.get(gene_symbol_str):
            key = f"{term_id}/{gene_symbol_str}"
            relation_obj, _ = OntologyTermGeneRelation.objects.update_or_create(
                term=self.add_term(term_id=term_id, stub=True),
                gene_symbol=gene_symbol,
                defaults={
                    "from_import": self.ontology_import,
                    "relation": relation
                }
            )
            return relation_obj

        print(f"Could not find gene symbol {gene_symbol_str}")
        return None

    def add_ontology_relation(self, source_term_id: str, dest_term_id: str, relation: str):
        relation, _ = OntologyTermRelation.objects.get_or_create(
            source_term=self.add_term(term_id=source_term_id, stub=True),
            dest_term=self.add_term(term_id=dest_term_id, stub=True),
            relation=relation,
            defaults={
                "from_import": self.ontology_import
            }
        )
        return relation


@transaction.atomic
def load_mondo(filename: str):
    print("This may take a few minutes")

    data_file = None
    with open(filename, 'r') as json_file:
        data_file = json.load(json_file)

    hash = file_md5sum(filename)

    ontology_builder = OntologyBuilder(filename=filename, ontology_set=OntologySet.MONDO, hash=hash)
    node_to_gene_symbol: [str, str] = dict()
    node_to_mondo: [str, str] = dict()
    count = 0

    for graph in data_file.get("graphs", []):
        print("Reviewing Nodes")
        for node in graph.get("nodes", []):
            if node.get("type") == "CLASS":

                count += 1
                if count % 1000 == 0:
                    print(f"Processing node {count}")

                if node_id_full := node.get("id"):
                    if match := ID_EXTRACT_P.match(node_id_full):
                        type = match[1]
                        index = match[2]
                        full_id = f"{type}:{index}"
                        if type != "MONDO":
                            continue

                        if meta := node.get("meta"):
                            label = node.get("lbl")
                            extra = dict()
                            synonyms = list()

                            defn = meta.get("definition", {}).get("val")
                            # if defn:
                            #    genes_mentioned = genes_mentioned.union(
                            #        set(GENE_SYMBOL_SEARCH.findall(defn)))
                            if synonyms := meta.get("synonyms"):
                                extra["synonyms"] = synonyms

                            # for synonym in meta.get("synonyms", []):
                            #    synonym_value = synonym.get("val")
                            # if synonym_value:
                            #    synonyms.append(synonym_value)
                            #    genes_mentioned = genes_mentioned.union(set(GENE_SYMBOL_SEARCH.findall(synonym_value)))

                            ontology_builder.add_term(
                                term_id=full_id,
                                name=label,
                                definition=defn,
                                extra=extra if extra else None
                            )
                            node_to_mondo[node_id_full] = full_id

                            for bp in meta.get("basicPropertyValues", []):
                                val = bp.get("val")
                                if omim_match := OMIM_URL_P.match(val):
                                    omim = f"OMIM:{omim_match[1]}"
                                    pred = bp.get("pred")
                                    pred = MATCH_TYPES.get(pred, pred)

                                    ontology_builder.add_ontology_relation(
                                        source_term_id=full_id,
                                        dest_term_id=omim,
                                        relation=pred
                                    )

                    # copy of id for gene symbol to gene symbol
                    elif match := HGNC_EXTRACT_P.match(node_id_full):
                        gene_symbol = node.get("lbl")
                        # hgnc_id = match[1]
                        node_to_gene_symbol[node_id_full] = gene_symbol

        print("Reviewing edges")
        for edge in graph.get("edges", []):
            subject_id: str = edge.get("sub")
            obj_id: str = edge.get("obj")
            relationship = edge.get("pred")

            if mondo_subject_id := node_to_mondo.get(subject_id):
                if gene_symbol := node_to_gene_symbol.get(obj_id):
                    relationship = GENE_RELATIONS.get(relationship, relationship)
                    ontology_builder.add_gene_relation(
                        term_id=mondo_subject_id,
                        gene_symbol_str=gene_symbol,
                        relation=relationship
                    )

                # TODO support other relationships
                elif mondo_obj_id := node_to_mondo.get(obj_id):
                    if relationship == "is_a":
                        ontology_builder.add_ontology_relation(
                            source_term_id=mondo_subject_id,
                            dest_term_id=mondo_obj_id,
                            relation=OntologyRelation.IS_A
                        )


def load_panelapp_au(gene_symbols: Iterable[GeneSymbol]):
    # TODO merge this with the other panel app loading code
    # can't pre-emptively do all of this - way too many gene symbols
    # so then how do we track what we've already asked for that's given us zero results
    ontology_builder = OntologyBuilder(filename="https://panelapp.agha.umccr.org/api/v1/genes/", ontology_set=OntologySet.PANEL_APP_AU, hash="-")

    PANEL_APP_OMIM = re.compile(r"([0-9]{5,})")

    for gene_symbol in gene_symbols:
        results = requests.get(f'https://panelapp.agha.umccr.org/api/v1/genes/{gene_symbol.symbol}/').json().get(
            "results")

        for panel_app_result in results:
            phenotype_row: str
            for phenotype_row in panel_app_result.get("phenotypes", []):
                # TODO look for terms other than OMIM
                for omim_match in PANEL_APP_OMIM.finditer(phenotype_row):
                    omim_id = int(omim_match[1])  # not always a valid id

                    # do we duplicate this for MONDO?
                    ontology_builder.add_gene_relation(
                        term_id=OntologySet.index_to_id(OntologySet.OMIM, omim_id),
                        gene_symbol_str=gene_symbol.symbol,
                        relation=OntologySet.PANEL_APP_AU,
                        extra={
                            "phenotype_row": phenotype_row,
                            "evidence": panel_app_result.get('evidence') or []
                        }
                    )
        print(gene_symbol.symbol)

class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--mondo_json', required=True)

    @transaction.atomic
    def handle(self, *args, **options):
        if filename := options.get("mondo_json"):
            load_mondo(filename)
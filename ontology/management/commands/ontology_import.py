import json
import re

from django.core.management import BaseCommand
from django.db import transaction

from library.file_utils import file_md5sum
from ontology.models import OntologySet, OntologyRelation
from ontology.ontology_builder import OntologyBuilder, OntologyBuilderDataUpToDateException

ID_EXTRACT_P = re.compile(r"^.*\/([A-Z]+)_([0-9]+)$")
HGNC_EXTRACT_P = re.compile(r"http://identifiers.org/hgnc/([0-9]+)")
OMIM_URL_P = re.compile(r"http://identifiers.org/omim/([0-9]+)")
GENE_SYMBOL_SEARCH = re.compile(r"([A-Z][A-Z,0-9]{2,})")
OMIM_P = re.compile("OMIM:[0-9]+")

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

@transaction.atomic
def load_mondo(filename: str, force: bool):

    data_file = None
    with open(filename, 'r') as json_file:
        data_file = json.load(json_file)

    file_hash = file_md5sum(filename)

    ontology_builder = OntologyBuilder(filename=filename, context="mondo_file", ontology_set=OntologySet.MONDO)
    if not force:
        ontology_builder.ensure_hash_changed(data_hash=file_hash)  # don't re-import if hash hasn't changed
    else:
        ontology_builder.data_hash=file_hash

    print("This may take a few minutes")
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
                            if defn:
                                defn = defn.replace(".nn", ".\n")

                            # if defn:
                            #    genes_mentioned = genes_mentioned.union(
                            #        set(GENE_SYMBOL_SEARCH.findall(defn)))
                            if synonyms := meta.get("synonyms"):
                                extra["synonyms"] = synonyms
                                for synonym in synonyms:
                                    pred = synonym.get("pred")
                                    if pred == "hasExactSynonym":
                                        val = synonym.get("val")
                                        for xref in synonym.get("xrefs", []):
                                            if OMIM_P.match(xref):
                                                ontology_builder.add_term(
                                                    term_id=xref,
                                                    name=label,
                                                    definition=f"Name copied from synonym {full_id}",
                                                    stub=True
                                                )
                                                ontology_builder.add_ontology_relation(
                                                    source_term_id=full_id,
                                                    dest_term_id=xref,
                                                    relation=OntologyRelation.EXACT
                                                )

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
        count = 0
        for edge in graph.get("edges", []):
            count += 1
            if count % 1000 == 0:
                print(f"Processing edge {count}")

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
    print("Purging old data")
    ontology_builder.complete()

    print(f"Records inserted: {ontology_builder.insert_count}")
    print(f"Records updated: {ontology_builder.update_count}")
    print(f"Records deleted: {ontology_builder.delete_count}")

    print("Committing...")


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--force', action="store_true")
        parser.add_argument('--mondo_json', required=True)

    @transaction.atomic
    def handle(self, *args, **options):
        if filename := options.get("mondo_json"):
            try:
                load_mondo(filename, force=options.get("force"))
            except OntologyBuilderDataUpToDateException:
                print("File hash is the same as last import")
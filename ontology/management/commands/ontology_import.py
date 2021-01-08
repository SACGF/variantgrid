import enum
import json
import re
from datetime import timedelta, datetime, timezone
from typing import Optional, Dict, List, Iterable
from django.core.management import BaseCommand
from django.db import transaction
from lazy import lazy
from model_utils.models import now
from django.utils import timezone
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

class OntologyBuilderDataUpToDateException(Exception):
    pass


class OntologyBuilder:
    """
    Not the most efficient way of getting data in (no bulk inserts are used)
    But simplest way to update existing values and create new ones
    """
    class BuildState(enum.IntEnum):
        """
        Bit overkill but returning a boolean from a few methods results in unclear
        """
        NO_UPDATE_REQUIRED = 0
        UPDATE_REQUIRED = 1

    def __init__(self, filename: str, context: str, ontology_set: OntologySet):
        self.filename = filename.split("/")[-1]
        self.ontology_set = ontology_set
        self.context = context
        self.data_hash = None
        self.previous_import: Optional[OntologyImport] = OntologyImport.objects.filter(ontology_set=ontology_set, context=context, completed=True).order_by('-modified').first()
        self.insert_count = 0
        self.update_count = 0
        self.delete_count = 0

    def ensure_old(self, max_age: timedelta):
        """
        Raises OntologyBuilderDataUpToDateException if data is up to date
        """
        if self.previous_import:
            if self.previous_import.processed_date > datetime.now(tz=timezone.get_default_timezone()) - max_age:
                raise OntologyBuilderDataUpToDateException()

    def ensure_hash_changed(self, data_hash: str):
        self.data_hash = data_hash
        if self.previous_import and self.previous_import.hash == data_hash:
            self.previous_import.processed_date = now
            self.previous_import.save()
            raise OntologyBuilderDataUpToDateException()

    def complete(self, purge_old=True):
        """
        Returns how many old relationships were deleted from a previous run
        """
        if purge_old:
            for model in [OntologyTermGeneRelation, OntologyTermRelation, OntologyTerm]:
                count, _ = model.objects.filter(from_import__context=self.context, from_import__ontology_set=self.ontology_set).exclude(from_import=self._ontology_import).delete()
                self.delete_count += count

        self._ontology_import.completed = True
        self._ontology_import.save()

    @lazy
    def _ontology_import(self) -> OntologyImport:
        return OntologyImport.objects.create(ontology_set=self.ontology_set, context=self.context, filename=self.filename, processed_date=now, hash=self.data_hash)

    def _count(self, created: bool):
        if created:
            self.insert_count += 1
        else:
            self.update_count += 1

    def add_term(self, term_id: str, name: str = None, definition: str = None, extra: Dict = None, stub: bool = False) -> OntologyTerm:

        term_id = term_id.replace("_", ":")
        parts = term_id.split(":")
        ontology_set = parts[0]
        ontology_index = int(parts[1])

        if stub:
            if existing := OntologyTerm.objects.filter(id=term_id).first():
                existing_import = existing.from_import
                # don't override an existing record with a stub
                if existing_import.context != self.context or existing_import.ontology_set != self.ontology_set:
                    return existing
                # but if the old data is a stub from the same kind of import, update the from_import
                # so we don't purge it

        term, created = OntologyTerm.objects.update_or_create(
            id=term_id,
            defaults={
                "ontology_set": ontology_set,
                "index": ontology_index,
                "name": name,
                "definition": definition,
                "extra": extra,
                "from_import": self._ontology_import
            }
        )
        self._count(created)
        return term

    def add_gene_relation(self, term_id: str, gene_symbol_str: str, relation: str, extra = None) -> Optional[OntologyTermGeneRelation]:
        """
        Overrides any previous relationship - only one term -> gene relationship can exist regardless of the relation
        If the gene symbol doesn't exist, nothing will happen
        """
        if gene_symbol := GeneSymbol.objects.filter(symbol=gene_symbol_str).first():
            key = f"{term_id}/{gene_symbol_str}"
            relation_obj, created = OntologyTermGeneRelation.objects.update_or_create(
                term=self.add_term(term_id=term_id, stub=True),
                gene_symbol=gene_symbol,
                relation=relation,
                defaults={
                    "extra": extra,
                    "from_import": self._ontology_import
                }
            )
            self._count(created)
            return relation_obj

        print(f"Could not find gene symbol {gene_symbol_str}")
        return None

    def add_ontology_relation(self, source_term_id: str, dest_term_id: str, relation: str):
        relation, created = OntologyTermRelation.objects.get_or_create(
            source_term=self.add_term(term_id=source_term_id, stub=True),
            dest_term=self.add_term(term_id=dest_term_id, stub=True),
            relation=relation,
            defaults={
                "from_import": self._ontology_import
            }
        )
        self._count(created)
        return relation


@transaction.atomic
def load_mondo(filename: str, force: bool):

    data_file = None
    with open(filename, 'r') as json_file:
        data_file = json.load(json_file)

    file_hash = file_md5sum(filename)

    ontology_builder = OntologyBuilder(filename=filename, context="mondo_file", ontology_set=OntologySet.MONDO)
    if not force:
        ontology_builder.ensure_hash_changed(data_hash=file_hash)  # don't re-import if hash hasn't changed

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
    print("Purging old data")
    ontology_builder.complete()

    print(f"Records inserted: {ontology_builder.insert_count}")
    print(f"Records updated: {ontology_builder.update_count}")
    print(f"Records deleted: {ontology_builder.delete_count}")


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
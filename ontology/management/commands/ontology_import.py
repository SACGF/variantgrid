import itertools
import json
import re
from dataclasses import dataclass

import pandas as pd
import pronto
from django.core.management import BaseCommand
from django.db import transaction

from annotation.models import HPOSynonymScope
from library.file_utils import file_md5sum
from ontology.models import OntologyService, OntologyRelation, OntologyTerm, OntologyImportSource
from ontology.ontology_builder import OntologyBuilder, OntologyBuilderDataUpToDateException


"""
MONDO import file can be found http://www.obofoundry.org/ontology/mondo.html
Importing it will provide MONDO and OMIM terms
"""

GENE_SYMBOL_SEARCH = re.compile(r"([A-Z][A-Z,0-9]{2,})")

GENE_RELATIONS = {
    "http://purl.obolibrary.org/obo/RO_0004025": "disease causes dysfunction of",
    "http://purl.obolibrary.org/obo/RO_0004001": "has material basis in gain of function germline mutation in",
    "http://purl.obolibrary.org/obo/RO_0004021": "disease has basis in disruption of",
    "http://purl.obolibrary.org/obo/RO_0004020": "disease has basis in dysfunction of",
    "http://purl.obolibrary.org/obo/RO_0004026": "disease has location",
    "http://purl.obolibrary.org/obo/RO_0004027": "disease has inflammation site",
    "http://purl.obolibrary.org/obo/RO_0004030": "disease arises from structure"
}

MATCH_TYPES = {
    "http://www.w3.org/2004/02/skos/core#exactMatch": OntologyRelation.EXACT,
    "http://www.w3.org/2004/02/skos/core#closeMatch": OntologyRelation.CLOSE,
    "http://www.w3.org/2004/02/skos/core#broadMatch": OntologyRelation.BROAD,
    "http://www.w3.org/2004/02/skos/core#narrowMatch": OntologyRelation.NARROW,
    "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId": OntologyRelation.ALTERNATIVE,
    "http://purl.obolibrary.org/obo/IAO_0100001": OntologyRelation.REPLACED
}

"""
TODO - this runs pretty slow (due to many redundant update_or_create) calls.
Rework it so we can keep a cache of everything already updated or created this run
"""

ID_OBO = re.compile("^http://purl[.]obolibrary[.]org/obo/([A-Z]+)_([0-9]+)$")
ID_IDENTIFIERS = re.compile("http://identifiers[.]org/([A-Za-z]+)/([0-9]+)$")
ID_STRAIGHT = re.compile("^([A-Z]+):([0-9]+)$")

@dataclass
class TermId:
    type: str = None
    index: str = None

    def __init__(self, qualified_ref: str):
        for pattern in [ID_OBO, ID_IDENTIFIERS, ID_STRAIGHT]:
            if match := pattern.match(qualified_ref):
                self.type = match[1].upper()
                self.index = match[2]
                return

    @property
    def id(self) -> str:
        return f"{self.type}:{self.index}"


@transaction.atomic
def load_mondo(filename: str, force: bool):
    data_file = None

    file_hash = file_md5sum(filename)

    ontology_builder = OntologyBuilder(
        filename=filename,
        context="mondo_file",
        import_source=OntologyService.MONDO,
        force_update=force,
        processor_version=3)

    ontology_builder.ensure_hash_changed(data_hash=file_hash)  # don't re-import if hash hasn't changed

    with open(filename, 'r') as json_file:
        data_file = json.load(json_file)

    print("This may take a few minutes")
    node_to_hgnc_id: [str, str] = dict()
    node_to_mondo: [str, str] = dict()

    for graph in data_file.get("graphs", []):

        print("** Reviewing Nodes")
        for node in graph.get("nodes", []):

            if node_id_full := node.get("id"):
                term = TermId(node_id_full)
                if term.type == "MONDO":
                    full_id = term.id
                    node_to_mondo[node_id_full] = full_id

                    if meta := node.get("meta"):
                        label = node.get("lbl")
                        extra = dict()

                        defn = meta.get("definition", {}).get("val")
                        if defn:
                            defn = defn.replace(".nn", ".\n")

                        # if defn:
                        #    genes_mentioned = genes_mentioned.union(
                        #        set(GENE_SYMBOL_SEARCH.findall(defn)))
                        if synonyms := meta.get("synonyms"):
                            extra["synonyms"] = synonyms

                        # make the term early so we don't have to create stubs for it if we find relationships
                        ontology_builder.add_term(
                            term_id=full_id,
                            name=label,
                            definition=defn,
                            extra=extra if extra else None,
                            primary_source=True
                        )

                        synonym_set = set()
                        if synonyms := meta.get("synonyms"):

                            # only allow 1 relationship between any 2 terms (though DB does allow more)
                            # storing all of them would be more "accurate" but gets in the way of our usage
                            # prioritise relationships as EXACT, RELATED, related terms, XREF
                            for synonym in synonyms:
                                pred = synonym.get("pred")
                                if pred == "hasExactSynonym":
                                    # val = synonym.get("val")
                                    for xref in synonym.get("xrefs", []):
                                        xref_term = TermId(xref)
                                        if xref_term.type in {"HP", "OMIM"}:
                                            ontology_builder.add_term(
                                                term_id=xref,
                                                name=label,
                                                definition=f"Name copied from synonym {full_id}",
                                                primary_source=False
                                            )
                                            ontology_builder.add_ontology_relation(
                                                source_term_id=full_id,
                                                dest_term_id=xref,
                                                relation=OntologyRelation.EXACT
                                            )
                                            synonym_set.add(xref)

                            # look at related synonyms second, if we have RELATED, don't bother with any other relationships
                            for synonym in synonyms:
                                pred = synonym.get("pred")
                                if pred == "hasRelatedSynonym":
                                    # val = synonym.get("val")
                                    for xref in synonym.get("xrefs", []):
                                        xref_term = TermId(xref)
                                        if xref_term.type in {"HP", "OMIM"} and not xref_term.id in synonym_set:
                                            ontology_builder.add_term(
                                                term_id=xref_term.id,
                                                name=label,
                                                definition=f"Name copied from related synonym {full_id}",
                                                primary_source=False
                                            )
                                            ontology_builder.add_ontology_relation(
                                                source_term_id=full_id,
                                                dest_term_id=xref_term.id,
                                                relation=OntologyRelation.RELATED
                                            )
                                            synonym_set.add(xref_term.id)
                        #end synonymns

                        for bp in meta.get("basicPropertyValues", []):
                            val = TermId(bp.get("val"))
                            if val.type in {"HP", "OMIM"} and not val.id in synonym_set:
                                pred = bp.get("pred")
                                pred = MATCH_TYPES.get(pred, pred)

                                ontology_builder.add_term(
                                    term_id=val.id,
                                    name=label,
                                    definition=f"Name copied from {pred} synonym {full_id}",
                                    primary_source=False
                                )
                                ontology_builder.add_ontology_relation(
                                    source_term_id=full_id,
                                    dest_term_id=val.id,
                                    relation=pred
                                )
                            synonym_set.add(val.id)

                        if xrefs := meta.get("xrefs"):
                            for xref in xrefs:
                                val = TermId(xref.get("val"))
                                if val.type in {"HP", "OMIM"} and val.id not in synonym_set:
                                    ontology_builder.add_term(
                                        term_id=val.id,
                                        name=label,
                                        definition=f"Name copied from xref synonym {full_id}",
                                        primary_source=False
                                    )
                                    ontology_builder.add_ontology_relation(
                                        source_term_id=full_id,
                                        dest_term_id=val.id,
                                        relation=OntologyRelation.XREF
                                    )

                # copy of id for gene symbol to gene symbol
                elif term.type == "HGNC":
                    gene_symbol = node.get("lbl")
                    ontology_builder.add_term(
                        term_id=term.id,
                        name=gene_symbol,
                        definition=None,
                        primary_source=False
                    )
                    node_to_hgnc_id[node_id_full] = term.id

        print("** Reviewing axioms")
        if axioms := graph.get("logicalDefinitionAxioms"):
            for axiom in axioms:
                defined_class_id = TermId(axiom.get("definedClassId"))
                if defined_class_id.type != "MONDO":
                    continue

                genus_ids = [TermId(genus) for genus in axiom.get("genusIds")]
                mondo_genus = [term for term in genus_ids if term.type == "MONDO"]

                for restriction in axiom.get("restrictions"):

                    filler = TermId(restriction.get("fillerId"))
                    if filler.type != "HGNC":
                        continue

                    relation = GENE_RELATIONS.get(restriction.get("propertyId"))
                    if relation is None:
                        print("Unexpected relationship " + restriction.get("propertyId"))
                        continue

                    ontology_builder.add_ontology_relation(
                        source_term_id=defined_class_id.id,
                        dest_term_id=filler.id,
                        extra={"via": ", ".join([m.id for m in mondo_genus])},
                        relation=relation
                    )
                    for term in mondo_genus:
                        ontology_builder.add_ontology_relation(
                            source_term_id=term.id,
                            dest_term_id=filler.id,
                            relation=relation
                        )

        print("** Reviewing edges")

        for edge in graph.get("edges", []):

            subject_id: str = edge.get("sub")
            obj_id: str = edge.get("obj")
            relationship = edge.get("pred")

            if mondo_subject_id := node_to_mondo.get(subject_id):
                if hgnc_id := node_to_hgnc_id.get(obj_id):
                    relationship = GENE_RELATIONS.get(relationship, relationship)
                    ontology_builder.add_ontology_relation(
                        source_term_id=mondo_subject_id,
                        dest_term_id=hgnc_id,
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
    ontology_builder.report()
    print("Committing...")


@transaction.atomic
def load_hpo(filename: str, force: bool):
    ontology_builder = OntologyBuilder(
        filename=filename,
        context="hpo_file",
        import_source=OntologyImportSource.HPO,
        force_update=force)

    file_hash = file_md5sum(filename)
    ontology_builder.ensure_hash_changed(data_hash=file_hash)  # don't re-import if hash hasn't changed
    print("About to pronto the file")
    ot = pronto.Ontology(filename)
    print("Pronto complete")
    scope_lookup = {v.upper(): k for k, v in HPOSynonymScope.choices}

    for term in ot.terms():
        term_id = str(term.id)

        extra = None
        synonyms = []
        for synonym in term.synonyms:
            scope = scope_lookup[synonym.scope]
            synonyms.append({
                "name": synonym.description,
                "scope": scope
            })
        if synonyms:
            extra = {"synonyms": synonyms}

        if not term_id.startswith(OntologyService.HPO):
            continue

        ontology_builder.add_term(
            term_id=term_id,
            name=term.name,
            extra=extra,
            definition=term.definition,
            primary_source=True
        )

        children = itertools.islice(term.subclasses(), 1, None)
        for kid_term in children:

            if not kid_term.id.startswith(OntologyService.HPO):
                continue

            ontology_builder.add_ontology_relation(
                dest_term_id=term.id,
                source_term_id=kid_term.id,
                relation=OntologyRelation.IS_A
            )
    ontology_builder.complete()
    ontology_builder.report()
    print("Committing...")


@transaction.atomic
def load_hpo_disease(filename: str, force: bool):
    ontology_builder = OntologyBuilder(
        filename=filename,
        context="hpo_disease",
        import_source=OntologyImportSource.HPO,
        processor_version=4,
        force_update=force)
    file_hash = file_md5sum(filename)
    ontology_builder.ensure_hash_changed(data_hash=file_hash)  # don't re-import if hash hasn't changed
    df = pd.read_csv(filename, index_col=None, comment='#', sep='\t',
                     names=['disease_id', 'gene_symbol', 'gene_id', 'hpo_id', 'hpo_name'],
                     dtype={"gene_id": int})

    by_hpo = df.groupby(["hpo_id"])
    for hpo_id, by_hpo_data in by_hpo:
        # make HPO stubs in case HPO import hasn't happened yet
        hpo_name = by_hpo_data["hpo_name"].iloc[0]
        ontology_builder.add_term(
            term_id=hpo_id,
            name=hpo_name,
            definition=None,
            primary_source=False
        )

        by_omim = by_hpo_data.groupby(["disease_id"])
        for omim_id, by_omim_data in by_omim:
            # link HPO -> OMIM
            ontology_builder.add_ontology_relation(
                source_term_id=hpo_id,
                dest_term_id=omim_id,
                relation=OntologyRelation.ALL_FREQUENCY
            )

    by_gene = df.groupby(["gene_symbol"])
    for gene_symbol, gene_data in by_gene:
        # IMPORTANT, the gene_id is the entrez gene_id, not the HGNC gene id
        try:
            hgnc_term = OntologyTerm.get_gene_symbol(gene_symbol)

            by_omim = gene_data.groupby(["disease_id"])
            for omim_id, by_omim_data in by_omim:
                ontology_builder.add_ontology_relation(
                    source_term_id=omim_id,
                    dest_term_id=hgnc_term.id,
                    relation=OntologyRelation.ENTREZ_ASSOCIATION
                )
        except ValueError:
            print(f"Could not resolve gene symbol {gene_symbol} to HGNC ID")

    ontology_builder.complete()
    ontology_builder.report()


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--force', action="store_true")
        parser.add_argument('--mondo_json', required=False)
        parser.add_argument('--hpo_owl', required=False)
        parser.add_argument('--omim_frequencies', required=False)

    @transaction.atomic
    def handle(self, *args, **options):
        force = options.get("force")
        if filename := options.get("mondo_json"):
            try:
                load_mondo(filename, force)
            except OntologyBuilderDataUpToDateException:
                print("MONDO File hash is the same as last import")
        if filename := options.get("hpo_owl"):
            try:
                load_hpo(filename, force)
            except OntologyBuilderDataUpToDateException:
                print("HPO File hash is the same as last import")
        if filename := options.get("omim_frequencies"):
            try:
                load_hpo_disease(filename, force)
            except OntologyBuilderDataUpToDateException:
                print("HPO Disease hash is the same as last import")

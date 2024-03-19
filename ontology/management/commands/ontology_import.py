import csv
import itertools
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

import pandas as pd
import pronto
from django.core.management import BaseCommand
from model_utils.models import now

from annotation.models.models_enums import HPOSynonymScope
from genes.models import HGNC
from library.utils.file_utils import file_md5sum
from ontology.gencc import load_gencc
from ontology.models import OntologyService, OntologyRelation, OntologyTerm, OntologyImportSource, OntologyImport, \
    OntologyTermRelation, OntologyVersion, OntologyTermStatus
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
    "http://purl.obolibrary.org/obo/RO_0004030": "disease arises from structure",
    "http://purl.obolibrary.org/obo/RO_0004003": "has material basis in germline mutation in",
    "http://purl.obolibrary.org/obo/RO_0004004": "has material basis in somatic mutation in",
    "http://purl.obolibrary.org/obo/RO_0004028": "realized in response to stimulus",

}

MATCH_TYPES = {
    "http://www.w3.org/2004/02/skos/core#exactMatch": OntologyRelation.EXACT,
    "http://www.w3.org/2004/02/skos/core#closeMatch": OntologyRelation.CLOSE,
    "http://www.w3.org/2004/02/skos/core#broadMatch": OntologyRelation.BROAD,
    "http://www.w3.org/2004/02/skos/core#narrowMatch": OntologyRelation.NARROW,
    "http://www.geneontology.org/formats/oboInOwl#hasAlternativeId": OntologyRelation.ALTERNATIVE,
    "http://www.geneontology.org/formats/oboInOwl#consider": OntologyRelation.CONSIDER,
    "http://purl.obolibrary.org/obo/IAO_0100001": OntologyRelation.REPLACED
}

"""
TODO - this runs pretty slow (due to many redundant update_or_create) calls.
Rework it so we can keep a cache of everything already updated or created this run
"""

ID_OBO = re.compile("^http://purl[.]obolibrary[.]org/obo/([A-Z]+)_([0-9]+)$")
ID_IDENTIFIERS = re.compile("http://identifiers[.]org/([A-Za-z]+)/([0-9]+)$")
ID_STRAIGHT = re.compile("^([A-Z]+):([0-9]+)$")
OMIM_URL = re.compile("^https://(omim)[.]org/entry/([0-9]+)$")


@dataclass
class TermId:
    type: str = None
    index: str = None

    def __init__(self, qualified_ref: str):
        for pattern in [ID_OBO, ID_IDENTIFIERS, ID_STRAIGHT, OMIM_URL]:
            if match := pattern.match(qualified_ref):
                self.type = match[1].upper()
                self.index = match[2]
                return

    @property
    def id(self) -> str:
        return f"{self.type}:{self.index}"


def load_mondo(filename: str, force: bool):
    file_hash = file_md5sum(filename)

    ontology_builder = OntologyBuilder(
        filename=filename,
        context="mondo_file",
        import_source=OntologyService.MONDO,
        force_update=force,
        processor_version=17)

    ontology_builder.ensure_hash_changed(data_hash=file_hash)  # don't re-import if hash hasn't changed
    ontology_builder.cache_everything()

    data_file: dict
    with open(filename, 'r') as json_file:
        data_file = json.load(json_file)

    node_to_hgnc_id: [str, str] = {}
    node_to_mondo: [str, str] = {}

    for graph in data_file.get("graphs", []):

        print("** Reviewing Nodes")
        for node in graph.get("nodes", []):

            if node_id_full := node.get("id"):
                term = TermId(node_id_full)
                if term.type == "MONDO":
                    full_id = term.id
                    node_to_mondo[node_id_full] = full_id

                    if meta := node.get("meta"):
                        label = node.get("lbl") or ""
                        deprecated = meta.get("deprecated")
                        extra = {}

                        defn = meta.get("definition", {}).get("val")
                        if defn:
                            defn = defn.replace(".nn", ".\n")

                        # if defn:
                        #    genes_mentioned = genes_mentioned.union(
                        #        set(GENE_SYMBOL_SEARCH.findall(defn)))
                        if synonyms := meta.get("synonyms"):
                            extra["synonyms"] = synonyms

                        aliases = []
                        term_relation_types: dict[str, list[str]] = defaultdict(list)

                        for bp in meta.get("basicPropertyValues", []):
                            val = TermId(bp.get("val"))
                            if val.type in {"HP", "OMIM", "MONDO"}:
                                pred = bp.get("pred")
                                if pred := MATCH_TYPES.get(pred):
                                    term_relation_types[val.id].append(pred)

                        if xrefs := meta.get("xrefs"):
                            for xref in xrefs:
                                val = TermId(xref.get("val"))
                                if val.type in {"HP", "OMIM"}:
                                    term_relation_types[val.id].append("xref")

                        if synonyms := meta.get("synonyms"):
                            # only allow 1 relationship between any 2 terms (though DB does allow more)
                            # storing all of them would be more "accurate" but gets in the way of our usage
                            # prioritise relationships as EXACT, RELATED, related terms, XREF

                            for pred_type in ["hasExactSynonym", "hasRelatedSynonym"]:
                                relation = {
                                    "hasExactSynonym": OntologyRelation.EXACT_SYNONYM,
                                    "hasRelatedSynonym": OntologyRelation.RELATED_SYNONYM
                                }[pred_type]
                                for synonym in synonyms:
                                    pred = synonym.get("pred")
                                    if pred == pred_type:
                                        # could potentially get name here from exact synonym
                                        # but there can be multiple
                                        aliases.append(synonym.get("val"))
                                        for xref in synonym.get("xrefs", []):
                                            xref_term = TermId(xref)
                                            if xref_term.type in {"HP", "OMIM"}:
                                                term_relation_types[xref_term.id].append(relation)

                        for key, relations in term_relation_types.items():
                            unique_relations = sorted(set(relations))

                            ontology_builder.add_term(
                                term_id=key,
                                name="",
                                definition="",
                                primary_source=False
                            )
                            ontology_builder.add_ontology_relation(
                                source_term_id=full_id,
                                dest_term_id=key,
                                relation=relations[0],
                                extra={"all_relations": unique_relations}
                            )

                        # end synonyms
                        # occasionally split up MONDO name into different aliases
                        if label and "MONDO" in full_id:
                            if ";" in label:
                                parts = [part.strip() for part in label.split(";")]
                                if len(parts) == 2:
                                    if parts[1].isupper():
                                        for part in parts:
                                            if part not in aliases:
                                                aliases.append(part)

                        ontology_builder.add_term(
                            term_id=full_id,
                            name=label,
                            definition=defn,
                            extra=meta,
                            aliases=aliases,
                            primary_source=True,
                            status=OntologyTermStatus.DEPRECATED if deprecated else None
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

                for mondo_genu in mondo_genus:
                    ontology_builder.add_ontology_relation(
                        source_term_id=defined_class_id.id,
                        dest_term_id=mondo_genu.id,
                        relation=OntologyRelation.IS_A,
                        extra={"subtype": "genus"}
                    )

                if restrictions := axiom.get("restrictions"):
                    for restriction in restrictions:
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

    ontology_builder.complete(verbose=True)


def load_hpo(filename: str, force: bool):
    ontology_builder = OntologyBuilder(
        filename=filename,
        context="hpo_file",
        import_source=OntologyImportSource.HPO,
        processor_version=2,
        force_update=force)

    file_hash = file_md5sum(filename)
    ontology_builder.ensure_hash_changed(data_hash=file_hash)  # don't re-import if hash hasn't changed
    ontology_builder.cache_everything()
    print("About to pronto the file")
    ot = pronto.Ontology(filename)
    print("Pronto complete")
    scope_lookup = {v.upper(): k for k, v in HPOSynonymScope.choices}

    for term in ot.terms():
        term_id = str(term.id)

        extra = None
        detailed_aliases = []
        aliases = set()
        for synonym in term.synonyms:
            aliases.add(synonym.description)
            scope = scope_lookup[synonym.scope]
            detailed_aliases.append({
                "name": synonym.description,
                "scope": scope
            })
        if detailed_aliases:
            detailed_aliases.sort(key=lambda x: x.get("name"))
            extra = {"synonyms": detailed_aliases}

        if not term_id.startswith(OntologyService.HPO):
            continue

        status = None
        if not term.name:
            # HP file seems to just wipe names if a term has been moved
            status = OntologyTermStatus.DEPRECATED

        ontology_builder.add_term(
            term_id=term_id,
            name=term.name,
            extra=extra,
            definition=term.definition,
            primary_source=True,
            status=status,
            aliases=list(sorted(aliases))
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
    ontology_builder.complete(verbose=True)
    print("Committing...")


def load_phenotype_to_genes(filename: str, force: bool):
    ontology_builder = OntologyBuilder(
        filename=filename,
        context="phenotype_to_genes",
        import_source=OntologyImportSource.HPO,
        processor_version=1,
        force_update=force)
    file_hash = file_md5sum(filename)
    ontology_builder.ensure_hash_changed(data_hash=file_hash)  # don't re-import if hash hasn't changed
    ontology_builder.cache_everything()
    df = pd.read_csv(filename, index_col=None, comment='#', sep='\t',
                     names=['hpo_id', 'hpo_name', 'entrez_gene_id', 'entrez_gene_symbol', 'status', 'source', 'omim_id'])

    df = df[df.source.isin(["mim2gene"])]

    # first kill of old data
    old_imports = OntologyImport.objects.filter(filename="OMIM_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt")
    if old_imports.exists():
        print("Removing any data from OMIM_ALL_FREQUENCIES")
        deleted = OntologyTermRelation.objects.filter(from_import__in=old_imports).delete()
        print(f"{deleted} relationships removed")

    by_hpo = df.groupby("hpo_id")
    for hpo_id, by_hpo_data in by_hpo:
        # make HPO stubs in case HPO import hasn't happened yet
        hpo_name = by_hpo_data["hpo_name"].iloc[0]
        ontology_builder.add_term(
            term_id=hpo_id,
            name=hpo_name,
            definition=None,
            primary_source=False
        )

        by_omim = by_hpo_data.groupby("omim_id")
        for omim_id, by_omim_data in by_omim:
            # link HPO -> OMIM
            ontology_builder.add_ontology_relation(
                source_term_id=hpo_id,
                dest_term_id=omim_id,
                relation=OntologyRelation.ASSOCIATED
            )

    by_gene = df.groupby("entrez_gene_symbol")
    for gene_symbol, gene_data in by_gene:
        # IMPORTANT, the gene_id is the entrez gene_id, not the HGNC gene id
        try:
            hgnc_term = OntologyTerm.get_gene_symbol(gene_symbol)

            by_omim = gene_data.groupby("omim_id")
            for omim_id, by_omim_data in by_omim:
                ontology_builder.add_ontology_relation(
                    source_term_id=omim_id,
                    dest_term_id=hgnc_term.id,
                    relation=OntologyRelation.MIM_2_GENE
                )
        except ValueError:
            print(f"Could not resolve gene symbol {gene_symbol} to HGNC ID")

    ontology_builder.complete(verbose=True)


def load_biomart(filename: str, force: bool):
    MIM_DESCRIPTION = "MIM morbid description"
    MIM_ACCESSION = "MIM morbid accession"

    ontology_builder = OntologyBuilder(
        filename=filename,
        context="biomart_omim_aliases",
        import_source="biomart",
        processor_version=3,
        force_update=force)
    file_hash = file_md5sum(filename)
    ontology_builder.ensure_hash_changed(data_hash=file_hash)
    ontology_builder.cache_everything()
    # Create MIMMorbid from BioMart file
    mim_biomart_df = pd.read_csv(filename, sep='\t').dropna().astype({"MIM morbid accession": int})
    for expected_col in [MIM_DESCRIPTION, MIM_ACCESSION]:
        if expected_col not in mim_biomart_df.columns:
            msg = f"file {filename} missing column: '{expected_col}': columns: '{mim_biomart_df.columns}'"
            raise ValueError(msg)

    mim_biomart_df = mim_biomart_df.set_index(MIM_ACCESSION)
    description_series = mim_biomart_df[MIM_DESCRIPTION]

    for mim_accession_id, description in description_series.items():
        descriptions_list = [x for x in str(description).split(";;")]
        name = descriptions_list[0]
        # aliases = [name] + [term for term in [term.strip() for term in str(description).split(";")] if term]
        aliases = [term for term in [term.strip() for term in str(description).split(";")] if term and term != name]
        ontology_builder.add_term(
            term_id=f"OMIM:{mim_accession_id}",
            name=name,
            definition=None,
            primary_source=False,  # primary source is now the OMIM file if it's available
            aliases=aliases
        )

    ontology_builder.complete(verbose=True)


def load_omim(filename: str, force: bool):

    ontology_builder = OntologyBuilder(
        filename=filename,
        context="omim_file",
        import_source=OntologyImportSource.OMIM,
        processor_version=7,
        force_update=force)

    file_hash = file_md5sum(filename)
    ontology_builder.ensure_hash_changed(data_hash=file_hash)  # don't re-import if hash hasn't changed
    ontology_builder.cache_everything()

    with open(filename, "r") as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        next(csv_reader)  # title row
        next(csv_reader)  # date row (worth reading e.g. "Generated: 20201-02-04")
        header = next(csv_reader)
        OMIM_EXPECTED_HEADER = ["# Prefix", "MIM Number", "Preferred Title; symbol", "Alternative Title(s); symbol(s)", "Included Title(s); symbols"]
        MOVED_TO = re.compile("MOVED TO ([0-9]+)")
        if header != OMIM_EXPECTED_HEADER:
            raise ValueError(f"Header not as expected, got {header}")

        PREFIXES = {
            "Asterisk": "Gene",
            "Plus": "Gene and phenotype, combined",
            "Number Sign": "Phenotype, molecular basis known",
            "Percent": "Phenotype or locus, molecular basis unknown",
            "NULL": "Other, mainly phenotypes with suspected mendelian basis",
            "Caret": "Entry has been removed from the database or moved to another entry"
        }
        VALID_PREFIXES = {
            "Number Sign",
            "Percent",
            "NULL"
        }

        for row in csv_reader:
            prefix = row[0]
            omim_type = PREFIXES.get(prefix)
            if len(row) <= 1 or not omim_type:
                continue
            mim_number = row[1]
            preferred_title = row[2]
            alternative_terms = row[3]
            included_titles = row[4]

            moved_to: Optional[str] = None
            aliases = []

            if match := MOVED_TO.match(preferred_title):
                # This will only happen if you uncomment Caret
                moved_to = match.group(1)
                # preferred_title = f"obsolete, see OMIM:{moved_to}"
            else:
                # aliases.append(preferred_title)  # other tables don't have the name copied into aliases
                # so don't do it here
                aliases += [term for term in [term.strip() for term in (preferred_title + ";" + alternative_terms).split(";")] if term and term != preferred_title]

            extras = {"type": omim_type}
            if included_titles:
                extras["included_titles"] = included_titles

            ontology_builder.add_term(
                term_id=f"OMIM:{mim_number}",
                name=preferred_title,
                definition=None,
                aliases=aliases,
                extra=extras,
                primary_source=True,
                status=None if prefix in VALID_PREFIXES else OntologyTermStatus.NON_CONDITION
            )
            if moved_to:
                ontology_builder.add_ontology_relation(
                    source_term_id=f"OMIM:{moved_to}",
                    dest_term_id=f"OMIM:{mim_number}",
                    relation=OntologyRelation.REPLACED
                )
    ontology_builder.complete(purge_old_terms=True)


def sync_hgnc():
    uploads: list[OntologyTerm] = []

    o_import = OntologyImport.objects.create(
        import_source="HGNC Sync",
        filename="HGNC Aliases",
        context="hgnc_sync",
        hash="N/A",
        processor_version=1,
        processed_date=now,
        completed=True)

    for hgnc in HGNC.objects.all():
        uploads.append(OntologyTerm(
            id=f"HGNC:{hgnc.id}",
            ontology_service=OntologyService.HGNC,
            name=hgnc.gene_symbol_id,
            index=hgnc.id,
            definition=hgnc.approved_name,
            from_import=o_import
        ))

    old_count = OntologyTerm.objects.filter(ontology_service=OntologyService.HGNC).count()
    OntologyTerm.objects.bulk_create(uploads, ignore_conflicts=True)
    updated_count = OntologyTerm.objects.filter(ontology_service=OntologyService.HGNC).count()
    delta = updated_count - old_count
    print(f"Inserted {delta} records from HGNCGeneNames into OntologyTerm")


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--force', action="store_true")
        parser.add_argument('--mondo', required=False, help="mondo json file")
        parser.add_argument('--hpo', required=False, help="hpo owl file")
        parser.add_argument('--omim_frequencies', required=False)  # note this is deprecated
        parser.add_argument('--phenotype_to_genes', required=False)
        parser.add_argument('--hgnc_sync', action="store_true", required=False)
        parser.add_argument('--biomart', required=False)
        parser.add_argument('--omim', required=False)
        parser.add_argument('--gencc', required=False)

    def handle(self, *args, **options):
        force = options.get("force")

        if options.get("hgnc_sync"):
            print("Syncing HGNC")
            sync_hgnc()

        if filename := options.get("gencc"):
            try:
                file_hash = file_md5sum(filename)
                load_gencc(filename, file_hash, force=force)
            except OntologyBuilderDataUpToDateException:
                print("GenCC File hash is the same as last import")

        if filename := options.get("omim"):
            try:
                load_omim(filename, force=force)
            except OntologyBuilderDataUpToDateException:
                print("OMIM File hash is the same as last import")

        if filename := options.get("biomart"):
            try:
                load_biomart(filename, force=force)
            except OntologyBuilderDataUpToDateException:
                print("BioMart File hash is the same as last import")

        if filename := options.get("mondo"):
            try:
                load_mondo(filename, force)
            except OntologyBuilderDataUpToDateException:
                print("MONDO File hash is the same as last import")

        if filename := options.get("hpo"):
            try:
                load_hpo(filename, force)
            except OntologyBuilderDataUpToDateException:
                print("HPO File hash is the same as last import")

        if filename := options.get("phenotype_to_genes"):
            try:
                load_phenotype_to_genes(filename, force)
            except OntologyBuilderDataUpToDateException:
                print("Phenotype to Genes File hash is the same as last import")

        if options.get("omim_frequencies"):
            print("THIS FILE IS DEPRECATED, please use phenotype_to_genes.txt instead")

        print("*** If your instance uses Condition Text Matching, you might want to run:")
        print("python3 manage.py sync_condition_text_matches --obsolete")
        print("*** To make sure no matched terms have become obsolete")

        # Create a new OntologyVersion with all the new imports
        OntologyVersion.latest()

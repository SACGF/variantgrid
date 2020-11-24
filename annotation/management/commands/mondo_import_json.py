import json
import re
from typing import Optional

from django.core.management import BaseCommand

from annotation.models import MonarchDiseaseOntology, MonarchDiseaseOntologyGeneRelationship, MIMMorbid, \
    MonarchDiseaseOntologyMIMMorbid
from genes.models import GeneSymbol, GeneSymbolAlias

ID_EXTRACT_P = re.compile(r"^.*\/([A-Z]+)_([0-9]+)$")
HGNC_EXTRACT_P = re.compile(r"http://identifiers.org/hgnc/([0-9]+)")
OMIM_URL_P = re.compile(r"http://identifiers.org/omim/([0-9]+)")
GENE_SYMBOL_SEARCH = re.compile(r"([A-Z][A-Z,0-9]{2,})")

RELATIONS = {
    "http://purl.obolibrary.org/obo/RO_0004025": "disease causes dysfunction of",
    "http://purl.obolibrary.org/obo/RO_0004001": "has material basis in gain of function germline mutation in",
    "http://purl.obolibrary.org/obo/RO_0004021": "disease has basis in disruption of",
    "http://purl.obolibrary.org/obo/RO_0004020": "disease has basis in dysfunction of"
}
MATCH_TYPES = {
    "http://www.w3.org/2004/02/skos/core#exactMatch": "exact",
    "http://www.w3.org/2004/02/skos/core#closeMatch": "close",
    "http://www.w3.org/2004/02/skos/core#broadMatch": "broad",
    "http://www.w3.org/2004/02/skos/core#narrowMatch": "narrow"
}

class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--file', required=True)

    def handle(self, *args, **options):
        data_file = None

        with open(options["file"], 'r') as json_file:
            data_file = json.load(json_file)

        mondo_dict = dict()
        gene_dict = dict()
        relation_counts = dict()
        gene_relations = 0
        for graph in data_file.get("graphs", []):
            for node in graph.get("nodes", []):
                raw_id = None
                type = None
                if node.get("type") == "CLASS":
                    if node_id_full := node.get("id"):
                        if match := ID_EXTRACT_P.match(node_id_full):
                            type = match[1]
                            raw_id = match[2]
                            if type != "MONDO":
                                continue

                            label = node.get("lbl")
                            defn = None
                            synonyms = []
                            omim_relationships = []
                            genes_mentioned = set()
                            if meta := node.get("meta"):
                                defn = meta.get("definition", {}).get("val")
                                if defn:
                                    genes_mentioned = genes_mentioned.union(
                                        set(GENE_SYMBOL_SEARCH.findall(defn)))

                                for synonym in meta.get("synonyms", []):
                                    synonym_value = synonym.get("val")
                                    if synonym_value:
                                        synonyms.append(synonym_value)
                                        genes_mentioned = genes_mentioned.union(set(GENE_SYMBOL_SEARCH.findall(synonym_value)))

                                for bp in meta.get("basicPropertyValues", []):
                                    val = bp.get("val")
                                    if omim_match := OMIM_URL_P.match(val):
                                        omim = omim_match[1]
                                        pred = bp.get("pred")
                                        pred = MATCH_TYPES.get(pred)
                                        omim_relationships.append({
                                            "omim": omim,
                                            "pred": pred
                                        })

                            gene_relationships = dict()
                            for gs in genes_mentioned:
                                gene_relationships[gs] = "Mondo definition text includes gene symbol"

                            mondo_dict[node_id_full] = {
                                "id": f"{type}:{raw_id}",
                                "type": type,
                                "type_id": int(raw_id),
                                "label": label,
                                "description": defn,
                                # "synonyms": synonyms,
                                "gene_relationships": gene_relationships,
                                "omim_relationships": omim_relationships
                            }

                        # copy of id for gene symbol to gene symbol
                        elif match := HGNC_EXTRACT_P.match(node_id_full):
                            gene_symbol = node.get("lbl")
                            gene_dict[node_id_full] = {"hgnc_id": match[1], "gene_symbol": gene_symbol}

            for edge in graph.get("edges", []):
                if mondo_sub := mondo_dict.get(edge.get("sub")):
                    if gene_obj := gene_dict.get(edge.get("obj")):
                        relationship = edge.get("pred")
                        relationship = RELATIONS.get(relationship, relationship)
                        mondo_sub.get("gene_relationships")[gene_obj.get("gene_symbol")] = relationship
                        """
                        mondo_sub.get("gene_relationships").append({
                            "type": relationship,
                            "gene_symbol": gene_obj.get("gene_symbol")
                        })
                        """

        mondo_list = list()
        mondo_gene_list = list()
        mondo_omim_list = list()

        text_only_mondo_count = 0
        text_only_gene_count = 0

        for mondo in mondo_dict.values():

            mondo_id = mondo.get("type_id")

            md = MonarchDiseaseOntology(
                pk=mondo_id,
                name=mondo.get("label") or "",
                definition=mondo.get("description")
            )
            mondo_list.append(md)
            already_exist = set()

            for gene_symbol_str, relationship_type in mondo.get("gene_relationships").items():

                use_gene_symbol: Optional[GeneSymbol] = None

                use_gene_symbol = GeneSymbol.objects.filter(symbol=gene_symbol_str).first()
                if not use_gene_symbol:
                    alias: Optional[GeneSymbolAlias]
                    if alias := GeneSymbolAlias.objects.filter(alias=gene_symbol_str).first():
                        use_gene_symbol = alias.gene_symbol

                        if use_gene_symbol.symbol in mondo.get("gene_relationships"):
                            # we already have a direct link, ignore the dodgy via alias one
                            continue

                        if relationship_type == "Mondo definition text includes gene symbol":
                            relationship_type = f"CAUTION: Mondo definition text includes gene symbol ({gene_symbol_str}) which is an alias to "

                if use_gene_symbol:

                    key = f"{mondo_id}*{use_gene_symbol.symbol}"
                    if key in already_exist:
                        continue
                    already_exist.add(key)

                    mdgr = MonarchDiseaseOntologyGeneRelationship(
                        mondo_id=mondo_id,
                        relationship=relationship_type,
                        gene_symbol=use_gene_symbol
                    )

                    mondo_gene_list.append(mdgr)
                else:
                    pass
                    # print(f"Gene symbol {gene_relation.get('gene_symbol')} doesn't exist")

            for omim_relation in mondo.get("omim_relationships"):
                moim = MonarchDiseaseOntologyMIMMorbid(
                    mondo_id=mondo.get("type_id"),
                    relationship=omim_relation.get("pred"),
                    omim_id=omim_relation.get("omim")
                )
                mondo_omim_list.append(moim)

        print(f"Found {text_only_mondo_count} mondo records with {text_only_gene_count} genes in description and no formal relationship")
        print("About to update database")

        MonarchDiseaseOntologyGeneRelationship.objects.all().delete()
        MonarchDiseaseOntology.objects.all().delete()
        MonarchDiseaseOntologyMIMMorbid.objects.all().delete()

        MonarchDiseaseOntology.objects.bulk_create(mondo_list)
        MonarchDiseaseOntologyGeneRelationship.objects.bulk_create(mondo_gene_list)
        MonarchDiseaseOntologyMIMMorbid.objects.bulk_create(mondo_omim_list)

        print("Update complete")

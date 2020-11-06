import json
import re

from django.core.management import BaseCommand

from annotation.models import MonarchDiseaseOntology, MonarchDiseaseOntologyGeneRelationship
from genes.models import GeneSymbol

ID_EXTRACT_P = re.compile(r"^.*\/([A-Z]+)_([0-9]+)$")
HGNC_EXTRACT_P = re.compile(r"http://identifiers.org/hgnc/([0-9]+)")
RELATIONS = {
    "http://purl.obolibrary.org/obo/RO_0004025": "disease causes dysfunction of",
    "http://purl.obolibrary.org/obo/RO_0004001": "has material basis in gain of function germline mutation in",
    "http://purl.obolibrary.org/obo/RO_0004021": "disease has basis in disruption of",
    "http://purl.obolibrary.org/obo/RO_0004020": "disease has basis in dysfunction of"
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
                            if meta := node.get("meta"):
                                defn = meta.get("definition", {}).get("val")
                                for synonym in meta.get("synonyms", []):
                                    synonym_valu = synonym.get("val")
                                    if synonym_valu:
                                        synonyms.append(synonym_valu)
                            mondo_dict[node_id_full] = {
                                "id": f"{type}:{raw_id}",
                                "type": type,
                                "type_id": int(raw_id),
                                "label": label,
                                "description": defn,
                                # "synonyms": synonyms,
                                "relationships": []
                            }

                        elif match := HGNC_EXTRACT_P.match(node_id_full):
                            gene_symbol = node.get("lbl")
                            gene_dict[node_id_full] = {"hgnc_id": match[1], "gene_symbol": gene_symbol}

            for edge in graph.get("edges", []):
                if mondo_sub := mondo_dict.get(edge.get("sub")):
                    if gene_obj := gene_dict.get(edge.get("obj")):
                        relationship = edge.get("pred")
                        relationship = RELATIONS.get(relationship, relationship)
                        mondo_sub.get("relationships").append({
                            "type": relationship,
                            "gene_symbol": gene_obj.get("gene_symbol")
                        })

        mondo_list = list()
        mondo_relationship_list = list()

        for mondo in mondo_dict.values():
            md = MonarchDiseaseOntology(
                pk=mondo.get("type_id"),
                name=mondo.get("label") or "",
                definition=mondo.get("description")
            )
            mondo_list.append(md)
            for relation in mondo.get("relationships"):
                if gene_symbol := GeneSymbol.objects.filter(symbol=relation.get("gene_symbol")).first():
                    mdgr = MonarchDiseaseOntologyGeneRelationship(
                        mondo_id=mondo.get("type_id"),
                        relationship=relation.get("type"),
                        gene_symbol=gene_symbol
                    )
                    mondo_relationship_list.append(mdgr)
                else:
                    print(f"Gene symbol {relation.get('gene_symbol')} doesn't exist")

        print("About to update database")

        MonarchDiseaseOntologyGeneRelationship.objects.all().delete()
        MonarchDiseaseOntology.objects.all().delete()

        MonarchDiseaseOntology.objects.bulk_create(mondo_list)
        MonarchDiseaseOntologyGeneRelationship.objects.bulk_create(mondo_relationship_list)

        print("Update complete")

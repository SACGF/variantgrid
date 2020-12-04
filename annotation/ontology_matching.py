import urllib
from typing import Dict, Optional, List, Any
import re

import requests
from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK
from rest_framework.views import APIView

from annotation.models import MonarchDiseaseOntologyMIMMorbid, MonarchDiseaseOntologyGeneRelationship, \
    MonarchDiseaseOntology, MIMMorbid, HumanPhenotypeOntology


class OntologyContext:

    @property
    def name(self):
        raise NotImplementedError("Ontology context must implement name")

    def merge(self, merge_me):
        pass

    def as_json(self) -> Dict:
        raise NotImplementedError("Ontology context must implement as_json")


class OntologyContextSimilarMatch(OntologyContext):

    def __init__(self, matches_gene: bool, vce_id: int):
        self.matches_gene = matches_gene
        self.vce_id = vce_id

    def merge(self, merge_me):
        if self < merge_me:
            self.matches_gene = merge_me.matches_gene
            self.vce_id = merge_me.vce_id

    def as_json(self) -> Dict:
        return {
            "matches_gene": self.matches_gene,
            "vce_id": self.vce_id
        }

    def __lt__(self, other):
        if not self.matches_gene and other.matches_gene:
            return True
        return False

    @property
    def name(self):
        return "similar_match"


class OntologyContextSearched(OntologyContext):

    @property
    def name(self):
        return "searched"

    # TODO might include search score in future
    def as_json(self) -> Any:
        return True


class OntologyContextDirectReference(OntologyContext):

    @property
    def name(self):
        return "referenced"

    def as_json(self) -> Any:
        return True


class OntologyContextSelected(OntologyContext):

    @property
    def name(self):
        return "selected"

    def as_json(self) -> Any:
        return True


class OntologyContextMonarchLink(OntologyContext):

    def __init__(self, gene_symbol: str, relation: str):
        self.gene_symbol = gene_symbol
        self.relation = relation

    @property
    def name(self):
        return "Monarch"

    def as_json(self) -> Dict:
        return {
            "gene_symbol": self.gene_symbol,
            "relation": self.relation
        }


class OntologyContextPanelApp(OntologyContext):

    def __init__(self, gene_symbol: str, omim_id: int, phenotype_row: str, evidence: List[str]):
        self.gene_symbol = gene_symbol
        self.omim_id = omim_id
        self.phenotype_row = phenotype_row
        self.evidence = evidence

    def merge(self, merge_me):
        self.evidence = self.evidence + merge_me.evidence

    @property
    def name(self):
        return "PanelApp AU"

    def as_json(self):
        sorted_evidence = list(set(self.evidence))
        sorted_evidence.sort()
        return {
            "gene_symbol": self.gene_symbol,
            "omim_id": self.omim_id,
            "phenotype_row": self.phenotype_row,
            "evidence": sorted_evidence
        }


class OntologyMeta:

    # TODO support non-mondo
    # TODO if taken from a 3rd party they might have a name or description that we don't
    # could be done in meta if had to be?
    def __init__(self, term_id: str):
        self.term_id = term_id
        self.name = None
        self.definition = None

        if term_id.startswith(MonarchDiseaseOntology.PREFIX):
            if mondo := MonarchDiseaseOntology.objects.filter(pk=MonarchDiseaseOntology.id_as_int(term_id)).first():
                self.name = mondo.name
                self.definition = mondo.definition
        elif term_id.startswith(MIMMorbid.PREFIX):
            if omim := MIMMorbid.objects.filter(pk=MIMMorbid.id_as_int(term_id)).first():
                self.name = omim.description
            if mondo_mim := MonarchDiseaseOntologyMIMMorbid.objects.filter(omim_id=MIMMorbid.id_as_int(term_id)).select_related("mondo").first():
                mondo: MonarchDiseaseOntology = mondo_mim.mondo
                if not self.name:
                    self.name = mondo.name
                self.definition = f"(Taken from synonym {mondo.id_str})\n{mondo.definition}"
        elif term_id.startswith(HumanPhenotypeOntology.PREFIX):
            if hpo := HumanPhenotypeOntology.objects.filter(pk=HumanPhenotypeOntology.id_as_int(term_id)).first():
                self.name = hpo.name
                self.definition = hpo.definition

        # FIXME can we do this in a method on Mondo?
        if self.definition:
            self.definition = self.definition.replace(".nn", ".\n")

        self.contexts: Dict[str, OntologyContext] = dict()

    def add_context(self, context: OntologyContext):
        if existing := self.contexts.get(context.name):
            existing.merge(context)
        else:
            self.contexts[context.name] = context

    @property
    def url(self):
        # don't want to use class definitions as an instance might not exist
        # FIXME handle other types
        if self.term_id.startswith("MONDO"):
            return f'https://ontology.dev.data.humancellatlas.org/ontologies/mondo/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2{self.term_id.replace(":", "_")}'
        return None

    def as_json(self) -> Dict:
        sorted_names = list(self.contexts.keys())
        sorted_names.sort()
        contexts = {key: self.contexts[key].as_json() for key in sorted_names}

        return {
            "id": self.term_id,
            "title": self.name,
            "url": self.url,
            "definition": self.definition,
            "contexts": contexts
        }


class OntologyMatching:

    PANEL_APP_OMIM = re.compile(r"([0-9]{5,})")

    def __init__(self):
        self.term_map: Dict[str, OntologyMeta] = dict()

    def find_or_create(self, term_id: str) -> OntologyMeta:
        mondo = self.term_map.get(term_id)
        if not mondo:
            mondo = OntologyMeta(term_id=term_id)
            self.term_map[term_id] = mondo
        return mondo

    def populate_panel_app_remote(self, gene_symbol: str):
        results = requests.get(f'https://panelapp.agha.umccr.org/api/v1/genes/{gene_symbol}/').json().get(
            "results")

        for panel_app_result in results:
            phenotype_row: str
            for phenotype_row in panel_app_result.get("phenotypes", []):
                for omim_match in OntologyMatching.PANEL_APP_OMIM.finditer(phenotype_row):
                    omim_id = int(omim_match[1])
                    if mondo_id := MonarchDiseaseOntologyMIMMorbid.objects.filter(omim_id=omim_id).values_list(
                            "mondo_id", flat=True).first():
                        mondo_meta = self.find_or_create(MonarchDiseaseOntology.int_as_id(mondo_id))
                        mondo_meta.add_context(OntologyContextPanelApp(
                            gene_symbol=gene_symbol,
                            omim_id=omim_id,
                            phenotype_row=phenotype_row,
                            evidence=panel_app_result.get('evidence') or []
                        ))

    def populate_monarch_local(self, gene_symbol: str):
        gene_relationships = MonarchDiseaseOntologyGeneRelationship.objects.filter(
            gene_symbol=gene_symbol)

        mdgr: MonarchDiseaseOntologyGeneRelationship
        for mdgr in gene_relationships:
            self.find_or_create(term_id=MonarchDiseaseOntology.int_as_id(mdgr.mondo_id)).add_context(
                OntologyContextMonarchLink(
                    gene_symbol=gene_symbol,
                    relation=mdgr.relationship
                )
            )

    def select_term(self, term):
        self.find_or_create(term).add_context(OntologyContextSelected())

    def reference_term(self, term):
        self.find_or_create(term).add_context(OntologyContextDirectReference())

    def as_json(self) -> Dict:
        return {value.term_id: value.as_json() for key, value in self.term_map.items()}


class SearchMondoText(APIView):

    def get(self, request, **kwargs) -> Response:

        ontologyMatches = OntologyMatching()

        search_term = request.GET.get('search_term')
        # gene_symbol = request.GET.get('gene_symbol')
        # a regular escape / gets confused for a URL divider
        urllib.parse.quote(search_term).replace('/', '%252F')

        results = requests.get(f'https://api.monarchinitiative.org/api/search/entity/autocomplete/{search_term}', {
            "prefix": "MONDO",
            "rows": 6,
            "minimal_tokenizer": "false",
            "category": "disease"
        }).json().get("docs")

        for result in results:
            o_id = result.get('id')
            label = result.get('label')
            if label:
                label = label[0]
            # match = result.get('match')
            # if extracted := ID_EXTRACT_MINI_P.match(o_id):
                # id_part = extracted[1]
                # defn = mondo_defns.get(int(id_part))
            # highlight = result.get('highlight')

            onto = ontologyMatches.find_or_create(term_id=o_id)
            if not onto.name:
                onto.name = label
            onto.add_context(OntologyContextSearched())

        return Response(status=HTTP_200_OK, data=ontologyMatches.as_json())

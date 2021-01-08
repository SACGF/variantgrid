import urllib
from typing import Dict, Optional, List, Any
import re

import requests
from rest_framework.response import Response
from rest_framework.status import HTTP_200_OK
from rest_framework.views import APIView
from ontology.models import OntologyTerm, OntologyTermRelation, OntologySet
from ontology.panel_app_ontology import get_or_fetch_gene_relations


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


    # FIXME, maybe we consolidate all ontologys into a single table
    # TODO if taken from a 3rd party they might have a name or description that we don't
    # could be done in meta if had to be?
    def __init__(self, term_id: str):
        self.term_id = term_id
        self.term: Optional[OntologyTerm] = None
        self.parent: Optional[OntologyTerm] = None
        if term := OntologyTerm.objects.filter(id=term_id).first():
            self.object = term
            self.parent = OntologyTermRelation.parent_of(term)
        self.contexts: Dict[str, OntologyContext] = dict()

    def add_context(self, context: OntologyContext):
        if existing := self.contexts.get(context.name):
            existing.merge(context)
        else:
            self.contexts[context.name] = context

    @property
    def url(self):
        # don't want to use OntologyTerm definitions as entry might not exist
        # FIXME handle other types
        # FIXME can we share the URLs between db_regexes and here?
        if self.term_id.startswith("MONDO"):
            return f'https://vm-monitor.monarchinitiative.org/disease/{self.term_id}'
        return None

    def as_json(self) -> Dict:
        sorted_names = list(self.contexts.keys())
        sorted_names.sort()
        contexts = {key: self.contexts[key].as_json() for key in sorted_names}

        data = {
            "id": self.term_id,
            "url": self.url,
            "contexts": contexts
        }
        if self.term:
            data["title"] = self.term.name
            data["definition"] = self.term.definition
        else:
            data["title"] = "Not stored locally"

        if self.parent:
            data["parent_id"] = self.parent.id
            data["parent_title"] = self.parent.name
            """
            sibling_count = MonarchDiseaseOntologyRelationship.objects.filter(object=self.parent, relationship="is_a").count() - 1
            child_count = MonarchDiseaseOntologyRelationship.objects.filter(object=self.object, relationship="is_a").count()

            data["sibling_count"] = sibling_count
            data["children_count"] = child_count
            """
        return data


class OntologyMatching:

    # FIXME use the auto_matching!
    # FIXME use the auto_matching!
    # FIXME use the auto_matching!


    PANEL_APP_OMIM = re.compile(r"([0-9]{5,})")

    def __init__(self):
        self.term_map: Dict[str, OntologyMeta] = dict()

    def find_or_create(self, term_id: str) -> OntologyMeta:
        mondo = self.term_map.get(term_id)
        if not mondo:
            mondo = OntologyMeta(term_id=term_id)
            self.term_map[term_id] = mondo
        return mondo

    def populate_relationships(self, gene_symbol: str):
        relationships = get_or_fetch_gene_relations(gene_symbol=gene_symbol)
        for relationship in relationships:
            # always convert to MONDO for now
            if mondo_term := OntologyTermRelation.mondo_version_of(relationship.term):
                mondo_meta = self.find_or_create(mondo_term.id)
                if relationship.relation == OntologySet.PANEL_APP_AU:
                    mondo_meta.add_context(OntologyContextPanelApp(
                        gene_symbol=gene_symbol,
                        omim_id=relationship.term_id,
                        phenotype_row=relationship.extra.get("phenotype_row"),
                        evidence=relationship.extra.get("evidence")
                    ))
                else:
                    mondo_meta.add_context(OntologyContextMonarchLink(
                        gene_symbol=gene_symbol,
                        relation=relationship.relation
                    ))

    def select_term(self, term):
        self.find_or_create(term).add_context(OntologyContextSelected())

    def reference_term(self, term):
        self.find_or_create(term).add_context(OntologyContextDirectReference())

    def as_json(self) -> Dict:
        return [value.as_json() for value in self.term_map.values()]

    def __iter__(self):
        # TODO probably should sort as well
        return iter(self.term_map.values())


class SearchMondoText(APIView):

    MONDO_PATTERN = re.compile("MONDO:[0-9]+")

    def get(self, request, **kwargs) -> Response:

        ontologyMatches = OntologyMatching()

        search_term = request.GET.get('search_term') or ''
        # gene_symbol = request.GET.get('gene_symbol')
        # a regular escape / gets confused for a URL divider
        urllib.parse.quote(search_term).replace('/', '%252F')
        selected = [term.strip() for term in (request.GET.get('selected') or '').split(",") if term.strip()]
        for term in selected:
            ontologyMatches.select_term(term)

        if gene_symbol := request.GET.get('gene_symbol'):
            ontologyMatches.populate_relationships(gene_symbol)

        if SearchMondoText.MONDO_PATTERN.match(search_term):
            ontologyMatches.find_or_create(search_term)

            # TODO add children
        else:
            # Call for Text Matches
            try:
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

                    onto = ontologyMatches.find_or_create(term_id=o_id)
                    if not onto.name:
                        onto.name = label
                    onto.add_context(OntologyContextSearched())
            except:
                pass
                # TODO communicate to the user couldn't search mondo text search

        return Response(status=HTTP_200_OK, data=ontologyMatches.as_json())

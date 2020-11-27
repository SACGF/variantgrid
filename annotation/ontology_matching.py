from typing import Dict, Optional, List, Any
import re

import requests

from annotation.models import MonarchDiseaseOntologyMIMMorbid, MonarchDiseaseOntologyGeneRelationship, \
    MonarchDiseaseOntology


class OntologyContext:

    @property
    def name(self):
        raise NotImplementedError("Ontology context must implement name")

    def merge(self, merge_me):
        raise NotImplementedError("Ontology context must implement merge")

    def as_json(self) -> Dict:
        raise NotImplementedError("Ontology context must implement as_json")


class OntologyContextSelected(OntologyContext):

    def __init__(self):
        pass

    def merge(self, merge_me):
        pass

    @property
    def name(self):
        return "selected"

    def as_json(self) -> Dict:
        return {
            "selected": True
        }


class OntologyContextMonarchLink(OntologyContext):

    def __init__(self, gene_symbol: str, relation: str):
        self.gene_symbol = gene_symbol
        self.relation = relation

    def merge(self, merge_me):
        pass

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

        if term_id.startswith("MONDO"):
            self.mondo: Optional[MonarchDiseaseOntology] = None
            if mondo := MonarchDiseaseOntology.objects.filter(id=MonarchDiseaseOntology.mondo_id_as_int(term_id)).first():
                self.name = mondo.name
                self.definition = mondo.definition

        self.contexts: Dict[str, OntologyContext] = dict()

    def add_context(self, context: OntologyContext):
        if existing := self.contexts.get(context.name):
            existing.merge(context)
        else:
            self.contexts[context.name] = context

    @property
    def url(self):
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

    def _find_or_create(self, term_id: str) -> OntologyMeta:
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
                        mondo_meta = self._find_or_create(MonarchDiseaseOntology.mondo_int_as_id(mondo_id))
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
            self._find_or_create(term_id=MonarchDiseaseOntology.mondo_int_as_id(mdgr.mondo_id)).add_context(
                OntologyContextMonarchLink(
                    gene_symbol=gene_symbol,
                    relation=mdgr.relationship
                )
            )

    def select_term(self, term):
        self._find_or_create(term).add_context(OntologyContextSelected())

    def as_json(self) -> Dict:
        return {value.term_id: value.as_json() for key, value in self.term_map.items()}

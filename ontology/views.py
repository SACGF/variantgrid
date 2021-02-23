from typing import List

from django.contrib import messages
from django.views.generic import TemplateView

from ontology.models import OntologyTerm, OntologyTermRelation, OntologyService, OntologySnake, OntologyRelation


class OntologyTermView(TemplateView):

    template_name = "ontology/ontology_term.html"

    def get_context_data(self, **kwargs):
        term_id = self.kwargs.get("term")
        term = OntologyTerm.get_or_stub(term_id)
        if not term.is_stub:
            gene_relationships = None
            if term.ontology_service != OntologyService.HGNC:
                gene_relationships = OntologySnake.snake_from(term=term, to_ontology=OntologyService.HGNC)

            all_relationships: List[OntologyTermRelation] = OntologyTermRelation.relations_of(term)
            family_relationships = [relationship for relationship in all_relationships if relationship.relation == OntologyRelation.IS_A]
            regular_relationshiops = [relationship for relationship in all_relationships if relationship.relation != OntologyRelation.IS_A]
            parent_relationships = [relationship for relationship in family_relationships if relationship.source_term == term]
            child_relationships = [relationship for relationship in family_relationships if relationship.dest_term == term]

            return {
                "term": term,
                "gene_relationships": gene_relationships,
                "parent_relationships": parent_relationships,
                "regular_relationships": regular_relationshiops,
                "child_relationships": child_relationships
            }
        messages.add_message(self.request, messages.ERROR, "This term is not stored in our database")
        return {"term": term}

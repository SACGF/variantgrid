from django.contrib import messages
from django.views.generic import TemplateView

from ontology.models import OntologyTerm, OntologyTermRelation, OntologyService, OntologySnake


class OntologyTermView(TemplateView):

    template_name = "ontology/ontology_term.html"

    def get_context_data(self, **kwargs):
        term_id = self.kwargs.get("term")
        term = OntologyTerm.get_or_stub(term_id)
        if not term.is_stub:
            gene_relationships = None
            if term.ontology_service != OntologyService.HGNC:
                gene_relationships = OntologySnake.snake_from(term=term, to_ontology=OntologyService.HGNC)
            return {
                "term": term,
                "gene_relationships": gene_relationships,
                "relationships": OntologyTermRelation.relations_of(term)
            }
        else:
            messages.add_message(self.request, messages.ERROR, "This term is not stored in our database")
            return {"term": term}

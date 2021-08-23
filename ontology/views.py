from typing import List

from django.contrib import messages
from django.views.generic import TemplateView

from annotation.models import patients_qs_for_ontology_term
from classification.views.classification_datatables import ClassificationColumns
from library.utils import LimitedCollection
from ontology.models import OntologyTerm, OntologyTermRelation, OntologyService, OntologySnake, OntologyRelation
from ontology.panel_app_ontology import update_gene_relations


class OntologyTermView(TemplateView):

    template_name = "ontology/ontology_term.html"

    def get_context_data(self, **kwargs):
        term_id = self.kwargs.get("term")
        term = OntologyTerm.get_or_stub(term_id)
        if not term.is_stub:
            gene_relationships = None
            is_gene = term.ontology_service == OntologyService.HGNC
            if is_gene:
                update_gene_relations(term.name)
            else:
                gene_relationships = LimitedCollection(OntologySnake.snake_from(term=term, to_ontology=OntologyService.HGNC), 250)

            all_relationships: List[OntologyTermRelation] = OntologyTermRelation.relations_of(term)
            regular_relationships = list()
            parent_relationships = list()
            child_relationships = list()
            for relationship in all_relationships:
                if relationship.relation == OntologyRelation.IS_A:
                    if relationship.source_term == term:
                        parent_relationships.append(relationship)
                    else:
                        child_relationships.append(relationship)
                else:
                    regular_relationships.append(relationship)

            patients_qs = patients_qs_for_ontology_term(self.request.user, term)
            return {
                "term": term,
                "is_ontology": not is_gene,
                # gene relationships can be double counted which is a bit misleading
                "relationship_count": (len(gene_relationships) if gene_relationships else 0) + (len(all_relationships) if all_relationships else 0),
                "gene_relationships": gene_relationships,
                "parent_relationships": LimitedCollection(parent_relationships, 250) if not is_gene else None,
                "regular_relationships": LimitedCollection(regular_relationships, 250),
                "child_relationships": LimitedCollection(child_relationships, 250) if not is_gene else None,
                "datatable_config": ClassificationColumns(self.request),
                "patients_qs": patients_qs,
            }

        messages.add_message(self.request, messages.ERROR, "This term is not stored in our database")
        return {"term": term}

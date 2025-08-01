from django.contrib import messages
from django.shortcuts import get_object_or_404, redirect
from django.views.generic import TemplateView

from annotation.models import patients_qs_for_ontology_term
from library.utils import LimitedCollection
from ontology.models import OntologyTerm, OntologyTermRelation, OntologyService, OntologySnake, OntologyRelation, \
    ONTOLOGY_RELATIONSHIP_MINIMUM_QUALITY_FILTER
from ontology.panel_app_ontology import update_gene_relations


def ontology_term_text(request, ontology_service, name):
    """ Occasionally we have service + name but not the ID - this is a way of building an URL for that """
    ontology_service = OntologyService(ontology_service)  # Ensure valid
    term = get_object_or_404(OntologyTerm, name=name, ontology_service=ontology_service)
    return redirect(term)


class OntologyTermView(TemplateView):

    template_name = "ontology/ontology_term.html"

    def get_context_data(self, **kwargs):
        term_id = self.kwargs.get("term")
        term = OntologyTerm.get_or_stub(term_id)
        if not term.is_stub:
            gene_relationships = None
            term_relationships = None
            is_gene = term.ontology_service == OntologyService.HGNC
            if is_gene:
                update_gene_relations(term.name)
                # the reverse of gene_relationships
                term_relationships = OntologySnake.snake_from(term=term, to_ontology=OntologyService.MONDO, quality_filter=ONTOLOGY_RELATIONSHIP_MINIMUM_QUALITY_FILTER).snakes + \
                                     OntologySnake.snake_from(term=term, to_ontology=OntologyService.OMIM,  quality_filter=ONTOLOGY_RELATIONSHIP_MINIMUM_QUALITY_FILTER).snakes
                term_relationships = list(sorted((snake.reverse() for snake in term_relationships), key=lambda snake: snake.leaf_term))
                term_relationships = LimitedCollection(term_relationships, 250)
            else:
                raw_gene_relationships = sorted(OntologySnake.snake_from(term=term, to_ontology=OntologyService.HGNC, quality_filter=ONTOLOGY_RELATIONSHIP_MINIMUM_QUALITY_FILTER), key=lambda snake: snake.leaf_relationship.dest_term.short)
                gene_relationships = LimitedCollection(raw_gene_relationships, 250)

            regular_relationships = None
            parent_relationships = []
            child_relationships = []
            relationship_count = 0
            term_link = None

            if is_gene:
                term_link = term.gene_symbol_url()
            else:
                regular_relationships = []
                all_relationships: list[OntologyTermRelation] = OntologyTermRelation.relations_of(term)
                relationship_count = len(all_relationships)
                for relationship in all_relationships:
                    if relationship.relation == OntologyRelation.IS_A:
                        if relationship.source_term == term:
                            parent_relationships.append(relationship)
                        else:
                            child_relationships.append(relationship)

                    elif is_gene or relationship.dest_term.ontology_service != OntologyService.HGNC:
                        # gene symbols go into gene_relationships, no need to list them again in direct relationships
                        # though currently gene symbols don't do reverse gene_relationships, so still show everything that links in gene_symbol
                        regular_relationships.append(relationship)

            patients_qs = patients_qs_for_ontology_term(self.request.user, term)
            has_hierarchy = term.ontology_service in {OntologyService.MONDO, OntologyService.HPO}
            return {
                "term": term,
                "term_link": term_link,
                "is_ontology": not is_gene,
                "gene_relationship_count": len(gene_relationships) if gene_relationships else 0,
                "gene_relationships": gene_relationships,
                "term_relationships": term_relationships,
                "relationship_count": relationship_count,
                "parent_relationships": LimitedCollection(parent_relationships, 250) if has_hierarchy else None,
                "regular_relationships": LimitedCollection(regular_relationships, 250),
                "child_relationships": LimitedCollection(child_relationships, 250) if has_hierarchy else None,
                "patients_qs": patients_qs,
            }

        messages.add_message(self.request, messages.ERROR, "This term is not stored in our database.")
        return {"term": term}

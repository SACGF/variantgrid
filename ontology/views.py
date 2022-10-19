from typing import List

from django.conf import settings
from django.contrib import messages
from django.views.generic import TemplateView

from annotation.models import patients_qs_for_ontology_term
from library.utils import LimitedCollection
from ontology.models import OntologyTerm, OntologyTermRelation, OntologyService, OntologySnake, OntologyRelation, \
    GeneDiseaseClassification
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
                raw_gene_relationships = sorted(OntologySnake.snake_from(term=term, to_ontology=OntologyService.HGNC), key=lambda snake: snake.leaf_relationship.dest_term.short)
                gene_relationships = LimitedCollection(raw_gene_relationships, 250)

            all_relationships: List[OntologyTermRelation] = OntologyTermRelation.relations_of(term)
            regular_relationships = list()
            parent_relationships = list()
            child_relationships = list()
            weak_relationships = list()
            for relationship in all_relationships:

                if relationship.relation == OntologyRelation.IS_A:
                    if relationship.source_term == term:
                        parent_relationships.append(relationship)
                    else:
                        child_relationships.append(relationship)

                else:
                    include_direct_relationship = False
                    if is_gene:
                        include_direct_relationship = True
                    elif relationship.dest_term.ontology_service != OntologyService.HGNC:
                        include_direct_relationship = True
                    else:
                        if extra := relationship.extra:
                            if strongest := extra.get('strongest_classification'):
                                allowed_set = GeneDiseaseClassification.get_above_min(GeneDiseaseClassification.STRONG)
                                if strongest not in allowed_set:
                                    # include_direct_relationship = True
                                    weak_relationships.append(relationship)

                    if include_direct_relationship:
                        # gene symbols go into gene_relationships, no need to list them again in direct relationships
                        # though currently gene symbols don't do reverse gene_relationships, so still show everyhting that links in gene_symbol
                        regular_relationships.append(relationship)

            patients_qs = patients_qs_for_ontology_term(self.request.user, term)
            has_hierarchy = term.ontology_service in {OntologyService.MONDO, OntologyService.HPO}
            return {
                "term": term,
                "is_ontology": not is_gene,
                "gene_relationship_count": len(gene_relationships) if gene_relationships else 0,
                "gene_relationships": gene_relationships,
                "weak_relationships": weak_relationships,
                "relationship_count": len(all_relationships) if all_relationships else 0,
                "parent_relationships": LimitedCollection(parent_relationships, 250) if has_hierarchy else None,
                "regular_relationships": LimitedCollection(regular_relationships, 250),
                "child_relationships": LimitedCollection(child_relationships, 250) if has_hierarchy else None,
                "patients_qs": patients_qs,
            }

        messages.add_message(self.request, messages.ERROR, "This term is not stored in our database.")
        return {"term": term}

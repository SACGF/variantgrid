from django.core.management import BaseCommand
from django.db.models import Count

from classification.models import ClinVarExport, ClinVarAllele
from ontology.models import AncestorCalculator, OntologySnake


class Command(BaseCommand):

    def handle(self, *args, **options):
        print("ClinVar\tAllele\tCount\tSCVs\tConditions\tUmbrella\tMax steps to Umbrella\tJSON")
        duplicates = ClinVarExport.objects.filter(classification_based_on__isnull=False).values('clinvar_allele').annotate(total=Count('clinvar_allele')).order_by('-total')
        for duplicate in duplicates:
            total = duplicate.get("total")
            if total > 1:
                clinvar_allele_id = duplicate.get('clinvar_allele')
                clinvar_allele = ClinVarAllele.objects.get(pk=clinvar_allele_id)
                clinvar_exports: list[ClinVarExport] = list(ClinVarExport.objects.filter(clinvar_allele=clinvar_allele, classification_based_on__isnull=False).all())
                all_conditions = set()
                all_scvs = set()
                includes_multi_condition = False
                for clinvar_export in clinvar_exports:
                    resolved_condition = clinvar_export.condition_resolved
                    if len(resolved_condition.terms) != 1:
                        includes_multi_condition = True
                    else:
                        all_conditions.add(resolved_condition.terms[0])
                    if scv := clinvar_export.scv:
                        all_scvs.add(scv)

                common_ansestor = AncestorCalculator.common_ancestor(all_conditions)
                all_distances: list[int] = []
                conditions_with_distance: list[str] = []
                for condition in all_conditions:
                    snakes = OntologySnake.check_if_ancestor(condition, common_ansestor)
                    min_length = 100
                    for snake in snakes:
                        snake_length = len(snake.show_steps())
                        if snake_length < min_length:
                            min_length = snake_length
                    all_distances.append(min_length)
                    conditions_with_distance.append(f"{condition} ({min_length})")

                all_conditions_string = ", ".join(str(condition) for condition in sorted(conditions_with_distance))
                max_distance = max(all_distances)
                all_scv_str = ", ".join(sorted(all_scvs))

                print(f"{clinvar_allele.clinvar_key.name}\t{clinvar_allele.allele:CA}\t{total}\t{all_scv_str}\t{all_conditions_string}\t{common_ansestor}\t{max_distance}")

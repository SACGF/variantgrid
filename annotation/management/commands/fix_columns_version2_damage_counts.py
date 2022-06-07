import operator
from functools import reduce

from django.conf import settings
from django.core.management.base import BaseCommand
from django.db.models import Case, Value, IntegerField, When, Q, F

from annotation.models import VariantAnnotationVersion, VariantAnnotation


class Command(BaseCommand):
    def handle(self, *args, **options):
        pathogenic_rankscore = settings.ANNOTATION_MIN_PATHOGENIC_RANKSCORE
        patho_kwargs = {
            'bayesdel_pathogenic': Case(
                When(bayesdel_noaf_rankscore__gte=pathogenic_rankscore, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            ),
            'cadd_raw_pathogenic': Case(
                When(cadd_raw_rankscore__gte=pathogenic_rankscore, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            ),
            'clinpred_pathogenic': Case(
                When(clinpred_rankscore__gte=pathogenic_rankscore, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            ),
            'metalr_pathogenic': Case(
                When(metalr_rankscore__gte=pathogenic_rankscore, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            ),
            'revel_pathogenic': Case(
                When(revel_rankscore__gte=pathogenic_rankscore, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            ),
            'vest4_pathogenic': Case(
                When(vest4_rankscore__gte=pathogenic_rankscore, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            ),
        }

        benign_kwargs = {
            'bayesdel_benign': Case(
                When(bayesdel_noaf_rankscore__lt=pathogenic_rankscore, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            ),
            'cadd_raw_benign': Case(
                When(cadd_raw_rankscore__lt=pathogenic_rankscore, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            ),
            'clinpred_benign': Case(
                When(clinpred_rankscore__lt=pathogenic_rankscore, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            ),
            'metalr_benign': Case(
                When(metalr_rankscore__lt=pathogenic_rankscore, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            ),
            'revel_benign': Case(
                When(revel_rankscore__lt=pathogenic_rankscore, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            ),
            'vest4_benign': Case(
                When(vest4_rankscore__lt=pathogenic_rankscore, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            ),
        }

        for vav in VariantAnnotationVersion.objects.filter(columns_version=2):
            print(f"Updating {vav}...")
            columns = []
            for col in vav.get_pathogenic_prediction_funcs().keys():
                columns.append(Q(**{f"{col}__isnull": False}))
            q_any_non_null = reduce(operator.or_, columns)

            va_qs = VariantAnnotation.objects.filter(q_any_non_null, version=vav)
            va_qs = va_qs.annotate(**patho_kwargs, **benign_kwargs)
            va_qs.update(
                predictions_num_pathogenic=reduce(operator.add, [F(p) for p in patho_kwargs]),
                predictions_num_benign=reduce(operator.add, [F(b) for b in benign_kwargs]),
            )

import operator
from functools import reduce

from django.core.management.base import BaseCommand
from django.db.models import Case, Value, IntegerField, When, Q, F

from annotation.models import VariantAnnotationVersion, VariantAnnotation
from annotation.pathogenicity_predictions import TOOLS


class Command(BaseCommand):
    def handle(self, *args, **options):
        patho_kwargs = {}
        benign_kwargs = {}
        non_null_q = []
        for t in TOOLS:
            if t.raw_field and t.raw_pathogenic_threshold is not None:
                threshold = t.raw_pathogenic_threshold
                patho_kwargs[f"{t.raw_field}_pathogenic"] = Case(
                    When(**{f"{t.raw_field}__gte": threshold}, then=Value(1)),
                    default=Value(0), output_field=IntegerField())
                benign_kwargs[f"{t.raw_field}_benign"] = Case(
                    When(**{f"{t.raw_field}__lt": threshold}, then=Value(1)),
                    default=Value(0), output_field=IntegerField())
                non_null_q.append(Q(**{f"{t.raw_field}__isnull": False}))
            if t.pred_field and t.pred_pathogenic_values:
                values = list(t.pred_pathogenic_values)
                patho_kwargs[f"{t.pred_field}_pathogenic"] = Case(
                    When(**{f"{t.pred_field}__in": values}, then=Value(1)),
                    default=Value(0), output_field=IntegerField())
                benign_kwargs[f"{t.pred_field}_benign"] = Case(
                    When(Q(**{f"{t.pred_field}__isnull": False})
                         & ~Q(**{f"{t.pred_field}__in": values}), then=Value(1)),
                    default=Value(0), output_field=IntegerField())
                non_null_q.append(Q(**{f"{t.pred_field}__isnull": False}))

        q_any_non_null = reduce(operator.or_, non_null_q)
        for vav in VariantAnnotationVersion.objects.filter(columns_version=4):
            print(f"Updating {vav}...")
            va_qs = VariantAnnotation.objects.filter(q_any_non_null, version=vav)
            va_qs = va_qs.annotate(**patho_kwargs, **benign_kwargs)
            va_qs.update(
                predictions_num_pathogenic=reduce(operator.add, [F(p) for p in patho_kwargs]),
                predictions_num_benign=reduce(operator.add, [F(b) for b in benign_kwargs]),
            )

import operator
from functools import reduce
from typing import List, Tuple

import numpy as np
from django.core.management.base import BaseCommand
from django.db.models import Q

from analysis.grids import VariantGrid
from analysis.models import AnalysisNode
from annotation.models import AnnotationVersion
from library.guardian_utils import admin_bot
from snpdb.models import CustomColumnsCollection, GenomeBuild


class Command(BaseCommand):
    """ Occasionally we run into a JSON serialization issue because some data was inserted without
        converting np.NaN to None - this looks for it: """
    def handle(self, *args, **options):
        genome_build = GenomeBuild.grch38()
        annotation_version = AnnotationVersion.latest(genome_build)
        node, float_paths = self._get_float_paths(annotation_version)
        q_list = []
        for fp in float_paths:
            q_list.append(Q(**{fp: np.NaN}))
        q = reduce(operator.or_, q_list)
        qs = node._get_model_queryset()
        if values := qs.filter(q).values(*float_paths).first():
            print(f"Example variant: {values['id']}")
            for k, v in values:
                if v == np.NaN:
                    print(f"{k}={v}")
        else:
            print(f"No NaN found in {len(float_paths)} float paths")
            print(f"Checking size of Variant query against {annotation_version=}.....")
            print(f"Query was {qs.count()} records...")

    @staticmethod
    def _get_float_paths(annotation_version) -> Tuple[AnalysisNode, List[float]]:
        user = admin_bot()
        all_columns = CustomColumnsCollection.objects.get(name='All columns')
        node = AnalysisNode.objects.filter(analysis__annotation_version=annotation_version).select_subclasses().first()
        if not node:
            raise AnalysisNode.DoesNotExist(f"No nodes for {annotation_version=}")
        node.analysis.custom_columns_collection = all_columns  # Don't save this!
        grid = VariantGrid(user, node)
        float_fields = []
        for cm in grid.get_colmodels():
            if cm.get("sorttype") == "float":
                float_fields.append(cm["index"])
        return node, float_fields

from dataclasses import dataclass
from typing import Dict, List

from django.utils.safestring import SafeString

from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from annotation.models import EnsemblGeneAnnotation
from snpdb.models import CustomColumnsCollection


def get_columns_qs(custom_columns_collection=None):
    if custom_columns_collection is None:
        custom_columns_collection = CustomColumnsCollection.get_system_default()
    return custom_columns_collection.customcolumn_set.filter(column__model_field=True).order_by('sort_order')


def get_variant_annotation_data(variant, annotation_version, columns_qs):
    column_data = list(columns_qs.values_list("column__grid_column_name", 'column__variant_column', 'column__description'))
    all_columns = [i[1] for i in column_data]
    qs = get_variant_queryset_for_annotation_version(annotation_version)
    variant_columns_dict = qs.filter(pk=variant.pk).values(*all_columns)[0]

    variant_data = []
    for (grid_column_name, variant_column, description) in column_data:
        row = (grid_column_name, description, variant_columns_dict[variant_column])
        variant_data.append(row)
    return variant_data


GENE_ANNOTATION_PREFIX = 'variantannotation__gene__ensemblgeneannotation__'

@dataclass
class GeneAnnotationColumnData:
    id: str
    name: str
    description: str
    values: List[str]

    @property
    def value(self) -> SafeString:
        if self.values:
            return '<br/>'.join(self.values)
        else:
            return '<span class="no-value">-</span>'


def get_gene_annotation_column_data(ensembl_gene_annotation) -> List[GeneAnnotationColumnData]:
    variant_columns_qs = get_columns_qs().filter(column__variant_column__startswith=GENE_ANNOTATION_PREFIX)
    column_data = list(variant_columns_qs.values_list("column__grid_column_name", 'column__variant_column', 'column__description'))

    def variant_to_gene_annotation_column(variant_column):
        return variant_column.replace(GENE_ANNOTATION_PREFIX, '')

    # Convert from variant -> gene columns
    gene_columns = []
    for cd in column_data:
        gene_columns.append(variant_to_gene_annotation_column(cd[1]))

    gene_columns_dict = EnsemblGeneAnnotation.objects.filter(pk=ensembl_gene_annotation.pk).values(*gene_columns)[0]

    variant_data: Dict[str, GeneAnnotationColumnData] = dict()

    description: str
    for (grid_column_name, variant_column, description) in column_data:
        gene_annotation_column = variant_to_gene_annotation_column(variant_column)
        value = gene_columns_dict[gene_annotation_column]

        gac = variant_data.get(description.lower())
        if not gac:
            gac = GeneAnnotationColumnData(id=grid_column_name, name=grid_column_name, description=description, values=list())
            variant_data[description.lower()] = gac

        if value is not None:
            gac.values.append(value)

    return variant_data.values()

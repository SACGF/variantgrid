import random

from analysis.models.nodes.filters.filter_node import FilterNode, FilterNodeItem
from analysis.models.nodes.filters.gene_list_node import GeneListNode
from analysis.models.nodes.node_utils import update_analysis
from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.models import GeneListCategory, CustomTextGeneList
from snpdb.models import VariantGridColumn


def create_filter_node(analysis, column_name, column_filter):
    filter_node = FilterNode.objects.create(analysis=analysis)

    vg_column = VariantGridColumn.objects.get(variant_column=column_name)

    name = str(vg_column.grid_column_name)
    if column_filter == 'null':
        operation = 'nu'
        name += " is null"
        column_filter = ''
    else:
        operation = 'eq'
        name += f" = {column_filter}"
    filter_item = FilterNodeItem(filter_node=filter_node,
                                 sort_order=0,
                                 operation=operation,
                                 field=vg_column.variant_column,
                                 data=column_filter)
    filter_item.save()
    filter_node.name = name[:50]  # column limit

    return filter_node


def create_gene_list_node(analysis, gene_string):
    custom_text_gene_list = CustomTextGeneList(text=gene_string)
    create_custom_text_gene_list(custom_text_gene_list, analysis.user.username,
                                 GeneListCategory.NODE_CUSTOM_TEXT, hidden=True)

    gene_list_node = GeneListNode.objects.create(analysis=analysis,
                                                 accordion_panel=GeneListNode.CUSTOM_GENE_LIST,
                                                 custom_text_gene_list=custom_text_gene_list)
    return gene_list_node


def create_filter_child_node(node, column_name, column_filter):
    if column_name == 'gene_symbol' and column_filter != 'null':
        child_node = create_gene_list_node(node.analysis, column_filter)
    else:
        child_node = create_filter_node(node.analysis, column_name, column_filter)

    child_node.x = node.x + 50 + random.randrange(-10, 10)
    child_node.y = node.y + 100 + random.randrange(-10, 10)
    child_node.add_parent(node)
    child_node.ready = False
    child_node.save()

    update_analysis(child_node.analysis_id)

    return child_node

from django.db import models
from django.db.models.aggregates import Count
from django.db.models.deletion import CASCADE

from analysis.models.nodes.analysis_node import NodeVersion
from annotation.models.models import VariantAnnotation
from genes.models import Gene


# TODO: Refactor classes below into using annotation.gene_counts
class NodeGenesCountCollection(models.Model):
    node_version = models.ForeignKey(NodeVersion, on_delete=CASCADE)

    @staticmethod
    def get_or_create_gene_counts_qs_for_node(node, queryset):
        ncscc, created = NodeGenesCountCollection.objects.get_or_create(node_version=node.node_version)
        if created:
            queryset = queryset.filter(**{VariantAnnotation.GENE_COLUMN + "__isnull": False})
            count_qs = queryset.values_list(VariantAnnotation.GENE_COLUMN).distinct().annotate(Count('id'))
            data_list = []
            for gene_id, count in count_qs:
                data = NodeGenesCount(collection=ncscc,
                                      gene_id=gene_id,
                                      count=count)
                data_list.append(data)

            if data_list:
                NodeGenesCount.objects.bulk_create(data_list)

        return ncscc.nodegenescount_set.all()


class NodeGenesCount(models.Model):
    collection = models.ForeignKey(NodeGenesCountCollection, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    count = models.IntegerField(null=False)

import json
import operator
from collections import defaultdict
from functools import reduce

import numpy as np
from django.db.models.query_utils import Q
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from rest_framework.generics import get_object_or_404, RetrieveAPIView
from rest_framework.response import Response
from rest_framework.views import APIView

from genes.models import GeneVersion
from genes.views import get_coverage_stats
from library.constants import WEEK_SECS
from library.utils import defaultdict_to_dict
from seqauto.models import GoldCoverageSummary, EnrichmentKit
from seqauto.serializers import EnrichmentKitSerializer, \
    GoldCoverageSummarySerializer, EnrichmentKitSummarySerializer


class EnrichmentKitSummaryView(RetrieveAPIView):
    """ Doesn't return gene list so is much faster """
    serializer_class = EnrichmentKitSummarySerializer
    lookup_field = 'pk'

    def get_queryset(self):
        return EnrichmentKit.objects.all()


class EnrichmentKitView(RetrieveAPIView):
    serializer_class = EnrichmentKitSerializer
    lookup_field = 'pk'

    def get_queryset(self):
        return EnrichmentKit.objects.all()


class EnrichmentKitGeneCoverageView(APIView):

    def get_coverage_q(self, enrichment_kit):
        return Q(gene_coverage_collection__qcgenecoverage__qc__bam_file__unaligned_reads__sequencing_sample__enrichment_kit=enrichment_kit)

    @method_decorator(cache_page(WEEK_SECS))
    def get(self, request, *args, **kwargs):
        enrichment_kit_id = self.kwargs['enrichment_kit_id']
        gene_symbol = self.kwargs['gene_symbol']

        enrichment_kit = get_object_or_404(EnrichmentKit, pk=enrichment_kit_id)

        FIELDS = ("mean", "percent_20x")
        field_enrichment_kit_genes = defaultdict(lambda: defaultdict(dict))
        all_gene_values = {}
        filter_q = self.get_coverage_q(enrichment_kit)

        for gene in GeneVersion.objects.filter(gene_symbol=gene_symbol):
            base_gene_coverage_qs = gene.genecoveragecanonicaltranscript_set.all()
            enrichment_kit_data = get_coverage_stats(base_gene_coverage_qs, filter_q, FIELDS)
            for field_name in FIELDS:
                field_data = enrichment_kit_data.get(field_name, [])
                field_enrichment_kit_genes[field_name][gene.pk] = field_data
                all_gene_values[field_name] = field_data

        means = defaultdict()
        for k, v in all_gene_values.items():
            mean_value = np.mean(v)
            if np.isnan(mean_value):
                mean_value = -1
            means[k] = mean_value

        # Need to convert each
        field_enrichment_kit_genes = defaultdict_to_dict(field_enrichment_kit_genes)
        data = {"fields": field_enrichment_kit_genes, "mean": means}
        return Response(data)


class EnrichmentKitGeneGoldCoverageView(EnrichmentKitGeneCoverageView):

    def get_coverage_q(self, enrichment_kit):
        gold_q = Q(gene_coverage_collection__qcgenecoverage__qc__bam_file__unaligned_reads__sequencing_sample__sample_sheet__sequencing_run__gold_standard=True)
        return reduce(operator.and_, [super().get_coverage_q(enrichment_kit), gold_q])


class GoldCoverageSummaryView(APIView):

    def get(self, request, *args, **kwargs):
        enrichment_kit_id = self.kwargs['enrichment_kit_id']
        gene_symbol = self.kwargs['gene_symbol']

        enrichment_kit = get_object_or_404(EnrichmentKit, pk=enrichment_kit_id)
        # TODO: Handle multi-transcripts??
        try:
            gcs = GoldCoverageSummary.filter_for_kit_and_gene_symbols(enrichment_kit, [gene_symbol]).get()
            serializer = GoldCoverageSummarySerializer(gcs)
            data = serializer.data
        except GoldCoverageSummary.MultipleObjectsReturned:
            raise  # Should never happen!
        except:
            data = {}  # No coverage
        return Response(data)


class BatchGoldCoverageSummaryView(APIView):
    """ Needs to be a post as we can send a large number of genes """

    def post(self, request, enrichment_kit_id):
        gene_symbols_json = request.data["gene_symbols_json"]
        gene_symbols = json.loads(gene_symbols_json)

        enrichment_kit = get_object_or_404(EnrichmentKit, pk=enrichment_kit_id)
        gcs_qs = GoldCoverageSummary.filter_for_kit_and_gene_symbols(enrichment_kit, gene_symbols)
        serializer = GoldCoverageSummarySerializer(gcs_qs, many=True)
        data = serializer.data
        return Response(data)

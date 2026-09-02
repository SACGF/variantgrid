
from django.shortcuts import redirect

from library.utils import full_class_name
from snpdb.graphs import graphcache
from snpdb.graphs.allele_frequency_graph import AlleleFrequencyHistogramGraph
from snpdb.graphs.chromosome_density_graph import SampleChromosomeDensityGraph
from snpdb.graphs.chromosome_intervals_graph import ChromosomeIntervalsGraph
from snpdb.graphs.homozygosity_percent_graph import HomozygosityPercentGraph


def genomic_intervals_graph(request, genomic_intervals_collection_id):
    graph_class_name = full_class_name(ChromosomeIntervalsGraph)
    cached_graph = graphcache.async_graph(graph_class_name, genomic_intervals_collection_id)
    return redirect(cached_graph)


def chrom_density_graph(request, sample_id, cmap):
    graph_class_name = full_class_name(SampleChromosomeDensityGraph)
    cached_graph = graphcache.async_graph(graph_class_name, cmap, sample_id)
    return redirect(cached_graph)


def homozygosity_graph(request, sample_id, cmap):
    graph_class_name = full_class_name(HomozygosityPercentGraph)
    cached_graph = graphcache.async_graph(graph_class_name, cmap, sample_id)
    return redirect(cached_graph)


def sample_allele_frequency_histogram_graph(request, sample_id, min_read_depth):
    graph_class_name = full_class_name(AlleleFrequencyHistogramGraph)
    cached_graph = graphcache.async_graph(graph_class_name, sample_id, min_read_depth)
    return redirect(cached_graph)

import re

from django.conf import settings

from snpdb.models import GenomeBuild, Contig
from snpdb.search import search_receiver, SearchInputInstance, SearchExample


@search_receiver(
    search_type=GenomeBuild,
    pattern=re.compile("(GRCh37|GRCh38|T2T|CHM13)", re.IGNORECASE),
    sub_name="GenomeBuild",
    example=SearchExample(
        note="Genome Builds",
        examples=["GRCh37", "GRCh38"]
    ),
    admin_only=settings.SEARCH_CONTIG_GENOME_BUILD_ADMIN_ONLY
)
def genome_build_search(search_input: SearchInputInstance):
    yield GenomeBuild.objects.filter(name__iexact=search_input.search_string)


@search_receiver(
    search_type=Contig,
    pattern=re.compile(r"(chr(\d+|X|Y|M|MT)|NC_\d+.\d)", re.IGNORECASE),
    example=SearchExample(
        note="Contigs or Chromosomes",
        examples=["chrX", "NC_000007.13"]
    ),
    admin_only=settings.SEARCH_CONTIG_GENOME_BUILD_ADMIN_ONLY
)
def contig_search(search_input: SearchInputInstance):
    q = Contig.get_q(search_input.search_string, case_sensitive=False)
    yield Contig.objects.filter(q, genomebuildcontig__genome_build__in=search_input.genome_builds)

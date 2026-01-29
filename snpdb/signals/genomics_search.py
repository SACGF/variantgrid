import re

from django.conf import settings
from django.db.models import Q

from library.enums.log_level import LogLevel
from snpdb.models import GenomeBuild, Contig
from snpdb.search import search_receiver, SearchInputInstance, SearchExample, SearchMessageOverall


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
    all_contigs = Contig.objects.filter(q)
    q_search_build = Q(genomebuildcontig__genome_build__in=search_input.genome_builds)
    if contigs := list(all_contigs.filter(q_search_build)):
        yield contigs
    else:
        if other_contig := all_contigs.exclude(q_search_build).first():
            contig_builds = other_contig.genomebuildcontig_set.all()
            build_names = ", ".join(contig_builds.order_by("genome_build__name").values_list("genome_build__name", flat=True))
            msg = f"contig={other_contig.refseq_accession} is from builds {build_names} which are not enabled/annotated."
            yield SearchMessageOverall(msg, severity=LogLevel.WARNING)
        elif re.match(r"^NC_\d+.\d$", search_input.search_string):
            enabled_builds = GenomeBuild.get_enabled_builds_comma_separated_string()
            msg = f"contig={search_input.search_string} not in enabled builds: {enabled_builds}"
            yield SearchMessageOverall(msg, severity=LogLevel.WARNING)

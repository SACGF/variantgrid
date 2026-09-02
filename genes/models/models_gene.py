import json
import logging
import os
import re
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from functools import cached_property, total_ordering
from typing import Any, Optional, Union
from urllib.error import URLError

from cache_memoize import cache_memoize
from django.conf import settings
from django.contrib.auth.models import User
from django.core.cache import cache
from django.db import models
from django.db.models import QuerySet, TextField
from django.db.models.deletion import CASCADE, SET_NULL
from django.db.models.functions import Collate, Upper
from django.db.models.query_utils import Q
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel
from requests import RequestException

from genes.models_enums import AnnotationConsortium, GeneSymbolAliasSource, HGNCStatus, MANEStatus

# Re-exported: callers have long imported these from genes.models
from genes.transcript_errors import BadTranscript, MissingTranscript, NoTranscript
from genes.transcript_parts import TranscriptParts, get_transcript_id_and_version
from genes.transcript_sequence_retrieval import FetchedTranscriptSequence, TranscriptSequenceFetcher
from genes.transcripts_utils import get_lrg_and_t
from library.cache import timed_cache
from library.constants import DAY_SECS, HOUR_SECS, WEEK_SECS
from library.django_utils import SortByPKMixin
from library.django_utils.django_object_managers import ObjectManagerCachingRequest
from library.preview_request import PreviewData, PreviewModelMixin
from snpdb.models import GenomicCoordinates, Wiki
from snpdb.models.models_genome import Contig, GenomeBuild


class HGNCImport(TimeStampedModel):
    pass


class HGNC(models.Model):
    # pk = HGNC id with HGNC: stripped out
    alias_symbols = models.TextField()
    approved_name = models.TextField()
    ccds_ids = models.TextField(null=True, blank=True)
    ensembl_gene_id = models.TextField(null=True, blank=True)
    gene_group_ids = models.TextField(null=True, blank=True)
    gene_groups = models.TextField(null=True, blank=True)
    # Believe it or not, gene_symbol is not unique - e.g. MMP21 has multiple entries
    gene_symbol = models.ForeignKey('GeneSymbol', on_delete=CASCADE)
    hgnc_import = models.ForeignKey(HGNCImport, on_delete=CASCADE)
    location = models.TextField(null=True, blank=True)
    mgd_ids = models.TextField(null=True, blank=True)
    omim_ids = models.TextField(null=True, blank=True)
    previous_symbols = models.TextField(null=True, blank=True)
    refseq_ids = models.TextField(null=True, blank=True)
    rgd_ids = models.TextField(null=True, blank=True)
    status = models.CharField(max_length=1, choices=HGNCStatus.choices)
    ucsc_ids = models.TextField(null=True, blank=True)
    uniprot_ids = models.TextField(null=True, blank=True)
    uniprot = models.ForeignKey('UniProt', null=True, on_delete=SET_NULL)

    def __str__(self):
        return f"HGNC:{self.pk} approved symbol: {self.gene_symbol}, " \
               f"previous symbols: {self.previous_symbols}, alias_symbols: {self.alias_symbols}"

    @property
    def hgnc_id(self) -> str:
        return f"HGNC:{self.pk}"

    def get_absolute_url(self):
        safe_hgnc = f"HGNC_{self.pk}"
        return reverse('ontology_term', kwargs={"term": safe_hgnc})

    def url(self):
        return f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:{self.pk}"

    @cached_property
    def ccds_list(self):
        return (self.ccds_ids or '').split(",")

    @cached_property
    def gene_group_id_list(self):
        return (self.gene_group_ids or '').split(",")

    @cached_property
    def mgd_list(self):
        return (self.mgd_ids or '').split(",")

    @cached_property
    def rgd_list(self):
        return (self.rgd_ids or '').split(",")

    @cached_property
    def ucsc_list(self):
        return (self.ucsc_ids or '').split(",")

    @cached_property
    def uniprot_list(self) -> list['UniProt']:
        ulist = []
        if self.uniprot_ids:
            uniprot_ids = self.uniprot_ids.split(",")
            ulist = list(UniProt.objects.filter(pk__in=uniprot_ids))
        return ulist

    @staticmethod
    def id_by_accession(hgnc_prefix=True) -> dict:
        pk_qs = HGNC.objects.all().values_list("pk", flat=True)
        if hgnc_prefix:
            return {f"HGNC:{pk}": pk for pk in pk_qs}
        return {pk: pk for pk in pk_qs}


class UniProt(models.Model):
    # accession = Primary (citable) accession number (1st element in SwissProt record)
    accession = models.TextField(primary_key=True)
    cached_web_resource = models.ForeignKey('annotation.CachedWebResource', on_delete=CASCADE)
    function = models.TextField(null=True, blank=True)
    pathway = models.TextField(null=True, blank=True)
    pathway_interaction_db = models.TextField(null=True, blank=True)
    reactome = models.TextField(null=True, blank=True)
    tissue_specificity = models.TextField(null=True, blank=True)

    def __str__(self):
        return self.accession

class GeneSymbol(models.Model, PreviewModelMixin):
    """
        If you want to do a like on symbol and get an error:

            GeneSymbol.objects.filter(symbol__startswith='GATA')
            django.db.utils.NotSupportedError: nondeterministic collations are not supported for operator class "text_pattern_ops"

        Instead you need to do:

        GeneSymbol.get_deterministic_queryset().filter(symbol_deterministic__startswith='GATA')
    """
    symbol = TextField(primary_key=True, db_collation='case_insensitive')  # See note above if 'like' breaks

    objects = ObjectManagerCachingRequest()

    class Meta:
        base_manager_name = 'objects'

    @classmethod
    def get_deterministic_queryset(cls) -> QuerySet['GeneSymbol']:
        """ Adds 'symbol_deterministic' you can do like queries on """
        qs = cls.objects.all()
        return qs.annotate(symbol_deterministic=Collate("symbol", "und-x-icu"))

    @staticmethod
    def cast(symbol: Union[str, 'GeneSymbol']) -> Optional['GeneSymbol']:
        if isinstance(symbol, str):
            return GeneSymbol._cast(symbol)
        return symbol

    @staticmethod
    @timed_cache(ttl=3600)
    def _cast(symbol_str: str) -> Optional['GeneSymbol']:
        return GeneSymbol.objects.filter(symbol=symbol_str).first()

    @property
    def metrics_logging_key(self) -> tuple[str, Any]:
        return "gene_symbol", self.symbol

    @property
    def name(self):
        """ For use by TextPhenotypeMatch """
        return self.symbol

    def get_genes(self) -> QuerySet:
        # To match HPO/OMIM so it can be used interchangeably during phenotype matching
        return Gene.objects.filter(~Q(identifier__startswith="unknown_"), geneversion__gene_symbol=self).distinct()

    @cached_property
    def genes(self) -> list['Gene']:
        # returns cached set of genes associated with this symbol
        # use over get_genes when possible
        return list(self.get_genes().all())

    def latest_gene_version(self, genome_build: GenomeBuild):
        return self.geneversion_set.filter(genome_build=genome_build).order_by("-version").first()

    def get_absolute_url(self):
        return reverse("view_gene_symbol", kwargs={"gene_symbol": self.symbol})

    @cached_property
    def alias_meta(self) -> 'GeneSymbolAliasesMeta':
        return GeneSymbolAliasesMeta(self)

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-dna"

    @property
    def preview(self):
        return self.preview_with(title="")

    def has_different_genes(self, other: 'GeneSymbol') -> bool:
        """
        Tries to work out if genes are equivilent, not that sometimes RefSeq or ensembl assign gene ids to both the
        symbol and the alias, but the other consortium only assigns to one. In that case we'd still like to treat them
        as the "same"
        """
        my_genes = set(self.genes)
        other_genes = set(other.genes)

        all_genes = my_genes.union(other_genes)
        source_has_extra = False
        other_has_extra = False
        for g in all_genes:
            if g in my_genes and g not in other_genes:
                source_has_extra = True
            elif g in other_genes and g not in my_genes:
                other_has_extra = True

        return source_has_extra and other_has_extra

    @staticmethod
    def overlapping_variant(variant, variant_annotation_version) -> QuerySet['GeneSymbol']:
        vta_qs = variant.varianttranscriptannotation_set.filter(version=variant_annotation_version)
        symbol_names = list(vta_qs.values_list("transcript_version__gene_version__gene_symbol", flat=True).distinct())
        return GeneSymbol.objects.filter(pk__in=symbol_names)

    def __lt__(self, other):
        return self.symbol < other.symbol

    def __str__(self):
        return self.symbol

    @staticmethod
    def get_upper_case_lookup():
        return dict(GeneSymbol.objects.annotate(uc_symbol=Upper("symbol")).values_list("uc_symbol", "symbol"))


class GeneSymbolAlias(TimeStampedModel):
    """ Gene Aliases record keep track of "source" and are from:
        NCBI:
         * Source: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
         * Code: python3 manage.py import_ncbi_gene_info <file>

        HGNC:
         * Source: https://www.genenames.org/cgi-bin/download
         * Code: python3 manage.py hgnc_gene_symbols_import <file>

        UCSC: We no longer use UCSC aliases, they will only exist upgraded legacy systems
         * Source: https://genome.ucsc.edu/cgi-bin/hgTables?command=start export kgAlias table
         * Code: N/A - obsolete
    """
    alias = TextField(db_collation='case_insensitive')
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    source = models.CharField(max_length=1, choices=GeneSymbolAliasSource.choices)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    description = models.TextField(null=True)

    class Meta:
        unique_together = ('alias', 'gene_symbol')

    @property
    def match_info(self) -> str:
        return f"{self.alias} is an alias for {self.gene_symbol_id} ({self.get_source_display()})"

    def __str__(self):
        return f"{self.gene_symbol_id} : {self.match_info}"

    def get_absolute_url(self):
        """ So search sends it to the symbol """
        return reverse("view_gene_symbol", kwargs={"gene_symbol": self.gene_symbol_id})

    @staticmethod
    def get_upper_case_lookup():
        return {a: (gs, alias_id) for a, gs, alias_id in GeneSymbolAlias.objects.values_list("alias", "gene_symbol", "id")}


@dataclass
@total_ordering
class GeneSymbolAliasSummary:
    other_obj: GeneSymbol
    other_symbol: str
    source: str  # HGNC etc
    my_symbol_is_main: bool  # true if the other symbol is an alias for this symbol, false if this symbol is an alias for the other
    different_genes: bool  # if true, then this should only be considered an alias with a priviso, and not used in automatic alias calculations

    def __lt__(self, other):
        return self.other_symbol < other.other_symbol

    @property
    def other_symbol_in_database(self) -> bool:
        return self.other_obj is not None


class GeneSymbolAliasesMeta:

    def __init__(self, gene_symbol: GeneSymbol):
        self.gene_symbol = gene_symbol
        self.alias_list: list[GeneSymbolAliasSummary] = []

        symbol = self.gene_symbol.symbol

        for alias in GeneSymbolAlias.objects.filter(alias=symbol):
            self.alias_list.append(
                GeneSymbolAliasSummary(
                    other_obj=alias.gene_symbol,
                    other_symbol=alias.gene_symbol.symbol,
                    source=alias.get_source_display(),
                    my_symbol_is_main=False,
                    different_genes=self.gene_symbol.has_different_genes(alias.gene_symbol)
                )
            )
        for alias in GeneSymbolAlias.objects.filter(gene_symbol=self.gene_symbol):
            other_gene_symbol = GeneSymbol.objects.filter(symbol=alias.alias).first()
            different_genes = False
            if other_gene_symbol:
                different_genes = self.gene_symbol.has_different_genes(other_gene_symbol)
            self.alias_list.append(
                GeneSymbolAliasSummary(
                    other_obj=other_gene_symbol,
                    other_symbol=alias.alias,
                    source=alias.get_source_display(),
                    my_symbol_is_main=True,
                    different_genes=different_genes
                )
            )

    @cached_property
    def genes(self) -> set['Gene']:
        """
        Returns a set of genes associated with all safe aliases to/from the primary Gene Symbol.
        (Even though we only look at "safe" aliases, e.g. ones where each symbol must be a subset of the other,
        looking through these aliases still catch where Refseq assigned a Gene ID to both but Ensembl only assigned
        their Gene ID to one and ignore the other)
        """
        gene_set: set[Gene] = set(self.gene_symbol.genes)
        for alias_summary in self.alias_list:
            if not alias_summary.different_genes and alias_summary.other_obj:
                gene_set = gene_set.union(alias_summary.other_obj.genes)
        return gene_set

    @cached_property
    def alias_symbol_strs(self) -> list[str]:
        gene_symbol_strs: set[str] = {self.gene_symbol.symbol}
        for alias_summary in self.alias_list:
            if not alias_summary.different_genes:
                gene_symbol_strs.add(alias_summary.other_symbol)
        return sorted(gene_symbol_strs)

    @cached_property
    def alias_symbols_in_db(self) -> list[GeneSymbolAlias]:
        return sorted([alias for alias in self.alias_list if not alias.different_genes and alias.other_symbol_in_database])

    @cached_property
    def aliases_out(self) -> list[GeneSymbolAliasSummary]:
        return sorted([alias for alias in self.alias_list if not alias.my_symbol_is_main])

    @cached_property
    def aliases_in(self) -> list[GeneSymbolAliasSummary]:
        return sorted([alias for alias in self.alias_list if alias.my_symbol_is_main])


class GeneAnnotationImport(TimeStampedModel):
    """ A GTF file imported via 'python3 manage import_gene_annotation'

        Many gene/transcript versions are shared among GTF annotations, so a GeneVersion/TranscriptVersion is only
        created the first time it's seen (linked back to input which created it via 'import_source') """
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    url = models.TextField()

    def __str__(self):
        return self.url


class Gene(PreviewModelMixin, models.Model):
    """ A stable identifier - build independent - has build specific versions with gene details """
    FAKE_GENE_ID_PREFIX = "unknown_"  # Legacy from when we allowed inserting GenePred w/o GFF3
    identifier = models.TextField(primary_key=True)
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)
    summary = models.TextField(null=True, blank=True)  # Only used by RefSeq

    @property
    def prefixed_identifier(self) -> str:
        if self.annotation_consortium == AnnotationConsortium.REFSEQ:
            return f"GeneID:{self.identifier}"
        else:
            return self.identifier

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-dna"

    @property
    def preview(self) -> 'PreviewData':
        return self.preview_with(
            identifier=self.prefixed_identifier,
            summary=self.summary
        )

    @property
    def is_legacy(self):
        """ Required internally, but probably shouldn't be shown to the user """
        return self.identifier.startswith(Gene.FAKE_GENE_ID_PREFIX)

    def get_external_url(self):
        if self.annotation_consortium == AnnotationConsortium.REFSEQ:
            return f"https://www.ncbi.nlm.nih.gov/gene/{self.identifier}"
        if self.annotation_consortium == AnnotationConsortium.ENSEMBL:
            return f"https://ensembl.org/Homo_sapiens/Gene/Summary?g={self.identifier}"
        raise ValueError(f"Unknown external url for {self}")

    def latest_gene_version(self, genome_build: GenomeBuild):
        return self.geneversion_set.filter(genome_build=genome_build).order_by("-version").first()

    def get_gene_symbol(self, genome_build: GenomeBuild) -> GeneSymbol:
        return self.latest_gene_version(genome_build).gene_symbol

    def get_symbols(self) -> QuerySet:
        """ This can change over time as versions are assigned different symbols... """
        return GeneSymbol.objects.filter(geneversion__gene=self).distinct()

    def has_versions(self):
        # RefSeq doesn't have gene versions
        return self.annotation_consortium == AnnotationConsortium.ENSEMBL

    def get_absolute_url(self):
        return reverse("view_gene", kwargs={"gene_id": self.pk})

    @staticmethod
    def known_gene_ids(annotation_consortium=None):
        qs = Gene.objects.all()
        if annotation_consortium:
            qs = qs.filter(annotation_consortium=annotation_consortium)
        return set(qs.values_list("identifier", flat=True))

    @staticmethod
    def delete_orphaned_fake_genes():
        used_genes = TranscriptVersion.objects.filter(gene_version__gene__identifier__startswith=Gene.FAKE_GENE_ID_PREFIX).values_list("gene_version__gene")
        qs = Gene.objects.filter(identifier__startswith=Gene.FAKE_GENE_ID_PREFIX).exclude(identifier__in=used_genes)
        ret = qs.delete()
        if ret:
            logging.info("Deleted orphaned %s records: %s", Gene.FAKE_GENE_ID_PREFIX, ret)

    def get_vep_canonical_transcript(self, variant_annotation_version: 'VariantAnnotationVersion') -> Optional['Transcript']:
        """ This may be slow. It requires an annotated (non-ref) variant in the gene """
        vta = self.varianttranscriptannotation_set.filter(version=variant_annotation_version, canonical=True).first()
        transcript = None
        if vta:
            transcript = vta.transcript
        return transcript

    def __lt__(self, other):
        if self.annotation_consortium == AnnotationConsortium.REFSEQ:
            try:
                return int(self.identifier) < int(other.identifier)
            except ValueError:
                pass

        return self.identifier < other.identifier

    def __str__(self):
        return f"{self.prefixed_identifier} ({self.get_annotation_consortium_display()})"


class GeneVersion(models.Model):
    """ A specific version of a Gene for a particular version/genome build
        Genes/TranscriptVersion needs to be able to represent both RefSeq and Ensembl """
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    version = models.IntegerField()  # RefSeq GeneIDs are always 0 (not versioned) need non-null for unique_together
    # symbol can be null as Ensembl has genes w/o symbols, e.g. ENSG00000238009 (lncRNA)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    # HGNC assignment - hgnc_identifier is set from the annotation file
    #                   hgnc (ForeignKey) is linked to our HGNC models
    hgnc_identifier = models.IntegerField(null=True)
    hgnc = models.ForeignKey(HGNC, null=True, on_delete=SET_NULL)
    description = models.TextField(null=True)
    biotype = models.TextField(null=True)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    import_source = models.ForeignKey(GeneAnnotationImport, on_delete=CASCADE)

    class Meta:
        unique_together = ("gene", "version", "genome_build")

    @cached_property
    def accession(self):
        # RefSeq has no versions, so is always 0
        # Ensembl had some old ones with no version provided, they are 0 as well
        if self.version:
            acc = f"{self.gene_id}.{self.version}"
        else:
            acc = self.gene_id
        return acc

    @property
    def name(self):
        return self.gene_symbol_id

    @property
    def chrom(self):
        return self._transcript_extents["chrom"]

    @property
    def start(self):
        return int(self._transcript_extents["start"])

    @property
    def end(self):
        return int(self._transcript_extents["end"])

    @property
    def strand(self):
        return self._transcript_extents["strand"]

    @cached_property
    def coordinate(self) -> str:
        """ 1-based for humans """
        try:
            return f"{self.chrom}:{self.start + 1}-{self.end} ({self.strand})"
        except Exception:
            return ""

    @cached_property
    def _transcript_extents(self):
        """ Stores chroms/min_start/max_end/strands - calculated from all linked TranscriptVersions """

        FIELDS = ["contig", "strand"]
        all_transcripts_data = defaultdict(set)
        for tv in self.transcriptversion_set.all():
            for f in FIELDS:
                all_transcripts_data[f].add(tv.genome_build_data[f])
            all_transcripts_data["start"].add(tv.start)
            all_transcripts_data["end"].add(tv.end)

        # Sometimes chrom is a contig so we'll end up with chrom as "3,NC_000003.11" - check only 1 and convert to name
        chrom_set = all_transcripts_data["contig"]
        contigs = {self.genome_build.chrom_contig_mappings[chrom] for chrom in chrom_set}
        if len(contigs) != 1:
            raise ValueError(f"{self}: 'chrom' ({all_transcripts_data['contig']}) didn't map to exactly 1 contig: '{contigs}'")

        strand_set = all_transcripts_data["strand"]
        if len(strand_set) != 1:
            raise ValueError(f"{self}: Not exactly 1 value for 'strand', was: {strand_set}")

        return {
            "chrom": contigs.pop().name,
            "start": min(all_transcripts_data["start"]),
            "end": max(all_transcripts_data["end"]),
            "strand": strand_set.pop(),
        }

    @staticmethod
    def get_gene_id_and_version(gene_accession: str) -> tuple[str, Optional[int]]:
        parts = gene_accession.split(".")
        if len(parts) == 2:
            identifier = str(parts[0])
            version = int(parts[1])
        else:
            identifier, version = gene_accession, None
        return identifier, version

    @staticmethod
    def id_by_accession(genome_build: GenomeBuild = None, annotation_consortium=None) -> dict[str, int]:
        filter_kwargs = {}
        if genome_build:
            filter_kwargs["genome_build"] = genome_build
        if annotation_consortium:
            filter_kwargs["gene__annotation_consortium"] = annotation_consortium
        gene_version_qs = GeneVersion.objects.filter(**filter_kwargs)
        ids_by_accession = {}  # Uses version if non-zero
        for (pk, gene_id, version) in gene_version_qs.values_list("pk", "gene_id", "version"):
            if version:
                gene_accession = f"{gene_id}.{version}"
            else:
                gene_accession = gene_id
            ids_by_accession[gene_accession] = pk
        return ids_by_accession

    def __str__(self):
        return f"{self.accession} ({self.gene_symbol}/{self.genome_build})"


class Transcript(models.Model, PreviewModelMixin):
    """ A stable identifier - has versions with actual transcript details """
    identifier = models.TextField(primary_key=True)
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-timeline"

    def get_absolute_url(self):
        kwargs = {"transcript_id": self.identifier}
        return reverse("view_transcript", kwargs=kwargs)

    def get_external_url(self, genome_build: GenomeBuild):
        if self.annotation_consortium == AnnotationConsortium.REFSEQ:
            return f"https://www.ncbi.nlm.nih.gov/nuccore/{self.identifier}"
        if self.annotation_consortium == AnnotationConsortium.ENSEMBL:
            if genome_build.name == "GRCh37":
                return f"https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?t={self.identifier}"
            return f"https://ensembl.org/Homo_sapiens/Transcript/Summary?t={self.identifier}"
        raise ValueError(f"Unknown external url for {self}")

    def latest_version(self, genome_build: GenomeBuild):
        return self.transcriptversion_set.filter(genome_build=genome_build).order_by("version").last()

    @staticmethod
    def known_transcript_ids(genome_build=None, annotation_consortium=None) -> set[str]:
        filter_kwargs = {}
        require_distinct = False
        if genome_build:
            filter_kwargs["transcriptversion__genome_build"] = genome_build
            require_distinct = True
        if annotation_consortium:
            filter_kwargs["annotation_consortium"] = annotation_consortium

        qs = Transcript.objects.all()
        if filter_kwargs:
            qs = qs.filter(**filter_kwargs)
        if require_distinct:
            qs = qs.distinct()
        return set(qs.values_list("identifier", flat=True))

    def __str__(self):
        return self.identifier


class TranscriptVersion(SortByPKMixin, models.Model, PreviewModelMixin):
    """ We store the ID and version separately, ie:
        ENST00000284274.4 => transcript=ENST00000284274, version=4

        Ensembl ID info: https://m.ensembl.org/Help/Faq?id=488

        There's currently multiple TranscriptVersion per genome build, this should probably be changed to only having
        1, merging in TranscriptVersionSequenceInfo and moving the data (which contains exons etc) into a related object

        A useful query to get the latest version for each transcript is:
        qs.order_by("transcript_id", "-version").distinct("transcript_id")
    """
    transcript = models.ForeignKey(Transcript, on_delete=CASCADE)
    version = models.IntegerField()
    gene_version = models.ForeignKey(GeneVersion, on_delete=CASCADE)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    contig = models.ForeignKey(Contig, on_delete=CASCADE)  # Optimisation to restrict Variant queries
    import_source = models.ForeignKey(GeneAnnotationImport, on_delete=CASCADE)
    biotype = models.TextField(null=True)  # Ensembl has gene + transcript biotypes
    data = models.JSONField(null=False, blank=True, default=dict)  # for cdot data

    @classmethod
    def preview_icon(cls) -> str:
        return Transcript.preview_icon()

    @property
    def preview(self) -> PreviewData:
        return self.preview_with(
            identifier=f"{self.transcript.identifier}.{self.version}"
        )

    # These are in data.tags
    CANONICAL_SCORES = {
        "MANE Select": 2,
        "MANE_Select": 2,
        "RefSeq Select": 1,
        # Some way to find canonical in Ensembl GRCh37?
    }

    class Meta:
        unique_together = ("transcript", "version", "genome_build")

    @cached_property
    def _transcript_regions(self) -> tuple[list, list, list]:
        """ Returns 5'UTR, CDS, 3'UTR """
        cds_start = self.genome_build_data.get("cds_start")
        cds_end = self.genome_build_data.get("cds_end")
        left_utr = []
        cds = []
        right_utr = []
        for exon_start, exon_end, *_ in self.genome_build_data["exons"]:
            if cds_start:
                if exon_start < cds_start:
                    left_end = min(cds_start, exon_end)
                    left_utr.append((exon_start, left_end))

                if exon_end > cds_start and exon_start < cds_end:
                    start = max(exon_start, cds_start)
                    end = min(exon_end, cds_end)
                    cds.append((start, end))

                if exon_end > cds_end:
                    right_start = max(exon_start, cds_end)
                    right_utr.append((right_start, exon_end))
            else:
                left_utr.append((exon_start, exon_end))

        if self.genome_build_data["strand"] == '+':
            return left_utr, cds, right_utr
        else:
            return right_utr, cds, left_utr

    @property
    def is_coding(self) -> bool:
        return "start_codon" in self.data

    @cached_property
    def fivep_utr(self):
        return self._transcript_regions[0]

    @cached_property
    def cds(self):
        return self._transcript_regions[1]

    @cached_property
    def threep_utr(self):
        return self._transcript_regions[2]

    @property
    def as_parts(self):
        return TranscriptParts(self.transcript.identifier, self.version)

    @staticmethod
    def transcript_parts(identifier: str) -> TranscriptParts:
        # TODO - is this redundant to HGVSMatcher's get_transcript_parts?
        t_regex = re.compile(r"^([_A-Z0-9]+)(?:[.]([0-9]+))?$", re.RegexFlag.IGNORECASE)
        if m := t_regex.match(identifier):
            version = m.group(2)
            if version:
                version = int(version)
            return TranscriptParts(m.group(1), version)
        raise ValueError(f'Invalid transcript identifier {identifier}')

    @staticmethod
    def get(accession: str, genome_build: GenomeBuild, annotation_consortium=None):
        """
        @param accession transcript_id w/optional version
        """
        transcript_id, version = TranscriptVersion.transcript_parts(accession)
        kwargs = {}
        if annotation_consortium is not None:
            kwargs["transcript__annotation_consortium"] = annotation_consortium
        if version:
            kwargs["version"] = version
        t_qs = TranscriptVersion.objects.filter(transcript_id=transcript_id,
                                                genome_build=genome_build, **kwargs)
        if version:
            transcript_version = t_qs.get()  # Should only be 1
        else:
            transcript_version = t_qs.order_by("version").last()
            if not transcript_version:
                t_qs.get()  # Throw does not exist
        return transcript_version

    @staticmethod
    def get_ensembl(identifier, genome_build: GenomeBuild):
        return TranscriptVersion.get(identifier, genome_build, AnnotationConsortium.ENSEMBL)

    @staticmethod
    def get_refseq(identifier, genome_build: GenomeBuild):
        return TranscriptVersion.get(identifier, genome_build, AnnotationConsortium.REFSEQ)

    @staticmethod
    def get_accession(transcript_id, version):
        if version is not None:
            acc = f"{transcript_id}.{version}"
        else:
            acc = transcript_id
        return acc

    @cached_property
    def accession(self):
        return TranscriptVersion.get_accession(self.transcript_id, self.version)

    @property
    def annotation_consortium(self):
        """ This is used in search results """
        return self.transcript.annotation_consortium

    @property
    def gene(self):
        return self.gene_version.gene

    @cached_property
    def gene_symbol(self):
        """ Returns HGNC symbol if available (to keep consistency between builds) or GeneVersion symbol (from GFF)
            GeneVersion symbol from GFF can diverge e.g. Entrez GeneID: 6901 - TAZ(37) and TAFAZZIN(38) """
        if hgnc := self.gene_version.hgnc:
            gene_symbol = hgnc.gene_symbol
        else:
            gene_symbol = self.gene_version.gene_symbol
        return gene_symbol

    @cached_property
    def hgvs_ok(self) -> bool:
        if self.has_valid_data:
            if settings.HGVS_VALIDATE_REFSEQ_TRANSCRIPT_LENGTH:
                if self.transcript.annotation_consortium == AnnotationConsortium.REFSEQ:
                    return bool(self.sequence_length_matches_exon_length_ignoring_poly_a_tail)
            return True
        return False

    @property
    def sequence_info(self):
        return TranscriptVersionSequenceInfo.get(self.accession)

    @property
    def sequence_poly_a_tail(self) -> int:
        """ Returns length of polyA tail if ALL bases after sum of exon lengths are A """
        if self.sequence_info.length > self.length:
            seq_end = self.sequence_info.sequence[self.length:]
            if not seq_end.upper().replace("A", ""):
                return len(seq_end)
        return 0

    @property
    def cdna_match_diff(self) -> str:
        """ Human readable """
        match_summary = ""
        if cdna_errors := self._validate_cdna_match():
            match_summary = ", ".join(cdna_errors)
        elif exons := self.genome_build_data.get("exons"):
            gap_operations = Counter()
            for _, _, _, _, _, gap in exons:
                if gap:
                    for gap_op in gap.split():
                        code = gap_op[0]
                        length = int(gap_op[1:])
                        gap_operations[code] += length

            if gap_operations:
                gap_summary = []
                for code, label in {"I": "Insertion", "D": "Deletion"}.items():
                    if value := gap_operations.get(code):
                        gap_summary.append(f"{value}bp {label}")
                match_summary = ", ".join(gap_summary)
                if match_summary:
                    match_summary = f"Transcript had {match_summary} vs genome reference"

        return match_summary

    @property
    def sequence_length_matches_exon_length_ignoring_poly_a_tail(self) -> bool:
        # We can't know exactly how long a polyA tail is (only that subtracting it from length is all A's)
        return self.sequence_info.length == self.length or self.sequence_poly_a_tail

    @cached_property
    def alignment_gap(self) -> bool:
        if self.transcript.annotation_consortium == AnnotationConsortium.REFSEQ:
            # Sometimes RefSeq transcripts have gaps when aligning to the genome
            for ex in self.genome_build_data["exons"]:
                if ex[-1] is not None:
                    return True
            return not self.sequence_length_matches_exon_length_ignoring_poly_a_tail

        # Ensembl transcripts use genomic sequence so there is never any gap
        return False

    @staticmethod
    def get_transcript_id_and_version(transcript_accession: str) -> TranscriptParts:
        return get_transcript_id_and_version(transcript_accession)

    @staticmethod
    def transcript_versions_by_id(genome_build: GenomeBuild = None, annotation_consortium=None) -> \
            dict[str, dict[str, int]]:
        """ {transcript_id: {1: PK of TranscriptVersion.1, 2: PK of TranscriptVersion.2} """
        filter_kwargs = {}
        if genome_build:
            filter_kwargs["genome_build"] = genome_build
        if annotation_consortium:
            filter_kwargs["transcript__annotation_consortium"] = annotation_consortium

        qs = TranscriptVersion.objects.filter(**filter_kwargs)
        tv_by_id = defaultdict(dict)
        for pk, transcript_id, version in qs.values_list("pk", "transcript_id", "version"):
            tv_by_id[transcript_id][version] = pk
        return tv_by_id

    @staticmethod
    def id_by_accession(genome_build: GenomeBuild = None, annotation_consortium=None) -> dict[str, int]:
        filter_kwargs = {}
        if genome_build:
            filter_kwargs["genome_build"] = genome_build
        if annotation_consortium:
            filter_kwargs["transcript__annotation_consortium"] = annotation_consortium

        qs = TranscriptVersion.objects.filter(**filter_kwargs)
        tv_values = qs.values_list("pk", "transcript_id", "version")
        return {f"{transcript_id}.{version}": pk for (pk, transcript_id, version) in tv_values}

    @staticmethod
    def get_for_parts(genome_build: GenomeBuild, transcript_parts: TranscriptParts) -> Optional['TranscriptVersion']:
        """ Exact lookup for an already split accession - None if the version wasn't given, or we don't have it """
        if transcript_parts.identifier and transcript_parts.version:
            return TranscriptVersion.objects.filter(genome_build=genome_build,
                                                    transcript_id=transcript_parts.identifier,
                                                    version=transcript_parts.version).first()
        return None

    @staticmethod
    def filter_by_accession(accession, genome_build=None):
        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(accession)
        kwargs = {"transcript_id": transcript_id}
        if version is not None:
            kwargs["version"] = version
        if genome_build:
            kwargs["genome_build"] = genome_build
        return TranscriptVersion.objects.filter(**kwargs)

    @staticmethod
    def raise_bad_or_missing_transcript(transcript_accession):
        """ Checks whether a transcript we can't match is wrong (their fault) or we don't have it (our fault) """

        annotation_consortium = AnnotationConsortium.get_from_transcript_accession(transcript_accession).label
        key_base = f"transcript_exists:{transcript_accession}"
        transcript_connection_error_key = key_base + ":ERROR"
        if not cache.get(transcript_connection_error_key):
            bad_transcript_key = key_base + "BAD"
            if message := cache.get(bad_transcript_key):
                raise BadTranscript(message)

            no_transcript_key = key_base + "NO"
            if message := cache.get(no_transcript_key):
                raise NoTranscript(message)
            try:
                TranscriptVersionSequenceInfo.get(transcript_accession)  # Throws BadTranscript
                raise MissingTranscript(f"Transcript '{transcript_accession}' valid but missing from our database.")
            except BadTranscript as bt:
                # Only cache if we don't have it (DB will have it if we do)
                cache.set(bad_transcript_key, str(bt), timeout=WEEK_SECS)
                raise
            except NoTranscript as nt:
                # Only cache if we don't have it (DB will have it if we do)
                cache.set(no_transcript_key, str(nt), timeout=WEEK_SECS)
                raise
            except (RequestException, URLError):
                cache.set(transcript_connection_error_key, True, timeout=HOUR_SECS)

        raise NoTranscript(f"Transcript '{transcript_accession}' missing from our DB - validity with {annotation_consortium} unknown")

    @staticmethod
    def get_transcript_version(genome_build: GenomeBuild, transcript_name, best_attempt=True) -> Optional['TranscriptVersion']:
        """ @param best_attempt if we don't have the exact version specified in transcript_name, grab the closest
            one that is larger (or the largest one if we don't have any larger than the requested) """
        transcript_version = None
        transcript_versions_qs = TranscriptVersion.objects.filter(genome_build=genome_build)
        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_name)
        transcript_versions_qs = transcript_versions_qs.filter(transcript_id=transcript_id).order_by("version")
        if version is not None:
            try:
                transcript_version = transcript_versions_qs.get(version=version)
            except TranscriptVersion.DoesNotExist as exc:
                possible_versions = set(transcript_versions_qs.values_list('version', flat=True))
                possible_versions = [int(p) for p in possible_versions if p is not None]
                possible_versions.sort()
                if possible_versions:
                    if best_attempt:
                        use_version = None
                        for possible_ver in possible_versions:
                            use_version = possible_ver
                            if use_version > int(version):
                                break
                        transcript_version = transcript_versions_qs.filter(version=use_version).last()
                    else:
                        version_list = ', '.join(str(v) for v in possible_versions)
                        raise MissingTranscript(f"No Transcript for '{transcript_name}' (build: {genome_build}) - but there are entries for versions {version_list}") from exc
        else:
            transcript_version = transcript_versions_qs.last()

        if transcript_version is None:
            TranscriptVersion.raise_bad_or_missing_transcript(transcript_name)

        if not transcript_version.has_valid_data:
            # only going to happen if we have legacy data in the database, transcripts that use the default for data {}
            data_str = json.dumps(transcript_version.data)
            raise MissingTranscript(f"Transcript for '{transcript_name}' (build: {genome_build}),"
                                    f" but did not have complete data {data_str}")

        return transcript_version

    @staticmethod
    def get_for_lrg(genome_build: GenomeBuild, lrg_identifier: str) -> Optional['TranscriptVersion']:
        """ Attempts to load RefSeq TranscriptVersion we have from LRG identifier """
        transcript_version: Optional[TranscriptVersion] = None
        if lrg_identifier == "LRG_199t1":
            transcript_version = TranscriptVersion.get_transcript_version(genome_build, "NM_004006.2")
        return transcript_version

    @staticmethod
    def _sum_intervals(intervals: list[tuple]):
        return sum(b - a for a, b in intervals)

    @classmethod
    def data_is_current_cdot_format(cls, genome_build: Optional['GenomeBuild'] = None) -> bool:
        """ Current cdot nests per-build data under a 'genome_builds' key; pre-cdot imports stored a flat
            structure with no 'genome_builds', which makes genome_build_data (and thus HGVS resolution)
            raise KeyError on every record. Returns True once any TranscriptVersion is in the current
            format - ie import_cdot_latest / import_gene_annotation has been run. """
        qs = cls.objects.all()
        if genome_build:
            qs = qs.filter(genome_build=genome_build)
        return qs.filter(data__has_key="genome_builds").exists()

    @cached_property
    def genome_build_data(self) -> dict:
        return self.data["genome_builds"][self.genome_build.name]

    @cached_property
    def length(self) -> Optional[int]:
        exons = self.genome_build_data["exons"]
        strand = self.genome_build_data["strand"]
        if strand == '-':
            transcript_end_match = exons[0]
        else:
            transcript_end_match = exons[-1]
        return transcript_end_match[4]

    @cached_property
    def fivep_utr_length(self) -> int:
        return self._sum_intervals(self.fivep_utr)

    @cached_property
    def coding_length(self) -> int:
        return self._sum_intervals(self.cds)

    @cached_property
    def threep_utr_length(self) -> int:
        return self._sum_intervals(self.threep_utr)

    @property
    def num_codons(self) -> int:
        return self.coding_length // 3

    @property
    def protein_length(self) -> int:
        return self.num_codons - 1  # stop codon doesn't count

    def get_coordinates(self) -> Optional[GenomicCoordinates]:
        """ Where this transcript sits on the genome for its build, or None if the cdot data lacks it
            (older imports can be sparse - @see genome_build_data). Reuses the start/end/chrom/strand
            properties so there's one source of truth for reading the cdot coordinates. """
        build_data = self.data.get("genome_builds", {}).get(self.genome_build.name)
        if not (build_data and build_data.get("contig") is not None
                and build_data.get("exons") and build_data.get("strand") is not None):
            return None
        return GenomicCoordinates(contig=self.contig, chrom=self.chrom,
                                  start=self.start, end=self.end, strand=self.strand)

    @property
    def chrom(self):
        raw_contig = self.genome_build_data["contig"]
        return self.genome_build.chrom_contig_mappings[raw_contig].name

    @property
    def start(self) -> int:
        exons = self.genome_build_data["exons"]
        return exons[0][0]

    @property
    def end(self) -> int:
        exons = self.genome_build_data["exons"]
        return exons[-1][1]

    @property
    def strand(self) -> str:
        return self.genome_build_data["strand"]

    @property
    def coordinates(self):
        """ 1-based for humans """
        return f"{self.chrom}:{self.start + 1}-{self.end} ({self.strand})"

    @cached_property
    def tags(self) -> list[str]:
        """ 'tag' has been in cdot since 0.2.12. A record with no data for this build has no tags -
            that's a missing annotation rather than a broken record, so read it without raising
            (@see genome_build_data) """
        REMOVE_TAGS = {"basic"}  # This is on pretty much every Ensembl transcript
        tag_list = []
        # 'tag' was in the transcript in versions 0.2.12 - 0.2.13
        # It is inside genome build data after 0.2.14
        build_data = self.data.get("genome_builds", {}).get(self.genome_build.name) or {}
        if tag_list_str := build_data.get("tag") or self.data.get("tag"):
            tag_list = sorted(tag.strip() for tag in tag_list_str.split(",") if tag.strip() not in REMOVE_TAGS)
        return tag_list

    @cached_property
    def canonical_tag(self) -> str:
        """ We only want to return the most important one """
        tags = set(self.tags)
        for ct in self.CANONICAL_SCORES:
            if ct in tags:
                return ct
        return ""

    @property
    def canonical_score(self) -> int:
        """ This is so you can sort multiple transcripts and find 'most canonical' """
        return self.CANONICAL_SCORES.get(self.canonical_tag, 0)

    @property
    def is_canonical(self) -> bool:
        return bool(self.canonical_score)

    def get_contigs(self) -> set[Contig]:
        contigs = {self.genome_build_data["contig"]}
        if other_contigs := self.genome_build_data.get("other_chroms"):
            contigs.update(other_contigs)
        return {self.genome_build.chrom_contig_mappings[c] for c in contigs}

    def get_chromosomes(self) -> set[str]:
        return {c.name for c in self.get_contigs()}

    def _validate_cdna_match(self) -> list[str]:
        cdna_match_errors = []
        if exons := self.genome_build_data.get('exons'):
            # cdna_match = (genomic start, genomic end, cDNA start, cDNA end, gap) (genomic=0 based, transcript=1)
            if self.genome_build_data["strand"] == '-':
                exons = list(reversed(exons))

            last_end = None
            for _, _, exon_id, cdna_start, cdna_end, _ in exons:
                if exon_id == 0:
                    if cdna_start != 1:
                        cdna_match_errors.append(f"cDNA match starts at {cdna_start} not 1")

                if last_end:
                    missing = cdna_start - (last_end + 1)
                    if missing:
                        msg = f"cDNA match missing transcript: cDNA start: {cdna_start} last cDNA end {last_end}" \
                              f" (missing {missing} bp)"
                        cdna_match_errors.append(msg)
                last_end = cdna_end
        return cdna_match_errors

    def hgvs_data_errors(self) -> dict[str, str]:
        # Legacy pre-cdot rows have no 'genome_builds' at all - they report as missing every field
        build_data = self.data.get("genome_builds", {}).get(self.genome_build.name, {})
        data_errors = {}
        for key in ["contig", "strand", "exons"]:
            if key not in build_data:
                data_errors[key] = "Field missing"

        if error := (self.data.get("error") or build_data.get("error")):
            data_errors["error"] = error

        if "exons" not in data_errors:
            if cdna_errors := self._validate_cdna_match():
                data_errors["cdna_match"] = ", ".join(cdna_errors)

        return data_errors

    @property
    def has_valid_data(self) -> bool:
        return not self.hgvs_data_errors()

    @property
    def hgvs_error_tooltip(self) -> str:
        field_errors = []
        for k, v in self.hgvs_data_errors().items():
            field_errors.append(f"{k}: {v}")
        return ", ".join(field_errors)

    def get_differences(self, transcript_version):
        """ Used to inform while HGVS may resolve differently """
        differences = {}
        FIELDS = ["transcript_id", "version", "length"]
        for f in FIELDS:
            mine = getattr(self, f)
            other = getattr(transcript_version, f)
            if mine != other:
                differences[f] = (mine, other)

        if self.data and transcript_version.data:
            my_chrom = self.genome_build_data["contig"]
            other_chrom = transcript_version.genome_build_data["contig"]
            if my_chrom != other_chrom:
                try:
                    # Could be different but map to the same thing - try resolving it to contig name
                    other_cleaned_chrom = transcript_version.genome_build.chrom_contig_mappings[other_chrom].name
                    if self.chrom != other_cleaned_chrom:
                        differences["contig"] = (f"{my_chrom} (contig name: {self.chrom})",
                                                 f"{other_chrom} (contig name: {other_cleaned_chrom})")
                except Exception:
                    # Can't convert - just show differences
                    differences["contig"] = (my_chrom, other_chrom)

            my_exon_count = len(self.genome_build_data["exons"])
            other_exon_count = len(transcript_version.genome_build_data["exons"])
            if my_exon_count != other_exon_count:
                differences["exon count"] = (my_exon_count, other_exon_count)
            else:
                exon = 1
                my_exons = self.genome_build_data["exons"]
                other_exons = transcript_version.genome_build_data["exons"]
                if self.genome_build_data["strand"] == "-":
                    my_exons = reversed(my_exons)
                    other_exons = reversed(other_exons)

                for my_exon, other_exon in zip(my_exons, other_exons):
                    my_len = my_exon[1] - my_exon[0]
                    other_len = other_exon[1] - other_exon[0]
                    if my_len != other_len:
                        differences[f"exon {exon}"] = (my_len, other_len)
                    exon += 1
        else:
            my_keys = set(self.genome_build_data.keys())
            other_keys = set(transcript_version.genome_build_data.keys())
            if my_keys ^ other_keys:
                differences["data"] = (my_keys, other_keys)

        return differences

    @staticmethod
    def get_preferred_transcript(data: dict[str, 'TranscriptVersion']) -> Optional['TranscriptVersion']:
        for transcript_key in settings.VARIANT_ANNOTATION_TRANSCRIPT_PREFERENCES:
            if transcript := data.get(transcript_key):
                return transcript
        return None

    @cached_property
    def protein_domains_and_accession(self) -> tuple[QuerySet, str]:
        """ Gets custom ProteinDomain if available, falling back on Pfam """
        PD_ARGS = ("protein_domain__name", "protein_domain__description", "start", "end")
        protein_domains = self.proteindomaintranscriptversion_set.all().order_by("start").values_list(*PD_ARGS)
        if protein_domains.exists():
            used_transcript_version = self.accession
        else:
            pfam_qs = self.transcript.pfamsequenceidentifier_set.all()
            pfam = pfam_qs.filter(version=self.version).first()  # Try our version
            if not pfam:
                pfam = pfam_qs.order_by("version").last()

            if pfam:
                PFAM_ARGS = ("pfam__pfam_id", "pfam__description", "start", "end")
                protein_domains = pfam.pfam_sequence.pfamdomains_set.order_by("start").values_list(*PFAM_ARGS)
                used_transcript_version = pfam.accession
            else:
                used_transcript_version = None
        return protein_domains, used_transcript_version

    def get_absolute_url(self):
        kwargs = {"transcript_id": self.transcript_id, "version": self.version}
        return reverse("view_transcript_version", kwargs=kwargs)

    def get_external_url(self):
        return self.transcript.get_external_url(self.genome_build) + f".{self.version}"

    @staticmethod
    def update_accessions(accessions: Iterable[str], genome_build=None, **update_kwargs):
        """ A way to quickly update lots of accessions, relying on the fact that there aren't that many versions """

        tv_by_version = defaultdict(set)
        for transcript_accession in accessions:
            transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_accession)
            tv_by_version[version].add(transcript_id)

        for version, transcripts in tv_by_version.items():
            filter_kwargs = {"version": version, "transcript_id__in": transcripts}
            if genome_build:
                filter_kwargs["genome_build"] = genome_build
            TranscriptVersion.objects.filter(**filter_kwargs).update(**update_kwargs)

    def __str__(self):
        return f"{self.accession} ({self.gene_version.gene_symbol}/{self.genome_build_id})"


class TranscriptVersionSequenceInfoFastaFileImport(TimeStampedModel):
    sha256_hash = models.TextField(unique=True)
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)
    filename = models.TextField()

    def __str__(self):
        basename = os.path.basename(self.filename)
        return f"{basename} ({self.get_annotation_consortium_display()})"


class TranscriptVersionSequenceInfo(TimeStampedModel):
    """ Current main use of this is to download transcript version lengths from the web
        and check if TranscriptVersion exons sum to the same length

        Lengths not matching means there is a gap, but we can't be sure there isn't a gap """
    transcript = models.ForeignKey(Transcript, on_delete=CASCADE)
    version = models.IntegerField()
    # Data from Fasta file will have this set, from API will be populated with api_response
    fasta_import = models.ForeignKey(TranscriptVersionSequenceInfoFastaFileImport, null=True, on_delete=CASCADE)
    api_response = models.TextField(null=True)  # null if loaded from a file
    sequence = models.TextField()
    length = models.IntegerField()

    class Meta:
        unique_together = ("transcript", "version")

    def __str__(self):
        return f"{self.accession} ({self.length}bp)"

    @cached_property
    def accession(self):
        return TranscriptVersion.get_accession(self.transcript_id, self.version)

    def raise_any_errors(self):
        # Some will be inserted from files, and thus will have no API resposne
        if self.api_response:
            # Ensembl is JSON, while RefSeq use their own text based format
            if self.transcript.annotation_consortium == AnnotationConsortium.ENSEMBL:
                data = json.loads(self.api_response)
                if self.version != data["version"]:
                    raise NoTranscript(f"Only latest version: (v{data['version']}) can be retrieved via API")

    @staticmethod
    def get(transcript_accession: str, retrieve=True) -> Optional['TranscriptVersionSequenceInfo']:
        """ Returns DB copy if we have it, or retrieves + stores from API """

        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_accession)
        if tvi := TranscriptVersionSequenceInfo.objects.filter(transcript_id=transcript_id, version=version).first():
            tvi.raise_any_errors()
            return tvi

        if retrieve:
            fetched = TranscriptSequenceFetcher.instance().fetch(transcript_accession)
            tvi = TranscriptVersionSequenceInfo.store_fetched(fetched)
        return tvi

    @staticmethod
    def store_fetched(fetched: FetchedTranscriptSequence) -> 'TranscriptVersionSequenceInfo':
        Transcript.objects.get_or_create(pk=fetched.transcript_id,
                                         annotation_consortium=fetched.annotation_consortium)
        defaults = {"sequence": fetched.sequence, "length": fetched.length, "api_response": fetched.api_response}
        requested_tvi = None
        for version in fetched.versions:
            tvi, _ = TranscriptVersionSequenceInfo.objects.get_or_create(transcript_id=fetched.transcript_id,
                                                                         version=version, defaults=defaults)
            if version == fetched.version:
                requested_tvi = tvi
        requested_tvi.raise_any_errors()
        return requested_tvi

    @staticmethod
    def get_refseq_transcript_versions(transcript_accessions: Iterable[str], entrez_batch_size: int = 100, fail_on_error=True) -> dict[str, 'TranscriptVersionSequenceInfo']:
        """ Batch method - returns DB copies if we have it, retrieves + stores from API """
        # Find the ones we already have so we don't need to re-retrieve
        all_transcript_accessions = set(transcript_accessions)
        transcript_ids = [TranscriptVersion.get_transcript_id_and_version(a)[0] for a in transcript_accessions]
        tvi_by_id = {}
        for tvi in TranscriptVersionSequenceInfo.objects.filter(transcript_id__in=transcript_ids):
            if tvi.accession in all_transcript_accessions:
                tvi_by_id[tvi.accession] = tvi

        unknown_accessions = all_transcript_accessions - set(tvi_by_id)
        if unknown_accessions:
            fetched_list = TranscriptSequenceFetcher.instance().fetch_refseq_batch(
                unknown_accessions, entrez_batch_size=entrez_batch_size, fail_on_error=fail_on_error)
            tvi_by_id.update(TranscriptVersionSequenceInfo.bulk_store_fetched(fetched_list))
        return tvi_by_id

    @staticmethod
    def bulk_store_fetched(fetched_list: list[FetchedTranscriptSequence]) -> dict[str, 'TranscriptVersionSequenceInfo']:
        new_records = [TranscriptVersionSequenceInfo(transcript_id=f.transcript_id, version=f.version,
                                                     sequence=f.sequence, length=f.length,
                                                     api_response=f.api_response)
                       for f in fetched_list]
        # Write them as we go so any failure only loses some
        tvi_by_id = {}
        if new_records:
            # Create Transcript objects in case they don't exist
            transcript_records = {Transcript(pk=f.transcript_id, annotation_consortium=f.annotation_consortium)
                                  for f in fetched_list}
            Transcript.objects.bulk_create(transcript_records, ignore_conflicts=True, batch_size=2000)
            TranscriptVersionSequenceInfo.objects.bulk_create(new_records, ignore_conflicts=True, batch_size=2000)
            for tvi in new_records:
                tvi_by_id[tvi.accession] = tvi
        return tvi_by_id


class LRGRefSeqGene(models.Model):
    """ A Locus Reference Genomic (LRG) is a manually curated record that contains stable and thus, un-versioned
        reference sequences designed specifically for reporting sequence variants with clinical implications.

        These map to RefSeq, see https://www.ncbi.nlm.nih.gov/refseq/rsg/lrg/

        We don't link directly to TranscriptVersion as they are per-build """
    cached_web_resource = models.ForeignKey('annotation.CachedWebResource', on_delete=CASCADE)
    lrg = models.TextField()
    rna = models.TextField()  # RefSeq transcript accession
    t = models.TextField(null=True)
    category = models.CharField(max_length=1, choices=GeneSymbolAliasSource.choices)

    class Meta:
        # LRG_994 has multiple entries for LRG/t as there are multiple protein versions for t1 and t2
        # But we don't need protein so will ignore that. RNA/t is unique per LRG
        unique_together = ("lrg", "t")

    @staticmethod
    def get_lrg_and_t(lrg_identifier: str) -> tuple[Optional[str], Optional[str]]:
        return get_lrg_and_t(lrg_identifier)

    @staticmethod
    def get_transcript_version(genome_build: GenomeBuild, lrg_identifier: str) -> Optional[TranscriptVersion]:
        """ Attempts to load RefSeq TranscriptVersion we have from LRG identifier """
        transcript_version: Optional[TranscriptVersion] = None
        lrg, t = LRGRefSeqGene.get_lrg_and_t(lrg_identifier)
        if lrg:
            if t is None:
                raise ValueError(f"We require a t version with LRG '{lrg_identifier}'")
            if lrg_ref_seq_gene := LRGRefSeqGene.objects.filter(lrg=lrg, t=t).first():
                transcript_version = TranscriptVersion.get_transcript_version(genome_build, lrg_ref_seq_gene.rna)
        return transcript_version


class GeneSymbolWiki(Wiki):
    gene_symbol = models.OneToOneField(GeneSymbol, on_delete=CASCADE)

    def save(self, *args, **kwargs):
        super().save(*args, **kwargs)
        self.update_citations()

    def update_citations(self):
        existing_citations = {}
        for gc in self.gene_symbol.genesymbolcitation_set.all():
            existing_citations[gc.citation.pk] = gc

        gene_citation_ids_to_keep = []
        from annotation.models import CitationFetchRequest
        for citation in CitationFetchRequest.get_unfetched_citations(self.markdown):
            gene_citation = existing_citations.get(citation.pk)
            if not gene_citation:
                gene_citation = self.gene_symbol.genesymbolcitation_set.create(citation=citation)
            gene_citation_ids_to_keep.append(gene_citation.pk)

        # Delete any non-used ones
        num_deleted = self.gene_symbol.genesymbolcitation_set.exclude(pk__in=gene_citation_ids_to_keep).delete()[0]
        logging.info("Kept %s, deleted %d", gene_citation_ids_to_keep, num_deleted)


class ProteinDomain(models.Model):
    """ For custom domains, can be used to override PFam - used in TranscriptVersion.protein_domains_and_accession """
    name = models.TextField(unique=True)
    description = models.TextField()


class ProteinDomainTranscriptVersion(models.Model):
    """ No constraints as there can be multiple of same domains within a transcript """
    transcript_version = models.ForeignKey(TranscriptVersion, on_delete=CASCADE)
    protein_domain = models.ForeignKey(ProteinDomain, on_delete=CASCADE)
    start = models.IntegerField()
    end = models.IntegerField()


class Pfam(models.Model):
    PFAM_ACCESSION_PATTERN = re.compile(r"PF(\d{5})")
    # Use accession (with PF stripped off) as PK, ie "PF00017" => 17
    pfam_id = models.TextField(unique=True)
    description = models.TextField()

    @staticmethod
    def get_pk_from_accession(accession: str) -> int:
        if m := Pfam.PFAM_ACCESSION_PATTERN.match(accession):
            return int(m.group(1))
        raise ValueError(f"'{accession}' didn't match Pfam pattern")


class PfamSequence(models.Model):
    seq_id = models.TextField(primary_key=True)
    # Set once InterPro has been asked for this sequence's domains - a protein with no Pfam match is
    # a legitimate answer, so "no PfamDomains rows" can't be used to mean "not fetched yet"
    domains_imported = models.DateTimeField(null=True)


class PfamDomains(models.Model):
    """ No constraints as there can be multiple of same domains within a sequence """
    pfam_sequence = models.ForeignKey(PfamSequence, on_delete=CASCADE)
    pfam = models.ForeignKey(Pfam, on_delete=CASCADE)
    start = models.IntegerField()
    end = models.IntegerField()


class PfamSequenceIdentifier(models.Model):
    """ From HUMAN_9606_idmapping - Used to map Transcripts to PFam sequences """
    pfam_sequence = models.ForeignKey(PfamSequence, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, on_delete=CASCADE)
    # PFam provides transcript versions, but is build independent, while our transcript versions have builds
    # So we'll just store the version number
    version = models.IntegerField(null=True)

    @property
    def accession(self) -> str:
        acc = f"{self.transcript_id}"
        if self.version:
            acc += f".{self.version}"
        return acc

    def __str__(self):
        return f"{self.pfam_sequence}: {self.transcript}.{self.version}"


class MANE(models.Model):
    """ Matched Annotation from NCBI and EMBL-EBI (MANE)
        @see https://www.ncbi.nlm.nih.gov/refseq/MANE/ """
    ncbi_gene_version = models.ForeignKey(GeneVersion, related_name="mane_ncbi", null=True, on_delete=CASCADE)
    ensembl_gene_version = models.ForeignKey(GeneVersion, related_name="mane_ensembl", null=True, on_delete=CASCADE)
    hgnc = models.ForeignKey(HGNC, null=True, on_delete=CASCADE)
    symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    refseq_transcript_version = models.ForeignKey(TranscriptVersion, related_name="mane_refseq",
                                                  null=True, on_delete=CASCADE)
    ensembl_transcript_version = models.ForeignKey(TranscriptVersion, related_name="mane_ensembl",
                                                   null=True, on_delete=CASCADE)
    status = models.CharField(max_length=1, choices=MANEStatus.choices)

    @cache_memoize(DAY_SECS)
    @staticmethod
    def has_mane_transcripts() -> bool:
        return MANE.objects.exists()

    @staticmethod
    def get_mane_and_aliases_list_from_symbol(gene_symbol_str: str) -> list[tuple['MANE', Optional[GeneSymbolAlias]]]:
        if not MANE.has_mane_transcripts():
            raise ValueError("MANE transcripts are not loaded")

        mane_and_aliases = []

        for mane in MANE.objects.filter(symbol=gene_symbol_str):
            mane_and_aliases.append((mane, None))

        for alias in GeneSymbolAlias.objects.filter(alias=gene_symbol_str, gene_symbol__mane__isnull=False):
            if mane := MANE.objects.filter(symbol=alias.gene_symbol).first():
                mane_and_aliases.append((mane, alias))

        return mane_and_aliases

    def get_transcript_version(self, annotation_consortium: AnnotationConsortium) -> TranscriptVersion:
        transcript_version = None
        if annotation_consortium == AnnotationConsortium.REFSEQ:
            transcript_version = self.refseq_transcript_version
        elif annotation_consortium == AnnotationConsortium.ENSEMBL:
            transcript_version = self.ensembl_transcript_version
        return transcript_version

import json
import logging
import os
import re
import shutil
import types
from collections import defaultdict, Counter
from dataclasses import dataclass
from datetime import timedelta
from functools import cached_property, total_ordering
from io import StringIO
from typing import Optional, Union, Iterable, Any
from urllib.error import URLError, HTTPError

import requests
from Bio import Entrez, SeqIO
from cache_memoize import cache_memoize
from cdot.pyhgvs.pyhgvs_transcript import PyHGVSTranscriptFactory
from django.conf import settings
from django.contrib.auth.models import User, Group
from django.contrib.postgres.aggregates import StringAgg
from django.core.cache import cache
from django.core.exceptions import PermissionDenied, ObjectDoesNotExist, MultipleObjectsReturned
from django.db import models, IntegrityError, transaction
from django.db.models import QuerySet, TextField
from django.db.models.deletion import CASCADE, SET_NULL, PROTECT
from django.db.models.functions import Upper, Collate
from django.db.models.query_utils import Q
from django.db.models.signals import post_save, pre_delete
from django.dispatch import receiver
from django.shortcuts import get_object_or_404
from django.urls.base import reverse
from django.utils import timezone
from django.utils.timezone import localtime
from django_extensions.db.models import TimeStampedModel
from guardian.shortcuts import get_objects_for_user
from requests import RequestException

from genes.gene_coverage import load_gene_coverage_df
from genes.models_enums import AnnotationConsortium, HGNCStatus, GeneSymbolAliasSource, MANEStatus
from library.cache import timed_cache
from library.constants import HOUR_SECS, WEEK_SECS, MINUTE_SECS, DAY_SECS
from library.django_utils import SortByPKMixin
from library.django_utils.django_object_managers import ObjectManagerCachingRequest
from library.django_utils.django_partition import RelatedModelsPartitionModel
from library.guardian_utils import assign_permission_to_user_and_groups, DjangoPermission, admin_bot, \
    add_public_group_read_permission
from library.log_utils import log_traceback
from library.preview_request import PreviewData, PreviewModelMixin
from library.utils import get_single_element, iter_fixed_chunks, FormerTuple
from library.utils.file_utils import mk_path
from snpdb.models import Wiki, Company, Sample, DataState
from snpdb.models.models_enums import ImportStatus
from snpdb.models.models_genome import GenomeBuild, Contig
from upload.vcf.sql_copy_files import write_sql_copy_csv, gene_coverage_canonical_transcript_sql_copy_csv, \
    gene_coverage_sql_copy_csv, GENE_COVERAGE_HEADER


class HGNCImport(TimeStampedModel):
    pass


class NoTranscript(ValueError):
    """
    Extends ValueError for backwards compatibility.
    Indicates the transcript we are looking for is not in our database
    """

class NoTranscriptVersion(NoTranscript):
    """ """


class MissingTranscript(NoTranscript):
    """
    Transcript exists in RefSeq/Ensembl, so c.hgvs (or otherwise) might be okay.
    """


class BadTranscript(NoTranscript):
    """
    Transcript not found in Ensembl or RefSeq (User error)
    """


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
        return list(sorted(gene_symbol_strs))

    @cached_property
    def alias_symbols_in_db(self) -> list[GeneSymbolAlias]:
        return list(sorted([alias for alias in self.alias_list if not alias.different_genes and alias.other_symbol_in_database]))

    @cached_property
    def aliases_out(self) -> list[GeneSymbolAliasSummary]:
        return list(sorted([alias for alias in self.alias_list if not alias.my_symbol_is_main]))

    @cached_property
    def aliases_in(self) -> list[GeneSymbolAliasSummary]:
        return list(sorted([alias for alias in self.alias_list if alias.my_symbol_is_main]))


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
            print(f"Deleted orphaned {Gene.FAKE_GENE_ID_PREFIX} records:")
            print(ret)

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
        except:
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
            "end": min(all_transcripts_data["end"]),
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


@dataclass
class TranscriptParts(FormerTuple):
    identifier: str
    version: Optional[int]

    @property
    def as_tuple(self) -> tuple:
        return self.identifier, self.version

    def __repr__(self):
        if self.version:
            return f"{self.identifier}.{self.version}"
        return self.identifier


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
        """ """
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
            # We've modified PyHGVS to be able to handle this
            for ex in self.genome_build_data["exons"]:
                if ex[-1] is not None:
                    return True
            return not self.sequence_length_matches_exon_length_ignoring_poly_a_tail

        # Ensembl transcripts use genomic sequence so there is never any gap
        return False

    @staticmethod
    def get_transcript_id_and_version(transcript_accession: str) -> TranscriptParts:
        parts = transcript_accession.split(".")
        if len(parts) == 2:
            identifier = str(parts[0])
            version = int(parts[1])
        else:
            identifier, version = transcript_accession, None
        return TranscriptParts(identifier, version)

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
            except TranscriptVersion.DoesNotExist:
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
                        version_list = ', '.join((str(v) for v in possible_versions))
                        raise MissingTranscript(f"No Transcript for '{transcript_name}' (build: {genome_build}) - but there are entries for versions {version_list}")
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

    @cached_property
    def genome_build_data(self) -> dict:
        return self.data["genome_builds"][self.genome_build.name]

    @cached_property
    def pyhgvs_data(self):
        transcript_json = self.data.copy()
        # Legacy data stored gene_name in JSON, but that could lead to diverging values vs TranscriptVersion relations
        # so use DB as source of truth and replace into PyHGVS at last minute
        if self.gene_symbol:
            transcript_json["gene_name"] = str(self.gene_symbol)
        tf = PyHGVSTranscriptFactory(transcripts={self.accession: transcript_json})
        return tf.get_pyhgvs_data(self.accession, self.genome_build.name, sacgf_pyhgvs_fork=True)

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
        """ 'tag' has been in cdot since 0.2.12 """
        REMOVE_TAGS = {"basic"}  # This is on pretty much every Ensembl transcript
        tag_list = []
        # 'tag' was in the transcript in versions 0.2.12 - 0.2.13
        # It is inside genome build data after 0.2.14
        if tag_list_str := self.genome_build_data.get("tag") or self.data.get("tag"):
            tag_list = sorted(tag for tag in tag_list_str.split(",") if tag not in REMOVE_TAGS)
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
        data_errors = {}
        for key in ["contig", "strand", "exons"]:
            try:
                _ = self.genome_build_data[key]
            except KeyError:
                data_errors[key] = f"Field missing"

        if error := (self.data.get("error") or self.genome_build_data.get("error")):
            data_errors["error"] = error

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
                except:
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
                    raise NoTranscript(f"Only latest version: (v{data['version']}) can be retrieved via API - asked for TVSI({self.pk}) {self.transcript.identifier} v{self.version}")

    @staticmethod
    def get(transcript_accession: str, retrieve=True) -> Optional['TranscriptVersionSequenceInfo']:
        """ Returns DB copy if we have it, or retrieves + stores from API """

        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_accession)
        if tvi := TranscriptVersionSequenceInfo.objects.filter(transcript_id=transcript_id, version=version).first():
            tvi.raise_any_errors()
            return tvi

        if retrieve:
            annotation_consortium = AnnotationConsortium.get_from_transcript_accession(transcript_accession)
            if annotation_consortium == AnnotationConsortium.REFSEQ:
                tvi = TranscriptVersionSequenceInfo._get_and_store_from_refseq_api(transcript_accession)
            else:
                tvi = TranscriptVersionSequenceInfo._get_and_store_from_ensembl_api(transcript_accession)
        return tvi

    @staticmethod
    def _get_kwargs_from_genbank_record(record):
        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(record.id)
        sequence = str(record.seq)
        length = len(sequence)
        return {"transcript_id": transcript_id,
                "version": version,
                "sequence": sequence,
                "length": length}

    @staticmethod
    def _get_and_store_from_refseq_api(transcript_accession):
        try:
            data = Entrez.efetch(db='nuccore', id=transcript_accession, rettype='gb', retmode='text')
        except HTTPError as e:
            if e.code == 400:
                raise BadTranscript(f"Bad Transcript: Entrez API reports \"{transcript_accession}\" not found")
            raise e
        api_response = data.read()
        with StringIO(api_response) as f:
            records = list(SeqIO.parse(f, "genbank"))
            record = get_single_element(records)
            kwargs = TranscriptVersionSequenceInfo._get_kwargs_from_genbank_record(record)
            transcript_id = kwargs.pop("transcript_id")
            version = kwargs.pop("version")
            Transcript.objects.get_or_create(pk=transcript_id,
                                             annotation_consortium=AnnotationConsortium.REFSEQ)
            defaults = kwargs
            defaults["api_response"] = api_response
            return TranscriptVersionSequenceInfo.objects.get_or_create(transcript_id=transcript_id, version=version,
                                                                       defaults=defaults)[0]

    @staticmethod
    def _get_and_store_from_ensembl_api(transcript_accession, genome_build: GenomeBuild = None):
        if genome_build is None:
            genome_build = GenomeBuild.grch38()
        ENSEMBL_REST_BASE_URLS = {
            "GRCh37": "https://grch37.rest.ensembl.org",
            "GRCh38": "https://rest.ensembl.org",
        }
        base_url = ENSEMBL_REST_BASE_URLS[genome_build.name]
        transcript_id, requested_version = TranscriptVersion.get_transcript_id_and_version(transcript_accession)
        url = f"{base_url}/sequence/id/{transcript_id}?type=cdna"
        r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=MINUTE_SECS)
        data = r.json()

        if r.ok:
            transcript, _ = Transcript.objects.get_or_create(identifier=data["id"],
                                                             annotation_consortium=AnnotationConsortium.ENSEMBL)

            # Ensembl only returns the latest version via API - so requested/returned version may be different
            returned_version = data["version"]
            kwargs = {
                "transcript": transcript,
                "version": requested_version,
                "defaults": {
                    "api_response": r.text,
                    "sequence": data["seq"],
                    "length": len(data["seq"])
                },
            }
            requested_tvsi = TranscriptVersionSequenceInfo.objects.get_or_create(**kwargs)[0]
            if requested_version != returned_version:
                # Store what was actually returned
                kwargs["version"] = returned_version
                TranscriptVersionSequenceInfo.objects.get_or_create(**kwargs)

            requested_tvsi.raise_any_errors()
        else:
            error = data.get("error")
            if error:
                if "not found" in error:
                    if genome_build != GenomeBuild.grch37():
                        # Try 37 this time
                        return TranscriptVersionSequenceInfo._get_and_store_from_ensembl_api(transcript_accession,
                                                                                             GenomeBuild.grch37())
                    raise BadTranscript(f"Ensembl API reports '{transcript_id}' not found")
            raise NoTranscript(f"Unable to understand Ensembl API response: {data}")

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

        for id_list in iter_fixed_chunks(unknown_accessions, entrez_batch_size):
            id_param = ",".join(id_list)
            try:
                search_results = Entrez.read(Entrez.epost("nuccore", id=id_param))
                fetch_handle = Entrez.efetch(
                    db="nuccore",
                    rettype="gb",
                    retmode="text",
                    webenv=search_results["WebEnv"],
                    query_key=search_results["QueryKey"],
                    idtype="acc",
                )
                tvi_by_id.update(TranscriptVersionSequenceInfo._insert_from_genbank_handle(fetch_handle))
            except RuntimeError as e:
                print(f"Entrez failed w/params: {id_param}")
                if fail_on_error:
                    raise e

        return tvi_by_id

    @staticmethod
    def _insert_from_genbank_handle(handle) -> dict[str, 'TranscriptVersionSequenceInfo']:
        new_records = []
        for record in SeqIO.parse(handle, "genbank"):
            # Store raw data so that we can retrieve more stuff from it later
            s = StringIO()
            SeqIO.write(record, s, "genbank")
            s.seek(0)
            api_response = s.read()
            kwargs = TranscriptVersionSequenceInfo._get_kwargs_from_genbank_record(record)
            tvi = TranscriptVersionSequenceInfo(**kwargs, api_response=api_response)
            new_records.append(tvi)

        # Write them as we go so any failure only loses some
        tvi_by_id = {}
        if new_records:
            # Create Transcript objects in case they don't exist
            transcript_records = {Transcript(pk=tvsi.transcript_id, annotation_consortium=AnnotationConsortium.REFSEQ)
                                  for tvsi in new_records}
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
    def get_lrg_and_t(lrg_identifier: str) -> tuple[str, Optional[str]]:
        if m := re.match(r"(LRG_\d+)(t\d+)?", lrg_identifier, flags=re.IGNORECASE):
            lrg, t = m.groups()
            lrg = lrg.upper()
            if t:
                t = t.lower()
        else:
            lrg = t = None
        return lrg, t

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


class GeneAnnotationRelease(models.Model):
    """
        import_gene_annotation management command now supports --release which is used to link all
        GeneVersion/TranscriptVersion from a particular GTF file.

        This release can be set on a VariantAnnotationVersion to be able to get genes/transcripts from a VEP build
    """
    version = models.TextField()  # Needs to support e.g. "109.20190607"
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    gene_annotation_import = models.ForeignKey(GeneAnnotationImport, on_delete=CASCADE)

    class Meta:
        unique_together = ('version', 'annotation_consortium', 'genome_build')

    def get_genes(self):
        return Gene.objects.filter(annotation_consortium=self.annotation_consortium,
                                   geneversion__releasegeneversion__release=self)

    @staticmethod
    def get_for_latest_annotation_versions_for_builds() -> list['GeneAnnotationRelease']:
        from annotation.models import VariantAnnotationVersion
        gene_annotation_releases = []
        for genome_build in GenomeBuild.builds_with_annotation().order_by("name"):
            if vav := VariantAnnotationVersion.latest(genome_build):
                if vav.gene_annotation_release:
                    gene_annotation_releases.append(vav.gene_annotation_release)
        return gene_annotation_releases

    def genes_for_symbols(self, gene_symbols) -> QuerySet:
        rgsg_qs = ReleaseGeneSymbolGene.objects.filter(release_gene_symbol__release=self,
                                                       release_gene_symbol__gene_symbol__in=gene_symbols)
        return Gene.objects.filter(pk__in=rgsg_qs.values_list("gene_id", flat=True))

    def genes_for_symbol(self, gene_symbol) -> QuerySet:
        return self.genes_for_symbols([gene_symbol])

    def transcript_versions_for_transcript(self, transcript) -> QuerySet[TranscriptVersion]:
        return TranscriptVersion.objects.filter(releasetranscriptversion__release=self,
                                                transcript=transcript)

    def transcript_versions_for_gene(self, gene) -> QuerySet:
        return TranscriptVersion.objects.filter(releasetranscriptversion__release=self,
                                                gene_version__gene=gene)

    def transcript_versions_for_symbol(self, gene_symbol) -> QuerySet:
        return TranscriptVersion.objects.filter(releasetranscriptversion__release=self,
                                                gene_version__gene__in=self.genes_for_symbol(gene_symbol))

    def __str__(self):
        return f"{self.genome_build_id}/{self.get_annotation_consortium_display()} - v{self.version}"


class ReleaseGeneVersion(models.Model):
    """ Used to be able to reconstruct what versions were in an import """
    release = models.ForeignKey(GeneAnnotationRelease, on_delete=CASCADE)
    gene_version = models.ForeignKey(GeneVersion, on_delete=CASCADE)

    class Meta:
        unique_together = ('release', 'gene_version')


class ReleaseTranscriptVersion(models.Model):
    """ Used to be able to reconstruct what versions were in an import """
    release = models.ForeignKey(GeneAnnotationRelease, on_delete=CASCADE)
    transcript_version = models.ForeignKey(TranscriptVersion, on_delete=CASCADE)

    class Meta:
        unique_together = ('release', 'transcript_version')


class ReleaseGeneSymbol(TimeStampedModel):
    """ Collects (via GeneSymbolReleaseGene) how gene symbol matches genes for a particular release
        Works as global cache, created as symbols appear in gene lists """
    release = models.ForeignKey(GeneAnnotationRelease, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)

    class Meta:
        unique_together = ('release', 'gene_symbol')

    def __str__(self):
        return f"{self.release}: {self.gene_symbol}"


class ReleaseGeneSymbolGene(TimeStampedModel):
    release_gene_symbol = models.ForeignKey(ReleaseGeneSymbol, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    match_info = models.TextField(null=True)

    class Meta:
        unique_together = ('release_gene_symbol', 'gene')


class GeneSymbolWiki(Wiki):
    gene_symbol = models.OneToOneField(GeneSymbol, on_delete=CASCADE)

    def save(self, **kwargs):
        super().save(**kwargs)
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


class GeneListCategory(models.Model):
    # This was a separate table rather than hard coded choices so that we can use dependent fields / chaining forms.
    # Used as CSS classes so have to have no spaces
    NODE_CUSTOM_TEXT = 'NodeCustomText'
    SAMPLE_GENE_LIST = 'SampleGeneList'
    QC_COVERAGE_CUSTOM_TEXT = 'QCCoverageCustomText'
    PATHOLOGY_TEST = 'PathologyTest'
    PATHOLOGY_TEST_ORDER = 'PathologyTestOrder'
    PANEL_APP_CACHE = 'PanelAppCache'
    GENE_INFO = "GeneInfo"

    name = models.TextField()
    company = models.OneToOneField(Company, null=True, blank=True, on_delete=CASCADE)
    icon_css_class = models.TextField(blank=True)
    hidden = models.BooleanField(default=False)  # for special ones
    public = models.BooleanField(default=False)
    description = models.TextField()

    @staticmethod
    def get_or_create_category(category_name, hidden=False):
        category, created = GeneListCategory.objects.get_or_create(name=category_name)
        if created:
            category.hidden = hidden
            category.save()
        return category

    @staticmethod
    def _get_pathology_test_gene_category(category_name, set_company=False):
        """ Requires settings.COMPANY to be set """
        category = None
        company = Company.get_our_company()
        if company:
            category, created = GeneListCategory.objects.get_or_create(name=category_name)
            if created:
                category.hidden = True
                if set_company:  # 1-to-1 field so can't always do it
                    category.company = company
                category.icon_css_class = "company-" + company.name.lower() + "-icon"
                category.save()
        return category

    @staticmethod
    def get_pathology_test_gene_category():
        return GeneListCategory._get_pathology_test_gene_category(GeneListCategory.PATHOLOGY_TEST, set_company=True)

    @staticmethod
    def get_pathology_test_order_gene_category():
        return GeneListCategory._get_pathology_test_gene_category(GeneListCategory.PATHOLOGY_TEST_ORDER)

    @staticmethod
    def get_gene_list_categories(extra_filters=None):
        categories = [{"name": "User", "pk": None}]  # Blank / no category
        qs = GeneListCategory.objects.exclude(genelist__isnull=True).exclude(hidden=True)
        if extra_filters:
            qs = qs.filter(extra_filters)
        for glc in qs.order_by("name"):
            categories.append({"name": glc.name,
                               "pk": glc.pk,
                               "icon_css_class": glc.icon_css_class,
                               "instance": glc})
        return categories

    def __str__(self):
        return self.name


class GeneList(TimeStampedModel):
    """ Stores a gene/transcript list (to be used as a filter) """

    category = models.ForeignKey(GeneListCategory, null=True, blank=True, on_delete=CASCADE)
    name = models.TextField()
    user = models.ForeignKey(User, on_delete=CASCADE)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    error_message = models.TextField(null=True, blank=True)
    locked = models.BooleanField(default=False)
    url = models.TextField(null=True, blank=True)

    def get_q(self, variant_annotation_version):
        """ For a Variant queryset """
        from annotation.models.models import VariantTranscriptAnnotation
        genes_qs = self.get_genes(variant_annotation_version.gene_annotation_release)
        return VariantTranscriptAnnotation.get_overlapping_genes_q(variant_annotation_version,
                                                                   genes_qs)

    @staticmethod
    def get_gene_ids_for_gene_lists(release: GeneAnnotationRelease, gene_lists: list['GeneList']):
        """ For GeneList node, we need to get query for multiple lists - and it's much faster to build the merged
            query here, than via joining separate queryset """
        rgs_qs = ReleaseGeneSymbol.objects.filter(release=release, gene_symbol__genelistgenesymbol__gene_list__in=gene_lists)
        return rgs_qs.values_list("releasegenesymbolgene__gene", flat=True)

    def get_genes(self, release: GeneAnnotationRelease):
        """ Get Genes (from a release) for symbols in this gene list """
        rgs_qs = ReleaseGeneSymbol.objects.filter(release=release, gene_symbol__genelistgenesymbol__gene_list=self)
        return Gene.objects.filter(releasegenesymbolgene__release_gene_symbol__in=rgs_qs)

    def get_gene_names(self):
        return self.genelistgenesymbol_set.filter(gene_symbol__isnull=False).values_list("gene_symbol", flat=True).distinct()

    def clone(self):
        """ Needed to clone analysis node GeneListNode """
        genelistgenesymbol_set = set(self.genelistgenesymbol_set.all())
        copy = self
        copy.pk = None
        copy.save()

        genes = 0
        for glgs in genelistgenesymbol_set:
            glgs.pk = None
            glgs.gene_list = copy
            glgs.save()
            genes += 1

        return copy

    def save(self, **kwargs):
        assign_permissions = kwargs.pop("assign_permissions", False)
        initial_save = not self.pk
        super().save(**kwargs)
        if initial_save or assign_permissions:
            # logging.info("GeneList: assign_permission_to_user_and_groups")
            assign_permission_to_user_and_groups(self.user, self)

    def can_view(self, user_or_group: Union[User, Group]) -> bool:
        read_perm = DjangoPermission.perm(self, DjangoPermission.READ)
        return user_or_group.has_perm(read_perm, self)

    def can_write(self, user_or_group: Union[User, Group]) -> bool:
        write_perm = DjangoPermission.perm(self, DjangoPermission.WRITE)
        return user_or_group.has_perm(write_perm, self) and not self.locked

    def get_warnings(self, release: GeneAnnotationRelease) -> list[str]:
        counts = {"unmatched symbols": self.unmatched_gene_symbols.count(),
                  "aliased": self.aliased_genes.count(),
                  "unmatched genes": self.unmatched_genes(release).count()}

        warnings = []
        for name, num in counts.items():
            if num:
                warnings.append(f"{name} x {num}")
        return warnings

    @property
    def unmatched_gene_symbols(self):
        return self.genelistgenesymbol_set.filter(gene_symbol__isnull=True).order_by("original_name")

    @property
    def aliased_genes(self):
        return self.genelistgenesymbol_set.filter(gene_symbol_alias__isnull=False)

    def unmatched_genes(self, release: GeneAnnotationRelease):
        string_agg = GeneListGeneSymbol.get_joined_genes_qs_annotation_for_release(release)
        return self.genelistgenesymbol_set.annotate(matched_genes=string_agg).filter(matched_genes__isnull=True)

    def add_and_remove_gene_symbols(self, gene_symbol_additions, gene_symbol_deletions,
                                    gene_additions_modification_info=None):
        """ Either works for all, or fails and makes no modification
            returns (num_added, num_deleted) """
        if gene_additions_modification_info is None:
            gene_additions_modification_info = {}

        if self.locked:
            msg = "Can't modify a locked GeneList!"
            raise PermissionDenied(msg)

        num_added = 0

        for gene_symbol in gene_symbol_additions:
            modification_info = gene_additions_modification_info.get(gene_symbol)
            defaults = {"original_name": gene_symbol, "modification_info": modification_info}
            _, created = GeneListGeneSymbol.objects.get_or_create(gene_list=self, gene_symbol_id=gene_symbol,
                                                                  defaults=defaults)
            if created:
                num_added += 1

        # Match original name and symbol so we can delete non-matched
        name_or_symbol_match = Q(original_name__in=gene_symbol_deletions) | Q(gene_symbol__in=gene_symbol_deletions)
        qs = self.genelistgenesymbol_set.filter(name_or_symbol_match)
        num_deleted = qs.delete()

        if num_added or num_deleted:
            self.set_modified_to_now()

        return num_added, num_deleted

    @staticmethod
    def get_for_user(user, gene_list_id, success_only=True):
        try:
            return GeneList.filter_for_user(user, success_only).get(pk=gene_list_id)
        except GeneList.DoesNotExist:
            # Need to distinguish between does not exist and no permission
            get_object_or_404(GeneList, pk=gene_list_id)  # potentially throws GeneList.DoesNotExist
            # If we're here, object exists but we have a permission error
            msg = f"You don't have permission to access gene_list {gene_list_id}"
            raise PermissionDenied(msg)

    @staticmethod
    def filter_for_user(user, success_only=True):
        perm = DjangoPermission.perm(GeneList, DjangoPermission.READ)
        user_qs = get_objects_for_user(user, perm, klass=GeneList, accept_global_perms=False)

        # Sample gene lists for samples we have permission to see
        samples_qs = Sample.filter_for_user(user)
        sample_gene_list_qs = GeneList.objects.filter(category__name=GeneListCategory.SAMPLE_GENE_LIST,
                                                      customtextgenelist__qcgenelist__qc__bam_file__unaligned_reads__sequencing_sample__samplefromsequencingsample__sample__in=samples_qs)

        qs = user_qs | sample_gene_list_qs
        if success_only:
            qs = qs.filter(import_status=ImportStatus.SUCCESS)
        return qs.distinct()  # Sometimes gene list can be reached multiple ways via 'or' above

    @staticmethod
    def _visible_gene_lists(gene_lists_qs):
        hidden_categories_qs = GeneListCategory.objects.filter(hidden=True)
        return gene_lists_qs.exclude(category__in=hidden_categories_qs)  # Allow through no category

    @classmethod
    def visible_gene_lists_containing_gene_symbol(cls, gene_lists_qs, gene_symbol):
        qs = cls._visible_gene_lists(gene_lists_qs)
        gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol)
        return qs.filter(genelistgenesymbol__gene_symbol=gene_symbol)

    def set_modified_to_now(self):
        GeneList.objects.filter(pk=self.pk).update(modified=timezone.now())

    def __str__(self):
        return f"{self.name} ({self.genelistgenesymbol_set.count()} x genes)"

    def get_absolute_url(self):
        url = self.url
        if url is None:
            url = reverse('view_gene_list', args=[str(self.id)])
        return url

    @classmethod
    def get_listing_url(cls):
        return reverse('gene_lists')


def create_fake_gene_list(*args, **kwargs):
    """  Originally FakeGeneList had abstract=True but with Django 3.2 got "Abstract models cannot be instantiated" """

    def get_absolute_url(self):
        return None

    gene_list = GeneList(*args, **kwargs)
    gene_list.get_absolute_url = types.MethodType(get_absolute_url, gene_list)
    return gene_list


class GeneListGeneSymbol(models.Model):
    gene_list = models.ForeignKey(GeneList, on_delete=CASCADE)
    original_name = models.TextField(null=True, blank=True)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    gene_symbol_alias = models.ForeignKey(GeneSymbolAlias, null=True, on_delete=CASCADE)
    modification_info = models.TextField(null=True)

    class Meta:
        unique_together = ('gene_list', 'original_name')

    @staticmethod
    def get_joined_genes_qs_annotation_for_release(release: GeneAnnotationRelease):
        """ Used to annotate GeneListGeneSymbol queryset """
        return StringAgg("gene_symbol__releasegenesymbol__releasegenesymbolgene__gene",
                         delimiter=',', distinct=True, output_field=TextField(),
                         filter=Q(gene_symbol__releasegenesymbol__release=release))

    def __str__(self):
        if self.gene_symbol:
            description = f"{self.gene_symbol_id}"
            if self.gene_symbol_alias:
                description += f" (via {self.gene_symbol_alias})"
        else:
            description = f"'{self.original_name}' (unrecognised symbol)"
        return description


class CustomTextGeneList(models.Model):
    """' Some human entered text which gets pulled apart to create a gene list """
    sha256_hash = models.TextField()
    name = models.TextField()
    text = models.TextField()
    gene_list = models.OneToOneField(GeneList, null=True, on_delete=SET_NULL)

    def clone(self):
        copy = self
        old_gene_list = self.gene_list
        if self.gene_list:
            copy.gene_list = self.gene_list.clone()
        copy.pk = None
        copy.save()
        self.gene_list = old_gene_list
        return copy

    def __str__(self):
        text = self.text
        if len(text) > 50:
            text = text[:50] + " ..."
        return f"CustomTextGeneList for '{text}': {self.gene_list}"


class SampleGeneList(TimeStampedModel):
    """ There can be multiple SampleGeneLists per sample, but only 1 active one.
        If multiple exist, the active one must be set manually """
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    gene_list = models.ForeignKey(GeneList, on_delete=CASCADE)
    visible = models.BooleanField(default=True, blank=False)

    class Meta:
        unique_together = ('sample', 'gene_list')
        ordering = ['created']

    def __str__(self):
        s = f"{self.sample}: {localtime(self.modified)}"
        if not self.visible:
            s += " (hidden)"
        return s


@receiver(post_save, sender=SampleGeneList)
def sample_gene_list_created(sender, instance, created, **kwargs):  # pylint: disable=unused-argument
    if created:
        sample = instance.sample
        if SampleGeneList.objects.filter(sample=sample).count() > 1:
            # Multiple exist, so need to set manually
            ActiveSampleGeneList.objects.filter(sample=sample).delete()
        else:
            try:
                with transaction.atomic():
                    # There can only be 1 - if this works it's active
                    ActiveSampleGeneList.objects.create(sample=sample, sample_gene_list=instance)
            except IntegrityError:
                ActiveSampleGeneList.objects.filter(sample=sample).delete()


class ActiveSampleGeneList(TimeStampedModel):
    """ Use 1-to-1 to enforce there's only 1 in DB
        (as compared to an "active" flag on SampleGeneList) """
    sample = models.OneToOneField(Sample, on_delete=CASCADE)
    sample_gene_list = models.ForeignKey(SampleGeneList, on_delete=CASCADE)


class GeneListWiki(Wiki):
    gene_list = models.OneToOneField(GeneList, on_delete=CASCADE)

    def _get_restricted_object(self):
        return self.gene_list


class PanelAppServer(models.Model):
    name = models.TextField(unique=True)
    url = models.TextField(unique=True)
    icon_css_class = models.TextField()

    def __str__(self):
        return self.name

    @staticmethod
    def australia_instance() -> 'PanelAppServer':
        return PanelAppServer.objects.get(name="PanelApp Australia")

    @staticmethod
    def england_instance() -> 'PanelAppServer':
        return PanelAppServer.objects.get(name="Genomics England PanelApp")


class PanelAppPanel(TimeStampedModel):
    """ Populated from PanelApp cached web resource task, updated to latest version
        It's not clear how PanelApp removes panels, so we don't do that. """
    server = models.ForeignKey(PanelAppServer, on_delete=CASCADE)
    panel_id = models.IntegerField()
    disease_group = models.TextField()
    disease_sub_group = models.TextField()
    name = models.TextField()
    status = models.TextField()
    current_version = models.TextField()

    class Meta:
        unique_together = ('server', 'panel_id')

    @property
    def url(self) -> str:
        return f"{self.server.url}/api/v1/panels/{self.panel_id}"

    @property
    def cache_valid(self) -> bool:
        # Attempt to use cache if recent and present, otherwise fall through and do a query
        max_age = timedelta(days=settings.PANEL_APP_CACHE_DAYS)
        return timezone.now() < self.modified + max_age

    def __str__(self):
        return self.name


class PanelAppPanelRelevantDisorders(models.Model):
    panel_app_panel = models.ForeignKey(PanelAppPanel, on_delete=CASCADE)
    name = models.TextField()

    class Meta:
        unique_together = ('panel_app_panel', 'name')


class PanelAppPanelLocalCache(TimeStampedModel):
    panel_app_panel = models.ForeignKey(PanelAppPanel, on_delete=CASCADE)
    version = models.TextField()

    class Meta:
        unique_together = ("panel_app_panel", "version")

    def get_gene_list(self, panel_app_confidence) -> GeneList:
        from genes.gene_matching import GeneSymbolMatcher

        min_level = int(panel_app_confidence)
        name = f"{self.panel_app_panel.name} v.{self.version}.min_{min_level}"

        # We'll try to re-use gene lists - but it's possible due to race conditions we may occasionally make
        # a duplicate, but this should work fine and won't affect much
        category = GeneListCategory.get_or_create_category(GeneListCategory.PANEL_APP_CACHE, hidden=True)
        gene_list_kwargs = {
            "category": category,
            "name": name,
            "user": admin_bot(),
            "import_status": ImportStatus.SUCCESS,
            "url": self.panel_app_panel.url,
        }
        if gene_list := GeneList.objects.filter(**gene_list_kwargs).order_by("pk").first():
            print(f"Reused existing gene list: {gene_list.pk}")
        else:
            gene_list = GeneList.objects.create(**gene_list_kwargs)
            print(f"Created gene list: {gene_list.pk}")
            gene_names_list = []
            for pap_lc_gs in self.panelapppanellocalcachegenesymbol_set.all():
                confidence_level = int(pap_lc_gs.data["confidence_level"])
                if confidence_level >= min_level:
                    gene_symbol = pap_lc_gs.data["gene_data"]["gene_symbol"]
                    gene_names_list.append(gene_symbol)

            print(f"Creating symbols: {gene_names_list}")
            gene_matcher = GeneSymbolMatcher()
            gene_matcher.create_gene_list_gene_symbols(gene_list, gene_names_list)
            gene_list.import_status = ImportStatus.SUCCESS
            gene_list.save()

            # PanelApp gene list should be public
            add_public_group_read_permission(gene_list)

        print(f"Returning gene list: {gene_list.pk}")
        return gene_list

    def __str__(self):
        return f"PanelApp cache for {self.panel_app_panel} v{self.version} (mod: {self.modified})"


class PanelAppPanelLocalCacheGeneSymbol(models.Model):
    panel_app_local_cache = models.ForeignKey(PanelAppPanelLocalCache, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    data = models.JSONField(null=False, blank=True, default=dict)  # API response


class CachedThirdPartyGeneList(models.Model):
    cached_web_resource = models.ForeignKey('annotation.CachedWebResource', on_delete=CASCADE)
    company = models.ForeignKey(Company, on_delete=CASCADE)


class GeneInfo(models.Model):
    """ These represent some kind of special tagging on genes """
    name = models.TextField(primary_key=True)
    description = models.TextField(blank=True)
    icon_css_class = models.TextField()
    gene_list = models.OneToOneField(GeneList, null=True, on_delete=PROTECT)

    @staticmethod
    def get_for_gene_symbol(gene_symbol):
        return GeneInfo.objects.filter(gene_list__genelistgenesymbol__gene_symbol=gene_symbol).distinct()


class CanonicalTranscriptCollection(TimeStampedModel):
    description = models.TextField(blank=True)
    filename = models.TextField(blank=True)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)
    file_sha256sum = models.TextField()

    @staticmethod
    def get_default():
        ctc_id = settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID
        ctc = None
        if ctc_id:
            try:
                ctc = CanonicalTranscriptCollection.objects.get(pk=ctc_id)
            except CanonicalTranscriptCollection.DoesNotExist:
                logging.error("setting.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID=%s - record not found", ctc_id)
        return ctc

    def get_absolute_url(self):
        return reverse('view_canonical_transcript_collection', kwargs={"pk": self.pk})

    def __str__(self):
        count = self.canonicaltranscript_set.count()
        name = f"{count} transcripts"
        description = self.description or os.path.basename(self.filename)
        if description:
            name = f"'{description}' ({name})"
        return name


class CanonicalTranscript(models.Model):
    collection = models.ForeignKey(CanonicalTranscriptCollection, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, null=True, blank=True, on_delete=SET_NULL)
    transcript_version = models.ForeignKey(TranscriptVersion, null=True, blank=True, on_delete=SET_NULL)
    original_gene_symbol = models.TextField()
    original_transcript = models.TextField()


class GeneCoverageCollection(RelatedModelsPartitionModel):
    """ Note both GeneCoverage and GeneCoverageCanonicalTranscript point off same collection object """
    RECORDS_BASE_TABLE_NAMES = ["genes_genecoverage", "genes_genecoveragecanonicaltranscript"]
    RECORDS_FK_FIELD_TO_THIS_MODEL = "gene_coverage_collection_id"
    PARTITION_LABEL_TEXT = "collection"

    path = models.TextField()
    data_state = models.CharField(max_length=1, choices=DataState.choices)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)

    @staticmethod
    def get_gene_coverage_for_sample(sample):
        """ returns coverage (whether via seqauto or standalone uploaded coverage """

        # TODO: Try and get from UploadedGeneCoverage

        gene_coverage = None
        try:
            sequencing_sample = sample.samplefromsequencingsample.sequencing_sample
            bam_file = sequencing_sample.get_single_bam()
            try:
                qc = bam_file.qc_set.get()
                gene_coverage = qc.qcgenecoverage.gene_coverage_collection
            except (ObjectDoesNotExist, MultipleObjectsReturned):
                logging.error("Wasn't exactly 1 qc for bam_file %s", bam_file)
        except ObjectDoesNotExist:
            pass
        except:
            log_traceback()

        return gene_coverage

    def load_from_file(self, enrichment_kit, **kwargs):
        logging.debug("GeneCoverageCollection.load_for_qc()")
        try:
            gene_matcher = kwargs["gene_matcher"]
            canonical_transcript_manager = kwargs.get("canonical_transcript_manager")
            transcript_versions_by_id = kwargs.get("transcript_versions_by_id")
        except KeyError as ke:
            missing_key = ke.args[0]
            logging.error("You need to pass '%s' to kwargs", missing_key)
            raise ke

        missing_genes = 0
        missing_transcripts = 0

        canonical_transcript_collection = canonical_transcript_manager.get_canonical_collection_for_enrichment_kit(enrichment_kit)
        _, original_canonical_transcript_accessions = canonical_transcript_manager.get_canonical_transcripts(canonical_transcript_collection)

        gene_coverage_df = load_gene_coverage_df(self.path)

        gene_coverage_tuples = []
        gene_coverage_canonical_transcript_tuples = []
        warnings = []

        for _, row in gene_coverage_df.iterrows():
            original_gene_symbol = row["original_gene_symbol"]
            original_transcript = row["original_transcript"]

            gene_symbol_id = gene_matcher.get_gene_symbol_id(original_gene_symbol)
            transcript_id, version = TranscriptVersion.get_transcript_id_and_version(original_transcript)
            if transcript_versions := transcript_versions_by_id.get(transcript_id):
                transcript_version_id = transcript_versions.get(version)
            else:
                transcript_id = None
                transcript_version_id = None

            gc_dict = {
                "gene_coverage_collection_id": self.pk,
                "transcript_id": transcript_id,
                "transcript_version_id": transcript_version_id,
                "gene_symbol_id": gene_symbol_id,
            }
            gc_dict.update(row.to_dict())
            gc_tup = tuple(gc_dict.get(f, '') for f in GENE_COVERAGE_HEADER)

            if settings.SEQAUTO_QC_GENE_COVERAGE_STORE_ALL:
                gene_coverage_tuples.append(gc_tup)

            if settings.SEQAUTO_QC_GENE_COVERAGE_STORE_CANONICAL:
                if original_transcript in original_canonical_transcript_accessions:
                    gene_coverage_canonical_transcript_tuples.append((canonical_transcript_collection.pk,) + gc_tup)

        processing_dir = os.path.join(settings.IMPORT_PROCESSING_DIR, "gene_coverage",
                                      f"gene_coverage_collection_{self.pk}")
        mk_path(processing_dir)
        if gene_coverage_tuples:
            csv_filename = os.path.join(processing_dir, f"gene_coverage_{self.pk}.csv")
            write_sql_copy_csv(gene_coverage_tuples, csv_filename)
            partition_table = self.get_partition_table(base_table_name="genes_genecoverage")
            gene_coverage_sql_copy_csv(csv_filename, partition_table)

        if settings.SEQAUTO_QC_GENE_COVERAGE_STORE_CANONICAL:
            if gene_coverage_canonical_transcript_tuples:
                csv_filename = os.path.join(processing_dir, f"gene_coverage_canonical_transcript_{self.pk}.csv")
                write_sql_copy_csv(gene_coverage_canonical_transcript_tuples, csv_filename)
                partition_table = self.get_partition_table(base_table_name="genes_genecoveragecanonicaltranscript")
                gene_coverage_canonical_transcript_sql_copy_csv(csv_filename, partition_table)
            else:
                logging.warning("GeneCoverage had no canonical transcripts")

        if not settings.DEBUG:
            shutil.rmtree(processing_dir)

        logging.info("%d missing genes, %d missing transcripts", missing_genes, missing_transcripts)
        return warnings

    def get_uncovered_gene_symbols(self, gene_symbols, min_coverage):
        # Do as inner query to ensure we restrict to gene coverage partition
        covered_qs = GeneCoverageCanonicalTranscript.objects.filter(gene_coverage_collection=self,
                                                                    min__gte=min_coverage)
        covered_gene_symbols = covered_qs.values_list("gene_symbol", flat=True)
        return gene_symbols.exclude(pk__in=covered_gene_symbols)

    def __str__(self):
        return "GeneCoverageCollection " + os.path.basename(self.path)


@receiver(pre_delete, sender=GeneCoverageCollection)
def gene_coverage_collection_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    try:
        instance.delete_related_objects()
    except:
        pass


class AbstractGeneCoverage(models.Model):
    gene_coverage_collection = models.ForeignKey(GeneCoverageCollection, on_delete=CASCADE)  # rename to "coverage"?
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, null=True, blank=True, on_delete=SET_NULL)
    transcript_version = models.ForeignKey(TranscriptVersion, null=True, blank=True, on_delete=SET_NULL)
    original_gene_symbol = models.TextField()  # raw input string
    original_transcript = models.TextField()  # raw input string
    size = models.IntegerField(null=True)
    min = models.IntegerField()
    max = models.IntegerField(null=True)
    mean = models.FloatField()
    std_dev = models.FloatField()
    # As per QCExecSummary - different enrichment kits store different coverages
    # So some will be None
    percent_1x = models.FloatField(null=True)
    percent_2x = models.FloatField(null=True)
    percent_5x = models.FloatField(null=True)
    percent_10x = models.FloatField(null=True)
    percent_15x = models.FloatField(null=True)
    percent_20x = models.FloatField(null=True)
    percent_25x = models.FloatField(null=True)
    percent_30x = models.FloatField(null=True)
    percent_40x = models.FloatField(null=True)
    percent_50x = models.FloatField(null=True)
    percent_60x = models.FloatField(null=True)
    percent_80x = models.FloatField(null=True)
    percent_100x = models.FloatField(null=True)
    percent_150x = models.FloatField(null=True)
    percent_200x = models.FloatField(null=True)
    percent_250x = models.FloatField(null=True)

    class Meta:
        abstract = True

    @classmethod
    def get_for_symbol(cls, genome_build, gene_symbol):
        return cls.objects.filter(gene_symbol=gene_symbol, gene_coverage_collection__genome_build=genome_build)


class GeneCoverage(AbstractGeneCoverage):
    """ We attempt to match transcript from refseq_transcript_id then if that fails,
        match on gene. So it's possible to have a gene match but no transcript """


class GeneCoverageCanonicalTranscript(AbstractGeneCoverage):
    """ A single transcript per gene to be used for some metrics """

    canonical_transcript_collection = models.ForeignKey(CanonicalTranscriptCollection, null=True, on_delete=SET_NULL)

    @staticmethod
    def filter_for_kit_and_gene_symbol(enrichment_kit, genome_build, gene_symbol):
        sequencing_sample = "gene_coverage_collection__qcgenecoverage__qc__bam_file__unaligned_reads__sequencing_sample"
        kwargs = {sequencing_sample + "__enrichment_kit": enrichment_kit,
                  # Ensure we only get current SampleSheet
                  sequencing_sample + "__sample_sheet__sequencingruncurrentsamplesheet__isnull": False}
        return GeneCoverageCanonicalTranscript.get_for_symbol(genome_build, gene_symbol).filter(**kwargs)


class GnomADGeneConstraint(models.Model):
    """ @see https://gnomad.broadinstitute.org/downloads#v4-constraint """
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, null=True, on_delete=CASCADE)
    transcript_version = models.ForeignKey(TranscriptVersion, null=True, on_delete=CASCADE)
    cached_web_resource = models.ForeignKey('annotation.CachedWebResource', on_delete=CASCADE)

    # Boolean indicator for MANE Select transcript
    mane_select = models.BooleanField(null=True, blank=True)

    # Number of observed high and low confidence loss-of-function (pLoF) variants
    lof_hc_lc_obs = models.IntegerField(null=True, blank=True)

    # Number of expected high and low confidence pLoF variants
    lof_hc_lc_exp = models.IntegerField(null=True, blank=True)

    # Number of possible high and low confidence pLoF variants
    lof_hc_lc_possible = models.IntegerField(null=True, blank=True)

    # Observed over expected ratio for high and low confidence pLoF variants
    lof_hc_lc_oe = models.FloatField(null=True, blank=True)

    # Mutation rate for high and low confidence pLoF variants
    lof_hc_lc_mu = models.FloatField(null=True, blank=True)

    # Probability of loss-of-function intolerance
    lof_hc_lc_pli = models.FloatField(null=True, blank=True)

    # Probability for recessive genes
    lof_hc_lc_prec = models.FloatField(null=True, blank=True)

    # Probability for unconstrained genes
    lof_hc_lc_pnull = models.FloatField(null=True, blank=True)

    # Number of observed high confidence pLoF variants
    lof_obs = models.IntegerField(null=True, blank=True)

    # Number of expected high confidence pLoF variants
    lof_exp = models.FloatField(null=True, blank=True)

    # Number of possible high confidence pLoF variants
    lof_possible = models.IntegerField(null=True, blank=True)

    # Observed over expected ratio for high confidence pLoF variants
    lof_oe = models.FloatField(null=True)

    # Mutation rate for high confidence pLoF variants
    lof_mu = models.FloatField(null=True, blank=True)

    # Probability of loss-of-function intolerance (high confidence)
    lof_pli = models.FloatField(null=True, blank=True)

    # Probability for recessive genes (high confidence)
    lof_prec = models.FloatField(null=True, blank=True)

    # Probability for unconstrained genes (high confidence)
    lof_pnull = models.FloatField(null=True, blank=True)

    # Lower bound of 90% CI for o/e ratio (high confidence pLoF variants)
    lof_oe_ci_lower = models.FloatField(null=True)

    # Upper bound of 90% CI for o/e ratio (high confidence pLoF variants)
    lof_oe_ci_upper = models.FloatField(null=True)

    # Raw Z-score for pLoF variants
    lof_z_raw = models.FloatField(null=True, blank=True)

    # Z-score for pLoF variants
    lof_z_score = models.FloatField(null=True, blank=True)

    # Number of observed missense variants
    mis_obs = models.IntegerField(null=True, blank=True)

    # Number of expected missense variants
    mis_exp = models.FloatField(null=True, blank=True)

    # Number of possible missense variants
    mis_possible = models.IntegerField(null=True, blank=True)

    # Observed over expected ratio for missense variants
    mis_oe = models.FloatField(null=True, blank=True)

    # Mutation rate for missense variants
    mis_mu = models.FloatField(null=True, blank=True)

    # Lower bound of 90% CI for o/e ratio (missense variants)
    mis_oe_ci_lower = models.FloatField(null=True, blank=True)

    # Upper bound of 90% CI for o/e ratio (missense variants)
    mis_oe_ci_upper = models.FloatField(null=True, blank=True)

    # Raw Z-score for missense variants
    mis_z_raw = models.FloatField(null=True, blank=True)

    # Z-score for missense variants
    mis_z_score = models.FloatField(null=True, blank=True)

    # Number of observed "probably damaging" missense variants (PolyPhen-2)
    mis_pphen_obs = models.IntegerField(null=True, blank=True)

    # Number of expected "probably damaging" missense variants (PolyPhen-2)
    mis_pphen_exp = models.IntegerField(null=True, blank=True)

    # Number of possible "probably damaging" missense variants (PolyPhen-2)
    mis_pphen_possible = models.IntegerField(null=True, blank=True)

    # O/E ratio for "probably damaging" missense variants (PolyPhen-2)
    mis_pphen_oe = models.FloatField(null=True, blank=True)

    # Number of observed synonymous variants
    syn_obs = models.IntegerField(null=True, blank=True)

    # Number of expected synonymous variants
    syn_exp = models.FloatField(null=True, blank=True)

    # Number of possible synonymous variants
    syn_possible = models.IntegerField(null=True, blank=True)

    # Observed over expected ratio for synonymous variants
    syn_oe = models.FloatField(null=True, blank=True)

    # Mutation rate for synonymous variants
    syn_mu = models.FloatField(null=True, blank=True)

    # Lower bound of 90% CI for o/e ratio (synonymous variants)
    syn_oe_ci_lower = models.FloatField(null=True, blank=True)

    # Upper bound of 90% CI for o/e ratio (synonymous variants)
    syn_oe_ci_upper = models.FloatField(null=True, blank=True)

    # Raw Z-score for synonymous variants
    syn_z_raw = models.FloatField(null=True, blank=True)

    # Z-score for synonymous variants
    syn_z_score = models.FloatField(null=True, blank=True)

    # Reason transcript does not have constraint metrics
    constraint_flags = models.JSONField(default=list)

    GNOMAD_BASE_URL = "https://gnomad.broadinstitute.org"

    @property
    def transcript_url(self):
        if tv := self.transcript_version:
            transcript_id = tv.transcript_id
        else:
            transcript_id = self.transcript_id
        return f"{self.GNOMAD_BASE_URL}/transcript/{transcript_id}"

    @property
    def gene_url(self):
        return f"{self.GNOMAD_BASE_URL}/gene/{self.gene_id}"

    @property
    def gene_symbol_url(self):
        return f"{self.GNOMAD_BASE_URL}/gene/{self.gene_symbol_id}"

    def _get_oe_summary(self, prefix: str) -> str:
        oe = getattr(self, f"{prefix}_oe")
        summary = "n/a"
        if oe is not None:
            oe_ci_lower = getattr(self, f"{prefix}_oe_ci_lower")
            oe_ci_upper = getattr(self, f"{prefix}_oe_ci_upper")
            summary = f"{oe:.2f} ({oe_ci_lower:.2f} - {oe_ci_upper:.2f})"
        return summary

    @property
    def mis_oe_summary(self) -> str:
        return self._get_oe_summary("mis")

    @property
    def syn_oe_summary(self) -> str:
        return self._get_oe_summary("syn")

    @property
    def lof_oe_summary(self) -> str:
        return self._get_oe_summary("lof")

    @property
    def transcript_accession(self) -> str:
        if self.transcript_version:
            return str(self.transcript_version)
        return str(self.transcript)

    @staticmethod
    def get_for_transcript_version(transcript_version: TranscriptVersion) -> Optional['GnomADGeneConstraint']:
        """ GnomADGeneConstraint uses Ensembl gene/transcripts - so load the most specific
            possible (transcript, gene, then symbol) """
        return GnomADGeneConstraint.get_for_transcript_version_with_method_and_url(transcript_version)[0]

    @staticmethod
    def get_for_transcript_version_with_method_and_url(transcript_version: TranscriptVersion) \
            -> tuple[Optional['GnomADGeneConstraint'], Optional[str], Optional[str]]:
        """ GnomADGeneConstraint uses Ensembl gene/transcripts - so load the most specific
            possible (transcript version, transcript then symbol) """
        ggc = None
        method = None
        url = None
        qs = GnomADGeneConstraint.objects.all()
        # May be able to get via transcript version or just transcript
        if ggc := qs.filter(transcript_version=transcript_version).first():
            url = ggc.transcript_url
            method = "transcript version"
        elif ggc := qs.filter(transcript=transcript_version.transcript).first():
            url = ggc.transcript_url
            method = "transcript"

        if ggc is None:
            if ggc := qs.filter(gene_symbol=transcript_version.gene_version.gene_symbol).first():
                url = ggc.gene_symbol_url
                method = "gene symbol"
        return ggc, method, url


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

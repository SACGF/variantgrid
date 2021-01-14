import json
import logging
import os
import re
from collections import namedtuple
from dataclasses import dataclass
from functools import total_ordering
from typing import Tuple, Optional, Dict, List, Set

from django.conf import settings
from django.contrib.auth.models import User
from django.contrib.postgres.aggregates import StringAgg
from django.core.exceptions import PermissionDenied, ObjectDoesNotExist, MultipleObjectsReturned
from django.db import models, IntegrityError, transaction
from django.db.models import Min, Max, QuerySet
from django.db.models.deletion import CASCADE, SET_NULL, PROTECT
from django.db.models.fields.json import KeyTextTransform
from django.db.models.functions import Upper
from django.db.models.query_utils import Q
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.shortcuts import get_object_or_404
from django.urls.base import reverse
from django.utils.timezone import localtime
from django_extensions.db.models import TimeStampedModel
from guardian.shortcuts import get_objects_for_user
from lazy import lazy
from pyhgvs.utils import make_transcript

from genes.annotation_consortium_api import transcript_exists
from genes.gene_coverage import load_gene_coverage_df
from genes.models_enums import AnnotationConsortium, HGNCStatus, GeneSymbolAliasSource
from library.django_utils import SortByPKMixin
from library.guardian_utils import assign_permission_to_user_and_groups, DjangoPermission
from library.log_utils import log_traceback
from library.utils import empty_dict
from seqauto.models_enums import DataState
from snpdb.models import Wiki, Company, Sample
from snpdb.models.models_enums import ImportStatus
from snpdb.models.models_genome import GenomeBuild


class HGNCGeneNamesImport(TimeStampedModel):
    pass


class MissingTranscript(ValueError):
    """
    Extends ValueError for backwards compatibility.
    Indicates the transcript we are looking for is not in our database,
    but does exist in RefSeq/Ensembl - so the c.hgvs (or otherwise) might be okay.
    """
    pass


class BadTranscript(ValueError):
    """
    Extends ValueError for backwards compatibility.
    Indicates the transcript is not found in Ensembl or RefSeq (User error)
    """
    pass


class HGNCGeneNames(models.Model):
    # pk = HGNC id with HGNC: stripped out
    hgnc_import = models.ForeignKey(HGNCGeneNamesImport, on_delete=CASCADE)
    approved_symbol = models.TextField()
    approved_name = models.TextField()
    status = models.CharField(max_length=1, choices=HGNCStatus.choices)
    previous_symbols = models.TextField()
    synonyms = models.TextField()
    refseq_ids = models.TextField()

    def __str__(self):
        return f"HGNC:{self.pk} approved symbol: {self.approved_symbol}, " \
               f"previous symbols: {self.previous_symbols}, synonyms: {self.synonyms}"


gene_symbol_withdrawn_str = '~withdrawn'


class GeneSymbol(models.Model):
    symbol = models.TextField(primary_key=True)

    @property
    def accession(self):
        """ For use by TextPhenotypeMatch """
        return self.symbol

    @property
    def name(self):
        """ For use by TextPhenotypeMatch """
        return self.symbol

    def get_genes(self) -> QuerySet:
        # To match HPO/OMIM so it can be used interchangeably during phenotype matching
        return Gene.objects.filter(~Q(identifier__startswith="unknown_"), geneversion__gene_symbol=self).distinct()

    @lazy
    def genes(self) -> List['Gene']:
        # returns cached set of genes associated with this symbol
        # use over get_genes when possible
        return list(self.get_genes().all())

    def get_absolute_url(self):
        return reverse("view_gene_symbol", kwargs={"gene_symbol": self.symbol})

    @lazy
    def alias_meta(self) -> 'GeneSymbolAliasesMeta':
        return GeneSymbolAliasesMeta(self)

    def has_different_genes(self, other: 'GeneSymbol') -> bool:
        """
        Tries to work out if genes are equivilant, not that sometimes refseq or ensembl assign gene ids to both the
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
         * Code: python3.8 manage.py import_ncbi_gene_info <file>

        HGNC:
         * Source: https://www.genenames.org/cgi-bin/download
         * Code: python3.8 manage.py hgnc_gene_symbols_import <file>

        UCSC: We no longer use UCSC aliases, they will only exist upgraded legacy systems
         * Source: https://genome.ucsc.edu/cgi-bin/hgTables?command=start export kgAlias table
         * Code: N/A - obsolete
    """
    alias = models.TextField(unique=True)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    source = models.CharField(max_length=1, choices=GeneSymbolAliasSource.choices)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    description = models.TextField(null=True)

    def __str__(self):
        return f"{self.gene_symbol_id} : {self.alias} is an alias for {self.gene_symbol_id} ({self.get_source_display()})"

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
        self.alias_list: List[GeneSymbolAliasSummary] = list()

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

    @lazy
    def genes(self) -> Set['Gene']:
        """
        Returns a set of genes associated with all safe aliases to/from the primary Gene Symbol.
        (Even though we only look at "safe" aliases, e.g. ones where each symbol must be a subset of the other,
        looking through these aliases still catch where Refseq assigned a Gene ID to both but Ensembl only assigned
        their Gene ID to one and ignore the other)
        """
        gene_set: Set[Gene] = set(self.gene_symbol.genes)
        for alias_summary in self.alias_list:
            if not alias_summary.different_genes and alias_summary.other_obj:
                gene_set = gene_set.union(alias_summary.other_obj.genes)
        return gene_set

    @lazy
    def alias_symbol_strs(self) -> List[str]:
        gene_symbol_strs: Set[str] = {self.gene_symbol.symbol}
        for alias_summary in self.alias_list:
            if not alias_summary.different_genes:
                gene_symbol_strs.add(alias_summary.other_symbol)
        return list(sorted(gene_symbol_strs))

    @lazy
    def aliases_out(self) -> List[GeneSymbolAliasSummary]:
        return list(sorted([alias for alias in self.alias_list if not alias.my_symbol_is_main]))

    @lazy
    def aliases_in(self) -> List[GeneSymbolAliasSummary]:
        return list(sorted([alias for alias in self.alias_list if alias.my_symbol_is_main]))


class GeneAnnotationImport(TimeStampedModel):
    """ A GTF file imported via 'python3 manage import_gene_annotation'

        Many gene/transcript versions are shared among GTF annotations, so a GeneVersion/TranscriptVersion is only
        created the first time it's seen (linked back to input which created it via 'import_source') """
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    filename = models.TextField()

    def __str__(self):
        return os.path.basename(self.filename)


class Gene(models.Model):
    """ A stable identifier - build independent - has build specific versions with gene details """
    FAKE_GENE_ID_PREFIX = "unknown_"  # Legacy from when we allowed inserting GenePred w/o GFF3
    identifier = models.TextField(primary_key=True)
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)

    @property
    def is_legacy(self):
        """ Required internally, but probably shouldn't be shown to the user """
        return self.identifier.startswith('unknown_')

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

    def __str__(self):
        if self.annotation_consortium == AnnotationConsortium.REFSEQ:
            gene_id_summary = f"GeneID:{self.identifier}"
        else:
            gene_id_summary = self.identifier

        return f"{gene_id_summary} ({self.get_annotation_consortium_display()})"


class GeneVersion(models.Model):
    """ A specific version of a Gene for a particular version/genome build
        Genes/TranscriptVersion needs to be able to represent both RefSeq and Ensembl """
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    version = models.IntegerField()  # RefSeq GeneIDs are always 1 (not versioned)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    hgnc = models.ForeignKey(HGNCGeneNames, null=True, on_delete=CASCADE)
    description = models.TextField(null=True)
    biotype = models.TextField(null=True)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    import_source = models.ForeignKey(GeneAnnotationImport, on_delete=CASCADE)

    class Meta:
        unique_together = ("gene", "version", "genome_build")

    @lazy
    def accession(self):
        if self.version is not None and self.gene.has_versions():
            acc = f"{self.gene_id}.{self.version}"
        else:
            acc = self.gene_id
        return acc

    @property
    def name(self):
        return self.gene_symbol_id

    @property
    def chrom(self):
        return self._transcript_extents["chroms"]

    @property
    def start(self):
        return int(self._transcript_extents["min_start"])

    @property
    def end(self):
        return int(self._transcript_extents["max_end"])

    @property
    def strand(self):
        return self._transcript_extents["strands"]

    @lazy
    def _transcript_extents(self):
        """ Stores chroms/min_start/max_end/strands - which are aggregate of all linked TranscriptVersions """
        DELIMITER = ","
        qs = self.transcriptversion_set.filter(data__strand__isnull=False)
        qs = qs.annotate(**{k: KeyTextTransform(k, "data") for k in ["chrom", "start", "end", "strand"]})
        data = qs.aggregate(chroms=StringAgg("chrom", delimiter=DELIMITER, distinct=True),
                            min_start=Min("start"), max_end=Max("end"),
                            strands=StringAgg("strand", delimiter=DELIMITER, distinct=True))
        for k in ["chroms", "strands"]:
            v = data[k]
            if len(v.split(DELIMITER)) != 1:
                raise ValueError(f"{self}: Not exactly 1 value for {k}, was: {v}")
        return data

    @staticmethod
    def get_gene_id_and_version(gene_accession: str) -> Tuple[str, int]:
        parts = gene_accession.split(".")
        if len(parts) == 2:
            identifier = str(parts[0])
            version = int(parts[1])
        else:
            identifier, version = gene_accession, None
        return identifier, version

    def __str__(self):
        return f"{self.accession} ({self.gene_symbol}/{self.genome_build})"


TranscriptParts = namedtuple('TranscriptParts', ['identifier', 'version'])


class Transcript(models.Model):
    """ A stable identifier - has versions with actual transcript details """
    identifier = models.TextField(primary_key=True)
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)

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
    def known_transcript_ids(genome_build=None, annotation_consortium=None) -> Set[str]:
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


class TranscriptVersion(SortByPKMixin, models.Model):
    """ We store the ID and version separately, ie:
        ENST00000284274.4 => transcript=ENST00000284274, version=4

        Ensembl ID info: https://m.ensembl.org/Help/Faq?id=488

        A useful query to get the latest version for each transcript is:
        qs.order_by("transcript_id", "-version").distinct("transcript_id")
    """
    transcript = models.ForeignKey(Transcript, on_delete=CASCADE)
    version = models.IntegerField()
    gene_version = models.ForeignKey(GeneVersion, on_delete=CASCADE)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    import_source = models.ForeignKey(GeneAnnotationImport, on_delete=CASCADE)
    biotype = models.TextField(null=True)  # Ensembl has gene + transcript biotypes
    data = models.JSONField(null=False, blank=True, default=empty_dict)  # for pyHGVS

    class Meta:
        unique_together = ("transcript", "version", "genome_build")

    @lazy
    def cds(self):
        cds_start = self.data["cds_start"]
        cds_end = self.data["cds_end"]
        cds = []
        for exon_start, exon_end in self.data["exons"]:
            if exon_end > cds_start and exon_start < cds_end:
                start = max(exon_start, cds_start)
                end = min(exon_end, cds_end)
                cds.append((start, end))
        return cds

    @property
    def num_amino_acids(self):
        return sum([b-a for a, b in self.cds]) / 3

    @staticmethod
    def transcript_parts(identifier: str) -> TranscriptParts:
        t_regex = re.compile(r"^([_A-Z0-9]+)(?:[.]([0-9]+))?$", re.RegexFlag.IGNORECASE)
        if m := t_regex.match(identifier):
            version = m.group(2)
            if version:
                version = int(version)
            return TranscriptParts(m.group(1), version)
        raise ValueError(f'Invalid transcript identifier {identifier}')

    @staticmethod
    def get(accession: str, genome_build: GenomeBuild, annotation_consortium):
        """
        @param accession transcript_id w/optional version
        """
        parts = TranscriptVersion.transcript_parts(accession)
        t_qs = TranscriptVersion.objects.filter(transcript_id=parts[0],
                                                transcript__annotation_consortium=annotation_consortium,
                                                genome_build=genome_build)
        if parts[1]:
            t_qs = t_qs.filter(version=parts[1])
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

    @lazy
    def accession(self):
        return TranscriptVersion.get_accession(self.transcript_id, self.version)

    @property
    def gene(self):
        return self.gene_version.gene

    @staticmethod
    def get_transcript_id_and_version(transcript_name: str) -> Tuple[str, int]:
        parts = transcript_name.split(".")
        if len(parts) == 2:
            identifier = str(parts[0])
            version = int(parts[1])
        else:
            identifier, version = transcript_name, None
        return identifier, version

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
    def raise_bad_or_missing_transcript(genome_build: GenomeBuild, identifier, version):
        """ Checks whether a transcript we can't match is wrong (their fault) or we don't have it (our fault) """

        annotation_consortium, exists = transcript_exists(genome_build, identifier, version)
        if exists:
            raise MissingTranscript(f"Transcript '{identifier}' valid but missing from our (build: {genome_build}) database.")
        accession = TranscriptVersion.get_accession(identifier, version)
        raise BadTranscript(f"The transcript '{accession}' does not exist in the {genome_build} {annotation_consortium} database.")

    @staticmethod
    def get_refgene(genome_build: GenomeBuild, transcript_name, best_attempt=True):
        """ Return a pyhgvs RefGene object
        @param best_attempt if we don't have the exact version specified in transcript_name, grab the closest
        one that is larger (or the largest one if we don't have any larger than the requested)
        """
        transcript_version = None
        transcript_versions_qs = TranscriptVersion.objects.filter(genome_build=genome_build)
        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_name)
        transcript_versions_qs = transcript_versions_qs.filter(transcript_id=transcript_id).order_by("version")
        if version is not None:
            try:
                transcript_version = transcript_versions_qs.get(version=version)
            except TranscriptVersion.DoesNotExist:
                possible_versions = set(transcript_versions_qs.values_list('version', flat=True).all())
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
            TranscriptVersion.raise_bad_or_missing_transcript(genome_build, transcript_id, version)

        if 'id' not in transcript_version.data:
            # only going to happen if we have legacy data in the database, transcripts that use the default for data {}
            data_str = json.dumps(transcript_version.data)
            raise MissingTranscript(f"Transcript for '{transcript_name}' (build: {genome_build}),"
                                    f" but did not have complete data {data_str}")

        # Legacy data stored gene_name in JSON, but that could lead to diverging values vs TranscriptVersion relations
        # so replace it with the DB records
        transcript_version.data["gene_name"] = transcript_version.gene_version.gene_symbol_id
        return make_transcript(transcript_version.data)

    @property
    def length(self):
        if 'exons' in self.data:
            return sum([e - s for s, e in self.data["exons"]])
        return None

    @property
    def coordinates(self):
        coords = None
        if self.data:
            coords = f"{self.data['chrom']}:{self.data['start']}-{self.data['end']} ({self.data['strand']})"
        return coords

    def get_differences(self, transcript):
        """ Used to inform while HGVS may resolve differently """
        differences = {}
        FIELDS = ["transcript_id", "version", "length"]
        for f in FIELDS:
            mine = getattr(self, f)
            other = getattr(transcript, f)
            if mine != other:
                differences[f] = (mine, other)

        if self.data and transcript.data:
            my_chrom = self.data["chrom"]
            other_chrom = transcript.data["chrom"]
            if my_chrom != other_chrom:
                differences["contig"] = (my_chrom, other_chrom)

            my_exon_count = len(self.data["exons"])
            other_exon_count = len(transcript.data["exons"])
            if my_exon_count != other_exon_count:
                differences["exon count"] = (my_exon_count, other_exon_count)
            else:
                exon = 1
                for my_exon, other_exon in zip(self.data["exons"], transcript.data["exons"]):
                    my_len = my_exon[1] - my_exon[0]
                    other_len = other_exon[1] - other_exon[0]
                    if my_len != other_len:
                        differences[f"exon {exon}"] = (my_len, other_len)
                    exon += 1
        else:
            my_keys = set(self.data.keys())
            other_keys = set(transcript.data.keys())
            if my_keys ^ other_keys:
                differences["data"] = (my_keys, other_keys)

        return differences

    @staticmethod
    def get_preferred_transcript(data: Dict[str, 'TranscriptVersion']) -> Optional['TranscriptVersion']:
        for transcript_key in settings.VARIANT_ANNOTATION_TRANSCRIPT_PREFERENCES:
            transcript = data.get(transcript_key)
            if transcript:
                return transcript
        return None

    def get_domains(self) -> List[Tuple]:
        """ Gets ProteinDomain if available, falling back on Pfam """
        PD_ARGS = ("protein_domain__name", "protein_domain__description", "start", "end")
        values_qs = self.proteindomaintranscriptversion_set.all().order_by("start").values_list(*PD_ARGS)
        if not values_qs.exists():
            pfam = self.pfamsequenceidentifier_set.first()
            if not pfam:
                pfam = self.transcript.pfamsequenceidentifier_set.first()

            if pfam:
                PFAM_ARGS = ("pfam__pfam_id", "pfam__description", "start", "end")
                values_qs = pfam.pfam_sequence.pfamdomains_set.order_by("start").values_list(*PFAM_ARGS)
        return values_qs

    def get_absolute_url(self):
        kwargs = {"transcript_id": self.transcript_id, "version": self.version}
        return reverse("view_transcript_version", kwargs=kwargs)

    def get_external_url(self):
        return self.transcript.get_external_url(self.genome_build) + f".{self.version}"

    def __str__(self):
        return f"{self.accession} ({self.gene_version.gene_symbol}/{self.genome_build.name})"


class GeneAnnotationRelease(models.Model):
    """
        import_gene_annotation management command now supports --release which is used to link all
        GeneVersion/TranscriptVersion from a particular GTF file.

        This release can be set on a VariantAnnotationVersion to be able to get genes/transcripts from a VEP build
    """
    version = models.IntegerField()
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    gene_annotation_import = models.ForeignKey(GeneAnnotationImport, on_delete=CASCADE)

    class Meta:
        unique_together = ('version', 'annotation_consortium', 'genome_build')

    def get_genes(self):
        return Gene.objects.filter(annotation_consortium=self.annotation_consortium,
                                   geneversion__releasegeneversion__release=self)

    @staticmethod
    def get_for_latest_annotation_versions_for_builds() -> List['GeneAnnotationRelease']:
        """ """
        gene_annotation_releases = []
        for genome_build in GenomeBuild.builds_with_annotation().order_by("name"):
            vav = genome_build.latest_variant_annotation_version
            if vav.gene_annotation_release:
                gene_annotation_releases.append(vav.gene_annotation_release)
        return gene_annotation_releases

    def __str__(self):
        return f"{self.genome_build.name}/{self.get_annotation_consortium_display()} - v{self.version}"


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
        from annotation.models import Citation

        existing_citations = {}
        for gc in self.gene_symbol.genesymbolcitation_set.all():
            existing_citations[gc.citation.unique_code()] = gc

        gene_citation_ids_to_keep = []
        for citation in Citation.citations_from_text(self.markdown):
            gene_citation = existing_citations.get(citation.unique_code())
            if not gene_citation:
                citation.save()
                gene_citation = self.gene_symbol.genesymbolcitation_set.create(citation=citation)
            gene_citation_ids_to_keep.append(gene_citation.pk)

        # Delete any non-used ones
        num_deleted = self.gene_symbol.genesymbolcitation_set.exclude(pk__in=gene_citation_ids_to_keep).delete()
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


class GeneList(models.Model):
    """ Stores a gene/transcript list (to be used as a filter) """

    category = models.ForeignKey(GeneListCategory, null=True, blank=True, on_delete=CASCADE)
    name = models.TextField()
    user = models.ForeignKey(User, on_delete=CASCADE)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    error_message = models.TextField(null=True, blank=True)
    locked = models.BooleanField(default=False)
    url = models.TextField(null=True, blank=True)

    def get_q(self, release: GeneAnnotationRelease):
        """ For a Variant queryset """
        from annotation.models.models import VariantTranscriptAnnotation
        return VariantTranscriptAnnotation.get_overlapping_genes_q(self.get_genes(release))

    @staticmethod
    def get_gene_ids_for_gene_lists(release: GeneAnnotationRelease, gene_lists: List['GeneList']):
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
        genelistgenesymbol_set = self.genelistgenesymbol_set.all()
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

    def can_view(self, user):
        read_perm = DjangoPermission.perm(self, DjangoPermission.READ)
        return user.has_perm(read_perm, self)

    def can_write(self, user):
        write_perm = DjangoPermission.perm(self, DjangoPermission.WRITE)
        return user.has_perm(write_perm, self) and not self.locked

    @lazy
    def warning(self):
        warnings = {"unmatched": self.unmatched_genes.count(),
                    "aliased": self.aliased_genes.count()}

        warning_list = []
        prefix = None
        for w, num in warnings.items():
            if num:
                if prefix is None:
                    prefix = "There were " if num > 1 else "There was "

                msg = f"{num} {w} gene"
                if num > 1:
                    msg += "s"
                warning_list.append(msg)

        if warning_list:
            message = prefix + " and ".join(warning_list) + "."
        else:
            message = ""
        return message

    @property
    def unmatched_genes(self):
        return self.genelistgenesymbol_set.filter(gene_symbol__isnull=True).order_by("original_name")

    @property
    def aliased_genes(self):
        return self.genelistgenesymbol_set.filter(gene_symbol_alias__isnull=False)

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

        qs = self.genelistgenesymbol_set.filter(gene_symbol__in=gene_symbol_deletions)
        num_deleted = qs.delete()
        return num_added, num_deleted

    @staticmethod
    def get_for_user(user, gene_list_id, success_only=True):
        try:
            return GeneList.filter_for_user(user, success_only).get(pk=gene_list_id)
        except:
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
        return qs

    @staticmethod
    def _visible_gene_lists(gene_lists_qs):
        hidden_categories_qs = GeneListCategory.objects.filter(hidden=True)
        return gene_lists_qs.exclude(category__in=hidden_categories_qs)  # Allow through no category

    @classmethod
    def visible_gene_lists_containing_gene_symbol(cls, gene_lists_qs, gene_symbol):
        qs = cls._visible_gene_lists(gene_lists_qs)
        gene_symbol = get_object_or_404(GeneSymbol, pk=gene_symbol)
        return qs.filter(genelistgenesymbol__gene_symbol=gene_symbol)

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


class FakeGeneList(GeneList):
    """ Fake for serializing """

    class Meta:
        abstract = True

    def get_absolute_url(self):
        return None


class GeneListGeneSymbol(models.Model):
    """ Replacement for GeneListGene - will keep both in place during switchover """
    gene_list = models.ForeignKey(GeneList, on_delete=CASCADE)
    original_name = models.TextField(null=True, blank=True)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    gene_symbol_alias = models.ForeignKey(GeneSymbolAlias, null=True, on_delete=CASCADE)
    modification_info = models.TextField(null=True)

    class Meta:
        unique_together = ('gene_list', 'original_name')

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
    md5_hash = models.CharField(max_length=32)
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
def sample_gene_list_created(sender, instance, created, **kwargs):
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

    def __str__(self):
        return self.name


class PanelAppPanelRelevantDisorders(models.Model):
    panel_app_panel = models.ForeignKey(PanelAppPanel, on_delete=CASCADE)
    name = models.TextField()


class PanelAppPanelLocalCacheGeneList(TimeStampedModel):
    panel_app_panel = models.ForeignKey(PanelAppPanel, on_delete=CASCADE)
    version = models.TextField()
    gene_list = models.ForeignKey(GeneList, on_delete=CASCADE)

    class Meta:
        unique_together = ("panel_app_panel", "version")


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
    file_md5sum = models.TextField()

    @staticmethod
    def get_default():
        ctc_id = settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID
        ctc = None
        if ctc_id:
            ctc = CanonicalTranscriptCollection.objects.get(pk=ctc_id)
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
    transcript = models.ForeignKey(Transcript, null=True, blank=True, on_delete=CASCADE)
    original_gene_symbol = models.TextField()
    original_transcript_id = models.TextField()


class GeneCoverageCollection(models.Model):
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
        except KeyError as ke:
            missing_key = ke.args[0]
            logging.error("You need to pass '%s' to kwargs", missing_key)
            raise ke

        missing_genes = 0
        missing_transcripts = 0

        canonical_transcript_collection = canonical_transcript_manager.get_canonical_collection_for_enrichment_kit(enrichment_kit)
        _, original_canonical_transcript_ids = canonical_transcript_manager.get_canonical_transcripts(canonical_transcript_collection)

        gene_coverage_df = load_gene_coverage_df(self.path)

        gene_coverage_list = []
        gene_coverage_canonical_list = []
        warnings = []

        known_transcript_ids = Transcript.known_transcript_ids(canonical_transcript_collection.genome_build,
                                                               canonical_transcript_collection.annotation_consortium)
        for _, row in gene_coverage_df.iterrows():
            original_gene_symbol = row["gene"]
            original_transcript_id = row["transcript"]

            gene_symbol_id = gene_matcher.get_gene_symbol_id(original_gene_symbol)
            if original_transcript_id in known_transcript_ids:
                transcript_id = original_transcript_id
            else:
                transcript_id = None

                kwargs = {"gene_coverage_collection": self,
                          "transcript_id": transcript_id,
                          "gene_symbol_id": gene_symbol_id,
                          "original_gene_symbol": original_gene_symbol,
                          "original_transcript_id": original_transcript_id,
                          "min": row["min"],
                          "mean": row["mean"],
                          "std_dev": row["std_dev"],
                          "percent_0x": row["percent_0x"],
                          "percent_10x": row.get("percent_10x"),
                          "percent_20x": row["percent_20x"],
                          "percent_100x": row.get("percent_100x"),
                          "sensitivity": row["sensitivity"]}

                gene_coverage = GeneCoverage(**kwargs)
                gene_coverage_list.append(gene_coverage)

                if original_transcript_id in original_canonical_transcript_ids:
                    kwargs["canonical_transcript_collection"] = canonical_transcript_collection
                    gene_coverage_canonical = GeneCoverageCanonicalTranscript(**kwargs)
                    gene_coverage_canonical_list.append(gene_coverage_canonical)

        if gene_coverage_list:
            GeneCoverage.objects.bulk_create(gene_coverage_list)

        if gene_coverage_canonical_list:
            GeneCoverageCanonicalTranscript.objects.bulk_create(gene_coverage_canonical_list)
        else:
            logging.warning("GeneCoverage had no canonical transcripts")

        logging.info("%d missing genes, %d missing transcripts", missing_genes, missing_transcripts)
        return warnings

    def get_uncovered_gene_symbols(self, gene_symbols, min_coverage):
        return gene_symbols.exclude(genecoverage__gene_coverage_collection=self, genecoverage__min__gte=min_coverage)

    def __str__(self):
        return "GeneCoverageCollection " + os.path.basename(self.path)


class AbstractGeneCoverage(models.Model):
    gene_coverage_collection = models.ForeignKey(GeneCoverageCollection, on_delete=CASCADE)  # rename to "coverage"?
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, null=True, blank=True, on_delete=CASCADE)
    original_gene_symbol = models.TextField()  # raw input string
    original_transcript_id = models.TextField()  # raw input string
    min = models.IntegerField()
    mean = models.FloatField()
    std_dev = models.FloatField()
    # As per QCExecSummary - different enrichment kits store different coverages
    # So some will be None
    percent_0x = models.FloatField()
    percent_10x = models.FloatField(null=True)
    percent_20x = models.FloatField()
    percent_100x = models.FloatField(null=True)
    sensitivity = models.FloatField()  # Estimated Sensitivity of Detection

    class Meta:
        abstract = True

    @classmethod
    def get_for_symbol(cls, genome_build, gene_symbol):
        return cls.objects.filter(gene_symbol=gene_symbol, gene_coverage_collection__genome_build=genome_build)


class GeneCoverage(AbstractGeneCoverage):
    """ We attempt to match transcript from refseq_transcript_id then if that fails,
        match on gene. So it's possible to have a gene match but no transcript """
    pass


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
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, null=True, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, null=True, on_delete=CASCADE)
    cached_web_resource = models.ForeignKey('annotation.CachedWebResource', on_delete=CASCADE)

    # @see https://gnomad.broadinstitute.org/downloads#gene-constraint
    # gnomAD per-gene constraint scores, @see https://macarthurlab.org/2018/10/17/gnomad-v2-1/ ("Constraint" section)
    # oe = observed / expected score. LOW oe = strong intolerance. lower/upper = 90% CI
    oe_lof = models.FloatField(null=True)
    oe_lof_lower = models.FloatField(null=True)
    oe_lof_upper = models.FloatField(null=True)
    # There are dozens of other fields we could add...

    GNOMAD_BASE_URL = "https://gnomad.broadinstitute.org"

    @property
    def transcript_url(self):
        return f"{self.GNOMAD_BASE_URL}/transcript/{self.transcript_id}"

    @property
    def gene_url(self):
        return f"{self.GNOMAD_BASE_URL}/gene/{self.gene_id}"

    @property
    def gene_symbol_url(self):
        return f"{self.GNOMAD_BASE_URL}/gene/{self.gene_symbol_id}"

    @property
    def oe_lof_summary(self):
        summary = "n/a"
        if self.oe_lof is not None:
            summary = f"{self.oe_lof} ({self.oe_lof_lower} - {self.oe_lof_upper})"
        return summary

    @staticmethod
    def get_for_transcript_version(transcript_version: TranscriptVersion) -> Optional['GnomADGeneConstraint']:
        """ GnomADGeneConstraint uses Ensembl gene/transcripts - so load the most specific
            possible (transcript, gene, then symbol) """
        return GnomADGeneConstraint.get_for_transcript_version_with_method_and_url(transcript_version)[0]

    @staticmethod
    def get_for_transcript_version_with_method_and_url(transcript_version: TranscriptVersion) \
            -> Tuple[Optional['GnomADGeneConstraint'], Optional[str], Optional[str]]:
        """ GnomADGeneConstraint uses Ensembl gene/transcripts - so load the most specific
            possible (transcript, gene, then symbol) """
        ggc = None
        method = None
        url = None
        qs = GnomADGeneConstraint.objects.all()
        if transcript_version.transcript.annotation_consortium == AnnotationConsortium.ENSEMBL:
            # May be able to get via transcript or gene
            if ggc := qs.filter(transcript=transcript_version.transcript).first():
                url = ggc.transcript_url
                method = "transcript"

            if ggc is None:
                if ggc := qs.filter(gene=transcript_version.gene).first():
                    url = ggc.gene_url
                    method = "gene"

        if ggc is None:
            if ggc := qs.filter(gene_symbol=transcript_version.gene_version.gene_symbol).first():
                url = ggc.gene_symbol_url
                method = "gene symbol"
        return ggc, method, url


class ProteinDomain(models.Model):
    """ For custom domains, can be used to override PFam - used in TranscriptVersion.get_domains() """
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
    transcript_version = models.ForeignKey(TranscriptVersion, null=True, on_delete=CASCADE)

    @classmethod
    def get_transcript_for_gene(cls, gene: Gene, variant_annotation_version):
        return cls._get_transcript_for_genes([gene], variant_annotation_version)

    @classmethod
    def get_transcript_for_gene_symbol(cls, gene_symbol: GeneSymbol, variant_annotation_version):
        genes_qs = gene_symbol.get_genes().filter(annotation_consortium=variant_annotation_version.annotation_consortium)
        return cls._get_transcript_for_genes(list(genes_qs), variant_annotation_version)

    @classmethod
    def _get_transcript_for_genes(cls, genes: List[Gene], variant_annotation_version):
        """ Get best transcript (canonical transcript if available), then longest)
            This can be quite slow ie 1-2 seconds """
        pfam_qs = PfamSequenceIdentifier.objects.all()
        for gene in genes:
            if canonical_transcript := gene.get_vep_canonical_transcript(variant_annotation_version):
                canonical_pfam_q = Q(transcript=canonical_transcript) | Q(transcript_version__transcript=canonical_transcript)
                if pfam_qs.filter(canonical_pfam_q).exists():
                    return canonical_transcript

        q_has_pfam = Q(pfamsequenceidentifier__isnull=False) | Q(transcript__pfamsequenceidentifier__isnull=False)
        genome_build = variant_annotation_version.genome_build
        tv_qs = TranscriptVersion.objects.filter(q_has_pfam, genome_build=genome_build, gene_version__gene__in=genes)
        latest_transcript_versions = list(tv_qs.order_by("transcript_id", "-version").distinct("transcript_id"))

        transcript = None
        if latest_transcript_versions:
            longest_transcript_versions = sorted(latest_transcript_versions, key=lambda tv: tv.length or -1, reverse=True)
            transcript = longest_transcript_versions[0].transcript
        return transcript

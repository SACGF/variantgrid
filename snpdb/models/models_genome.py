import itertools
import operator
import os
import re
from functools import cached_property, reduce
from typing import Optional, Iterable

from django.conf import settings
from django.db import models
from django.db.models import QuerySet
from django.db.models.deletion import CASCADE
from django.db.models.query_utils import Q
from django.urls import reverse
from django.utils.text import slugify

from genes.models_enums import AnnotationConsortium
from library.cache import timed_cache
from library.django_utils import SortMetaOrderingMixin
from library.django_utils.django_object_managers import ObjectManagerCachingImmutable
from library.genomics.fasta_wrapper import FastaFileWrapper
from library.preview_request import PreviewModelMixin
from library.utils import invert_dict
from snpdb.genome.fasta_index import load_genome_fasta_index
from snpdb.models.models_enums import SequenceRole, AssemblyMoleculeType


class GenomeBuild(models.Model, SortMetaOrderingMixin, PreviewModelMixin):
    """ There is only 1 patch version of a genome build on the system

        enabled means a build can be loaded via get_name_or_alias()
        Whether there is any annotation (and thus will be used in an environment) is determined by settings, ie:
        settings.ANNOTATION["GRCh38"]["enabled"]

        Build & Contig are populated via migration snpdb 0006 """

    objects = ObjectManagerCachingImmutable()

    name = models.TextField(primary_key=True)
    accession = models.TextField(null=True, blank=True)
    alias = models.TextField(null=True, unique=True, blank=True)
    slug = models.TextField(null=True, unique=True, blank=True)  # We need to slugify T2T-CHMv2.0
    enabled = models.BooleanField(default=True, blank=True)
    igv_genome = models.TextField(null=True, blank=True)

    class Meta:
        ordering = ["name"]
        base_manager_name = 'objects'

    def get_absolute_url(self):
        return reverse("view_genome_build", kwargs={"genome_build_name": self.name})

    def is_version(self, version: int) -> bool:
        return str(version) in self.name

    @classmethod
    def grch37(cls) -> 'GenomeBuild':
        return cls.objects.get(pk='GRCh37')

    @classmethod
    def grch38(cls) -> 'GenomeBuild':
        return cls.objects.get(pk='GRCh38')

    @classmethod
    def t2tv2(cls) -> 'GenomeBuild':
        return cls.objects.get(pk='T2T-CHM13v2.0')

    @classmethod
    def legacy_build(cls) -> 'GenomeBuild':
        """ Use this for hacks - makes it easy to find / fix later """
        return cls.objects.get(pk='GRCh37')

    @cached_property
    def is_annotated(self):
        return GenomeBuild.builds_with_annotation().filter(pk=self.pk).exists()

    @staticmethod
    def get_from_fuzzy_string(genome_build_str: str):
        if "h37" in genome_build_str:
            return GenomeBuild.grch37()
        return GenomeBuild.grch38()

    @staticmethod
    @timed_cache(ttl=60)
    def get_name_or_alias(build_name) -> 'GenomeBuild':
        """ Get by insensitive name or alias """
        q_or_list = [
            Q(name__iexact=build_name),
            Q(alias__iexact=build_name),
            Q(slug=slugify(build_name)),
        ]
        try:
            build_no_patch = GenomeBuildPatchVersion.get_build_no_patch(build_name)
            q_or_list.extend([
                Q(name__iexact=build_no_patch),
                Q(alias__iexact=build_no_patch),
                Q(slug=slugify(build_no_patch)),
            ])
        except ValueError:
            pass
        q = reduce(operator.or_, q_or_list)
        return GenomeBuild.objects.get(q, enabled=True)

    @staticmethod
    def detect_from_filename(filename) -> Optional['GenomeBuild']:
        """ Attempt to detect from filename
            e.g. "foo_grch37.bed" (returns GRCh37) - must only contain name/alias from 1 build """

        lc_filename = filename.lower()
        possible_builds = set()
        for genome_build in GenomeBuild.builds_with_annotation():
            for n in [genome_build.name, genome_build.alias]:
                if n.lower() in lc_filename:
                    possible_builds.add(genome_build)
        if len(possible_builds) == 1:
            return possible_builds.pop()
        return None

    @staticmethod
    def builds_with_annotation() -> QuerySet['GenomeBuild']:
        enabled_annotation = []
        for build_name, values in settings.ANNOTATION.items():
            if values.get("enabled"):
                enabled_annotation.append(build_name)
        return GenomeBuild.objects.filter(enabled=True, name__in=enabled_annotation).order_by('name')

    @staticmethod
    @timed_cache(ttl=60)
    def builds_with_annotation_cached() -> list['GenomeBuild']:
        return list(GenomeBuild.builds_with_annotation())

    @staticmethod
    def builds_with_annotation_priority(priority: 'GenomeBuild') -> list['GenomeBuild']:
        return [priority] + list(GenomeBuild.builds_with_annotation().exclude(pk=priority.pk).all())

    @staticmethod
    def available_names_or_aliases():
        values_qs = GenomeBuild.builds_with_annotation().values_list("name", "alias")
        return ", ".join(itertools.chain.from_iterable(values_qs))

    @cached_property
    def contigs(self):
        qs = Contig.objects.filter(genomebuildcontig__genome_build=self)
        return qs.order_by("genomebuildcontig__order")

    @cached_property
    def standard_contigs(self):
        return self.contigs.filter(role=SequenceRole.ASSEMBLED_MOLECULE)

    @cached_property
    def mitochondria_accession(self) -> Optional[str]:
        mito_ac = None
        if mito := self.contigs.filter(molecule_type=AssemblyMoleculeType.MITOCHONDRION).first():
            mito_ac = mito.refseq_accession
        return mito_ac

    @property
    def genome_builds(self):
        """ For PreviewModelMixin - search knows it's a genome build """
        return [self]

    @staticmethod
    def get_choices():
        choices = []
        for genome_build in GenomeBuild.builds_with_annotation():
            choices.append((genome_build.pk, str(genome_build)))
        return choices

    @staticmethod
    def get_build_annotation_descriptions():
        build_consortia = {}
        for genome_build in GenomeBuild.builds_with_annotation():
            build_consortia[genome_build.name] = genome_build.settings["annotation_consortium"]

        unique_values = set(build_consortia.values())
        if len(unique_values) == 1:
            only_build = unique_values.pop()
            if len(build_consortia) == 1:
                build_annotation_descriptions = only_build
            else:
                build_annotation_descriptions = f"{only_build} for all builds."
        else:
            build_annotation_descriptions = ", ".join([f"{k}: {v}" for k, v in build_consortia.items()])
        return build_annotation_descriptions

    def _get_chrom_contig_mappings(self, contigs: Iterable['Contig']) -> dict[str, 'Contig']:
        chrom_contig_mappings = {}
        for contig in contigs:
            chrom_contig_mappings[contig.name] = contig
            chrom_contig_mappings[contig.ucsc_name] = contig
            chrom_contig_mappings[contig.genbank_accession] = contig
            chrom_contig_mappings[contig.refseq_accession] = contig

        # Map lowercase "mt" -> "MT"
        if mt := chrom_contig_mappings.get("MT"):
            chrom_contig_mappings["mt"] = mt
        return chrom_contig_mappings

    @cached_property
    def chrom_contig_mappings(self) -> dict[str, 'Contig']:
        return self._get_chrom_contig_mappings(self.contigs)

    @cached_property
    def chrom_standard_contig_mappings(self) -> dict[str, 'Contig']:
        return self._get_chrom_contig_mappings(self.standard_contigs)

    def get_chrom_contig_id_mappings(self, standard_contigs_only=False) -> dict[str, int]:
        if standard_contigs_only:
            chrom_contig_mappings = self.chrom_standard_contig_mappings
        else:
            chrom_contig_mappings = self.chrom_contig_mappings
        return {k: v.pk for k, v in chrom_contig_mappings.items()}

    def convert_chrom_to_contig_accession(self, chrom: str) -> str:
        """ chrom = ucsc_name/genbank_accession/refseq accession """
        contig = self.chrom_contig_mappings[chrom]
        return contig.refseq_accession

    @property
    def settings(self):
        return settings.ANNOTATION[self.name]

    def get_annotation_consortium_display(self):
        return self.settings["annotation_consortium"]

    @property
    def annotation_consortium(self):
        consortia_dict = invert_dict(dict(AnnotationConsortium.choices))
        try:
            ac_str = self.settings["annotation_consortium"]
            return consortia_dict[ac_str]
        except Exception as e:
            choices = ",".join(consortia_dict)
            msg = f"Annotation settings 'vep_config.annotation_consortium' must be present and one of {choices}: {e}"
            raise ValueError(msg)

    @property
    def reference_fasta(self):
        return self.get_settings_file("reference_fasta")

    @property
    def genome_fasta(self):
        """ This can't be cached as need to be there for unit tests """
        return GenomeFasta.get_for_genome_build(self)

    @property
    def reference_fasta_has_chr(self):
        """ Does reference contigs start with chr? """

        SETTING_NAME = "reference_fasta_has_chr"
        has_chr = self.settings.get(SETTING_NAME)
        if has_chr is not None:
            return has_chr

        # attempt to work it out...
        if self.name == 'GRCh38':
            return False

        ac = self.settings.get("annotation_consortium")
        if ac == 'Ensembl':
            return False
        if ac == 'RefSeq':
            if self.name in ['hg19', 'GRCh37']:
                return True
        raise ValueError(f"Please set settings.ANNOTATION[{self.name}]{SETTING_NAME}")

    def is_equivalent(self, other: 'GenomeBuild') -> bool:
        if self.pk == other.pk:
            return True
        if self.name == other.alias or self.alias == other.name:
            return True
        return False

    def get_settings_file(self, name, mandatory=True):
        filename = self.settings.get(name)
        if filename:
            if not os.path.exists(filename):
                msg = f"Can't access '{filename}' (from build '{self.name}' setting '{name}')"
                raise FileNotFoundError(msg)
        elif mandatory:
            msg = f"No entry for mandatory settings.ANNOTATIONS['{self.name}']['{name}']"
            raise KeyError(msg)

        return self.settings.get(name)

    def get_build_with_patch(self, annotation_version=None) -> str:
        """ annotation_version: defaults to latest """
        if annotation_version is None:
            annotation_version = self.annotationversion_set.order_by("pk").last()
            if annotation_version is None:
                raise ValueError(f"{self} has no annotation versions")

        return annotation_version.variant_annotation_version.assembly

    @staticmethod
    def get_known_builds_comma_separated_string() -> str:
        return ", ".join(GenomeBuild.objects.all().order_by("name").values_list("name", flat=True))

    def __str__(self):
        return self.name


class GenomeBuildPatchVersion(models.Model):
    """
    Stores a GenomeBuild along with a patch version.
    At the time of this comment, we don't do anything with different patch versions, but use this class to future proof
    """

    name = models.TextField(primary_key=True)
    """ The patch version as it's primarily described e.g. GRCh37.p13 """

    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    """ The patchless genome build """

    patch_version = models.IntegerField(blank=True, null=True)
    """ The version of the patch, or None if the version is unknown, note this is different to an explicit patch_version of 0 """

    # Standard eg GRCh37.p13 (build=GRCh37, patch=13) or GRCh38.p14 (build=GRCh38, patch=14
    GRCH_WITH_PATCH_REGEX = re.compile(r"(?P<genome_build>GRCh[0-9]+)(?:[.]p(?P<patch_version>[0-9]+))?", re.IGNORECASE)
    GENOME_BUILD_VERSION_REGEXES = [
        GRCH_WITH_PATCH_REGEX,
        # T2T-CHM13v2.0 is ALL BUILD, v2.0 is not a patch
        re.compile(r"(?P<genome_build>T2T-CHM[0-9]+v[0-9]+\.[0-9]+)", re.IGNORECASE),
    ]

    class Meta:
        unique_together = ('genome_build', 'patch_version')

    def __str__(self):
        return self.name

    def __lt__(self, other: 'GenomeBuildPatchVersion'):
        def sort_key(obj: GenomeBuildPatchVersion):
            return obj.genome_build.name, obj.patch_version or 0
        return sort_key(self) < sort_key(other)

    @staticmethod
    def get_build_no_patch(build_name: str) -> str:
        for regex in GenomeBuildPatchVersion.GENOME_BUILD_VERSION_REGEXES:
            if match := regex.match(build_name):
                return match.group("genome_build")
        raise ValueError(f"Invalid Genome Build (with optional patch) '{build_name}'")

    @staticmethod
    def get_or_create(name: str) -> 'GenomeBuildPatchVersion':
        for regex in GenomeBuildPatchVersion.GENOME_BUILD_VERSION_REGEXES:
            if match := regex.match(name):
                genome_build = GenomeBuild.get_name_or_alias(match.group('genome_build'))
                patch_version: Optional[int] = None
                normalized_name: str
                if patch_version_str := match.groupdict().get('patch_version'):
                    patch_version = int(patch_version_str)
                    normalized_name = f"{genome_build.name}.p{patch_version}"
                else:
                    normalized_name = genome_build.name

                gbpv, _ = GenomeBuildPatchVersion.objects.get_or_create(
                    name=normalized_name,
                    genome_build=genome_build,
                    patch_version=patch_version
                )
                return gbpv

        raise ValueError(f"Invalid Genome Build Patch Version '{name}'")

    @staticmethod
    def get_unspecified_patch_version_for(genome_build: GenomeBuild) -> 'GenomeBuildPatchVersion':
        """
        Returns a GenomePatchVersion object where the patch is not known (represented by None)
        """
        gbpv, _ = GenomeBuildPatchVersion.objects.get_or_create(
            name=genome_build.name,
            genome_build=genome_build,
            patch_version=None
        )
        return gbpv


class Contig(models.Model, PreviewModelMixin):
    """ Contigs don't directly link to genome build so builds can share contigs
        e.g. hg19 and GRCh37 are the same except for mitochondrion
        @see annotation.reference_contigs.get_assembly_report_df """

    name = models.TextField()
    role = models.CharField(max_length=3, choices=SequenceRole.choices)
    assigned_molecule = models.TextField(null=True, blank=True)  # unlocalised may know what chrom it's from
    molecule_type = models.CharField(max_length=1, choices=AssemblyMoleculeType.choices, null=True, blank=True)
    genbank_accession = models.TextField(null=True)
    refseq_accession = models.TextField(unique=True)
    ucsc_name = models.TextField(null=True)
    length = models.IntegerField()

    class ContigNotInBuildError(ValueError):
        pass

    def get_absolute_url(self):
        return reverse("view_contig", kwargs={"contig_accession": self.refseq_accession})

    def get_genome_builds(self, require_annotation=True) -> Iterable[GenomeBuild]:
        if require_annotation:
            builds = GenomeBuild.builds_with_annotation()
        else:
            builds = GenomeBuild.objects.all()
        return builds.filter(genomebuildcontig__contig=self, enabled=True).order_by("name")

    @property
    def genome_builds(self):
        return self.get_genome_builds(require_annotation=True)

    @staticmethod
    def get_q(accession: str, case_sensitive=True) -> Q:
        q_list = []
        for field in ["genbank_accession", "refseq_accession", "name", "ucsc_name"]:
            if case_sensitive:
                f = field
            else:
                f = f"{field}__iexact"
            q_list.append(Q(**{f: accession}))

        # MT doesn't have ucsc_name, so need special case to load M with chr prefix
        if accession.lower() in ("chrm", "chrmt"):
            q_list.append(Q(name="MT"))

        return reduce(operator.or_, q_list)

    def __str__(self):
        return f"{self.name} - {self.refseq_accession} ({self.length})"


class GenomeBuildContig(models.Model):
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    contig = models.ForeignKey(Contig, on_delete=CASCADE)
    order = models.IntegerField()

    class Meta:
        unique_together = ("genome_build", "contig")

    def __str__(self):
        return f"{self.contig} ({self.genome_build})"


class GenomeFasta(models.Model):
    """ Incoming VCF chroms are mapped to contigs from VCF header then converted to fasta names (for VT) """
    filename = models.TextField(unique=True)
    index_filename = models.TextField()
    index_sha256sum = models.TextField()
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)

    class ContigNotInFastaError(ValueError):
        pass

    class Meta:
        unique_together = ("genome_build", "annotation_consortium")

    def get_contig_id_to_name_mappings(self):
        contig_id_to_name_mappings = dict(self.genomefastacontig_set.all().values_list("contig_id", "name"))
        if not contig_id_to_name_mappings:
            msg = f"Contig/Fasta names empty for {self}"
            raise ValueError(msg)
        return contig_id_to_name_mappings

    def __str__(self):
        return f"{self.genome_build} ({self.get_annotation_consortium_display()})"

    @staticmethod
    def get_for_genome_build(genome_build: GenomeBuild):
        genome_fasta, created = GenomeFasta.objects.get_or_create(filename=genome_build.reference_fasta,
                                                                  genome_build=genome_build)
        if created:
            load_genome_fasta_index(genome_fasta, genome_build)
        return genome_fasta

    def convert_chrom_to_fasta_sequence(self, chrom: str) -> str:
        chrom_contig_id, contig_id_to_name = self._contig_lookups
        contig_id = chrom_contig_id.get(chrom)
        if contig_id is None:
            raise Contig.ContigNotInBuildError(f"'{chrom}' not in {self.genome_build}")

        sequence_name = contig_id_to_name.get(contig_id)
        if sequence_name is None:
            raise self.ContigNotInFastaError(f"'{chrom}' (contig_id: {contig_id}) not in {self} ({self.filename})")
        return sequence_name

    @cached_property
    def _contig_lookups(self):
        chrom_contig_id = self.genome_build.get_chrom_contig_id_mappings()
        contig_id_to_name = self.get_contig_id_to_name_mappings()
        return chrom_contig_id, contig_id_to_name

    @property
    def fasta(self) -> FastaFileWrapper:
        ffw = FastaFileWrapper(self.filename, convert_chrom_func=self.convert_chrom_to_fasta_sequence)
        contig_1 = self.convert_chrom_to_fasta_sequence("1")  # Should be in all human genome builds
        assert contig_1 in ffw, f"Genome fasta '{self.filename}' does not contain contig '{contig_1}'"
        return ffw


class GenomeFastaContig(models.Model):
    genome_fasta = models.ForeignKey(GenomeFasta, on_delete=CASCADE)
    name = models.TextField()
    length = models.IntegerField()
    contig = models.ForeignKey(Contig, on_delete=CASCADE)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['genome_fasta', 'name'], name='one_name_per_build'),
            models.UniqueConstraint(fields=['genome_fasta', 'contig'], name='one_contig_per_build'),
        ]

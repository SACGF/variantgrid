from django.conf import settings
from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query_utils import Q
import itertools
import os
from typing import Dict, Optional, List

from lazy import lazy

from genes.models_enums import AnnotationConsortium
from library.cache import timed_cache
from library.django_utils import SortMetaOrderingMixin
from library.genomics.fasta_wrapper import FastaFileWrapper
from library.utils import invert_dict
from snpdb.genome.fasta_index import load_genome_fasta_index
from snpdb.models.models_enums import SequenceRole, AssemblyMoleculeType


class GenomeBuild(models.Model, SortMetaOrderingMixin):
    """ There is only 1 patch version of a genome build on the system

        enabled means a build can be loaded via get_name_or_alias()
        Whether there is any annotation (and thus will be used in an environment) is determined by settings, ie:
        settings.ANNOTATION["GRCh38"]["enabled"]

        Build & Contig are populated via migration snpdb 0006 """

    name = models.TextField(primary_key=True)
    accession = models.TextField(null=True)
    alias = models.TextField(null=True, unique=True)
    enabled = models.BooleanField(default=True)
    igv_genome = models.TextField(null=True, blank=True)

    class Meta:
        ordering = ["name"]

    @classmethod
    @timed_cache(ttl=60)
    def grch37(cls) -> 'GenomeBuild':
        return cls.objects.get(pk='GRCh37')

    @classmethod
    @timed_cache(ttl=60)
    def grch38(cls) -> 'GenomeBuild':
        return cls.objects.get(pk='GRCh38')

    @classmethod
    @timed_cache(ttl=60)
    def legacy_build(cls) -> 'GenomeBuild':
        """ Use this for hacks - makes it easy to find / fix later """
        return cls.objects.get(pk='GRCh37')

    @lazy
    def is_annotated(self):
        return GenomeBuild.builds_with_annotation().filter(pk=self.pk).exists()

    @staticmethod
    @timed_cache(ttl=60)
    def default_build():
        """ Used for initialising user settings or worst case fallback """
        return GenomeBuild.objects.get(name=settings.DEFAULT_BUILD)

    @staticmethod
    @timed_cache(ttl=60)
    def get_name_or_alias(build_name):
        """ Get by insensitive name or alias """

        build_no_patch = build_name.split(".", 1)[0]
        q = Q(name__iexact=build_name) | Q(alias__iexact=build_name)
        q_no_patch = Q(name__iexact=build_no_patch) | Q(alias__iexact=build_no_patch)
        return GenomeBuild.objects.get(q | q_no_patch, enabled=True)

    @staticmethod
    def detect_from_filename(filename) -> Optional['GenomeBuild']:
        """ Attempt to detect from filename
            eg "foo_grch37.bed" (returns GRCh37) - must only contain name/alias from 1 build) """

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
    def builds_with_annotation():
        enabled_annotation = []
        for build_name, values in settings.ANNOTATION.items():
            if values.get("enabled"):
                enabled_annotation.append(build_name)
        return GenomeBuild.objects.filter(name__in=enabled_annotation)

    @staticmethod
    def builds_with_annotation_priority(priority: 'GenomeBuild') -> List['GenomeBuild']:
        return [priority] + list(GenomeBuild.builds_with_annotation().exclude(pk=priority.pk).all())

    @staticmethod
    def available_names_or_aliases():
        values_qs = GenomeBuild.builds_with_annotation().values_list("name", "alias")
        return ", ".join(itertools.chain.from_iterable(values_qs))

    @property
    def contigs(self):
        qs = Contig.objects.filter(genomebuildcontig__genome_build=self)
        return qs.order_by("genomebuildcontig__order")

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

    @lazy
    def chrom_contig_mappings(self) -> Dict[str, 'Contig']:
        chrom_contig_mappings = {}
        for contig in self.contigs:
            chrom_contig_mappings[contig.name] = contig
            chrom_contig_mappings[contig.ucsc_name] = contig
            chrom_contig_mappings[contig.genbank_accession] = contig
            chrom_contig_mappings[contig.refseq_accession] = contig
        # Map lowercase "mt" -> "MT"
        chrom_contig_mappings["mt"] = chrom_contig_mappings["MT"]
        return chrom_contig_mappings

    def get_chrom_contig_id_mappings(self) -> Dict[str, int]:
        return {k: v.pk for k, v in self.chrom_contig_mappings.items()}

    def convert_chrom_to_contig_name(self, chrom: str) -> str:
        """ chrom = ucsc_name/genbank_accession/refseq accession """
        contig = self.chrom_contig_mappings[chrom]
        return contig.name

    @property
    def settings(self):
        return settings.ANNOTATION[self.name]

    def get_annotation_consortium_display(self):
        return self.settings["annotation_consortium"]

    @property
    def annotation_consortium(self):
        consortia_dict = invert_dict(dict(AnnotationConsortium.CHOICES))
        try:
            ac_str = self.settings["annotation_consortium"]
            return consortia_dict[ac_str]
        except Exception as e:
            choices = ",".join(consortia_dict)
            msg = f"Annotation settings 'vep_config.annotation_consortium' must be present and one of {choices}: {e}"
            raise ValueError(msg)

    @lazy
    def latest_variant_annotation_version(self):
        return self.variantannotationversion_set.order_by("annotation_date").last()

    @property
    def reference_fasta(self):
        return self.get_settings_file("reference_fasta")

    @lazy
    def genome_fasta(self):
        from snpdb.models import GenomeFasta
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

        return annotation_version.variant_annotation_version.assembly

    def __str__(self):
        name = self.name
        if self.alias:
            name += f" (aka {self.alias})"
        return name


class Contig(models.Model):
    """ Contigs don't directly link to genome build so builds can share contigs
        eg hg19 and GRCh37 are the same except for mitochondrion
        @see annotation.reference_contigs.get_assembly_report_df """

    name = models.TextField()
    role = models.CharField(max_length=3, choices=SequenceRole.CHOICES)
    assigned_molecule = models.TextField(null=True, blank=True)  # unlocalised may know what chrom it's from
    molecule_type = models.CharField(max_length=1, choices=AssemblyMoleculeType.CHOICES, null=True, blank=True)
    genbank_accession = models.TextField(null=True)
    refseq_accession = models.TextField(unique=True)
    ucsc_name = models.TextField(null=True)
    length = models.IntegerField()

    class ContigNotInBuildError(ValueError):
        pass

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
    index_md5sum = models.TextField()
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.CHOICES)

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

    @lazy
    def _contig_lookups(self):
        chrom_contig_id = self.genome_build.get_chrom_contig_id_mappings()
        contig_id_to_name = self.get_contig_id_to_name_mappings()
        return chrom_contig_id, contig_id_to_name

    @property
    def fasta(self):
        return FastaFileWrapper(self.filename, convert_chrom_func=self.convert_chrom_to_fasta_sequence)


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

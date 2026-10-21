"""
Gene-level annotation: GeneAnnotationVersion (partitioned GeneAnnotation rows per GeneAnnotationRelease, built by
`manage.py gene_annotation`), the dbNSFP gene table it links (DBNSFPGeneAnnotationVersion, list-partitioned
DBNSFPGeneAnnotation) and Human Protein Atlas tissue expression (HumanProteinAtlasAnnotationVersion).
"""

import logging
from typing import Optional

from django.db import connection, models
from django.db.models.deletion import CASCADE, DO_NOTHING, PROTECT, SET_NULL
from django.db.models.signals import pre_delete
from django.dispatch.dispatcher import receiver
from model_utils.models import TimeStampedModel
from psqlextra.models import PostgresPartitionedModel
from psqlextra.types import PostgresPartitioningMethod

from annotation.enums import (
    EssentialGeneCRISPR,
    EssentialGeneCRISPR2,
    EssentialGeneGeneTrap,
    HumanProteinAtlasAbundance,
)
from annotation.models.models_version import SubVersionPartition
from genes.models import (
    Gene,
    GeneAnnotationRelease,
    GeneSymbol,
    Transcript,
)
from library.django_utils.data_archive_mixin import DataArchiveMixin
from ontology.models import OntologyVersion


class DBNSFPGeneAnnotationVersion(TimeStampedModel):
    """ @see https://sites.google.com/site/jpopgen/dbNSFP
        This isn't updated every release, so can have same hash across diff versions """
    version = models.TextField(primary_key=True)
    sha256_hash = models.TextField()

    class Meta:
        unique_together = ('version', 'sha256_hash')

    def save(self, *args, **kwargs):
        created = not self.pk
        super().save(*args, **kwargs)
        if created:
            logging.info("Creating new DBNSFPGeneAnnotation partition")
            version = self.pk
            connection.schema_editor().add_list_partition(
                model=DBNSFPGeneAnnotation,
                name=f"version_{version}",
                values=[version],
            )

    @staticmethod
    def latest() -> Optional['DBNSFPGeneAnnotationVersion']:
        return DBNSFPGeneAnnotationVersion.objects.order_by("created").last()


class DBNSFPGeneAnnotation(PostgresPartitionedModel, TimeStampedModel):
    """ @see https://sites.google.com/site/jpopgen/dbNSFP """
    version = models.ForeignKey(DBNSFPGeneAnnotationVersion, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    refseq_transcript = models.ForeignKey(Transcript, related_name="refseq_dbnsfp_gene",
                                          null=True, blank=True, on_delete=CASCADE)
    ensembl_transcript = models.ForeignKey(Transcript, related_name="ensembl_dbnsfp_gene",
                                           null=True, blank=True, on_delete=CASCADE)
    gene_damage_index_score = models.FloatField(null=True)
    gene_damage_index_phred = models.FloatField(null=True)
    phi = models.FloatField(null=True)
    ghis = models.FloatField(null=True)
    prec = models.FloatField(null=True)
    hipred_score = models.FloatField(null=True)
    gnomad_pli = models.FloatField(null=True, blank=True)
    gnomad_prec = models.FloatField(null=True, blank=True)
    gnomad_pnull = models.FloatField(null=True, blank=True)
    loftool = models.FloatField(null=True, blank=True)
    gene_indispensability_score = models.FloatField(null=True, blank=True)
    hipred_prediction = models.BooleanField(null=True)
    gene_indispensability_pred = models.BooleanField(null=True)
    pathway_biocarta_full = models.TextField(null=True, blank=True)
    pathway_consensus_pathdb = models.TextField(null=True, blank=True)
    pathway_kegg_id = models.TextField(null=True, blank=True)
    pathway_kegg_full = models.TextField(null=True, blank=True)
    gwas_trait_association = models.TextField(null=True, blank=True)
    go_biological_process = models.TextField(null=True, blank=True)
    go_cellular_component = models.TextField(null=True, blank=True)
    go_molecular_function = models.TextField(null=True, blank=True)
    interactions_biogrid = models.TextField(null=True, blank=True)
    interactions_consensus_pathdb = models.TextField(null=True, blank=True)
    expression_egenetics = models.TextField(null=True, blank=True)
    expression_gnf_atlas = models.TextField(null=True, blank=True)
    essential_gene_crispr = models.CharField(max_length=1, null=True, blank=True,
                                             choices=EssentialGeneCRISPR.choices)
    essential_gene_crispr2 = models.CharField(max_length=1, null=True, blank=True,
                                              choices=EssentialGeneCRISPR2.choices)
    essential_gene_gene_trap = models.CharField(max_length=1, null=True, blank=True,
                                                choices=EssentialGeneGeneTrap.choices)

    def get_gene_indispensability_pred_display(self) -> str:
        if self.gene_indispensability_pred is None:
            return "n/a"
        elif self.gene_indispensability_pred:
            return "Essential"
        return "LoF tolerant"

    class PartitioningMeta:
        method = PostgresPartitioningMethod.LIST
        key = ["version_id"]


class GeneAnnotationVersion(DataArchiveMixin, SubVersionPartition):
    RECORDS_BASE_TABLE_NAMES = ["annotation_geneannotation"]
    gene_annotation_release = models.ForeignKey(GeneAnnotationRelease, on_delete=CASCADE)
    ontology_version = models.ForeignKey(OntologyVersion, null=True, on_delete=PROTECT)
    dbnsfp_gene_version = models.ForeignKey(DBNSFPGeneAnnotationVersion, null=True, on_delete=CASCADE)
    gnomad_import_date = models.DateTimeField()

    @property
    def genome_build(self):
        return self.gene_annotation_release.genome_build

    def __str__(self):
        return f"GAR: {self.gene_annotation_release}, Ont: {self.ontology_version}, DBnsFP: {self.dbnsfp_gene_version}"


@receiver(pre_delete, sender=GeneAnnotationVersion)
def gene_annotation_version_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    instance.delete_related_objects()


class GeneAnnotation(models.Model):
    """ This is generated against genes found via GeneAnnotationRelease
        so that data matches up in analyses """
    # DO_NOTHING as records are deleted by dropping the partition (a CASCADE locks every other version's)
    version = models.ForeignKey(GeneAnnotationVersion, on_delete=DO_NOTHING)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    dbnsfp_gene = models.ForeignKey(DBNSFPGeneAnnotation, null=True, on_delete=SET_NULL)
    hpo_terms = models.TextField(null=True)
    omim_terms = models.TextField(null=True)
    mondo_terms = models.TextField(null=True)
    gene_disease_moderate_or_above = models.TextField(null=True)
    gene_disease_supportive_or_below = models.TextField(null=True)
    gnomad_oe_lof = models.FloatField(null=True)  # Copied from genes.GnomADGeneConstraint

    class Meta:
        unique_together = ("version", "gene")

    @staticmethod
    def get_for_gene_version(gene_version):
        pass


class HumanProteinAtlasAnnotationVersion(DataArchiveMixin, SubVersionPartition):
    RECORDS_BASE_TABLE_NAMES = ["annotation_humanproteinatlasannotation"]
    filename = models.TextField()
    sha256_hash = models.TextField()
    hpa_version = models.FloatField()
    unit = models.TextField()  # What unit "value" is in

    def get_minimum_for_abundance_level(self, abundance: HumanProteinAtlasAbundance):
        """ Attempt to remain backwards compatibile - see how this is calculated here:
            https://github.com/SACGF/variantgrid/issues/9#issuecomment-563510678 """
        abundance_mins = {
            HumanProteinAtlasAbundance.NOT_DETECTED: 0,
            HumanProteinAtlasAbundance.LOW: 1.1,
            HumanProteinAtlasAbundance.MEDIUM: 15.5,
            HumanProteinAtlasAbundance.HIGH: 83.6,
        }

        if self.hpa_version == 15:
            abundance_mins = {
                HumanProteinAtlasAbundance.NOT_DETECTED: 0,
                HumanProteinAtlasAbundance.LOW: 0.5,
                HumanProteinAtlasAbundance.MEDIUM: 10,
                HumanProteinAtlasAbundance.HIGH: 50,
            }
        return abundance_mins[abundance]


@receiver(pre_delete, sender=HumanProteinAtlasAnnotationVersion)
def human_protein_annotation_version_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    instance.delete_related_objects()


class HumanProteinAtlasTissueSample(models.Model):
    name = models.TextField(unique=True)

    def __str__(self):
        return self.name


class HumanProteinAtlasAnnotation(models.Model):
    # DO_NOTHING as records are deleted by dropping the partition (a CASCADE locks every other version's)
    version = models.ForeignKey(HumanProteinAtlasAnnotationVersion, on_delete=DO_NOTHING)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, null=True, on_delete=CASCADE)  # Always Ensembl
    tissue_sample = models.ForeignKey(HumanProteinAtlasTissueSample, on_delete=CASCADE)
    value = models.FloatField()

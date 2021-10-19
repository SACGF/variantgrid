import os

from django.conf import settings
from django.db import models
from django.db.models import CASCADE
from django.db.models.signals import post_delete
from django.dispatch import receiver

from analysis.models.models_variant_tag import VariantTagsImport
from annotation.models.models import ClinVarCitationsCollection, ManualVariantEntryCollection, ClinVarVersion
from classification.models import ClassificationImport
from expression.models import CuffDiffFile
from genes.models import GeneList, GeneCoverageCollection
from library.file_utils import name_from_filename
from patients.models import PatientRecords
from pedigree.models import PedFile
from snpdb.models import GenomicIntervalsCollection, AnnotationLevel, ImportStatus, Sample
from snpdb.models.models_variant import Liftover
from upload.bed_file_processing import process_bed_file
from upload.models import UploadedFile
from upload.models.models_enums import ExpressionType


class UploadedGeneList(models.Model):
    uploaded_file = models.OneToOneField(UploadedFile, on_delete=CASCADE)
    gene_list = models.OneToOneField(GeneList, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.gene_list


class UploadedExpressionFile(models.Model):
    uploaded_file = models.ForeignKey(UploadedFile, on_delete=CASCADE)
    format = models.CharField(max_length=1, choices=ExpressionType.choices)
    annotation_level = models.CharField(max_length=1, choices=AnnotationLevel.choices)
    cuff_diff_file = models.OneToOneField(CuffDiffFile, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.cuff_diff_file


class UploadedBed(models.Model):
    uploaded_file = models.ForeignKey(UploadedFile, on_delete=CASCADE)
    genomic_intervals_collection = models.OneToOneField(GenomicIntervalsCollection, null=True, on_delete=CASCADE)

    def process_bed_file(self):
        bed_file = self.uploaded_file.get_filename()
        has_chr = self.genomic_intervals_collection.genome_build.reference_fasta_has_chr
        if has_chr:
            chrom_description = "has_chr"
        else:
            chrom_description = "no_chr"
        name = name_from_filename(bed_file)
        processed_base_name = f"{self.genomic_intervals_collection.pk}.{name}.{chrom_description}.processed.bed"
        processed_file = os.path.join(settings.PROCESSED_BED_FILES_DIR, processed_base_name)
        if not os.path.exists(processed_file):
            num_records = process_bed_file(bed_file, processed_file, has_chr)
            self.genomic_intervals_collection.processed_records = num_records

        self.genomic_intervals_collection.processed_file = processed_file
        self.genomic_intervals_collection.import_status = ImportStatus.SUCCESS
        self.genomic_intervals_collection.save()

    def get_data(self):
        return self.genomic_intervals_collection


class UploadedPedFile(models.Model):
    uploaded_file = models.ForeignKey(UploadedFile, on_delete=CASCADE)
    ped_file = models.OneToOneField(PedFile, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.ped_file


class UploadedPatientRecords(models.Model):
    uploaded_file = models.ForeignKey(UploadedFile, on_delete=CASCADE)
    patient_records = models.OneToOneField(PatientRecords, on_delete=CASCADE)

    def get_data(self):
        return self.patient_records


class UploadedGeneCoverage(models.Model):
    uploaded_file = models.OneToOneField(UploadedFile, on_delete=CASCADE)
    gene_coverage_collection = models.OneToOneField(GeneCoverageCollection, null=True, on_delete=CASCADE)
    sample = models.OneToOneField(Sample, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.gene_coverage_collection


class UploadedManualVariantEntryCollection(models.Model):
    uploaded_file = models.OneToOneField(UploadedFile, on_delete=CASCADE)
    collection = models.OneToOneField(ManualVariantEntryCollection, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.collection


class UploadedClassificationImport(models.Model):
    uploaded_file = models.OneToOneField(UploadedFile, on_delete=CASCADE)
    classification_import = models.OneToOneField(ClassificationImport, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.classification_import


class UploadedLiftover(models.Model):
    uploaded_file = models.OneToOneField(UploadedFile, on_delete=CASCADE)
    liftover = models.OneToOneField(Liftover, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.liftover


@receiver(post_delete, sender=UploadedGeneCoverage)
def uploaded_gene_coverage_post_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    if instance.gene_coverage_collection:
        instance.gene_coverage_collection.delete()


class UploadedClinVarVersion(models.Model):
    uploaded_file = models.OneToOneField(UploadedFile, on_delete=CASCADE)
    clinvar_version = models.OneToOneField(ClinVarVersion, null=True, on_delete=CASCADE)


class UploadedVariantTags(models.Model):
    uploaded_file = models.OneToOneField(UploadedFile, on_delete=CASCADE)
    variant_tags_import = models.OneToOneField(VariantTagsImport, on_delete=CASCADE)

    def get_data(self):
        return self.variant_tags_import

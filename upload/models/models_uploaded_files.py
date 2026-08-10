import logging
import os
from typing import Optional

from django.conf import settings
from django.db import models
from django.db.models import CASCADE
from django.db.models.signals import post_delete
from django.dispatch import receiver

from analysis.models.models_variant_tag import VariantTagsImport
from annotation.models.models import ClinVarVersion, ManualVariantEntryCollection
from classification.models import ClassificationImport
from genes.models import GeneCoverageCollection, GeneList
from library.utils.file_utils import name_from_filename
from patients.models import PatientRecords
from pedigree.models import PedFile
from snpdb.models import (
    GenomeBuild,
    GenomicIntervalsCollection,
    ImportedWikiCollection,
    ImportStatus,
    Sample,
)
from snpdb.models.models_variant import LiftoverRun
from upload.bed_file_processing import process_bed_file
from upload.models import FileUpload, UploadData


class UploadedGeneList(UploadData):
    file_upload = models.OneToOneField(FileUpload, on_delete=CASCADE)
    gene_list = models.OneToOneField(GeneList, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.gene_list


class UploadedBed(UploadData):
    file_upload = models.ForeignKey(FileUpload, on_delete=CASCADE)
    genomic_intervals_collection = models.OneToOneField(GenomicIntervalsCollection, null=True, on_delete=CASCADE)

    def process_bed_file(self):
        bed_file = self.file_upload.get_filename()
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


class UploadedPedFile(UploadData):
    file_upload = models.ForeignKey(FileUpload, on_delete=CASCADE)
    ped_file = models.OneToOneField(PedFile, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.ped_file


class UploadedPatientRecords(UploadData):
    file_upload = models.ForeignKey(FileUpload, on_delete=CASCADE)
    patient_records = models.OneToOneField(PatientRecords, on_delete=CASCADE)

    def get_data(self):
        return self.patient_records


@receiver(post_delete, sender=UploadedPatientRecords)
def uploaded_patient_records_post_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    logging.info("UploadedPatientRecords deleted")
    if instance.patient_records:
        logging.info("Deleting linked PatientRecords")
        instance.patient_records.delete()


class UploadedGeneCoverage(UploadData):
    file_upload = models.OneToOneField(FileUpload, on_delete=CASCADE)
    gene_coverage_collection = models.OneToOneField(GeneCoverageCollection, null=True, on_delete=CASCADE)
    sample = models.OneToOneField(Sample, null=True, on_delete=CASCADE)

    def delete(self, *args, **kwargs):
        logging.info("UploadedGeneCoverage.delete()")
        if gcc := self.gene_coverage_collection:
            logging.info("Delegating delete to gene_coverage_collection - will cascade delete")
            ret = gcc.delete()  # Will also cascade delete this
        else:
            ret = super().delete(*args, **kwargs)
        return ret

    def get_data(self):
        return self.gene_coverage_collection


class UploadedManualVariantEntryCollection(UploadData):
    created_by_pipeline = False  # @see annotation.manual_variant_entry.create_manual_variants

    file_upload = models.OneToOneField(FileUpload, on_delete=CASCADE)
    collection = models.OneToOneField(ManualVariantEntryCollection, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.collection

    @property
    def genome_build(self) -> Optional[GenomeBuild]:
        if collection := self.collection:
            return collection.genome_build
        return None


class UploadedClassificationImport(UploadData):
    created_by_pipeline = False  # @see classification.classification_import.process_classification_import

    file_upload = models.OneToOneField(FileUpload, on_delete=CASCADE)
    classification_import = models.OneToOneField(ClassificationImport, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.classification_import

    @property
    def genome_build(self) -> Optional[GenomeBuild]:
        if classification_import := self.classification_import:
            return classification_import.genome_build
        return None


class UploadedLiftover(UploadData):
    created_by_pipeline = False  # @see snpdb.liftover.create_liftover_pipelines

    file_upload = models.OneToOneField(FileUpload, on_delete=CASCADE)
    liftover = models.OneToOneField(LiftoverRun, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.liftover

    @property
    def genome_build(self) -> Optional[GenomeBuild]:
        if liftover := self.liftover:
            return liftover.genome_build
        return None


class UploadedWikiCollection(UploadData):
    file_upload = models.OneToOneField(FileUpload, on_delete=CASCADE)
    wiki_collection = models.OneToOneField(ImportedWikiCollection, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.wiki_collection

    @property
    def genome_build(self) -> Optional[GenomeBuild]:
        # Gene wiki collections have no build - only variant ones do
        if wiki_collection := self.wiki_collection:
            return wiki_collection.genome_build
        return None


class UploadedClinVarVersion(UploadData):
    file_upload = models.OneToOneField(FileUpload, on_delete=CASCADE)
    # It's possible someone could upload the same clinvar VCF again
    clinvar_version = models.ForeignKey(ClinVarVersion, null=True, on_delete=CASCADE)

    def get_data(self):
        return self.clinvar_version

    @property
    def genome_build(self) -> Optional[GenomeBuild]:
        if clinvar_version := self.clinvar_version:
            return clinvar_version.genome_build
        return None


class UploadedVariantTags(UploadData):
    file_upload = models.OneToOneField(FileUpload, on_delete=CASCADE)
    variant_tags_import = models.OneToOneField(VariantTagsImport, on_delete=CASCADE)

    def get_data(self):
        return self.variant_tags_import

    @property
    def genome_build(self) -> Optional[GenomeBuild]:
        return self.variant_tags_import.genome_build

"""
The upload side of sequencing data: BackendVCF links an imported UploadedVCF to the SingleSampleVCF /
JointCalledVCF record it was sent for (seqauto/upload_links.py follows the link), and the UploadData rows that
mark a DRAGEN TSO500 run-level file processed (seqauto/import_task_factories.py). The tables keep their
upload_* names from when these models lived in upload. upload's models are named by string, as upload.models
imports the apps (classification, analysis) that import these models.
"""
from typing import Optional

from django.db import models
from django.db.models import CASCADE

from seqauto.models.models_seqauto import (
    DragenTSO500CombinedVariantOutput,
    JointCalledVCF,
    LibraryQC,
    SingleSampleVCF,
    get_samples_by_sequencing_sample,
)
from seqauto.models.models_software import VariantCaller
from upload.upload_data import UploadData


class BackendVCF(models.Model):
    """ Link between UploadedVCF (upload) and Filesystem VCF (SeqAuto) """
    uploaded_vcf = models.OneToOneField("upload.UploadedVCF", on_delete=CASCADE)
    single_sample_vcf = models.OneToOneField(SingleSampleVCF, null=True, on_delete=CASCADE)
    joint_called_vcf = models.OneToOneField(JointCalledVCF, null=True, on_delete=CASCADE)

    class Meta:
        db_table = "upload_backendvcf"

    @property
    def filesystem_vcf(self):
        vcfs = [self.joint_called_vcf, self.single_sample_vcf]
        if all(vcfs):
            msg = f"{self} has both single_sample_vcf and joint_called_vcf set"
            raise ValueError(msg)
        if self.single_sample_vcf:
            return self.single_sample_vcf
        if self.joint_called_vcf:
            return self.joint_called_vcf
        raise ValueError(f"{self} has neither 'single_sample_vcf' or 'joint_called_vcf' set")

    @property
    def sample_sheet(self):
        return self.filesystem_vcf.sample_sheet

    @property
    def user(self):
        return self.uploaded_vcf.file_upload.user

    @property
    def vcf(self):
        return self.uploaded_vcf.vcf

    @property
    def variant_caller(self) -> VariantCaller:
        record = self.single_sample_vcf or self.joint_called_vcf
        return record.variant_caller

    def get_samples_by_sequencing_sample(self):
        return get_samples_by_sequencing_sample(self.filesystem_vcf.get_sequencing_samples(), self.vcf)

    def __str__(self):
        backend = self.single_sample_vcf or self.joint_called_vcf or 'None'
        return f"BackendVCF: backend: {backend}, uploaded_vcf: {self.uploaded_vcf}"


class UploadedDragenTSO500MetricsOutput(UploadData):
    """ DRAGEN's MetricsOutput.tsv - the per-library QC rows it wrote (@see seqauto.models.LibraryQC).

        The rows are keyed on the run and the file's pair columns rather than on this record, since a
        re-analysis of the run replaces them; this is what makes the upload 'processed', which is what
        the API's sha256 de-duplication keys on """
    file_upload = models.OneToOneField("upload.FileUpload", on_delete=CASCADE)

    class Meta:
        db_table = "upload_uploadeddragentso500metricsoutput"

    def get_data(self):
        """ The first row that found its extraction - a file landing before its specimens are
            accessioned has rows but nothing to show them against yet """
        return LibraryQC.objects.filter(file_upload=self.file_upload,
                                        specimen__isnull=False).order_by("pk").first()

    def get_data_url(self) -> Optional[str]:
        """ The run the file describes, which is every row's; a run seqauto has not registered falls
            back to the first specimen the rows reached """
        if linked := LibraryQC.objects.filter(file_upload=self.file_upload,
                                              sequencing_run__isnull=False).first():
            return linked.sequencing_run.get_absolute_url()
        if library_qc := self.get_data():
            return library_qc.specimen.get_absolute_url()
        return None


class UploadedDragenTSO500CombinedVariantOutput(UploadData):
    """ DRAGEN's CombinedVariantOutput.tsv - the pair's analysis row it wrote
        (@see seqauto.models.DragenTSO500CombinedVariantOutput).

        The row is keyed on (run, pair) rather than on this record, since a re-analysis of the run
        replaces it; this is what makes the upload 'processed', which is what the API's sha256
        de-duplication keys on """
    file_upload = models.OneToOneField("upload.FileUpload", on_delete=CASCADE)

    class Meta:
        db_table = "upload_uploadeddragentso500combinedvariantoutput"

    def get_data(self):
        """ The row this file wrote - None where the file named no run to key one against """
        return DragenTSO500CombinedVariantOutput.objects.filter(file_upload=self.file_upload).first()

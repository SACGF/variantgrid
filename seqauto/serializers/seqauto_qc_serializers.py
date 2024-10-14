from Bio.Align import write
from rest_framework import serializers

from genes.models import CustomTextGeneList, ActiveSampleGeneList
from genes.serializers import SampleGeneListSerializer, GeneCoverageCollectionSerializer
from library.utils import Value
from seqauto.models import IlluminaFlowcellQC, QCGeneList, QC, QCGeneCoverage, QCExecSummary, FastQC, SequencingSample, \
    SampleSheet, SequencingRun, BamFile, VCFFile
from seqauto.serializers.sequencing_serializers import SampleSheetLookupSerializer, FastqSerializer, SeqAutoViewMixin, \
    BamFilePathSerializer, VCFFilePathSerializer, SequencingSampleLookupSerializer
from snpdb.models import DataState


class FastQCSerializer(SeqAutoViewMixin, serializers.ModelSerializer):
    fastq = FastqSerializer()

    class Meta:
        model = FastQC
        fields = "__all__"


class IlluminaFlowcellQCSerializer(SeqAutoViewMixin, serializers.ModelSerializer):
    sample_sheet = SampleSheetLookupSerializer()

    class Meta:
        model = IlluminaFlowcellQC
        #fields = "__all__"
        exclude = ("sequencing_run", )  # Already part of sample_sheet


class QCSerializer(SeqAutoViewMixin, serializers.ModelSerializer):
    # Instead of dealing with all the bam/vcf etc - we'll just deal with sequencing_sample and
    # assume we're using the latest ones associated with that
    sequencing_sample = SequencingSampleLookupSerializer()
    bam_file = BamFilePathSerializer()
    vcf_file = VCFFilePathSerializer()

    class Meta:
        model = QC
        fields = ("sequencing_sample", "bam_file", "vcf_file")

    @staticmethod
    def get_object(data):
        # We are passed "sequencing_sample" - which we can use to get what we really want
        sequencing_sample_data = data.pop("sequencing_sample")
        sheet_data = sequencing_sample_data["sample_sheet"]
        sequencing_run = SequencingRun.objects.get(pk=sheet_data["sequencing_run"])
        sample_sheet = SampleSheet.objects.get(hash=sheet_data["hash"], sequencing_run=sequencing_run)
        sample_name = sequencing_sample_data["sample_name"]
        sequencing_sample = SequencingSample.objects.get(sample_sheet=sample_sheet, sample_name=sample_name)
        bam_file_data = data.pop("bam_file")
        bam_file = BamFile.objects.get(path=bam_file_data["path"],
                                       sequencing_run=sequencing_run,
                                       unaligned_reads__sequencing_sample=sequencing_sample)

        vcf_file_data = data.pop("vcf_file")
        vcf_file = VCFFile.objects.get(path=vcf_file_data["path"],
                                       bam_file=bam_file)

        qc, _ = QC.objects.get_or_create(
            path=data["path"],
            sequencing_run=sequencing_run,
            bam_file=bam_file,
            vcf_file=vcf_file,
        )
        qc.data_state = DataState.COMPLETE
        qc.save()
        return qc


class QCGeneListSerializer(SeqAutoViewMixin, serializers.ModelSerializer):
    """ When we retrieve this, we want to see linked sample gene list """
    qc = QCSerializer()
    sample_gene_list = SampleGeneListSerializer()

    class Meta:
        model = QCGeneList
        fields = ("path", "qc", "sample_gene_list")


class QCGeneListCreateSerializer(SeqAutoViewMixin, serializers.ModelSerializer):
    """ When we create, we just want to send up gene list

        This also handles complexity of setting active gene list
    """
    qc = QCSerializer()
    gene_list = serializers.ListField(
        child=serializers.CharField(),
        write_only=True
    )

    class Meta:
        model = QCGeneList
        fields = ("path", "qc", "gene_list")

    def create(self, validated_data):
        qc_data = validated_data.pop("qc")
        qc = QCSerializer.get_object(qc_data)
        gene_list_data = validated_data.pop("gene_list")
        gene_list_text = ",".join(gene_list_data)
        custom_text_gene_list = QCGeneList.create_gene_list(gene_list_text,
                                                            sequencing_sample=qc.sequencing_sample)
        instance = QCGeneList.objects.create(qc=qc,
                                             custom_text_gene_list=custom_text_gene_list)

        # With API - whatever we sent is always the active one
        instance.link_samples_if_exist(force_active=True)
        return instance


class QCGeneCoverageSerializer(SeqAutoViewMixin, serializers.ModelSerializer):
    qc = QCSerializer()
    gene_coverage_collection = GeneCoverageCollectionSerializer()

    class Meta:
        model = QCGeneCoverage
        fields = ("path", "qc", "gene_coverage_collection")


class QCExecSummarySerializer(SeqAutoViewMixin, serializers.ModelSerializer):
    qc = QCSerializer()

    class Meta:
        model = QCExecSummary
        exclude = ('gene_list', )

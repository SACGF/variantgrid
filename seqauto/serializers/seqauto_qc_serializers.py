from rest_framework import serializers

from genes.serializers import SampleGeneListSerializer, GeneCoverageCollectionSerializer
from seqauto.models import IlluminaFlowcellQC, QCGeneList, QC, QCGeneCoverage, QCExecSummary, FastQC, SequencingSample
from seqauto.serializers.sequencing_serializers import SampleSheetLookupSerializer, FastqSerializer, SeqAutoViewMixin, \
    BamFilePathSerializer, VCFFilePathSerializer, SequencingSampleLookupSerializer


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
    sequencing_sample = serializers.SerializerMethodField(read_only=True)
    bam_file = BamFilePathSerializer()  # These end up really big, think we just want to pass PKs
    vcf_file = VCFFilePathSerializer()

    class Meta:
        model = QC
        fields = ("sequencing_sample", "bam_file", "vcf_file")

    def get_sequencing_sample(self, instance):
        try:
            ss = instance.bam_file.unaligned_reads.sequencing_sample
            serializer = SequencingSampleLookupSerializer(ss)
            return serializer.data
        except SequencingSample.DoesNotExist:
            return None


class QCGeneListSerializer(SeqAutoViewMixin, serializers.ModelSerializer):
    qc = QCSerializer()
    sample_gene_list = SampleGeneListSerializer()

    class Meta:
        model = QCGeneList
        fields = ("path", "qc", "sample_gene_list")


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

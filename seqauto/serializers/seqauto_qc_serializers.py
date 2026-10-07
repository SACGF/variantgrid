from rest_framework import serializers

from genes.serializers import GeneCoverageCollectionSerializer, SampleGeneListSerializer
from seqauto.models import (
    QC,
    FastQC,
    IlluminaFlowcellQC,
    QCExecSummary,
    QCGeneCoverage,
    QCGeneList,
    SingleSampleVCF,
)
from seqauto.serializers.sequencing_serializers import (
    AlignmentFilePathSerializer,
    FastqSerializer,
    SampleSheetLookupSerializer,
    SequencingSampleLookupSerializer,
    SingleSampleVCFPathSerializer,
)


class FastQCSerializer(serializers.ModelSerializer):
    fastq = FastqSerializer()

    class Meta:
        model = FastQC
        fields = "__all__"


class IlluminaFlowcellQCSerializer(serializers.ModelSerializer):
    sample_sheet = SampleSheetLookupSerializer()

    class Meta:
        model = IlluminaFlowcellQC
        #fields = "__all__"
        exclude = ("sequencing_run", )  # Already part of sample_sheet


class QCSerializer(serializers.ModelSerializer):
    """ Finds a QC from its sequencing sample and VCF. The QC hangs off the VCF's own alignment file (the one it was
        called from), so the alignment file sent only picks between VCFs sharing a path - it can be any of the
        sample's alignment files, eg a CRAM sent after the BAM the VCF was called from.
        'bam_file' is the pre 'alignment_file' name and still accepted """
    sequencing_sample = SequencingSampleLookupSerializer()
    alignment_file = AlignmentFilePathSerializer(required=False)
    bam_file = AlignmentFilePathSerializer(write_only=True, required=False)
    vcf_file = SingleSampleVCFPathSerializer()

    class Meta:
        model = QC
        fields = ("sequencing_sample", "alignment_file", "bam_file", "vcf_file")

    def validate(self, attrs):
        alignment_file = attrs.get("alignment_file")
        if bam_file := attrs.pop("bam_file", None):
            if alignment_file and alignment_file["path"] != bam_file["path"]:
                raise serializers.ValidationError("'alignment_file' and deprecated 'bam_file' have different paths "
                                                  f"('{alignment_file['path']}' / '{bam_file['path']}')")
            attrs["alignment_file"] = bam_file
        return attrs

    def to_representation(self, instance):
        data = super().to_representation(instance)
        data["bam_file"] = data["alignment_file"]  # Older clients read 'bam_file'
        return data

    @staticmethod
    def get_object(data):
        # We are passed "sequencing_sample" - which we can use to get what we really want
        sequencing_sample = SequencingSampleLookupSerializer.get_object(data.pop("sequencing_sample"))
        sequencing_run = sequencing_sample.sequencing_run
        alignment_file_data = data.pop("alignment_file", None)
        vcf_file_data = data.pop("vcf_file")
        vcf_file_kwargs = {
            "path": vcf_file_data["path"],
            "alignment_file__sequencing_run": sequencing_run,
            "alignment_file__sequencing_sample": sequencing_sample,
        }
        vcf_qs = SingleSampleVCF.objects.filter(**vcf_file_kwargs)
        if alignment_file_data:
            called_from_alignment_qs = vcf_qs.filter(alignment_file__path=alignment_file_data["path"])
            if called_from_alignment_qs.exists():
                vcf_qs = called_from_alignment_qs
        vcf_file = vcf_qs.order_by("pk").first()
        if not vcf_file:
            raise SingleSampleVCF.DoesNotExist(f"No vcf file for {vcf_file_kwargs=}")

        defaults = {}
        qc_path = data.get("path")
        if qc_path is None:
            # We currently require path to define QC (in sequencing scans)
            # @see QC class docs
            qc_path = QC.get_path_from_vcf(vcf_file)
        defaults["path"] = qc_path

        qc, _ = QC.objects.get_or_create(
            sequencing_run=sequencing_run,
            alignment_file=vcf_file.alignment_file,
            vcf_file=vcf_file,
            defaults=defaults
        )
        return qc


class QCGeneListSerializer(serializers.ModelSerializer):
    """ When we retrieve this, we want to see linked sample gene list """
    qc = QCSerializer()
    sample_gene_list = SampleGeneListSerializer()

    class Meta:
        model = QCGeneList
        fields = ("path", "qc", "sample_gene_list")


class QCGeneListCreateSerializer(serializers.ModelSerializer):
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
        path = validated_data["path"]
        qc_data = validated_data.pop("qc")
        qc = QCSerializer.get_object(qc_data)
        gene_list_data = validated_data.pop("gene_list")
        gene_list_text = ",".join(gene_list_data)
        custom_text_gene_list = QCGeneList.create_gene_list(gene_list_text,
                                                            sequencing_sample=qc.sequencing_sample)
        defaults = {
            "custom_text_gene_list": custom_text_gene_list,
        }
        instance, _created = QCGeneList.objects.update_or_create(qc=qc,
                                                                 path=path,
                                                                 defaults=defaults)

        # With API - whatever we sent is always the active one
        instance.link_samples_if_exist(force_active=True)
        return instance


class QCGeneListBulkCreateSerializer(serializers.Serializer):
    records = QCGeneListCreateSerializer(many=True)

    def create(self, validated_data):
        records = validated_data.get("records", [])
        qcgl_serializer = QCGeneListCreateSerializer()
        created_records = []
        for record in records:
            qcgl = qcgl_serializer.create(record)
            created_records.append(qcgl)
        return {
            "records": created_records,
        }


class QCGeneCoverageSerializer(serializers.ModelSerializer):
    """ The goal here is to just set the path - so that when we upload the file (and path)
        we can match paths in ImportGeneCoverageTask """
    qc = QCSerializer()
    gene_coverage_collection = GeneCoverageCollectionSerializer(read_only=True)

    class Meta:
        model = QCGeneCoverage
        fields = ("path", "qc", "gene_coverage_collection")

    def create(self, validated_data):
        qc_data = validated_data.pop("qc")
        qc = QCSerializer.get_object(qc_data)
        path = validated_data["path"]

        defaults = {
            "path": path,
        }
        instance, _created = QCGeneCoverage.objects.update_or_create(qc=qc,
                                                                     defaults=defaults)
        return instance


class QCExecSummarySerializer(serializers.ModelSerializer):
    qc = QCSerializer()

    class Meta:
        model = QCExecSummary
        exclude = ('gene_list', )

    def create(self, validated_data):
        qc_data = validated_data.pop("qc")
        qc = QCSerializer.get_object(qc_data)
        validated_data["sequencing_run"] = qc.sequencing_run
        instance, _created = QCExecSummary.objects.update_or_create(qc=qc,
                                                                    defaults=validated_data)
        return instance


class QCExecSummaryBulkCreateSerializer(serializers.Serializer):
    records = QCExecSummarySerializer(many=True)

    def create(self, validated_data):
        records = validated_data.get("records", [])
        qces_serializer = QCExecSummarySerializer()
        created_records = []
        for record in records:
            qcgl = qces_serializer.create(record)
            created_records.append(qcgl)
        return {
            "records": created_records,
        }


class QCGeneCoverageBulkCreateSerializer(serializers.Serializer):
    records = QCGeneCoverageSerializer(many=True)

    def create(self, validated_data):
        records = validated_data.get("records", [])
        qcgc_serializer = QCGeneCoverageSerializer()
        created_records = []
        for record in records:
            qcgl = qcgc_serializer.create(record)
            created_records.append(qcgl)
        return {
            "records": created_records,
        }

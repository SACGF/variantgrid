from rest_framework import serializers

from genes.serializers import SampleGeneListSerializer, GeneCoverageCollectionSerializer
from seqauto.models import IlluminaFlowcellQC, QCGeneList, QC, QCGeneCoverage, QCExecSummary, FastQC, SequencingSample, \
    SampleSheet, SequencingRun, BamFile, VCFFile
from seqauto.serializers.sequencing_serializers import SampleSheetLookupSerializer, FastqSerializer, \
    BamFilePathSerializer, VCFFilePathSerializer, SequencingSampleLookupSerializer
from snpdb.models import DataState


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
        # Occasionally we could have multiples in there, we don't really care so take 1st
        bam_file_kwargs = {
            "path": bam_file_data["path"],
            "sequencing_run": sequencing_run,
            "unaligned_reads__sequencing_sample": sequencing_sample
        }
        bam_file = BamFile.objects.filter(**bam_file_kwargs).first()
        if not bam_file:
            raise BamFile.DoesNotExist(f"No bam file for {bam_file_kwargs=}")

        vcf_file_data = data.pop("vcf_file")
        vcf_file_kwargs = {
            "path": vcf_file_data["path"],
            "bam_file": bam_file,
        }
        vcf_file = VCFFile.objects.filter(**vcf_file_kwargs).first()
        if not vcf_file:
            raise VCFFile.DoesNotExist(f"No vcf file for {vcf_file_kwargs=}")

        defaults = {}
        if qc_path := data.get("path"):
            defaults["path"] = qc_path

        qc, _ = QC.objects.get_or_create(
            sequencing_run=sequencing_run,
            bam_file=bam_file,
            vcf_file=vcf_file,
            defaults=defaults
        )
        qc.data_state = DataState.COMPLETE
        qc.save()
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
    data_state = serializers.CharField(read_only=True)

    class Meta:
        model = QCExecSummary
        exclude = ('gene_list', )

    def create(self, validated_data):
        qc_data = validated_data.pop("qc")
        qc = QCSerializer.get_object(qc_data)
        validated_data["data_state"] = DataState.COMPLETE
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


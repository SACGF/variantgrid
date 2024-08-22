import logging

from rest_framework import serializers

from seqauto.models import Sequencer, Experiment, VariantCaller, SequencingRun, SequencerModel, SampleSheet, \
    SequencingSampleData, SequencingSample, UnalignedReads, Flagstats, FastQC, SampleSheetCombinedVCFFile, VCFFile, \
    BamFile, Fastq, Aligner
from seqauto.serializers import EnrichmentKitSerializer
from snpdb.models import Manufacturer


class SequencerModelSerializer(serializers.ModelSerializer):
    manufacturer = serializers.StringRelatedField()
    data_naming_convention = serializers.SerializerMethodField()

    class Meta:
        model = SequencerModel
        fields = "__all__"

    def create(self, validated_data):
        model = validated_data.get('model')
        manufacturer = validated_data.get('manufacturer')
        manufacturer, _ = Manufacturer.objects.get_or_create(name=manufacturer)

        # Check if the object already exists
        instance, _created = SequencerModel.objects.get_or_create(
            model=model,
            defaults={
                "manufacturer": manufacturer,
                "data_naming_convention": model.data_naming_convention,
            }
        )
        return instance  # Return the existing or new instance

    def get_data_naming_convention(self, obj):
        return obj.get_data_naming_convention_display()


class SequencerSerializer(serializers.ModelSerializer):
    sequencer_model = SequencerModelSerializer()

    class Meta:
        model = Sequencer
        fields = "__all__"

    def create(self, validated_data):
        name = validated_data.get('name')
        sequencer_model = validated_data.get('sequencer_model')
        logging.info("sequencer_model=%s", sequencer_model)

        instance, _created = Sequencer.objects.get_or_create(
            name=name,
            defaults={
                "sequencer_model": sequencer_model,
            }
        )
        return instance


class ExperimentSerializer(serializers.ModelSerializer):
    class Meta:
        model = Experiment
        fields = ["name"]

    def create(self, validated_data):
        name = validated_data.get('name')
        instance, _created = Experiment.objects.get_or_create(
            name=name
        )
        return instance


class AlignerSerializer(serializers.ModelSerializer):
    class Meta:
        model = Aligner
        fields = "__all__"

    def create(self, validated_data):
        name = validated_data.get('name')
        version = validated_data.get('version')

        instance, _created = Aligner.objects.get_or_create(
            name=name,
            version=version
        )
        return instance


class VariantCallerSerializer(serializers.ModelSerializer):
    class Meta:
        model = VariantCaller
        fields = "__all__"

    def create(self, validated_data):
        name = validated_data.get('name')
        version = validated_data.get('version')
        run_params = validated_data.get('run_params')

        instance, _created = VariantCaller.objects.get_or_create(
            name=name,
            version=version,
            defaults={"run_params": run_params}
        )
        return instance


class SequencingRunSerializer(serializers.ModelSerializer):
    sequencer = SequencerSerializer()
    experiment = ExperimentSerializer()
    enrichment_kit = EnrichmentKitSerializer()

    class Meta:
        model = SequencingRun
        fields = ("name", "date", "sequencer", "gold_standard", "bad", "hidden", "experiment", "enrichment_kit", "has_basecalls", "has_interop")

    def create(self, validated_data):
        name = validated_data.get('name')
        instance, _created = SequencingRun.objects.get_or_create(
            name=name,
            defaults=validated_data
        )
        return instance


class SequencingSampleDataSerializer(serializers.ModelSerializer):
    # sequencing_run = SequencingRunSerializer()

    class Meta:
        model = SequencingSampleData
        fields = ("column", "value")


class SampleSheetLookupSerializer(serializers.Serializer):
    """ This is when we want to refer to it in related objects in a minimal way
        The samplesheet MUST have already been created """
    sequencing_run = serializers.CharField()
    hash = serializers.CharField()

    def validate(self, attrs):
        sequencing_run = attrs.get('sequencing_run')
        hash = attrs.get('hash')

        try:
            sample_sheet = SampleSheet.objects.get(
                sequencing_run__name=sequencing_run,
                hash=hash,
            )
        except SampleSheet.DoesNotExist:
            raise serializers.ValidationError("SampleSheet not found.")

        attrs['sample_sheet'] = sample_sheet
        return attrs


class SequencingSampleLookupSerializer(serializers.Serializer):
    """ This is when we want to refer to it in related objects in a minimal way """
    sample_sheet = SampleSheetLookupSerializer()
    sample_name = serializers.CharField()

    def validate(self, attrs):
        sample_sheet_data = attrs.pop('sample_sheet')
        validated_sample_sheet = SampleSheetLookupSerializer(data=sample_sheet_data)
        validated_sample_sheet.is_valid(raise_exception=True)
        attrs['sample_sheet'] = validated_sample_sheet.validated_data['sample_sheet']
        return attrs


class SequencingSampleSerializer(serializers.ModelSerializer):
    """ This is when we want the whole object as a dict """
    sequencingsampledata_set = SequencingSampleDataSerializer(many=True, required=False)
    enrichment_kit = EnrichmentKitSerializer()

    class Meta:
        model = SequencingSample
        fields = ['sample_id', 'sample_name', 'sample_project', 'sample_number', 'lane', 'barcode', 'enrichment_kit', 'is_control', 'failed', 'sequencingsampledata_set']


# TODO: Write ViewSet
class SampleSheetSerializer(serializers.ModelSerializer):
    sequencing_run = SequencingRunSerializer()
    sequencingsample_set = SequencingSampleSerializer(many=True)

    class Meta:
        model = SampleSheet
        fields = ("path", "sequencing_run", "file_last_modified", "hash", "sequencingsample_set")

    @staticmethod
    def _create_sequencing_samples(sample_sheet, sequencing_samples_data):
        for sample_data in sequencing_samples_data:
            sample_data_data = sample_data.pop('sample_data', [])
            sequencing_sample = SequencingSample.objects.create(sample_sheet=sample_sheet, **sample_data)
            for data in sample_data_data:
                SequencingSampleData.objects.create(sequencing_sample=sequencing_sample, **data)

    def create(self, validated_data):
        sequencing_samples_data = validated_data.pop('sequencing_samples')
        sample_sheet, _created = SampleSheet.objects.update_or_create(
            sequencing_run=validated_data["sequencing_run"],
            hash=validated_data["hash"],
            defaults=validated_data,
        )
        self._create_sequencing_samples(sample_sheet, sequencing_samples_data)
        return sample_sheet

    def update(self, instance, validated_data):
        sequencing_samples_data = validated_data.pop('sequencing_samples')

        instance = super().update(instance, validated_data)
        # Clear existing samples and add the new ones
        instance.sequencing_samples.all().delete()
        self._create_sequencing_samples(instance, sequencing_samples_data)
        return instance


class FastqSerializer(serializers.ModelSerializer):
    class Meta:
        model = Fastq
        fields = ("path", "name", "read")


class UnalignedReadsSerializer(serializers.ModelSerializer):
    sequencing_sample = SequencingSampleLookupSerializer()
    fastq_r1 = FastqSerializer()
    fastq_r2 = FastqSerializer()

    class Meta:
        model = UnalignedReads
        fields = "__all__"

    def create(self, validated_data):
        sequencing_sample = validated_data.pop('sequencing_sample')
        sequencing_run = sequencing_sample.sequencing_run
        unaligned_reads_kwargs = {}

        for field_name in ["fastq_r1", "fastq_r2"]:
            if fastq_data := validated_data.pop(field_name):
                fastq = Fastq.objects.update_or_create(sequencing_sample=sequencing_sample,
                                                       defaults=fastq_data)
                unaligned_reads_kwargs[field_name] = fastq

        instance, _created = UnalignedReads.objects.update_or_create(sequencing_run=sequencing_run,
                                                                     sequencing_sample=sequencing_sample,
                                                                     defaults=unaligned_reads_kwargs)
        return instance


class FlagstatsSerializer(serializers.ModelSerializer):
    class Meta:
        model = Flagstats
        fields = ("total", "read1", "read2", "mapped", "properly_paired")


class BamFileSerializer(serializers.ModelSerializer):
    unaligned_reads = UnalignedReadsSerializer()
    aligner = AlignerSerializer()
    flagstats = FlagstatsSerializer()  # 1-to-1 field

    class Meta:
        model = BamFile
        fields = ("path", "unaligned_reads", "name", "aligner", "flagstats")

    def create(self, validated_data):
        flagstats_data = validated_data.pop('flagstats', None)
        bam_file = BamFile.objects.create(**validated_data)

        if flagstats_data:
            Flagstats.objects.create(bam_file=bam_file, **flagstats_data)

        return bam_file

    def update(self, instance, validated_data):
        flagstats_data = validated_data.pop('flagstats', None)
        instance = super().update(instance, validated_data)
        if flagstats_data:
            flagstats_instance = getattr(instance, 'flagstats', None)
            flagstats_serializer = FlagstatsSerializer(instance=flagstats_instance, data=flagstats_data)
            flagstats_serializer.is_valid(raise_exception=True)
            flagstats_serializer.save(bam_file=instance)
        return instance


# TODO: Write ViewSet - can we do in bulk?
class VCFFileSerializer(serializers.ModelSerializer):
    bam_file = BamFileSerializer()
    variant_caller = VariantCallerSerializer()

    class Meta:
        model = VCFFile
        fields = ("path", "bam_file", "variant_caller")


# TODO: Write ViewSet
class SampleSheetCombinedVCFFileSerializer(serializers.ModelSerializer):
    sample_sheet = SampleSheetLookupSerializer()
    variant_caller = VariantCallerSerializer()

    class Meta:
        model = SampleSheetCombinedVCFFile
        fields = ("path", "sample_sheet", "variant_caller")


# TODO: Write ViewSet
class FastQCSerializer(serializers.ModelSerializer):
    fastq = FastqSerializer()

    class Meta:
        model = FastQC
        fields = ("fastq", "total_sequences", "filtered_sequences", "gc")

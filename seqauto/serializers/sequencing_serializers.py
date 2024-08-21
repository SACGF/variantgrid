import logging

from rest_framework import serializers

from seqauto.models import Sequencer, Experiment, VariantCaller, SequencingRun, SequencerModel
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
    sequencer_model = serializers.StringRelatedField()

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

# STEP 2

#SampleSheet
#SequencingSample
#SequencingSampleData

#UnalignedReads  # Needs sequencing sample
#Fastq

#BamFile  # needs unaligned_reads
#VCFFile
#SampleSheetCombinedVCFFile  # sample sheet / variant caller


#FastQC
#Flagstats # Needs BamFile

#SampleFromSequencingSample

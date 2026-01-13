import os.path

from rest_framework import serializers

from seqauto.models import Sequencer, Experiment, VariantCaller, SequencingRun, SequencerModel, SampleSheet, \
    SequencingSampleData, SequencingSample, UnalignedReads, Flagstats, SampleSheetCombinedVCFFile, VCFFile, \
    BamFile, Fastq, Aligner, PairedEnd
from seqauto.serializers import EnrichmentKitSerializer, EnrichmentKitSummarySerializer
from snpdb.models import Manufacturer, DataState


class ManufacturerSerializer(serializers.ModelSerializer):
    class Meta:
        model = Manufacturer
        fields = ["name"]

    def to_internal_value(self, data):
        return Manufacturer.objects.get_or_create(name=data["name"])[0]


class SequencerModelSerializer(serializers.ModelSerializer):
    manufacturer = ManufacturerSerializer()
    data_naming_convention_display = serializers.CharField(source='get_data_naming_convention_display', read_only=True)

    class Meta:
        model = SequencerModel
        extra_kwargs = {'model': {'validators': []}}  # turn off UniqueValidator
        fields = ["model", "manufacturer", "data_naming_convention", "data_naming_convention_display"]

    def create(self, validated_data):
        instance, _created = SequencerModel.objects.get_or_create(
            model=validated_data["model"],
            defaults={
                "manufacturer": validated_data.get("manufacturer"),
                "data_naming_convention": validated_data.get("data_naming_convention"),
            }
        )
        return instance

    @staticmethod
    def get_object(validated_data):
        return SequencerModel.objects.get(model=validated_data["model"])


class SequencerSerializer(serializers.ModelSerializer):
    name = serializers.CharField(validators=[])  # Disable UniqueValidator
    sequencer_model = SequencerModelSerializer()

    class Meta:
        model = Sequencer
        fields = ["name", "sequencer_model"]

    def create(self, validated_data):
        name = validated_data.get('name')
        sequencer_model_data = validated_data.get('sequencer_model')
        sequencer_model = SequencerModelSerializer.get_object(sequencer_model_data)

        instance, _created = Sequencer.objects.get_or_create(
            name=name,
            defaults={
                "sequencer_model": sequencer_model,
            }
        )
        return instance


class ExperimentSerializer(serializers.ModelSerializer):
    # TODO: duplicated logic below in internal/create - I think we can just use SlugRelatedField
    name = serializers.CharField(validators=[])  # Remove the UniqueValidator

    class Meta:
        model = Experiment
        fields = ["name"]

    def to_internal_value(self, data):
        # When POSTing, we expect only the `name` to be passed
        name = data.get('name')
        if experiment := Experiment.objects.filter(name=name).first():
            data = self.to_representation(experiment)
        return super().to_internal_value(data)

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
    name = serializers.CharField(validators=[])  # disable UniqueValidator
    sequencer = serializers.PrimaryKeyRelatedField(queryset=Sequencer.objects.all())
    experiment = serializers.PrimaryKeyRelatedField(queryset=Experiment.objects.all())
    enrichment_kit = EnrichmentKitSummarySerializer()
    vcf_set = serializers.SerializerMethodField(read_only=True)

    class Meta:
        model = SequencingRun
        fields = ("path", "name", "date", "sequencer", "gold_standard", "bad", "hidden", "experiment", "enrichment_kit", "has_basecalls", "has_interop", "vcf_set")

    def validate_name(self, value):
        if error := SequencingRun.get_name_validation_errors(value):
            raise serializers.ValidationError(error)
        return value

    def create(self, validated_data):
        name = validated_data.get('name')
        if ek_data := validated_data.pop('enrichment_kit', None):
            enrichment_kit = EnrichmentKitSerializer.get_from_data(ek_data)
            validated_data['enrichment_kit'] = enrichment_kit
        instance, _created = SequencingRun.objects.get_or_create(
            name=name,
            defaults=validated_data
        )
        return instance

    def get_vcf_set(self, obj):
        vcfs = []
        for vsr in obj.vcffromsequencingrun_set.all().order_by("pk"):
            vcf = vsr.vcf
            data = {
                "pk": vcf.pk,
                "name": vcf.name,
            }
            try:
                # This is set on ones sent up via API
                data["path"] = vcf.uploadedvcf.uploaded_file.path
            except:
                pass
            vcfs.append(data)
        return vcfs


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

    @staticmethod
    def get_object(validated_data):
        sequencing_run = validated_data['sequencing_run']
        hash = validated_data['hash']

        try:
            return SampleSheet.objects.get(
                sequencing_run__name=sequencing_run,
                hash=hash,
            )
        except SampleSheet.DoesNotExist:
            raise serializers.ValidationError("SampleSheet not found.")


class SequencingSampleLookupSerializer(serializers.Serializer):
    """ This is when we want to refer to it in related objects in a minimal way """
    sample_sheet = SampleSheetLookupSerializer()
    sample_name = serializers.CharField()


class SequencingSampleSerializer(serializers.ModelSerializer):
    """ This is when we want the whole object as a dict """
    sequencingsampledata_set = SequencingSampleDataSerializer(many=True, required=False)
    enrichment_kit = EnrichmentKitSerializer()

    class Meta:
        model = SequencingSample
        fields = ['sample_id', 'sample_name', 'sample_project', 'sample_number', 'lane', 'barcode', 'enrichment_kit', 'is_control', 'failed', 'sequencingsampledata_set']


class SampleSheetSerializer(serializers.ModelSerializer):
    sequencing_run = serializers.PrimaryKeyRelatedField(queryset=SequencingRun.objects.all())
    sequencingsample_set = SequencingSampleSerializer(many=True)

    class Meta:
        model = SampleSheet
        fields = ("path", "sequencing_run", "file_last_modified", "hash", "sequencingsample_set")

    @staticmethod
    def _create_sequencing_samples(sample_sheet, sequencing_samples_data):
        for sample_data in sequencing_samples_data:
            if not sample_data.get("sample_name"):
                sample_data["sample_name"] = sample_data["sample_id"]

            ss_data = sample_data.pop('sequencingsampledata_set', [])
            if ek_data := sample_data.pop("enrichment_kit", None):
                enrichment_kit = EnrichmentKitSerializer.get_from_data(ek_data)
                sample_data["enrichment_kit"] = enrichment_kit
            sample_id = sample_data.pop("sample_id")
            sequencing_sample, _ = SequencingSample.objects.update_or_create(sample_sheet=sample_sheet,
                                                                             sample_id=sample_id,
                                                                             defaults=sample_data)
            for data in ss_data:
                column = data.pop("column")
                SequencingSampleData.objects.update_or_create(sequencing_sample=sequencing_sample,
                                                              column=column,
                                                              defaults=data)

    def create(self, validated_data):
        sequencing_samples_data = validated_data.pop('sequencingsample_set')
        validated_data["data_state"] = DataState.COMPLETE
        sequencing_run = validated_data["sequencing_run"]
        sample_sheet, created = SampleSheet.objects.update_or_create(
            sequencing_run=sequencing_run,
            hash=validated_data["hash"],
            defaults=validated_data,
        )
        self._create_sequencing_samples(sample_sheet, sequencing_samples_data)
        # Whatever is last sent via API is the current sample sheet
        sample_sheet.set_as_current_sample_sheet(sequencing_run, created)
        return sample_sheet

    def update(self, instance, validated_data):
        sequencing_samples_data = validated_data.pop('sequencing_samples')

        instance = super().update(instance, validated_data)
        # Clear existing samples and add the new ones
        instance.sequencing_samples.all().delete()
        self._create_sequencing_samples(instance, sequencing_samples_data)
        return instance


class FastqSerializer(serializers.ModelSerializer):
    name = serializers.CharField(read_only=True)
    read = serializers.SerializerMethodField(read_only=True)

    class Meta:
        model = Fastq
        fields = ("path", "name", "read")

    def create(self, validated_data):
        path = validated_data["path"]
        validated_data["name"] = os.path.basename(path)
        fastq, _created = Fastq.objects.update_or_create(**validated_data,
                                                         defaults={"data_state": DataState.COMPLETE})
        return fastq


class UnalignedReadsSerializer(serializers.ModelSerializer):
    """ UnalignedReads is a joiner class - not represented by a file thus not SeqAutoRecord """
    sequencing_sample = SequencingSampleLookupSerializer(required=False)  # Might not be there during validation
    fastq_r1 = FastqSerializer()
    fastq_r2 = FastqSerializer(required=False)

    class Meta:
        model = UnalignedReads
        fields = "__all__"

    def create(self, validated_data):
        sequencing_sample = validated_data['sequencing_sample']
        unaligned_reads_kwargs = {}

        fastq_serializer = FastqSerializer()
        # Always have R1
        fastq_r1_data = validated_data['fastq_r1']
        fastq_r1_data["read"] = PairedEnd.R1
        fastq_r1_data["sequencing_sample"] = sequencing_sample
        unaligned_reads_kwargs["fastq_r1"] = fastq_serializer.create(fastq_r1_data)

        if fastq_r2_data := validated_data.get("fastq_r2"):
            fastq_r2_data["read"] = PairedEnd.R2
            fastq_r2_data["sequencing_sample"] = sequencing_sample
            unaligned_reads_kwargs["fastq_r2"] = fastq_serializer.create(fastq_r2_data)
        else:
            unaligned_reads_kwargs["fastq_r2"] = None  # Be able to blank it out

        # Unaligned reads isn't a file so doesn't have 'data_state'
        instance, _created = UnalignedReads.objects.update_or_create(sequencing_sample=sequencing_sample,
                                                                     defaults=unaligned_reads_kwargs)
        return instance


class FlagstatsSerializer(serializers.ModelSerializer):
    class Meta:
        model = Flagstats
        fields = ("total", "read1", "read2", "mapped", "properly_paired")


class BamFilePathSerializer(serializers.ModelSerializer):
    class Meta:
        model = BamFile
        fields = ("path", )


class BamFileSerializer(serializers.ModelSerializer):
    unaligned_reads = UnalignedReadsSerializer(required=False)
    aligner = AlignerSerializer(required=False)
    flagstats = FlagstatsSerializer(read_only=True, required=False)  # 1-to-1 field
    name = serializers.CharField(read_only=True)

    class Meta:
        model = BamFile
        fields = ("path", "unaligned_reads", "name", "aligner", "flagstats")

    def create(self, validated_data):
        unaligned_reads_data = validated_data["unaligned_reads"]
        aligner_data = validated_data["aligner"]
        flagstats_data = validated_data.get('flagstats')

        unaligned_reads = UnalignedReadsSerializer().create(unaligned_reads_data)
        aligner = AlignerSerializer().create(aligner_data)
        path = validated_data["path"]
        name = os.path.basename(path)
        bam_file, _ = BamFile.objects.update_or_create(path=path,
                                                       sequencing_run=unaligned_reads.sequencing_run,
                                                       unaligned_reads=unaligned_reads,
                                                       aligner=aligner,
                                                       name=name,
                                                       defaults={"data_state": DataState.COMPLETE})

        if flagstats_data:
            flagstats_data["data_state"] = DataState.COMPLETE
            Flagstats.objects.create(bam_file=bam_file, **flagstats_data)

        return bam_file

    def update(self, instance, validated_data):
        flagstats_data = validated_data.get('flagstats')
        instance = super().update(instance, validated_data)
        if flagstats_data:
            flagstats_instance = getattr(instance, 'flagstats', None)
            flagstats_serializer = FlagstatsSerializer(instance=flagstats_instance, data=flagstats_data)
            flagstats_serializer.is_valid(raise_exception=True)
            flagstats_serializer.save(bam_file=instance)
        return instance


class VCFFilePathSerializer(serializers.ModelSerializer):
    class Meta:
        model = VCFFile
        fields = ("path", )


def validate_unique_vcf_path(klass, path, **kwargs):
    """ We need to ensure there is only 1 backend/sequencing VCF model with path
        as we need it as a unique result lookup in upload.vcf.vcf_import.create_backend_vcf_links

        We want the API to be idempotent (can be called multiple times) so ok if the lookup is the same
        (hence exists call if passed in same class below)
    """
    for vcf_class in [VCFFile, SampleSheetCombinedVCFFile]:
        qs = vcf_class.objects.filter(path=path)
        if klass == vcf_class:
            qs = qs.exclude(**kwargs)  # Allow for same object/lookup
        if obj := qs.first():
            if klass == vcf_class:
                diff = []
                for k, v in kwargs.items():
                    o_v = getattr(obj, k)
                    if o_v != v:
                        diff.append(f"{k} = {v}/{o_v}")
                diff_str = ", ".join(diff)
                msg = f"Cannot make record for {klass.__name__} with same path but different fields: {diff_str}"
            else:
                msg = f"Cannot make record for {klass.__name__}, when other sequencing VCF file type {vcf_class.__name__} with {path=} exists"

            raise serializers.ValidationError({
                "path": f"Existing {vcf_class.__name__}(pk={obj.pk}) with {path=} exists: {msg}",
            })


class VCFFileSerializer(serializers.ModelSerializer):
    bam_file = BamFileSerializer(required=False)
    variant_caller = VariantCallerSerializer()

    class Meta:
        model = VCFFile
        fields = ("path", "bam_file", "variant_caller")

    def create(self, validated_data):
        bam_file_data = validated_data['bam_file']
        variant_caller_data = validated_data['variant_caller']
        bam_file = BamFileSerializer().create(bam_file_data)
        variant_caller = VariantCallerSerializer().create(variant_caller_data)
        path = validated_data["path"]
        kwargs = {
            "sequencing_run": bam_file.sequencing_run,
            "bam_file": bam_file,
            "variant_caller": variant_caller,
        }
        validate_unique_vcf_path(VCFFile, path, **kwargs)
        vcf_file, _ = VCFFile.objects.update_or_create(**kwargs,
                                                       defaults={
                                                           "path": path,
                                                           "data_state": DataState.COMPLETE
                                                       })
        return vcf_file

class SequencingFilesSerializer(serializers.Serializer):
    sample_name = serializers.CharField()
    unaligned_reads = UnalignedReadsSerializer()
    bam_file = BamFileSerializer()
    vcf_file = VCFFileSerializer()

    def __init__(self, *args, **kwargs):
        self.sequencing_sample = kwargs.pop("sequencing_sample", None)
        super().__init__(*args, **kwargs)

    def create(self, validated_data):
        unaligned_reads_data = validated_data.pop('unaligned_reads')
        bam_file_data = validated_data.pop('bam_file')
        vcf_file_data = validated_data.pop('vcf_file')

        if self.sequencing_sample is None:
            raise ValueError("SequencingSample is required for create()")
        unaligned_reads_data["sequencing_sample"] = self.sequencing_sample
        unaligned_reads = UnalignedReadsSerializer().create(unaligned_reads_data)

        bam_file_data['unaligned_reads'] = unaligned_reads_data
        bam_file = BamFileSerializer().create(bam_file_data)

        vcf_file_data['bam_file'] = bam_file_data
        vcf_file = VCFFileSerializer().create(vcf_file_data)

        # Return the full data structure (optional depending on your needs)
        return {
            'sample_name': validated_data['sample_name'],
            'unaligned_reads': unaligned_reads,
            'bam_file': bam_file,
            'vcf_file': vcf_file
        }


class SequencingFilesBulkCreateSerializer(serializers.Serializer):
    sample_sheet = SampleSheetLookupSerializer()
    records = SequencingFilesSerializer(many=True)

    def create(self, validated_data):
        sample_sheet = SampleSheetLookupSerializer.get_object(validated_data.pop('sample_sheet'))
        records_data = validated_data['records']

        ss_by_name = sample_sheet.get_sequencing_samples_by_name()
        sequencing_files = []
        for record in records_data:
            ss_name = record["sample_name"]
            sequencing_sample = ss_by_name[ss_name]
            sf = SequencingFilesSerializer(sequencing_sample=sequencing_sample).create(record)
            sequencing_files.append(sf)

        return {
            "sample_sheet": sample_sheet,
            "records": sequencing_files
        }


class SampleSheetCombinedVCFFileSerializer(serializers.ModelSerializer):
    sample_sheet = SampleSheetLookupSerializer()
    variant_caller = VariantCallerSerializer()

    class Meta:
        model = SampleSheetCombinedVCFFile
        fields = ("path", "sample_sheet", "variant_caller")

    def create(self, validated_data):
        sample_sheet = SampleSheetLookupSerializer.get_object(validated_data.pop('sample_sheet'))
        variant_caller = VariantCallerSerializer().create(validated_data.pop('variant_caller'))
        path = validated_data["path"]
        kwargs = {
            "sequencing_run": sample_sheet.sequencing_run,
            "sample_sheet": sample_sheet,
            "variant_caller": variant_caller,
        }
        validate_unique_vcf_path(SampleSheetCombinedVCFFile, path, **kwargs)
        sscvcf, _ = SampleSheetCombinedVCFFile.objects.update_or_create(**kwargs,
                                                                        defaults={
                                                                            "path": validated_data["path"],
                                                                            "data_state": DataState.COMPLETE,
                                                                        })
        return sscvcf

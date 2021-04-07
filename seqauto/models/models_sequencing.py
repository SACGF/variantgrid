import re

from django.db import models
from django.db.models import CASCADE, PROTECT
from django.urls import reverse

from genes.models import GeneList, CanonicalTranscriptCollection
from seqauto.illumina import illumina_sequencers
from seqauto.models.models_enums import DataGeneration, EnrichmentKitType
from snpdb.models import Manufacturer, GenomicIntervalsCollection, SET_NULL, LabProject


class SequencerModel(models.Model):
    model = models.TextField(primary_key=True)
    manufacturer = models.ForeignKey(Manufacturer, null=True, on_delete=CASCADE)
    data_naming_convention = models.CharField(max_length=1, choices=DataGeneration.choices)

    @property
    def css_class(self):
        return str(self.model).replace(' ', '-').lower()

    def __str__(self):
        return f"{self.manufacturer} {self.model}"


class Sequencer(models.Model):
    name = models.TextField(primary_key=True)
    sequencer_model = models.ForeignKey(SequencerModel, on_delete=CASCADE)

    @staticmethod
    def get_sequencers_dict_by_name():
        sequencers = {}
        for sequencer in Sequencer.objects.all():
            sequencers[sequencer.pk] = sequencer
        return sequencers

    @staticmethod
    def get_or_create_sequencer(name, existing_sequencers=None):
        if existing_sequencers is None:
            existing_sequencers = Sequencer.get_sequencers_dict_by_name()

        sequencer = existing_sequencers.get(name)
        if not sequencer:
            model_name = illumina_sequencers.get_sequencer_model_from_name(name)
            if model_name:
                illumina, _ = Manufacturer.objects.get_or_create(name='Illumina')
                try:
                    sequencer_model = SequencerModel.objects.get(model=model_name,
                                                                 manufacturer=illumina)
                except:
                    if 'HiSeq' in model_name:
                        data_naming_convention = DataGeneration.HISEQ
                    else:
                        data_naming_convention = DataGeneration.MISEQ
                    sequencer_model = SequencerModel.objects.create(model=model_name,
                                                                    data_naming_convention=data_naming_convention,
                                                                    manufacturer=illumina)

                sequencer = Sequencer.objects.create(name=name,
                                                     sequencer_model=sequencer_model)
                existing_sequencers[sequencer.name] = sequencer  # Add in case this is being used elsewhere
            else:
                msg = f"Can't recognise Illumina model from name '{name}' and don't know how to handle non-Illumina sequencers yet"
                raise ValueError(msg)

        return sequencer

    def get_absolute_url(self):
        return reverse('view_sequencer', kwargs={'pk': self.pk})

    def __str__(self):
        return f"{self.name} ({self.sequencer_model})"


class EnrichmentKit(models.Model):
    """ A lab method to enrich a sample (eg Capture Panel or Amplicon etc) """
    name = models.TextField()
    version = models.IntegerField(default=1)
    enrichment_kit_type = models.CharField(max_length=1, choices=EnrichmentKitType.choices, null=True)
    manufacturer = models.ForeignKey(Manufacturer, null=True, blank=True, on_delete=CASCADE)
    genomic_intervals = models.ForeignKey(GenomicIntervalsCollection, null=True, blank=True, on_delete=SET_NULL)
    gene_list = models.ForeignKey(GeneList, null=True, blank=True, on_delete=PROTECT)
    canonical_transcript_collection = models.ForeignKey(CanonicalTranscriptCollection, null=True, blank=True, on_delete=SET_NULL)
    obsolete = models.BooleanField(default=False)

    @staticmethod
    def get_latest_version(name):
        return EnrichmentKit.objects.filter(name=name).order_by("version").last()

    @staticmethod
    def get_enrichment_kits(enrichment_kit_kwargs_list):
        enrichment_kits = []
        for enrichment_kit_kwargs in enrichment_kit_kwargs_list:
            # There can be more than 1 per enrichment_kit name
            for enrichment_kit in EnrichmentKit.objects.filter(**enrichment_kit_kwargs):
                enrichment_kits.append(enrichment_kit)
        return enrichment_kits

    def get_gold_sequencing_runs_qs(self):
        return self.sequencingrun_set.filter(gold_standard=True)

    def get_absolute_url(self):
        return reverse('view_enrichment_kit', kwargs={"pk": self.pk})

    def __str__(self):
        name = self.name
        if self.version > 1:
            name += f" (version {self.version})"
        return name


class Library(models.Model):
    name = models.TextField(primary_key=True)
    description = models.TextField(null=True, blank=True)
    manufacturer = models.ForeignKey(Manufacturer, null=True, on_delete=CASCADE)

    def get_absolute_url(self):
        return reverse('view_library', kwargs={'pk': self.pk})

    def __str__(self):
        return ' '.join(map(str, filter(lambda x: x, [self.name, self.description, self.manufacturer])))


class Assay(models.Model):
    library = models.ForeignKey(Library, on_delete=CASCADE)
    sequencer = models.ForeignKey(Sequencer, on_delete=CASCADE)
    enrichment_kit = models.ForeignKey(EnrichmentKit, on_delete=CASCADE)  # Genes captured

    def get_absolute_url(self):
        return reverse('view_assay', kwargs={'pk': self.pk})

    def __str__(self):
        return f"Seq: {self.sequencer}, Lib: {self.library}, EnrichmentKit: {self.enrichment_kit}"


# ExperimentManager manages the Experiment objects
# to ensure that experiment names are cleaned before creation or updating
class ExperimentManager(models.Manager):

    @staticmethod
    def fix_kwargs(kwargs):
        try:
            old_name = kwargs["name"]
            name = Experiment.clean_experiment_name(old_name)
            kwargs["name"] = name
        except:
            pass
        return kwargs

    def get(self, *args, **kwargs):
        self.fix_kwargs(kwargs)
        return super().get(*args, **kwargs)

    def get_or_create(self, default=None, **kwargs):
        self.fix_kwargs(kwargs)
        return super().get_or_create(default, **kwargs)


class Experiment(models.Model):
    """ What was sequenced in the flowcell - ie you can resequence something and 2 SequencingRuns will have the
        same Experiment. Set from RunParameters.xml ExperimentName in the SequencingRun directory """
    name = models.TextField(primary_key=True)
    created = models.DateTimeField(auto_now_add=True)
    objects = ExperimentManager()

    @staticmethod
    def clean_experiment_name(experiment_name):
        experiment_name = experiment_name.upper()
        experiment_name = experiment_name.replace("-", "_")
        experiment_name = experiment_name.replace(" ", "_")

        # Remove RPT off the end...
        experiment_name = re.sub("_RPT$", "", experiment_name)
        experiment_name = re.sub("RPT$", "", experiment_name)
        return experiment_name

    def save(self, **kwargs):
        old_name = self.name
        self.name = Experiment.clean_experiment_name(old_name)
        return super().save(**kwargs)

    def can_write(self, user):
        """ can't delete once you've linked to SequencingRun """
        return user.is_superuser and not self.sequencingrun_set.exists()

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse('view_experiment', kwargs={"experiment_id": self.pk})


class SequencingInfo(models.Model):
    lab_project = models.ForeignKey(LabProject, on_delete=CASCADE)
    doi = models.TextField(blank=True, null=True)
    paper_name = models.TextField(blank=True, null=True)
    year_published = models.IntegerField(null=True)
    enrichment_kit = models.ForeignKey(EnrichmentKit, on_delete=CASCADE, blank=True, null=True)
    sequencer = models.ForeignKey(Sequencer, on_delete=CASCADE, blank=True, null=True)
    seq_details = models.TextField(blank=True, null=True)
    file_type = models.TextField(blank=True, null=True)
    file_count = models.IntegerField(blank=True, default=0)
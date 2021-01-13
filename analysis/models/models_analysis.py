from collections import defaultdict
from typing import Tuple, Dict, List

from django.apps import apps
from django.contrib.auth.models import User
from django.db import models
from django.db.models import Model, Q
from django.db.models.deletion import SET_NULL, CASCADE, SET_DEFAULT, PROTECT, ProtectedError
from django.db.models.signals import pre_delete
from django.dispatch import receiver
from django.urls.base import reverse
from django.utils import timezone
from django_extensions.db.models import TimeStampedModel
from lazy import lazy

from analysis.models.enums import AnalysisType, AnalysisTemplateType
from annotation.models import AnnotationVersion
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsAutoInitialSaveMixin
from library.guardian_utils import admin_bot, assign_permission_to_user_and_groups
from snpdb.models import CustomColumnsCollection, CustomColumn, \
    UserSettings, AbstractNodeCountSettings, Sample
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_enums import BuiltInFilters


class Analysis(GuardianPermissionsAutoInitialSaveMixin, TimeStampedModel):
    # Changing some analysis settings alters node editors/grids - and we need to increment version to expire node cache
    VERSION_BUMP_FIELDS = ["custom_columns_collection", "default_sort_by_column"]

    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    version = models.IntegerField(default=0)  # By bumping this we can invalidate node caches
    # TODO: Remove 'analysis_type' by creating legacy 'AnalysisTemplateSnapshot'
    analysis_type = models.CharField(max_length=1, choices=AnalysisType.choices, null=True, blank=True)
    user = models.ForeignKey(User, on_delete=CASCADE)
    name = models.TextField()
    description = models.TextField(null=True, blank=True)
    custom_columns_collection = models.ForeignKey(CustomColumnsCollection,
                                                  default=CustomColumnsCollection.get_system_default_id,
                                                  on_delete=SET_DEFAULT)
    default_sort_by_column = models.ForeignKey(CustomColumn, null=True, blank=True, on_delete=SET_NULL)
    show_igv_links = models.BooleanField(default=True)
    analysis_panel_fraction = models.FloatField(default=0.25)
    annotation_version = models.ForeignKey(AnnotationVersion, null=True, on_delete=SET_NULL)
    lock_input_sources = models.BooleanField(default=False)
    visible = models.BooleanField(default=True)
    template_type = models.CharField(max_length=1, choices=AnalysisTemplateType.choices, null=True, blank=True)

    def __str__(self):
        name = self.name or f"Analysis {self.pk}"
        if self.description:
            name = f"{name} ({self.description})"
        return name

    @lazy
    def template(self):
        """ Works for both a template and a snapshot """
        try:
            return self.analysistemplate
        except AnalysisTemplate.DoesNotExist:
            try:
                return self.analysistemplateversion.template
            except AnalysisTemplateVersion.DoesNotExist:
                return None

    @lazy
    def last_lock(self):
        return self.analysislock_set.order_by("pk").last()

    def is_locked(self):
        is_snapshot = AnalysisTemplateVersion.objects.filter(analysis_snapshot=self).exists()
        return is_snapshot or (self.last_lock and self.last_lock.locked)

    def can_unlock(self, user):
        """ Use parent to see if we have Guardian permissions to write """
        is_snapshot = AnalysisTemplateVersion.objects.filter(analysis_snapshot=self).exists()
        return (not is_snapshot) and super().can_write(user)

    def lock_history(self):
        return self.analysislock_set.order_by("pk")

    def can_write(self, user):
        """ Disable modification when locked """
        if super().can_write(user):
            return not self.is_locked()
        return False

    def get_absolute_url(self):
        return reverse('analysis', kwargs={"analysis_id": self.pk})

    @property
    def gene_annotation_release(self):
        return self.annotation_version.variant_annotation_version.gene_annotation_release

    @classmethod
    def get_listing_url(cls):
        return reverse('analyses')

    def check_valid(self):

        def check_field(field):
            value = getattr(self, field)
            if value is None:
                msg = f"Analysis setting '{field}' is not set. "
                msg += "<a href='javascript:analysisSettings()'>Open Analysis Settings</a>"
                raise ValueError(msg)

        check_field('custom_columns_collection')
        check_field('annotation_version')

    def is_valid(self):
        try:
            self.check_valid()
            return True
        except ValueError:
            return False

    def set_defaults_and_save(self, user):
        self.user = user
        self.annotation_version = AnnotationVersion.latest(self.genome_build)

        # Initial config from user settings
        user_settings = UserSettings.get_for_user(user)
        self.custom_columns_collection = user_settings.columns
        self.default_sort_by_column = user_settings.default_sort_by_column
        self.save()

        default_node_count_config = user_settings.get_node_count_settings_collection()
        if default_node_count_config:
            self.set_node_count_types(default_node_count_config.get_node_count_filters())
            self.save()

    def get_node_count_types(self):
        """ Return ordered array of tuples """

        node_count_labels = []
        try:
            node_count_config = self.analysisnodecountconfiguration
            for nc in node_count_config.analysisnodecountconfigrecord_set.all().order_by("sort_order"):
                node_count_labels.append(nc.built_in_filter)
        except AnalysisNodeCountConfiguration.DoesNotExist:
            node_count_labels = BuiltInFilters.DEFAULT_NODE_COUNT_FILTERS

        return AbstractNodeCountSettings.get_types_from_labels(node_count_labels)

    def set_node_count_types(self, node_counts_array):
        """ OK to pass in arrays of blank (as no entries does), will just delete everything """
        node_count_config, _ = AnalysisNodeCountConfiguration.objects.get_or_create(analysis=self)
        record_set = node_count_config.analysisnodecountconfigrecord_set

        AbstractNodeCountSettings.save_count_configs_from_array(record_set, node_counts_array)

    def get_samples(self) -> List[Sample]:
        samples = set()
        for node in self.analysisnode_set.filter(analysisnode_parent__isnull=True).select_subclasses():
            samples.update(node.get_samples_from_node_only_not_ancestors())
        return list(sorted(samples))

    def clone(self, user: User = None):
        """ user - if provided, assign new ownership """
        from analysis.models import AnalysisNode, AnalysisEdge
        from analysis.models.nodes.node_utils import get_toposorted_nodes

        # TODO: Import all the other funcs locally as deps screwed up

        now = timezone.now()
        analysis_id = self.pk
        analysis_copy = self
        analysis_copy.pk = None  # to save as new
        analysis_copy.name = f"Copy of {self.name}"
        analysis_copy.created = now
        analysis_copy.modified = now
        if user:
            analysis_copy.user = user
        analysis_copy.template_type = None
        analysis_copy.save()

        nodes_qs = AnalysisNode.objects.filter(analysis_id=analysis_id).select_subclasses()
        topo_sorted = get_toposorted_nodes(nodes_qs)

        # Still need to copy the edges. Probably need to build up a table of old vs new
        old_version_and_new_mapping: Dict[Tuple] = {}
        old_new_map = {}

        nodes = []
        for group in topo_sorted:
            for node in group:
                old_id = node.pk
                old_version = node.version
                node.analysis = analysis_copy
                try:
                    copy = node.save_clone()
                except NotImplementedError as nie:
                    raise NotImplementedError(f"Could not clone node {old_id}") from nie
                copy.adjust_cloned_parents(old_new_map)
                copy.save()

                nodes.append(copy)

                old_version_and_new_mapping[old_id] = (old_version, copy.pk)
                old_new_map[old_id] = copy

        for edge in AnalysisEdge.objects.filter(parent__analysis_id=analysis_id):
            edge.parent_id = old_version_and_new_mapping[edge.parent_id][1]
            edge.child_id = old_version_and_new_mapping[edge.child_id][1]
            edge.pk = None
            edge.save()

        # re-save all the nodes to re-validate them after all joined up
        for node in nodes:
            node.save()

        for av in AnalysisVariable.objects.filter(node__analysis_id=analysis_id):
            new_node = old_new_map[av.node.pk]
            av.pk = None
            av.node = new_node
            av.save()

        return analysis_copy

    def get_warnings(self, user: User) -> List[str]:
        warnings = []
        if self.lock_input_sources:
            warnings.append("INPUT LOCKED - cannot create new input source nodes.")
        if not self.can_write(user) and not self.is_locked():  # Locked has own icon, no need for warning
            warnings.append("This analysis is READ-ONLY - you do not have write permission to modify anything")
        return warnings


@receiver(pre_delete, sender=Analysis)
def pre_delete_analysis(sender, instance, *args, **kwargs):
    """ Delete analysis template if not used for a run, otherwise soft delete it """
    try:
        analysis_template = instance.analysistemplate
        analysis_template.delete_or_soft_delete()
    except AnalysisTemplate.DoesNotExist:
        pass


class AnalysisLock(models.Model):
    """ Users can lock/unlock however much they like - the last one set is the current status """
    analysis = models.ForeignKey(Analysis, on_delete=CASCADE)
    locked = models.BooleanField()
    user = models.ForeignKey(User, on_delete=CASCADE)
    date = models.DateTimeField()


class AnalysisNodeCountConfiguration(models.Model):
    analysis = models.OneToOneField(Analysis, on_delete=CASCADE)


class AnalysisNodeCountConfigRecord(AbstractNodeCountSettings):
    node_count_config = models.ForeignKey(AnalysisNodeCountConfiguration, on_delete=CASCADE)


class AnalysisVariable(models.Model):
    """ Node fields are extracted to an AnalysisVariable so they can be set """
    node = models.ForeignKey("AnalysisNode", on_delete=CASCADE)
    field = models.TextField()
    class_name = models.TextField()

    class Meta:
        unique_together = ('node', 'field')
        ordering = ['node__analysis', 'node', 'field']

    @staticmethod
    def get_node_field_class_name(node, field):
        field_value = getattr(node, field)
        if isinstance(field_value, Model):
            # Easy way - just ask for Python type
            class_name = field_value._meta.label
        else:
            # Hard way, try and determine from model meta
            PYTHON_CLASS = {"IntegerField": "int",
                            "TextField": "str"}
            f = node._meta.get_field(field)
            internal_type = f.get_internal_type()
            if internal_type == "ForeignKey":
                class_name = f.related_model._meta.label
            else:
                class_name = PYTHON_CLASS.get(internal_type)
                if class_name is None:
                    raise ValueError(f"Don't know how to handle {node}.{field} (internal type: {internal_type})")
        return class_name

    def __str__(self):
        return f"{self.node_id}/{self.field}"


class AnalysisTemplate(GuardianPermissionsAutoInitialSaveMixin, TimeStampedModel):
    """ A snapshot of an analysis - locked-down to be used as a template """
    name = models.TextField(unique=True)
    description = models.TextField(blank=True)
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)
    analysis = models.OneToOneField(Analysis, null=True, on_delete=SET_NULL)  # deleted or soft deleted in pre_delete
    deleted = models.BooleanField(default=False)

    @classmethod
    def get_permission_class(cls):
        # Delegate to Analysis for permissions
        return Analysis

    def get_permission_object(self):
        # Delegate to Analysis for permissions
        return self.analysis

    @classmethod
    def _filter_from_permission_object_qs(cls, queryset):
        return cls.objects.filter(analysis__in=queryset)

    def delete_or_soft_delete(self):
        try:
            self.delete()  # Will succeed if no analysis runs have been created
        except ProtectedError:
            self.deleted = True
            self.save()

    @property
    def active(self):
        return self.analysistemplateversion_set.filter(active=True).first()

    def latest_version(self):
        last = self.analysistemplateversion_set.order_by("-pk").first()
        if last:
            return last.version
        return 0

    @classmethod
    def filter_for_user(cls, user, queryset=None, **kwargs):
        """ Hides deleted objects """
        qs = super().filter_for_user(user, queryset=queryset, **kwargs)
        return qs.filter(deleted=False)

    @staticmethod
    def filter(user: User, requires_sample_somatic=None, requires_sample_gene_list=None, class_name=None, atv_kwargs=None):
        """ requires_sample_somatic/requires_sample_gene_list - leave None for all """
        if atv_kwargs is None:
            atv_kwargs = {}

        template_versions_qs = AnalysisTemplateVersion.objects.filter(active=True, **atv_kwargs)

        if requires_sample_somatic is not None:
            template_versions_qs = template_versions_qs.filter(requires_sample_somatic=requires_sample_somatic)

        if requires_sample_gene_list is not None:
            template_versions_qs = template_versions_qs.filter(requires_sample_gene_list=requires_sample_gene_list)

        if class_name:
            # Restrict to template versions who's variables are all of class_name
            # unsupported_variables = AnalysisVariable.objects.exclude(class_name=class_name)
            q_class_name = Q(analysis_snapshot__analysisnode__analysisvariable__class_name=class_name)
            template_versions_qs = template_versions_qs.exclude(~q_class_name)

        return AnalysisTemplate.filter_for_user(user).filter(analysistemplateversion__in=template_versions_qs)

    def default_name_template(self):
        """ The initial analysis_name_template in form for save version """
        analysis_name_template = "%(template)s for %(input)s"  # default
        if self.active:
            # Use last value if available
            analysis_name_template = self.active.analysis_name_template
        return analysis_name_template

    def __str__(self):
        s = self.name
        if self.deleted:
            s += " (deleted)"
        return s


class AnalysisTemplateVersion(TimeStampedModel):
    template = models.ForeignKey(AnalysisTemplate, on_delete=CASCADE)
    version = models.IntegerField()
    analysis_name_template = models.TextField(null=True)  # Python string template
    analysis_snapshot = models.OneToOneField(Analysis, null=True, on_delete=PROTECT)
    active = models.BooleanField(default=True)
    appears_in_autocomplete = models.BooleanField(default=True)
    appears_in_links = models.BooleanField(default=False)
    requires_sample_somatic = models.BooleanField(default=False)
    requires_sample_gene_list = models.BooleanField(default=False)

    class Meta:
        unique_together = ('template', 'version')

    def __str__(self):
        return f"{self.template} v.{self.version}"


class AnalysisTemplateRun(TimeStampedModel):
    template_version = models.ForeignKey(AnalysisTemplateVersion, null=True, on_delete=PROTECT)
    analysis = models.OneToOneField(Analysis, on_delete=CASCADE)  # Created new analysis

    @staticmethod
    def create(analysis_template: AnalysisTemplate, genome_build: GenomeBuild, user: User = None):
        if user is None:
            user = admin_bot()

        template_version = analysis_template.active
        analysis = template_version.analysis_snapshot.clone()
        analysis.user = user
        analysis.genome_build = genome_build
        analysis.template_type = None
        analysis.visible = True
        analysis.name = f"TemplateRun from {analysis_template.name}"  # Will be set in populate arguments
        analysis.save()

        assign_permission_to_user_and_groups(user, analysis)
        return AnalysisTemplateRun.objects.create(template_version=template_version, analysis=analysis)

    def populate_arguments(self, data: Dict):
        for av in AnalysisVariable.objects.filter(node__analysis=self.analysis):
            if obj := data.get(av.field):
                # If variable is a model instance - need to load object
                variable_class = apps.get_model(av.class_name)
                if isinstance(variable_class, type(Model)):
                    if isinstance(obj, (str, int)):
                        obj = variable_class.objects.get(pk=obj)
                        # Guardian security check (if available)
                        if check_can_view := getattr(obj, "check_can_view", None):
                            check_can_view(self.analysis.user)

                object_class_name = obj._meta.label
                object_pk = None
                error = None
                if object_class_name == av.class_name:
                    if isinstance(obj, Model):
                        object_pk = obj.pk
                else:
                    error = f"{av} type was: {object_class_name} expected: {av.class_name}"

                AnalysisTemplateRunArgument.objects.create(template_run=self,
                                                           variable=av,
                                                           object_pk=object_pk,
                                                           value=str(obj),
                                                           error=error)

    def _get_unset_or_errored_variables(self):
        av_qs = AnalysisVariable.objects.filter(node__analysis=self.analysis)
        return av_qs.exclude(analysistemplaterunargument__error__isnull=True)

    def is_populated(self) -> bool:
        return not self._get_unset_or_errored_variables().exists()

    def check_populated(self):
        """ Throws exception if all AnalysisVariables not populated w/o error """
        unpopulated = self._get_unset_or_errored_variables()
        if unpopulated_fields := list(unpopulated.values_list("field", flat=True)):
            raise ValueError(f"Missing/Errored analysis variables: {', '.join(unpopulated_fields)}")

    def populate_analysis_name(self):
        """ Populate analysis_name_template with params based on AnalysisVariable fields, and the magic values:
                * template - TemplateVersion string representation
                * input - 1st we find of "sample", "cohort", "trio" """

        params = {"template": str(self.template_version)}
        for arg in self.analysistemplaterunargument_set.all():
            params[arg.variable.field] = arg.value

        for field in ["sample", "cohort", "trio", "pedigree"]:
            if field in params:
                params["input"] = field
                break

        try:
            self.analysis.name = self.template_version.analysis_name_template % params
            self.analysis.save()
        except (TypeError, ValueError, KeyError):
            pass

    def get_node_field_values(self):
        self.check_populated()

        node_field_values = defaultdict(dict)

        for arg in self.analysistemplaterunargument_set.all():
            # Get the value
            if arg.object_pk:
                klass = apps.get_model(arg.variable.class_name)
                value = klass.objects.get(pk=arg.object_pk)
                if validate_func := getattr(value, "check_can_view", None):
                    validate_func(self.analysis.user)
            else:
                value = arg.value

            node_field_values[arg.variable.node_id][arg.variable.field] = value
        return node_field_values


class AnalysisTemplateRunArgument(models.Model):
    """ Arguments used to create a template copy """
    template_run = models.ForeignKey(AnalysisTemplateRun, on_delete=CASCADE)
    variable = models.OneToOneField(AnalysisVariable, on_delete=CASCADE)
    object_pk = models.TextField(null=True, blank=True)  # for variable.class_name object
    value = models.TextField()
    error = models.TextField(null=True)


class SampleAnalysisTemplateRun(models.Model):
    """ Analysis (settings.ANALYSIS_TEMPLATES_AUTO_SAMPLE) automatically run once against a sample  """
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    analysis_template_run = models.ForeignKey(AnalysisTemplateRun, on_delete=CASCADE)

    class Meta:
        unique_together = ('sample', 'analysis_template_run')

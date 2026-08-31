import operator
from functools import reduce
from typing import Optional

from auditlog.registry import auditlog
from cache_memoize import cache_memoize
from django.db import models
from django.db.models import CASCADE, SET_NULL
from django.db.models.query_utils import Q

from analysis.models import GeneCoverageMixin
from analysis.models.enums import SampleNodeSourceLevel
from analysis.models.nodes.analysis_node import (
    AnalysisNode,
    NodeAlleleFrequencyFilter,
    NodeAuditLogMixin,
    NodeVCFFilter,
    annotate_and_filter_queryset,
    queryset_to_pk_in_q,
)
from analysis.models.nodes.cohort_mixin import SampleMixin
from analysis.models.nodes.node_display import NodeChip, NodeIcon
from analysis.models.nodes.stats_cache import (
    get_cached_label_count_for_cohort,
    get_handler_for_node,
)
from genes.models import SampleGeneList
from library.constants import DAY_SECS
from patients.models import Extraction, Patient, Specimen
from patients.models_enums import NucleicAcid, Sex, Zygosity
from patients.sample_grouping import SampleGroup, get_extraction_sample_group
from snpdb.models import Sample


class SampleNode(SampleMixin, GeneCoverageMixin, AnalysisNode):
    """ One entry point for a single VCF sample, or every sample belonging to a grouping object
        (an Extraction - the DNA arm's small variant, CNV and exon CNV calls in one node).

        Group levels resolve their samples at query time and OR a per-sample subquery, so nothing
        is pre-built to go stale as samples are linked to an extraction.

        Use restrict_to_qc_gene_list to keep track of that as sample_gene_list is cleared when a sample
        changes including in an AnalysisTemplates """
    source_level = models.CharField(max_length=1, choices=SampleNodeSourceLevel.choices,
                                    default=SampleNodeSourceLevel.SAMPLE)
    sample = models.ForeignKey(Sample, null=True, on_delete=SET_NULL)
    extraction = models.ForeignKey(Extraction, null=True, blank=True, on_delete=SET_NULL)
    specimen = models.ForeignKey(Specimen, null=True, blank=True, on_delete=SET_NULL)
    patient = models.ForeignKey(Patient, null=True, blank=True, on_delete=SET_NULL)
    # When setting sample, if restrict_to_qc_gene_list = True, sample_gene_list is set to active sample gene list
    sample_gene_list = models.ForeignKey(SampleGeneList, null=True, blank=True, on_delete=SET_NULL)
    has_gene_coverage = models.BooleanField(null=True)
    min_ad = models.IntegerField(default=0)
    min_dp = models.IntegerField(default=0)
    min_gq = models.IntegerField(default=0)
    max_pl = models.IntegerField(null=True, blank=True)
    zygosity_ref = models.BooleanField(default=False)
    zygosity_het = models.BooleanField(default=True)
    zygosity_hom = models.BooleanField(default=True)
    zygosity_unk = models.BooleanField(default=False)
    restrict_to_qc_gene_list = models.BooleanField(default=False)

    FIELDS_THAT_CHANGE_QUERYSET = ("min_ad", "min_dp", "min_gq", "max_pl", "restrict_to_qc_gene_list")
    SAMPLE_FIELD_MAPPINGS = [("ad", "allele_depth"),
                             ("dp", "read_depth"),
                             ("gq", "genotype_quality")]
    # Levels that resolve to a set of samples. Specimen/Patient join here with their own resolvers
    SOURCE_LEVEL_FIELDS = {
        SampleNodeSourceLevel.SAMPLE: "sample",
        SampleNodeSourceLevel.EXTRACTION: "extraction",
        SampleNodeSourceLevel.SPECIMEN: "specimen",
        SampleNodeSourceLevel.PATIENT: "patient",
    }
    IMPLEMENTED_SOURCE_LEVELS = {SampleNodeSourceLevel.SAMPLE, SampleNodeSourceLevel.EXTRACTION}
    # Above this the individual VCF chips stop fitting on the card, so collapse them into a count
    MAX_VCF_CHIPS = 3
    min_inputs = 0
    max_inputs = 0

    # ── What samples this node reads ──────────────────────────────────────────

    @property
    def is_group_level(self) -> bool:
        return self.source_level != SampleNodeSourceLevel.SAMPLE

    @property
    def source_field(self) -> str:
        return self.SOURCE_LEVEL_FIELDS[self.source_level]

    def get_source_object(self):
        """ The object this node groups on (a Sample at sample level) """
        return getattr(self, self.source_field)

    def get_sample_group(self) -> SampleGroup:
        """ The one place samples are resolved - query, warnings, badge and preview all come from here """
        if self.source_level == SampleNodeSourceLevel.SAMPLE:
            return SampleGroup(samples=[self.sample] if self.sample else [])
        if self.source_level == SampleNodeSourceLevel.EXTRACTION:
            return get_extraction_sample_group(self.analysis.user, self.extraction, self.analysis.genome_build)
        return SampleGroup()

    def get_source_samples(self) -> list[Sample]:
        return self.get_sample_group().samples

    @cache_memoize(DAY_SECS, args_rewrite=lambda s: (s.pk, s.version))
    def get_source_vcf_names(self) -> list[str]:
        """ For the canvas chips - get_node_chips runs for every node on every analysis render """
        return [str(vcf) for vcf in self.get_sample_group().vcfs]

    def _get_sample(self) -> Optional[Sample]:
        # Overrides SampleMixin - a group node's genotype joins are per resolved sample
        if self.is_group_level:
            return None
        return self.sample

    def _get_proband_sample_for_node(self) -> Optional[Sample]:
        """ The study's sample where that is unambiguous - it's what lets downstream AncestorSampleMixin
            nodes (Zygosity, MOI, GeneList, AlleleFrequency) auto-populate rather than give up """
        samples = self.get_source_samples()
        if len(samples) == 1:
            return samples[0]
        dna_samples = [s for s in samples
                       if s.extraction and s.extraction.nucleic_acid_source == NucleicAcid.DNA]
        if len(dna_samples) == 1:
            return dna_samples[0]
        return None

    def _get_cohorts_and_sample_visibility_for_node(self):
        if not self.is_group_level:
            return super()._get_cohorts_and_sample_visibility_for_node()

        cohorts = []
        visibility = {}
        for sample in self.get_source_samples():
            cohort = sample.vcf.cohort
            if cohort not in cohorts:
                cohorts.append(cohort)
            visibility[sample] = sample.has_genotype
        return cohorts, visibility

    def _get_annotation_kwargs_for_node(self, **kwargs) -> dict:
        annotation_kwargs = super()._get_annotation_kwargs_for_node(**kwargs)
        if self.is_group_level:
            kwargs["override"] = False
            for sample in self.get_source_samples():
                annotation_kwargs.update(self._get_sample_annotation_kwargs(sample, **kwargs))
        return annotation_kwargs

    @staticmethod
    def _get_sample_annotation_kwargs(sample: Sample, **kwargs) -> dict:
        """ The genotype join for one sample's VCF, plus its zygosity alias """
        annotation_kwargs = dict(sample.cohort_genotype_collection.get_annotation_kwargs(**kwargs))
        annotation_kwargs.update(sample.get_annotation_kwargs(**kwargs))
        return annotation_kwargs

    def _get_cache_key(self) -> str:
        """ CohortMixin folds in a single CGC pk so a reloaded VCF invalidates the cached query.
            A group spans collections, so every sample's CGC and cohort version has to be in there """
        cache_key = super()._get_cache_key()
        if self.is_group_level:
            parts = [cache_key]
            for sample in self.get_source_samples():
                cohort = sample.vcf.cohort
                parts.extend([str(sample.pk), str(cohort.version), str(sample.cohort_genotype_collection.pk)])
            cache_key = "_".join(parts)
        return cache_key

    # ── Thresholds, per sample ────────────────────────────────────────────────

    def get_sample_thresholds(self, sample: Sample) -> dict:
        """ The node's own values, overridden where the user set a row for this sample. A sample
            attached to the extraction after the node was configured gets the node defaults """
        thresholds = {f: getattr(self, f) for f in ("min_ad", "min_dp", "min_gq", "max_pl")}
        if override := self.samplenodesamplethreshold_set.filter(sample=sample).first():
            thresholds.update({f: getattr(override, f) for f in thresholds})
        return thresholds

    # ── Query ─────────────────────────────────────────────────────────────────

    def _get_zygosities(self, sample: Optional[Sample] = None):
        # No genotype - all will be unknown
        if sample and sample.has_genotype is False:
            return [Zygosity.UNKNOWN_ZYGOSITY]

        ZYGOSITY_LOOKUP = [
            (self.zygosity_ref, Zygosity.HOM_REF),
            (self.zygosity_het, Zygosity.HET),
            (self.zygosity_hom, Zygosity.HOM_ALT),
            (self.zygosity_unk, Zygosity.UNKNOWN_ZYGOSITY),
        ]
        return [z for s, z in ZYGOSITY_LOOKUP if s]

    def get_zygosity_description(self):
        zygosities = self._get_zygosities(self.sample)
        if len(zygosities) != len(Zygosity.CHOICES):  # Subset
            zyg_dict = dict(Zygosity.CHOICES)
            zyg = ','.join([zyg_dict[z] for z in zygosities])
        else:
            zyg = ''  # Any
        return zyg

    def _get_sample_arg_q_dict(self, sample: Sample) -> dict[Optional[str], dict[str, Q]]:
        """ The genotype filters for one sample - zygosity, thresholds and allele frequency """
        arg_q_dict = {}
        # _get_zygosities handles no genotype (UNKNOWN)
        if zygosity := self._get_zygosities(sample):
            alias, field = sample.get_cohort_genotype_alias_and_field("zygosity")
            q = Q(**{f"{field}__in": zygosity})
            arg_q_dict[alias] = {str(q): q}

        if sample.has_genotype:
            thresholds = self.get_sample_thresholds(sample)
            for node_field, ov_field in self.SAMPLE_FIELD_MAPPINGS:
                if min_value := thresholds[f"min_{node_field}"]:
                    alias, ov_path = sample.get_cohort_genotype_alias_and_field(ov_field)
                    q = Q(**{f"{ov_path}__gte": min_value})
                    self.merge_arg_q_dicts(arg_q_dict, {alias: {str(q): q}})

            if (max_pl := thresholds["max_pl"]) is not None:
                alias, pl_path = sample.get_cohort_genotype_alias_and_field("phred_likelihood")
                q = Q(**{f"{pl_path}__lte": max_pl})
                self.merge_arg_q_dicts(arg_q_dict, {alias: {str(q): q}})

            if sample_arg_q_dict := NodeAlleleFrequencyFilter.get_sample_arg_q_dict(self, sample):
                self.merge_arg_q_dicts(arg_q_dict, sample_arg_q_dict)

        self.merge_arg_q_dicts(arg_q_dict, self.get_vcf_locus_filters_arg_q_dict_for_sample(sample))
        return arg_q_dict

    def _get_qc_gene_list_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        """ Sample level only - a QC gene list has no meaning over a set of samples """
        if not self.restrict_to_qc_gene_list:
            return {}
        if self.sample_gene_list:
            q_hash = f"gene_list_{self.sample_gene_list.gene_list_id}"
            q = self.sample_gene_list.gene_list.get_q(self.analysis.annotation_version.variant_annotation_version)
        else:
            q = self.q_none()  # Safety - don't show anything if missing
            q_hash = str(q)
        return {None: {q_hash: q}}

    def _get_sample_pk_q(self, sample: Sample) -> Q:
        """ One sample's variants as pk IN (subquery). Annotated with only its own VCF's genotype join,
            so the subquery doesn't drag the group's other outer joins through with it """
        qs = self._get_model_queryset()
        a_kwargs = self._get_sample_annotation_kwargs(sample)
        qs, q_list = annotate_and_filter_queryset(qs, a_kwargs, self._get_sample_arg_q_dict(sample))
        if q_list:
            qs = qs.filter(reduce(operator.and_, q_list))
        return queryset_to_pk_in_q(qs)

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        samples = self.get_source_samples()
        if not samples:
            q_none = self.q_none()
            return {None: {str(q_none): q_none}}

        if len(samples) == 1:
            # A single sample degenerates to the alias path - byte for byte the query a sample level
            # node produces today
            arg_q_dict = self._get_sample_arg_q_dict(samples[0])
            self.merge_arg_q_dicts(arg_q_dict, self._get_qc_gene_list_arg_q_dict())
            return arg_q_dict

        # pk__in subqueries don't fan out rows the way joins do, so no distinct() is needed
        q = reduce(operator.or_, [self._get_sample_pk_q(sample) for sample in samples])
        return {None: {self._get_node_q_hash(): q}}

    # ── VCF FILTER - one node level selection, resolved per VCF ───────────────

    def get_vcf_locus_filters_arg_q_dict_for_sample(self, sample: Sample) -> dict[Optional[str], dict[str, Q]]:
        if not self.has_filters:
            return {}
        alias = sample.cohort_genotype_collection.cohortgenotype_alias
        return self._get_vcf_locus_filters_arg_q_dict(sample.vcf, alias)

    def get_vcf_locus_filters_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        # Group levels apply these per sample inside each subquery
        if self.is_group_level:
            return {}
        return super().get_vcf_locus_filters_arg_q_dict()

    def get_vcf_locus_filter_vcfs(self) -> list:
        if self.is_group_level:
            return self.get_sample_group().vcfs
        return super().get_vcf_locus_filter_vcfs()

    @property
    def has_filters(self):
        if self.is_group_level:
            return any(vcf.vcffilter_set.exists() for vcf in self.get_sample_group().vcfs)
        return super().has_filters

    def get_filter_code(self):
        if self.is_group_level:
            # Filter ids are node level - they translate into each VCF's own codes at query time
            if filter_ids := NodeVCFFilter.get_filter_ids(self):
                if filter_ids == {None}:  # PASS only
                    return 1
                return 2
            return 0
        return super().get_filter_code()

    def _get_filters_cohort_genotype_collections(self) -> list:
        if self.is_group_level:
            if not self.has_filters:
                return []
            cgcs = []
            for sample in self.get_source_samples():
                cgc = sample.cohort_genotype_collection
                if cgc not in cgcs:
                    cgcs.append(cgc)
            return cgcs
        return super()._get_filters_cohort_genotype_collections()

    # ── Node display ──────────────────────────────────────────────────────────

    def _get_sample_method_summary(self, sample: Sample) -> str:
        sample_description = f"Sample {sample.name}"
        thresholds = self.get_sample_thresholds(sample)
        for node_field, _ in self.SAMPLE_FIELD_MAPPINGS:
            if min_value := thresholds[f"min_{node_field}"]:
                sample_description += f" {node_field.upper()}>={min_value}"
        if (max_pl := thresholds["max_pl"]) is not None:
            sample_description += f" PL<={max_pl}"

        if sample.has_genotype:
            if zyg := self.get_zygosity_description():
                sample_description += f" ({zyg})"
        return f"{sample_description}, from VCF {sample.vcf}"

    def _get_method_summary(self):
        samples = self.get_source_samples()
        if not samples:
            return 'No sample selected'

        summaries = [self._get_sample_method_summary(sample) for sample in samples]
        if self.is_group_level:
            header = f"{self.get_source_level_display()} {self.get_source_object()} - {len(samples)} sample(s):"
            return "<br>".join([header] + summaries)
        return summaries[0]

    def get_node_name(self):
        name_parts = []
        if self.is_group_level:
            if source_object := self.get_source_object():
                name_parts.append(str(source_object))
                if (num_samples := len(self.get_source_samples())) != 1:
                    name_parts.append(f"\n({num_samples} samples)")
        elif self.sample:
            name_parts.append(self.sample.name)
            if self.sample.has_genotype:
                if zyg := self.get_zygosity_description():
                    name_parts.append(f"\n({zyg})")

        if name_parts:
            if filter_description := self.get_filter_description():
                name_parts.append(f"\n({filter_description})")

        return "\n".join(name_parts)

    @staticmethod
    def get_help_text() -> str:
        return "Variants from a VCF sample, or every sample of an extraction (usually one genotype - " \
               "patient, cell or organism)"

    def get_gene_lists(self) -> list:
        """ Used for gene coverage """
        gene_lists = []
        if self.sample_gene_list and self.restrict_to_qc_gene_list:
            gene_lists.append(self.sample_gene_list.gene_list)
        return gene_lists

    def _has_filters_that_affect_label_counts(self):
        # This returns None if no filters to apply
        if NodeAlleleFrequencyFilter.get_sample_arg_q_dict(self, self.sample):
            return True

        return any(getattr(self, f) for f in self.FIELDS_THAT_CHANGE_QUERYSET)

    def _get_cached_label_count(self, label):
        """ Input counts can be static, so use cached CohortGenotype*Stats if we can. """
        if self._has_filters_that_affect_label_counts():
            return None  # Have to do counts
        if self.sample is None:
            return None

        filter_key = get_handler_for_node(self).filter_key_for_node(self)

        filter_code = self.get_filter_code()
        if filter_code not in (0, 1):
            return None  # Custom filters defeat the cache

        zygosities = [self.zygosity_ref, self.zygosity_het, self.zygosity_hom, self.zygosity_unk]
        if not self.sample.has_genotype:
            zygosities = [True] * len(zygosities)  # Show everything

        return get_cached_label_count_for_cohort(
            cohort=self.sample.vcf.cohort,
            sample=self.sample,
            filter_key=filter_key,
            annotation_version=self.analysis.annotation_version,
            passing_filter=bool(filter_code),
            zygosities=zygosities,
            label=label,
        )

    @staticmethod
    def get_node_class_label():
        return "Sample"

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(symbol="node-icon-sample")

    def get_node_icon(self) -> NodeIcon:
        """ Pedigree notation on the badge - square/circle for sex, struck through if deceased """
        patient = self.sample.patient if self.sample else None
        if patient is None:
            return self.get_node_class_icon()
        sex = "female" if patient.sex == Sex.FEMALE else "male"
        deceased = "-deceased" if patient.deceased else ""
        return NodeIcon(symbol=f"node-icon-sample-{sex}{deceased}")

    def get_node_chips(self) -> list[NodeChip]:
        """ Why a group node returns fewer rows than it used to - a VCF was archived - is a question
            only the canvas can answer, so show what the node is actually reading """
        chips = super().get_node_chips()
        if not self.is_group_level:
            return chips

        if self.source_level == SampleNodeSourceLevel.EXTRACTION and self.extraction:
            extraction_id = self.extraction.reference_id or self.extraction.external_pk or f"({self.extraction.pk})"
            chips.append(NodeChip(text=extraction_id, icon="fa-solid fa-vial", title=str(self.extraction)))

        vcf_names = self.get_source_vcf_names()
        if len(vcf_names) > self.MAX_VCF_CHIPS:
            chips.append(NodeChip(text=f"VCF x{len(vcf_names)}", icon="fa-solid fa-file-lines",
                                  title="\n".join(vcf_names)))
        else:
            for vcf_name in vcf_names:
                chips.append(NodeChip(text="VCF", icon="fa-solid fa-file-lines", title=vcf_name))
        return chips

    def _get_configuration_check_cohorts(self) -> list:
        if self.is_group_level:
            cohorts, _ = self._get_cohorts_and_sample_visibility_for_node()
            return cohorts
        return super()._get_configuration_check_cohorts()

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if self.source_level not in self.IMPLEMENTED_SOURCE_LEVELS:
            errors.append(f"{self.get_source_level_display()} level is not supported yet.")
        elif self.get_source_object() is None:
            errors.append(f"No {self.get_source_level_display().lower()} selected.")
        elif self.is_group_level:
            if not self.get_source_samples():
                errors.append(f"{self.get_source_object()} has no samples in this analysis' genome build.")
        else:
            errors.extend(self._get_genome_build_errors("sample", self.sample.genome_build))
        if self.restrict_to_qc_gene_list and self.sample_gene_list is None:
            errors.append("Restricted to Sample Gene List, but none specified!")
        return errors

    def get_warnings(self) -> list[str]:
        """ A sample the group silently drops - archived VCF, wrong build - is the failure mode that
            bites here, since the user didn't hand pick the samples """
        return super().get_warnings() + self.get_sample_group().warnings

    def save_clone(self):
        thresholds = list(self.samplenodesamplethreshold_set.all())

        copy = super().save_clone()
        for threshold in thresholds:
            threshold.pk = None
            threshold.node = copy
            threshold.save()
        return copy



class SampleNodeSampleThreshold(NodeAuditLogMixin, models.Model):
    """ Per sample threshold overrides, so an extraction can apply a different min_ad per caller
        (sapath#301). The node's own fields stay the value for any sample without a row - so a sample
        reconcile attaches after the node was configured gets the node default rather than nothing """
    node = models.ForeignKey(SampleNode, on_delete=CASCADE)
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    min_ad = models.IntegerField(default=0)
    min_dp = models.IntegerField(default=0)
    min_gq = models.IntegerField(default=0)
    max_pl = models.IntegerField(null=True, blank=True)

    class Meta:
        unique_together = ("node", "sample")

    def _get_node(self):
        return self.node

    def __str__(self):
        return f"{self.sample}: AD>={self.min_ad} DP>={self.min_dp} GQ>={self.min_gq} PL<={self.max_pl}"


auditlog.register(SampleNode)
auditlog.register(SampleNodeSampleThreshold)

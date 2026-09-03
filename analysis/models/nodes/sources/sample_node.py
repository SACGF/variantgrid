import operator
from functools import cached_property, reduce
from typing import Optional

from auditlog.registry import auditlog
from cache_memoize import cache_memoize
from django.conf import settings
from django.db import models
from django.db.models import CASCADE, SET_NULL
from django.db.models.query_utils import Q

from analysis.models import GeneCoverageMixin
from analysis.models.nodes.analysis_node import (
    AnalysisNode,
    NodeAlleleFrequencyFilter,
    NodeAuditLogMixin,
    NodeVCFFilter,
    annotate_and_filter_queryset,
    queryset_to_pk_in_q,
)
from analysis.models.nodes.cohort_mixin import SampleMixin
from analysis.models.nodes.node_display import NodeChip, NodeIcon, NodeMenuEntry
from analysis.models.nodes.stats_cache import (
    get_cached_label_count_for_cohort,
    get_handler_for_node,
)
from genes.models import SampleGeneList
from library.constants import DAY_SECS
from library.utils import remove_duplicates_from_list
from patients.models import Extraction, Patient, Specimen
from patients.models_enums import NucleicAcid, SampleSourceLevel, Sex, Zygosity
from patients.sample_grouping import (
    SampleGroup,
    get_extraction_label,
    get_patient_for_source,
    get_sample_group,
    get_specimen_label,
)
from snpdb.models import Sample


class SampleNode(SampleMixin, GeneCoverageMixin, AnalysisNode):
    """ One entry point for a single VCF sample, or every sample belonging to a grouping object
        (an Extraction - the DNA arm's small variant, CNV and exon CNV calls in one node).

        Group levels resolve their samples at query time and OR a per-sample subquery, so nothing
        is pre-built to go stale as samples are linked to an extraction.

        Use restrict_to_qc_gene_list to keep track of that as sample_gene_list is cleared when a sample
        changes including in an AnalysisTemplates """
    source_level = models.CharField(max_length=1, choices=SampleSourceLevel.choices,
                                    default=SampleSourceLevel.SAMPLE)
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
    THRESHOLD_FIELDS = ("min_ad", "min_dp", "min_gq", "max_pl")
    SAMPLE_FIELD_MAPPINGS = [("ad", "allele_depth"),
                             ("dp", "read_depth"),
                             ("gq", "genotype_quality")]
    # The node's zygosity checkboxes, and the code each stands for - a sample's override stores the
    # whole set as one string of these
    ZYGOSITY_FIELD_CODES = [("zygosity_ref", Zygosity.HOM_REF),
                            ("zygosity_het", Zygosity.HET),
                            ("zygosity_hom", Zygosity.HOM_ALT),
                            ("zygosity_unk", Zygosity.UNKNOWN_ZYGOSITY)]
    # Levels that resolve to a set of samples. Specimen/Patient join here with their own resolvers
    SOURCE_LEVEL_FIELDS = {
        SampleSourceLevel.SAMPLE: "sample",
        SampleSourceLevel.EXTRACTION: "extraction",
        SampleSourceLevel.SPECIMEN: "specimen",
        SampleSourceLevel.PATIENT: "patient",
    }
    # Chip icons follow the models' own preview_icon()s, so chips, hover cards and search agree
    SOURCE_LEVEL_ICONS = {
        SampleSourceLevel.SAMPLE: NodeIcon(symbol="node-icon-sample"),
        SampleSourceLevel.EXTRACTION: NodeIcon(fa=Extraction.preview_icon()),
        SampleSourceLevel.SPECIMEN: NodeIcon(fa=Specimen.preview_icon()),
        SampleSourceLevel.PATIENT: NodeIcon(fa=Patient.preview_icon()),
    }
    VCF_CHIP_ICON = "fa-solid fa-file-lines"
    # Above this the class strip's ~45px runs out at 7px, so the card asks for the smaller rule
    MAX_STRIP_LABEL_CHARS = 8
    # Above this many siblings in one chip group the card runs out of room, so they collapse to a count
    MAX_GROUP_CHIPS = 3
    min_inputs = 0
    max_inputs = 0

    @classmethod
    def implemented_source_levels(cls) -> list[str]:
        """ A deployment with no patient data can trim the levels it offers """
        return [level for level in cls.SOURCE_LEVEL_FIELDS if level in settings.ANALYSIS_SAMPLE_NODE_LEVELS]

    # ── What samples this node reads ──────────────────────────────────────────

    @property
    def is_group_level(self) -> bool:
        return self.source_level != SampleSourceLevel.SAMPLE

    @property
    def source_field(self) -> str:
        return self.SOURCE_LEVEL_FIELDS[self.source_level]

    def get_source_object(self):
        """ The object this node groups on (a Sample at sample level) """
        return getattr(self, self.source_field)

    @cached_property
    def _sample_group(self) -> SampleGroup:
        if self.source_level == SampleSourceLevel.SAMPLE:
            # A hand picked sample is read as chosen - an archived or other build one is a
            # configuration error rather than a group silently narrowing
            return SampleGroup(samples=[self.sample] if self.sample else [])
        return get_sample_group(self.analysis.user, self.source_level, self.get_source_object(),
                                self.analysis.genome_build)

    def get_sample_group(self) -> SampleGroup:
        """ The one place samples are resolved - query, warnings, badge and chips all come from here.

            Cached per instance: the editor alone asks four times (errors, warnings, context, form
            kwargs) and the query build asks once per sample """
        return self._sample_group

    def get_source_samples(self) -> list[Sample]:
        return self.get_sample_group().samples

    def get_patient(self) -> Optional[Patient]:
        """ Every level resolves to one - it's what the pedigree badge has always drawn """
        return get_patient_for_source(self.source_level, self.get_source_object())

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

    # ── Filters, per sample ───────────────────────────────────────────────────

    @cached_property
    def _sample_filters(self) -> dict[int, "SampleNodeSampleFilter"]:
        """ Overrides by sample pk - the query build asks once per sample, and so does the editor.
            A node being built in an analysis template has no pk to hang rows off yet """
        if not self.pk:
            return {}
        return {sf.sample_id: sf for sf in self.samplenodesamplefilter_set.all()}

    def get_sample_filter(self, sample: Sample) -> Optional["SampleNodeSampleFilter"]:
        return self._sample_filters.get(sample.pk)

    def get_sample_thresholds(self, sample: Sample) -> dict:
        """ The node's own values, overridden where the sample's row sets one. A sample attached to
            the extraction after the node was configured gets the node defaults """
        thresholds = {f: getattr(self, f) for f in self.THRESHOLD_FIELDS}
        if override := self.get_sample_filter(sample):
            thresholds.update({f: value for f in thresholds
                               if (value := getattr(override, f)) is not None})
        return thresholds

    @cached_property
    def _node_pass_only(self) -> bool:
        return NodeVCFFilter.has_pass(self)

    def get_sample_pass_only(self, sample: Sample) -> bool:
        """ PASS is the one FILTER value that means the same thing in every VCF, so it is the one a
            sample can override - the rest are picked against the VCF that declared them """
        if (override := self.get_sample_filter(sample)) and override.pass_only is not None:
            return override.pass_only
        return self._node_pass_only

    def get_node_zygosity(self) -> str:
        """ The node's own selection, in the code string form a sample's override holds """
        return ''.join(code for field, code in self.ZYGOSITY_FIELD_CODES if getattr(self, field))

    # ── Query ─────────────────────────────────────────────────────────────────

    def _get_zygosities(self, sample: Optional[Sample] = None):
        # No genotype - all will be unknown
        if sample and sample.has_genotype is False:
            return [Zygosity.UNKNOWN_ZYGOSITY]

        if sample and (override := self.get_sample_filter(sample)) and override.zygosity is not None:
            return list(override.zygosity)
        return list(self.get_node_zygosity())

    def get_zygosity_description(self, sample: Optional[Sample] = None):
        zygosities = self._get_zygosities(sample or self.sample)
        if len(zygosities) != len(Zygosity.CHOICES):  # Subset
            zyg_dict = dict(Zygosity.CHOICES)
            zyg = ','.join([zyg_dict[z] for z in zygosities])
        else:
            zyg = ''  # Any
        return zyg

    def _get_sample_allele_frequency_arg_q_dict(self, sample: Sample) -> dict[Optional[str], dict[str, Q]]:
        """ The node's ranges, unless the sample's row sets its own - a caller that reports allele
            frequency differently is what that override is for """
        if (override := self.get_sample_filter(sample)) and override.has_allele_frequency:
            alias, allele_frequency_path = sample.get_cohort_genotype_alias_and_field("allele_frequency")
            if q := override.get_allele_frequency_q(allele_frequency_path, sample.vcf.allele_frequency_percent):
                return {alias: {str(q): q}}
            return {}
        return NodeAlleleFrequencyFilter.get_sample_arg_q_dict(self, sample)

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

            if sample_arg_q_dict := self._get_sample_allele_frequency_arg_q_dict(sample):
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
        return self._get_vcf_locus_filters_arg_q_dict(sample.vcf, alias,
                                                      pass_only=self.get_sample_pass_only(sample))

    def get_vcf_locus_filters_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        # Group levels apply these per sample inside each subquery
        if self.is_group_level:
            return {}
        return super().get_vcf_locus_filters_arg_q_dict()

    def get_vcf_locus_filter_vcfs(self) -> list:
        if self.is_group_level:
            return self.get_sample_group().vcfs
        return super().get_vcf_locus_filter_vcfs()

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
        if self.has_filters and self.get_sample_pass_only(sample):
            sample_description += " PASS"
        if (override := self.get_sample_filter(sample)) and override.has_allele_frequency:
            sample_description += f" AF {override.get_allele_frequency_description()}"

        if sample.has_genotype:
            if zyg := self.get_zygosity_description(sample):
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
        return "Variants from a VCF sample, or every sample of a patient, specimen or extraction " \
               "(so one node reads every caller's VCF)"

    def get_gene_lists(self) -> list:
        """ Used for gene coverage """
        gene_lists = []
        if self.sample_gene_list and self.restrict_to_qc_gene_list:
            gene_lists.append(self.sample_gene_list.gene_list)
        return gene_lists

    def _has_filters_that_affect_label_counts(self):
        # A sample filtered differently to the node is off the cached stats no matter what it sets
        if self._sample_filters:
            return True

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

    @classmethod
    def get_menu_entries(cls) -> list[NodeMenuEntry]:
        """ One row per level, so a user looking for "Patient" finds it where they look for it.
            The menu shows each level's own model icon; the card still draws the pedigree badge. """
        return [NodeMenuEntry(key=f"{cls.__name__}:{level}",
                              label=SampleSourceLevel(level).label,
                              icon=cls.SOURCE_LEVEL_ICONS[level],
                              initial_kwargs={"source_level": level})
                for level in cls.implemented_source_levels()]

    def get_node_icon(self) -> NodeIcon:
        """ Pedigree notation on the badge - square/circle for sex, struck through if deceased.

            The badge has always drawn the patient, and every level resolves to one, so what the
            node actually is goes on the class strip instead """
        if self.get_source_object() is None:
            return self.get_node_class_icon()
        if not (badge := self._get_card_data()["patient_badge"]):
            return self.get_node_class_icon()
        return NodeIcon(symbol=badge)

    def get_node_strip_label(self) -> str:
        return SampleSourceLevel(self.source_level).label

    def get_css_classes(self):
        css_classes = super().get_css_classes()
        if len(self.get_node_strip_label()) > self.MAX_STRIP_LABEL_CHARS:
            css_classes.append("node-klass-long")
        return css_classes

    def _get_patient_badge(self) -> Optional[str]:
        """ The sprite symbol for the patient's pedigree shape, or None where there is no patient """
        patient = self.get_patient()
        if patient is None:
            return None
        sex = "female" if patient.sex == Sex.FEMALE else "male"
        deceased = "-deceased" if patient.deceased else ""
        return f"node-icon-sample-{sex}{deceased}"

    @cache_memoize(DAY_SECS, args_rewrite=lambda s: (s.pk, s.version))
    def _get_card_data(self) -> dict:
        """ Everything the canvas card draws that needs the database - the pedigree badge, and the
            specimen -> extraction -> VCF names tree of what the node actually reads.

            Memoized together because get_rendering_dict runs for every node on every analysis
            render, and walking to the patient is up to three FK hops """
        group = self.get_sample_group()
        specimens = {}
        unlinked_vcfs = []
        for sample in group.samples:
            vcf_name = str(sample.vcf)
            if extraction := sample.extraction:
                specimen = extraction.specimen
                specimen_data = specimens.setdefault(specimen.pk, {
                    "label": get_specimen_label(specimen), "title": str(specimen), "extractions": {}})
                extraction_data = specimen_data["extractions"].setdefault(extraction.pk, {
                    "label": get_extraction_label(extraction), "title": str(extraction), "vcfs": []})
                if vcf_name not in extraction_data["vcfs"]:  # A caller can have several samples
                    extraction_data["vcfs"].append(vcf_name)
            else:
                unlinked_vcfs.append(vcf_name)

        return {
            "patient_badge": self._get_patient_badge(),
            "sample_vcf": str(self.sample.vcf) if self.sample else None,
            "specimens": [dict(sp, extractions=list(sp["extractions"].values())) for sp in specimens.values()],
            "unlinked_vcfs": remove_duplicates_from_list(unlinked_vcfs),
            "warnings": group.warnings,
        }

    @classmethod
    def _vcf_chip(cls, vcf_names: list[str]) -> NodeChip:
        """ Always one chip - a card has no room for a pill per caller, and the names are the hover """
        return NodeChip(text="VCF", icon=cls.VCF_CHIP_ICON, title="\n".join(vcf_names),
                        count=len(vcf_names) if len(vcf_names) > 1 else None)

    @classmethod
    def _group_chips(cls, items: list[dict], level: str, chip_func) -> tuple[NodeChip, ...]:
        """ Past MAX_GROUP_CHIPS siblings the card can't show them individually, so they become a count """
        if len(items) > cls.MAX_GROUP_CHIPS:
            return (NodeChip(text=SampleSourceLevel(level).label, icon=cls.SOURCE_LEVEL_ICONS[level].fa,
                             count=len(items), title="\n".join(i["title"] for i in items)),)
        return tuple(chip_func(item) for item in items)

    @classmethod
    def _extraction_chip(cls, extraction_data: dict) -> NodeChip:
        return NodeChip(text=extraction_data["label"], icon=cls.SOURCE_LEVEL_ICONS[SampleSourceLevel.EXTRACTION].fa,
                        title=extraction_data["title"], children=(cls._vcf_chip(extraction_data["vcfs"]),))

    @classmethod
    def _specimen_chip(cls, specimen_data: dict) -> NodeChip:
        children = cls._group_chips(specimen_data["extractions"], SampleSourceLevel.EXTRACTION, cls._extraction_chip)
        return NodeChip(text=specimen_data["label"], icon=cls.SOURCE_LEVEL_ICONS[SampleSourceLevel.SPECIMEN].fa,
                        title=specimen_data["title"], children=children)

    def get_node_chips(self) -> list[NodeChip]:
        """ The hierarchy from the picked object down, nested the way the relations are. Why a group
            node returns fewer rows than it used to - a VCF was archived - is a question only the
            canvas can answer, so what was left out gets an amber chip """
        chips = super().get_node_chips()
        if self.get_source_object() is None:
            return chips

        data = self._get_card_data()
        if self.source_level == SampleSourceLevel.SAMPLE:
            # A hand picked sample leaves nothing out, so there's no tree to walk
            chips.append(NodeChip(text="VCF", icon=self.VCF_CHIP_ICON, title=data["sample_vcf"]))
            return chips

        specimens = data["specimens"]
        if self.source_level == SampleSourceLevel.EXTRACTION:
            extractions = specimens[0]["extractions"] if specimens else []
            chips.extend(self._extraction_chip(e) for e in extractions)
        elif self.source_level == SampleSourceLevel.SPECIMEN:
            chips.extend(self._specimen_chip(sp) for sp in specimens)
        elif self.source_level == SampleSourceLevel.PATIENT:
            chips.extend(self._group_chips(specimens, SampleSourceLevel.SPECIMEN, self._specimen_chip))
            if unlinked_vcfs := data["unlinked_vcfs"]:
                chips.append(self._vcf_chip(unlinked_vcfs))

        if warnings := data["warnings"]:
            chips.append(NodeChip(text="", icon="fa-solid fa-triangle-exclamation", count=len(warnings),
                                  css_class="node-chip-warning", title="\n".join(warnings)))
        return chips

    def _get_configuration_check_cohorts(self) -> list:
        if self.is_group_level:
            cohorts, _ = self._get_cohorts_and_sample_visibility_for_node()
            return cohorts
        return super()._get_configuration_check_cohorts()

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if self.get_source_object() is None:
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
        sample_filters = list(self.samplenodesamplefilter_set.all())

        copy = super().save_clone()
        for sample_filter in sample_filters:
            sample_filter.pk = None
            sample_filter.node = copy
            sample_filter.save()
        return copy



class SampleNodeSampleFilter(NodeAuditLogMixin, models.Model):
    """ Per sample overrides of the node's genotype filters, so an extraction can apply a different
        min_ad per caller (sapath#301). Null is inherit: the node's own value is what any sample
        without a row gets, including one a reconcile attaches after the node was configured.

        Zygosity is the whole set in one column rather than a boolean each - it is one choice, and a
        row overriding HET while silently inheriting HOM is not a thing worth being able to store """
    node = models.ForeignKey(SampleNode, on_delete=CASCADE)
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    min_ad = models.IntegerField(null=True, blank=True)
    min_dp = models.IntegerField(null=True, blank=True)
    min_gq = models.IntegerField(null=True, blank=True)
    max_pl = models.IntegerField(null=True, blank=True)
    zygosity = models.CharField(max_length=4, null=True, blank=True)  # Zygosity codes to keep
    pass_only = models.BooleanField(null=True)
    # One range, where the node holds a list of them - a caller that reports allele frequency
    # differently is what this is for, not a different group operation
    af_min = models.FloatField(null=True, blank=True)
    af_max = models.FloatField(null=True, blank=True)

    class Meta:
        unique_together = ("node", "sample")

    def _get_node(self):
        return self.node

    @property
    def has_allele_frequency(self) -> bool:
        return self.af_min is not None or self.af_max is not None

    def get_allele_frequency_q(self, allele_frequency_path: str, allele_frequency_percent: bool) -> Optional[Q]:
        """ A full 0-1 range is no filter, the same as the node's own ranges """
        af_min = self.af_min if self.af_min is not None else 0
        af_max = self.af_max if self.af_max is not None else 1
        and_filters = []
        if af_min > 0:
            min_value = af_min * 100.0 if allele_frequency_percent else af_min
            and_filters.append(Q(**{f"{allele_frequency_path}__gte": min_value}))
        if af_max < 1:
            max_value = af_max * 100.0 if allele_frequency_percent else af_max
            and_filters.append(Q(**{f"{allele_frequency_path}__lte": max_value}))
        if and_filters:
            return reduce(operator.and_, and_filters)
        return None

    def get_allele_frequency_description(self) -> str:
        af_min = self.af_min if self.af_min is not None else 0
        af_max = self.af_max if self.af_max is not None else 1
        return f"{af_min}-{af_max}"

    def __str__(self):
        thresholds = {"AD>=": self.min_ad, "DP>=": self.min_dp, "GQ>=": self.min_gq, "PL<=": self.max_pl}
        parts = [f"{label}{value}" for label, value in thresholds.items() if value is not None]
        if self.zygosity is not None:
            parts.append(f"zygosity={self.zygosity}")
        if self.pass_only is not None:
            parts.append(f"pass_only={self.pass_only}")
        if self.has_allele_frequency:
            parts.append(f"AF {self.get_allele_frequency_description()}")
        return f"{self.sample}: {' '.join(parts) or 'no overrides'}"


auditlog.register(SampleNode)
auditlog.register(SampleNodeSampleFilter)

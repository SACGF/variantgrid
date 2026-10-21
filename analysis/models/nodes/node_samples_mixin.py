""" NodeSamplesMixin - the samples, cohorts and proband an AnalysisNode reads, its own plus every
    ancestor's. Entry points: get_samples / get_samples_with_genotype / get_cohorts_and_sample_visibility,
    get_proband (NodeProband - sample and patient resolved independently) and get_alignments_dict. """
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cached_property
from typing import Optional

from cache_memoize import cache_memoize

from library.constants import DAY_SECS
from patients.models import Patient
from patients.models_enums import SampleSourceLevel
from patients.sample_grouping import get_patient_for_source
from snpdb.enums import SampleFileType
from snpdb.models import Cohort, Sample, SampleFilePath


@dataclass(frozen=True)
class NodeProband:
    """ Who a node is about - the object of the study, resolved once per node up the DAG.
        Either half can be known without the other: an extraction level node whose DNA arm has two callers
        has no single sample but does have a patient (@see AnalysisNode.get_proband) """
    sample: Optional[Sample]
    patient: Optional[Patient]


class NodeSamplesMixin:
    """ Part of AnalysisNode. Hooks a node subclass overrides:
        - _get_cohorts_and_sample_visibility_for_node: the cohorts the node itself reads, and which of their
          samples are shown on the grid - ancestors are added by get_cohorts_and_sample_visibility
        - _get_proband_sample_for_node / _get_proband_patient_for_node: who the study is about, if this node
          knows - get_proband combines them with every ancestor's answer """

    @cached_property
    def num_samples_for_build(self) -> int:
        return Sample.objects.filter(vcf__genome_build=self.analysis.genome_build).count()

    def _get_cohorts_and_sample_visibility_for_node(self) -> tuple[Sequence[Cohort], dict]:
        """ Visibility = can see on grid """
        return [], {}

    @staticmethod
    def _get_visible_samples_from_cohort(cohorts, visibility):
        samples = set()
        for c in cohorts:
            for s in c.get_samples():
                if visibility.get(s):
                    samples.add(s)
        return sorted(samples)

    def get_cohorts_and_sample_visibility(self, sort=True) -> tuple[Sequence[Cohort], dict]:
        """ Returns all node + ancestor cohorts (and visibilities of their samples)

            The underlying data for all samples/cohorts/sub-cohorts/trios/pedigrees is Cohorts, so need to know which
            to retrieve from DB (and what sample info to extract from packed columns) to filter + show on grid """

        cohorts, visibility = self._get_cohorts_and_sample_visibility_for_node()
        cohorts = set(cohorts)
        if self.has_input():
            parents, _ = self.get_parent_subclasses_and_errors()
            for parent in parents:
                c, v = parent.get_cohorts_and_sample_visibility(sort=False)
                cohorts.update(c)
                visibility.update(v)

        # May have sub-cohorts, so get unique base cohorts
        cohorts = {c.get_base_cohort() for c in cohorts}
        if sort:
            cohorts = sorted(cohorts)
        return cohorts, visibility

    @cache_memoize(DAY_SECS, args_rewrite=lambda s: (s.pk, s.version))
    def get_sample_ids_with_genotype(self) -> list[int]:
        return [s.pk for s in self.get_samples_with_genotype()]

    def get_samples_from_node_only_not_ancestors(self):
        _, visibility = self._get_cohorts_and_sample_visibility_for_node()
        return sorted(visibility)

    def _get_proband_sample_for_node(self) -> Optional[Sample]:
        """ Sample of the object of a study, if known """
        return None

    def _get_proband_patient_for_node(self) -> Optional[Patient]:
        """ Patient of the object of a study, if known - the proband sample's patient, either way the
            sample is linked (@see patients/sample_grouping.py:get_patient_for_source) """
        return get_patient_for_source(SampleSourceLevel.SAMPLE, self._get_proband_sample_for_node())

    def get_proband(self, proband_by_node_id: Optional[dict[int, NodeProband]] = None) -> NodeProband:
        """ Who the node is about - sample and patient resolved independently by the same rule, so a node
            that cannot name a sample (two callers on the one extraction) can still name the patient.
            proband_by_node_id: answers already worked out, shared across a whole graph rather than kept on the
            node - a descendant asks each of its ancestors, so a diamond in the DAG asks the same node once per
            path to it (@see analysis/variant_tag_operations.py:get_proband_by_node_id) """
        if proband_by_node_id is not None and self.pk in proband_by_node_id:
            return proband_by_node_id[self.pk]

        proband_samples = set()
        proband_patients = set()
        if proband_sample := self._get_proband_sample_for_node():
            proband_samples.add(proband_sample)
        if proband_patient := self._get_proband_patient_for_node():
            proband_patients.add(proband_patient)

        if self.has_input():
            parents, _ = self.get_parent_subclasses_and_errors()
            for parent in parents:
                parent_proband = parent.get_proband(proband_by_node_id)
                if parent_proband.sample:
                    proband_samples.add(parent_proband.sample)
                if parent_proband.patient:
                    proband_patients.add(parent_proband.patient)

        # If ambiguous, then just give up
        proband = NodeProband(sample=proband_samples.pop() if len(proband_samples) == 1 else None,
                              patient=proband_patients.pop() if len(proband_patients) == 1 else None)
        if proband_by_node_id is not None:
            proband_by_node_id[self.pk] = proband
        return proband

    def get_proband_sample(self, proband_by_node_id: Optional[dict[int, NodeProband]] = None) -> Optional[Sample]:
        """ Sample of the object of a study if known """
        return self.get_proband(proband_by_node_id).sample

    def get_proband_patient(self, proband_by_node_id: Optional[dict[int, NodeProband]] = None) -> Optional[Patient]:
        """ Patient of the object of a study if known - a node above sample level has one where it has no
            single sample """
        return self.get_proband(proband_by_node_id).patient

    def get_samples_with_genotype(self) -> list[Sample]:
        """ Node + ancestor samples whose genotype we can show/filter on - ie variant-only VCFs
            (has_sample_columns=False) are left out. Use get_samples() for sample level data """
        cohorts, visibility = self.get_cohorts_and_sample_visibility(sort=False)
        return self._get_visible_samples_from_cohort(cohorts, visibility)

    def get_samples(self) -> list[Sample]:
        """ Every node + ancestor sample, including those from variant-only VCFs. Sample level data
            that doesn't come from the genotype - gene lists, coverage, alignment files, patients - and
            restricting sample fields to ancestors @see AncestorSampleMixin """
        _, visibility = self.get_cohorts_and_sample_visibility(sort=False)
        return sorted(visibility)  # Every sample the node knows about is a key, genotype or not

    def get_sample_ids(self) -> list[int]:
        return [s.pk for s in self.get_samples()]

    def get_alignments_dict(self):
        alignments_dict = defaultdict(set)
        sfp_qs = SampleFilePath.objects.filter(sample__in=self.get_samples(),
                                               file_type__in=SampleFileType.ALIGNMENT_TYPES)
        for sample_id, file_path in sfp_qs.values_list("sample_id", "file_path"):
            alignments_dict[sample_id].add(file_path)
        return {k: list(v) for k, v in alignments_dict.items()}

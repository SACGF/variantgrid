"""
Auto-launching analysis templates when VCFs import (#2898).

An AutoLaunchAnalysisTemplate rule matches a sample on its enrichment kit, name, the variant caller seqauto
registered its VCF under and the VCF's ##source - a blank criterion matches anything. Its launch_level says
what the run is for:

  * Sample - one run per matching sample, on import (and again once its QC gene list arrives).
  * Extraction / Specimen / Patient - one run per (template, group, sequencing run), the template's SampleNode
    taking the group as its variable so it reads every VCF the group has. A TSO 500 pair's DNA and RNA arms
    are two extractions of one specimen, so a Specimen rule gives one analysis per case.

A group launches once every file the pipeline registered for it on that run (the SingleSampleVCFs of its
SequencingSamples) has imported. Any of those in error stops the launch, with an event saying why - an
incomplete case analysis is worse than none. Group membership can arrive after the VCFs, so a sample linked
to its extraction / patient after import retries (analysis.signals.signal_handlers.handle_sample_saved). Either
way the run's trigger is ALL_LINKED_SEQUENCING_SAMPLES_IMPORTED. A VCF with no sequencing records has no
expected set: its group launches (trigger VCF_IMPORT) on the first matching VCF, and later imports reload it.

Entry points: auto_launch_analysis_templates_for_sample (sample level), auto_launch_group_analyses_for_sample
(group level), get_auto_launch_analysis_template_matches (the matching, for the auto launch page).
"""
import logging
import re
from dataclasses import dataclass
from typing import Optional

from django.contrib.auth.models import User
from django.db.models import Q, QuerySet

from analysis.analysis_templates import populate_analysis_from_template_run
from analysis.enums import AnalysisTemplateRunTrigger
from analysis.models.models_analysis import (
    AnalysisTemplate,
    AnalysisTemplateRun,
    AutoLaunchAnalysisTemplate,
    GroupAnalysisTemplateRun,
)
from analysis.models.nodes.node_utils import reload_analysis_nodes
from analysis.related_analyses import get_related_analysis_details_for_samples
from eventlog.models import create_event
from genes.models import ActiveSampleGeneList
from library.django_utils.database_utils import advisory_lock
from library.log_utils import LogLevel
from patients.enums import SampleSourceLevel
from patients.sample_grouping import get_source_level_q
from seqauto.models import (
    DragenTSO500CombinedVariantOutput,
    SequencingRun,
    SequencingSample,
    SingleSampleVCF,
    VCFFromSequencingRun,
)
from snpdb.enums import ImportStatus, ProcessingStatus
from snpdb.models import Sample

GROUP_LEVEL_FIELDS = {
    SampleSourceLevel.EXTRACTION: "extraction",
    SampleSourceLevel.SPECIMEN: "specimen",
    SampleSourceLevel.PATIENT: "patient",
}
# SequencingSample -> the group, through the extraction it was made from
SEQUENCING_SAMPLE_GROUP_PATHS = {
    SampleSourceLevel.EXTRACTION: "extraction",
    SampleSourceLevel.SPECIMEN: "extraction__specimen",
    SampleSourceLevel.PATIENT: "extraction__specimen__patient",
}
# DragenTSO500CombinedVariantOutput -> the group. A pair's CVO names both arms before either arm's samples
# are linked, so it is what keeps the RNA arm expected while only the DNA arm has arrived
CVO_GROUP_PATHS = {
    SampleSourceLevel.SPECIMEN: "specimen",
    SampleSourceLevel.PATIENT: "specimen__patient",
}
UPLOAD_FAILED_STATUSES = {ProcessingStatus.ERROR, ProcessingStatus.TERMINATED_EARLY, ProcessingStatus.TIMED_OUT}


@dataclass
class AutoLaunchSubject:
    """ What a rule is matched against - a sample, or the example on the auto launch page """
    enrichment_kit_name: Optional[str] = None
    sample_name: str = ""
    variant_caller_name: Optional[str] = None
    vcf_source: Optional[str] = None

    @staticmethod
    def from_sample(sample: Sample) -> 'AutoLaunchSubject':
        variant_caller_name = None
        if vfsr := VCFFromSequencingRun.objects.filter(vcf=sample.vcf).select_related("variant_caller").first():
            if vfsr.variant_caller:
                variant_caller_name = vfsr.variant_caller.name
        return AutoLaunchSubject(enrichment_kit_name=sample.enrichment_kit.name if sample.enrichment_kit else None,
                                 sample_name=sample.name,
                                 variant_caller_name=variant_caller_name,
                                 vcf_source=sample.vcf.source)


def _regex_match(regex: Optional[str], value: Optional[str]) -> bool:
    """ Blank means anything """
    if not regex:
        return True
    return value is not None and bool(re.match(regex, value))


@dataclass
class AutoLaunchAnalysisTemplateMatch:
    auto_launch: AutoLaunchAnalysisTemplate
    enrichment_kit_name: Optional[str]
    enrichment_kit_match: bool
    sample_regex_match: bool
    variant_caller_match: bool
    vcf_source_match: bool

    @property
    def analysis_template(self) -> AnalysisTemplate:
        return self.auto_launch.template

    @property
    def match(self) -> bool:
        return all([self.enrichment_kit_match, self.sample_regex_match, self.variant_caller_match,
                    self.vcf_source_match])


def get_auto_launch_analysis_template_matches(user, subject: AutoLaunchSubject) -> list[AutoLaunchAnalysisTemplateMatch]:
    matches = []
    templates_qs = AnalysisTemplate.filter_for_user(user)
    auto_launch_qs = AutoLaunchAnalysisTemplate.objects.filter(template__in=templates_qs)
    for auto_launch in auto_launch_qs.select_related("template", "enrichment_kit"):
        enrichment_kit_name = None
        if enrichment_kit := auto_launch.enrichment_kit:
            enrichment_kit_name = enrichment_kit.name

        matches.append(AutoLaunchAnalysisTemplateMatch(
            auto_launch=auto_launch,
            enrichment_kit_name=enrichment_kit_name,
            enrichment_kit_match=subject.enrichment_kit_name == enrichment_kit_name,
            sample_regex_match=_regex_match(auto_launch.sample_regex, subject.sample_name),
            variant_caller_match=_regex_match(auto_launch.variant_caller_regex, subject.variant_caller_name),
            vcf_source_match=_regex_match(auto_launch.vcf_source_regex, subject.vcf_source)))
    return matches


def _get_matching_auto_launches(user, sample: Sample, group_level: bool) -> list[AutoLaunchAnalysisTemplate]:
    matches = get_auto_launch_analysis_template_matches(user, AutoLaunchSubject.from_sample(sample))
    return [m.auto_launch for m in matches
            if m.match and (m.auto_launch.launch_level != SampleSourceLevel.SAMPLE) == group_level]


def auto_launch_analysis_templates_for_sample(user, sample, trigger: AnalysisTemplateRunTrigger,
                                              skip_already_analysed=False):
    if skip_already_analysed:
        if get_related_analysis_details_for_samples(user, [sample]):
            return

    for auto_launch in _get_matching_auto_launches(user, sample, group_level=False):
        analysis_template = auto_launch.template
        template_version = analysis_template.active
        if template_version is None:
            logging.warning("Skipping auto analysis '%s' for sample: %s. Template has no active version",
                            analysis_template, sample)
            continue

        template_arguments = {"sample": sample}
        if template_version.requires_sample_gene_list:
            try:
                template_arguments["sample_gene_list"] = sample.activesamplegenelist.sample_gene_list
            except ActiveSampleGeneList.DoesNotExist:
                logging.warning("Skipping auto analysis '%s' for sample: %s. Will try again if QC Gene Lists created", analysis_template, sample)
                continue
        template_run = AnalysisTemplateRun.create(analysis_template, sample.genome_build, trigger, user=user,
                                                  template_version=template_version)
        template_run.populate_arguments(template_arguments)
        populate_analysis_from_template_run(template_run)


def get_sample_group_object(sample: Sample, level: SampleSourceLevel):
    """ The extraction / specimen / patient a sample belongs to, None while it isn't linked to one """
    if level == SampleSourceLevel.PATIENT:
        return sample.patient  # Sample.save fills it from the extraction
    if extraction := sample.extraction:
        if level == SampleSourceLevel.EXTRACTION:
            return extraction
        return extraction.specimen
    return None


def get_group_expected_vcfs(level: SampleSourceLevel, group, sequencing_run: SequencingRun) -> QuerySet[SingleSampleVCF]:
    """ The files registered for the group on this run: the VCFs of the run's current SequencingSamples linked
        to the group by their extraction, by the Sample made from them, or as an arm of a pair the group claims """
    group_samples = Sample.objects.filter(get_source_level_q(level, group))
    sequencing_sample_q = Q(**{SEQUENCING_SAMPLE_GROUP_PATHS[level]: group})
    sequencing_sample_q |= Q(samplefromsequencingsample__sample__in=group_samples)
    if cvo_path := CVO_GROUP_PATHS.get(level):
        cvo_qs = DragenTSO500CombinedVariantOutput.objects.filter(sequencing_run=sequencing_run, **{cvo_path: group})
        sequencing_sample_q |= Q(pk__in=cvo_qs.values("dna_sequencing_sample"))
        sequencing_sample_q |= Q(pk__in=cvo_qs.values("rna_sequencing_sample"))
    sequencing_samples = SequencingSample.get_current().filter(sequencing_sample_q,
                                                               sample_sheet__sequencing_run=sequencing_run)
    return SingleSampleVCF.objects.filter(sequencing_sample__in=sequencing_samples)


@dataclass
class GroupFiles:
    waiting: list[SingleSampleVCF]
    errors: list[SingleSampleVCF]

    @staticmethod
    def check(expected_vcfs) -> 'GroupFiles':
        group_files = GroupFiles(waiting=[], errors=[])
        for single_sample_vcf in expected_vcfs:
            vcf = single_sample_vcf.vcf
            if vcf and vcf.import_status == ImportStatus.SUCCESS:
                continue
            upload_pipeline = single_sample_vcf.upload_pipeline
            if (vcf and vcf.import_status == ImportStatus.ERROR) or \
                    (upload_pipeline and upload_pipeline.status in UPLOAD_FAILED_STATUSES):
                group_files.errors.append(single_sample_vcf)
            else:
                group_files.waiting.append(single_sample_vcf)
        return group_files


def auto_launch_group_analyses_for_sample(user: User, sample: Sample, reload_existing: bool = True):
    """ reload_existing - a VCF joining a group that has already launched reloads its analysis, so node counts
        include it. Off for a retry that only exists to catch a late link, which a sample save also fires """
    if sample.import_status != ImportStatus.SUCCESS or sample.vcf.data_archived:
        return

    sequencing_run = None
    if vfsr := VCFFromSequencingRun.objects.filter(vcf=sample.vcf).first():
        sequencing_run = vfsr.sequencing_run

    launches = {(al.template, al.launch_level) for al in _get_matching_auto_launches(user, sample, group_level=True)}
    for analysis_template, level in launches:
        if group := get_sample_group_object(sample, level):
            _launch_or_reload_group_analysis(user, sample, analysis_template, level, group, sequencing_run,
                                             reload_existing)


def _launch_or_reload_group_analysis(user, sample: Sample, analysis_template: AnalysisTemplate,
                                     level: SampleSourceLevel, group, sequencing_run: Optional[SequencingRun],
                                     reload_existing: bool):
    group_field = GROUP_LEVEL_FIELDS[level]
    lock_name = f"{GroupAnalysisTemplateRun._meta.db_table}_{analysis_template.pk}_{level}_{sequencing_run and sequencing_run.pk}"
    with advisory_lock(lock_name, group.pk):
        existing = GroupAnalysisTemplateRun.objects.filter(template=analysis_template, source_level=level,
                                                           sequencing_run=sequencing_run,
                                                           **{group_field: group}).first()
        if existing:
            if reload_existing:
                reload_analysis_nodes(existing.analysis_template_run.analysis_id)
            return

        description = f"'{analysis_template}' for {group_field} '{group}'"
        if sequencing_run:
            group_files = GroupFiles.check(get_group_expected_vcfs(level, group, sequencing_run))
            if group_files.errors:
                failed = ", ".join(str(ssvcf) for ssvcf in group_files.errors)
                create_event(user, "auto_launch_analysis_template",
                             f"Not launching {description} (run {sequencing_run}): failed imports: {failed}",
                             severity=LogLevel.WARNING)
                return
            if group_files.waiting:
                logging.info("Auto analysis %s waiting on %d file(s)", description, len(group_files.waiting))
                return

        template_version = analysis_template.active
        if template_version is None:
            logging.warning("Skipping auto analysis %s. Template has no active version", description)
            return
        if template_version.requires_sample_gene_list:
            # A group has a gene list per sample, so there is no one to give it
            create_event(user, "auto_launch_analysis_template",
                         f"Not launching {description}: the template needs a sample gene list, which a "
                         f"{group_field} level launch can't supply", severity=LogLevel.WARNING)
            return

        # Whichever import or late link completed the set, the launch is the set being complete
        if sequencing_run:
            trigger = AnalysisTemplateRunTrigger.ALL_LINKED_SEQUENCING_SAMPLES_IMPORTED
        else:
            trigger = AnalysisTemplateRunTrigger.VCF_IMPORT
        template_run = AnalysisTemplateRun.create(analysis_template, sample.genome_build, trigger, user=user,
                                                  template_version=template_version)
        template_run.populate_arguments({group_field: group})
        GroupAnalysisTemplateRun.objects.create(analysis_template_run=template_run, template=analysis_template,
                                                source_level=level, sequencing_run=sequencing_run,
                                                **{group_field: group})
        populate_analysis_from_template_run(template_run)

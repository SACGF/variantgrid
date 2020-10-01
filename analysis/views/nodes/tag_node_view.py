from django.conf import settings

from analysis.forms import InputSamplesForm
from analysis.forms.forms_nodes import TagNodeForm
from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.filters.tag_node import TagNode
from analysis.views.nodes.node_view import NodeView
from annotation.transcripts_annotation_selections import VariantTranscriptSelections
from genes.hgvs import HGVSMatcher
from snpdb.models.models import Tag
from snpdb.models.models_vcf import Sample
from classification.models.classification import Classification


class TagNodeView(NodeView):
    model = TagNode
    form_class = TagNodeForm

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        if self.object.visible:
            requires_classification_data = []
        else:
            requires_classification_data = _get_requires_classification_data(self.object.analysis, self.request.user)
        context["requires_classification_data"] = requires_classification_data
        return context


def _get_requires_classification_data_for_tag(hgvs_matcher, analysis, samples, num_samples,
                                              analysis_single_sample, variant_tag):
    single_sample = analysis_single_sample
    if num_samples >= 2:  # could be more than 1
        if variant_tag.node:
            node = AnalysisNode.objects.filter(pk=variant_tag.node.pk).get_subclass()
            samples = node.get_samples()
            if len(samples) == 1:
                single_sample = samples[0]

    variant = variant_tag.variant
    variant_summary = str(variant)

    single_transcript_label = None
    single_transcript_id = None
    variant_transcript_selections = VariantTranscriptSelections(variant,
                                                                analysis.genome_build,
                                                                analysis.annotation_version)
    if variant_transcript_selections.transcript_data:
        if len(variant_transcript_selections.transcript_data) == 1:
            single_transcript_id = variant_transcript_selections.initial_transcript_id
            try:
                single_transcript_label = hgvs_matcher.variant_to_c_hgvs(variant, single_transcript_id)
            except:
                single_transcript_label = single_transcript_id
    else:
        single_transcript_id = "(intergenic)"
        single_transcript_label = single_transcript_id

    if single_sample:
        sample_form = None
    else:
        sample_form = InputSamplesForm(samples=samples)

    return variant_tag, variant_summary, single_transcript_label, single_transcript_id, variant_transcript_selections, single_sample, sample_form


def _get_requires_classification_data(analysis, user):
    """ "Special tab in Analysis wide Tags view accessed from top button
        returns (variant_tag, variant_summary, single_transcript_label, single_transcript_id, variant_transcript_selections, single_sample, sample_form)
    """
    requires_classification_data = []
    if Classification.can_create_via_web_form(user):
        tag = Tag.objects.filter(pk=settings.TAG_REQUIRES_CLASSIFICATION).first()
        if tag:
            hgvs_matcher = HGVSMatcher(analysis.genome_build)

            samples = analysis.get_samples()
            num_samples = len(samples)
            if num_samples == 1:
                analysis_single_sample = samples[0]
            else:
                analysis_single_sample = None

            for variant_tag in analysis.varianttag_set.filter(tag=tag).order_by("created"):
                data = _get_requires_classification_data_for_tag(hgvs_matcher, analysis,
                                                                 samples, num_samples,
                                                                 analysis_single_sample,
                                                                 variant_tag)
                requires_classification_data.append(data)

    return requires_classification_data

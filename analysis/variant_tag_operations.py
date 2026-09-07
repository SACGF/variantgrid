"""
Operations on taggings (VariantTag) - which sample a tagging is about, and resolving the to-do tags a
classification has satisfied.

A tag with Tag.requires_classification is a to-do item, and classifying the variant is what completes it.
The tagging is marked resolved and linked to the classification rather than deleted, so it stays as the record
of what was flagged and what it turned into, and withdrawing the classification puts the to-do back
(@see VariantTag.is_resolved). Resolved taggings still show everywhere a tagging shows.

VariantTag isn't registered with auditlog - taggings come and go all the time and we only want this one
deliberate resolution - so the LogEntry is written by hand. Putting analysis_id in additional_data is what
makes it show up in the analysis audit log (@see Analysis.log_entry_qs).

Tagging stays one click - the sample is the study's proband where the node knows it, left null otherwise and
never prompted for. Carrying the variant is not what makes a tagging someone's: a relative who is HET for the
proband's variant doesn't need their own classification.
"""
from collections import defaultdict
from collections.abc import Iterable
from typing import Optional

from auditlog.models import LogEntry
from django.contrib.auth.models import User
from django.db import transaction
from django.utils.timezone import now

from analysis.models import Analysis, AnalysisEdge, VariantTag
from analysis.models.nodes.node_utils import get_nodes_by_id
from classification.models import Classification
from patients.models_enums import Zygosity
from snpdb.models import Sample, SampleGenotype, Variant

VARIANT_TAG_CLASSIFIED = "classified"

# A sample carries the variant if it was called in it - reference and no-call don't make it this case's variant
CARRIER_ZYGOSITIES = {Zygosity.HET, Zygosity.HOM_ALT}


def _variant_for_sample(variant_tag: VariantTag, sample: Sample) -> Optional[Variant]:
    """ The tag's variant as the sample's VCF knows it - tags and samples can be in different builds """
    genome_build = sample.vcf.genome_build
    if variant_tag.genome_build_id == genome_build.pk:
        return variant_tag.variant
    if allele := variant_tag.allele:
        return allele.variant_for_build_optional(genome_build)
    return None


def get_sample_genotype_for_variant_tag(sample: Sample, variant_tag: VariantTag) -> Optional[SampleGenotype]:
    if variant := _variant_for_sample(variant_tag, sample):
        return sample.get_genotype(variant)
    return None


def sample_carries_variant(sample: Sample, variant_tag: VariantTag) -> bool:
    """ Shown on the classify dialog as zygosity - it is not what decides whose tagging this is """
    if sample_genotype := get_sample_genotype_for_variant_tag(sample, variant_tag):
        return sample_genotype.zygosity in CARRIER_ZYGOSITIES
    return False


def get_sample_for_variant_tag(variant_tag: VariantTag) -> Optional[Sample]:
    """ Which sample the tagging is about - the proband of the study the tagged node sits in, which is the same
        answer AncestorSampleMixin nodes auto-populate from. None when the node's ancestors disagree """
    if node := variant_tag.node:
        return node.get_subclass().get_proband_sample()
    return None


def get_proband_sample_by_node_id(analysis: Analysis) -> dict[int, Optional[Sample]]:
    """ Every node's answer to get_sample_for_variant_tag, from the analysis graph loaded once.
        Asking a tagging at a time walks the ancestors a subclass query at a time and re-walks them for
        the next tagging - one analysis has thousands of taggings across a handful of nodes """
    nodes_by_id = get_nodes_by_id(analysis.analysisnode_set.all().select_subclasses())
    parents = defaultdict(list)
    for parent_id, child_id in AnalysisEdge.objects.filter(parent__analysis=analysis).values_list("parent", "child"):
        parents[child_id].append(nodes_by_id[parent_id])
    for node_id, node in nodes_by_id.items():
        node._cached_parents = parents.get(node_id, [])

    proband_by_node_id = {}
    for node in nodes_by_id.values():
        node.get_proband_sample(proband_by_node_id)
    return proband_by_node_id


def classification_resolves_tag(variant_tag: VariantTag, classification: Classification) -> bool:
    """ Whether the classification is of the person the tagging is about. Where either side doesn't say, only a
        one-person analysis is unambiguous - otherwise it is the scientist's call (@see resolve_variant_tag) """
    if variant_tag.sample_id and classification.sample_id:
        return variant_tag.sample_id == classification.sample_id
    if analysis := variant_tag.analysis:
        return len(analysis.get_samples()) <= 1
    return True


def _log_variant_tag_classified(variant_tag: VariantTag, classification: Classification, user: User) -> LogEntry:
    return LogEntry.objects.log_create(
        variant_tag,
        force_log=True,
        action=LogEntry.Action.UPDATE,
        actor=user,
        additional_data={
            "operation": VARIANT_TAG_CLASSIFIED,
            "analysis_id": variant_tag.analysis_id,
            "node_id": variant_tag.node_id,
            "tag_id": variant_tag.tag_id,
            "variant_id": variant_tag.variant_id,
            "allele_id": variant_tag.allele_id,
            "classification_id": classification.pk,
        },
    )


def resolve_variant_tag(variant_tag: VariantTag, classification: Classification, user: User) -> VariantTag:
    """ Mark the to-do satisfied by this classification. Called automatically where the classification is of the
        tagging's own sample, and from the queue's button where the scientist says so """
    with transaction.atomic():
        variant_tag.resolved = now()
        variant_tag.resolved_by = user
        variant_tag.resolved_classification = classification
        variant_tag.save(update_fields=["resolved", "resolved_by", "resolved_classification", "modified"])
        _log_variant_tag_classified(variant_tag, classification, user)
    return variant_tag


def _resolve_unambiguous(variant_tags: Iterable[VariantTag], classification: Classification,
                         user: User) -> list[VariantTag]:
    resolved = []
    for variant_tag in variant_tags:
        if not variant_tag.is_resolved and classification_resolves_tag(variant_tag, classification):
            resolved.append(resolve_variant_tag(variant_tag, classification, user))
    return resolved


def resolve_requires_classification_tags(classification: Classification, analysis: Analysis,
                                         user: User) -> list[VariantTag]:
    """ Resolve the taggings the classification just satisfied - the variant is done for that person, so that's
        everyone's tagging of it in this analysis, not just the one that was clicked """
    variant_tags = VariantTag.objects.filter(variant=classification.variant, analysis=analysis,
                                             tag__requires_classification=True)
    return _resolve_unambiguous(variant_tags, classification, user)


def resolve_requires_classification_tags_for_samples(classification: Classification, samples: Iterable[Sample],
                                                     user: User) -> list[VariantTag]:
    """ Same as resolve_requires_classification_tags, for the sample/patient page where the case is a set of
        samples rather than an analysis - a tagging with no sample of its own belongs to the case when it was
        made in an analysis one of the case's samples is in. Someone with read-only access to the analysis can
        still classify, their tagging just stays as it is """
    variant = classification.variant
    if variant is None:
        return []

    sample_ids = {s.pk for s in samples}
    variant_tags = []
    for variant_tag in VariantTag.objects.filter(variant__in=variant.equivalent_variants,
                                                 tag__requires_classification=True):
        if variant_tag.sample_id:
            in_case = variant_tag.sample_id in sample_ids
        else:
            in_case = bool(variant_tag.analysis) and \
                bool(sample_ids.intersection(s.pk for s in variant_tag.analysis.get_samples()))
        if in_case and variant_tag.can_write(user):
            variant_tags.append(variant_tag)
    return _resolve_unambiguous(variant_tags, classification, user)

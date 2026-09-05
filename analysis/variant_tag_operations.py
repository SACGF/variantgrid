"""
Operations on taggings (VariantTag) - which sample a tagging is about, and retiring the to-do tags a
classification has satisfied.

The RequiresClassification tag is a to-do item - classifying the variant is what completes it, so the tagging
is deleted once a classification exists. The row going away takes with it any sign the variant was ever
flagged, so before it goes we log what it turned into.

VariantTag isn't registered with auditlog - taggings come and go all the time and we only want this one
deliberate retirement - so the LogEntry is written by hand. Putting analysis_id in additional_data is what
makes it show up in the analysis audit log (@see Analysis.log_entry_qs).

Tagging stays one click - the sample is worked out silently where the answer is free (a single sample node,
or a single carrier in the analysis) and left null otherwise, to be resolved at classification time where a
sample dropdown already exists.
"""
from collections.abc import Iterable
from typing import Optional

from auditlog.models import LogEntry
from django.conf import settings
from django.contrib.auth.models import User
from django.db import transaction

from analysis.models import Analysis, VariantTag
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
    if sample_genotype := get_sample_genotype_for_variant_tag(sample, variant_tag):
        return sample_genotype.zygosity in CARRIER_ZYGOSITIES
    return False


def get_sample_for_variant_tag(variant_tag: VariantTag) -> Optional[Sample]:
    """ Which sample the tagging is about, where that can be answered without asking the user:
        (a) the node it was tagged in has exactly one sample, or
        (b) exactly one of the analysis's samples carries the variant """
    if node := variant_tag.node:
        if samples := node.get_subclass().get_samples():
            if len(samples) == 1:
                return samples[0]

    if analysis := variant_tag.analysis:
        carriers = [s for s in analysis.get_samples() if sample_carries_variant(s, variant_tag)]
        if len(carriers) == 1:
            return carriers[0]
    return None


def _log_variant_tag_classified(variant_tag: VariantTag, classification: Classification, user: User) -> LogEntry:
    return LogEntry.objects.log_create(
        variant_tag,
        force_log=True,
        action=LogEntry.Action.DELETE,
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


def _retire_variant_tags(variant_tags: list[VariantTag], classification: Classification, user: User) -> int:
    with transaction.atomic():
        for variant_tag in variant_tags:
            _log_variant_tag_classified(variant_tag, classification, user)
        VariantTag.objects.filter(pk__in=[vt.pk for vt in variant_tags]).delete()
    return len(variant_tags)


def retire_requires_classification_tags(classification: Classification, analysis: Analysis, user: User) -> int:
    """ Retire the taggings the classification just satisfied - the whole variant is done, so that's everyone's
        tagging of it in this analysis, not just the one that was clicked """
    variant_tags = list(VariantTag.objects.filter(variant=classification.variant, analysis=analysis,
                                                 tag_id=settings.TAG_REQUIRES_CLASSIFICATION))
    return _retire_variant_tags(variant_tags, classification, user)


def retire_requires_classification_tags_for_samples(classification: Classification, samples: Iterable[Sample],
                                                    user: User) -> int:
    """ Same as retire_requires_classification_tags, for the sample/patient page where the case is a set of
        samples rather than an analysis - a tagging with no sample of its own is only retired when it was
        made in an analysis one of the case's samples is in. Someone with read-only access to the analysis can
        still classify, their tagging just stays where it is """
    variant = classification.variant
    if variant is None:
        return 0

    sample_ids = {s.pk for s in samples}
    variant_tags = []
    for variant_tag in VariantTag.objects.filter(variant__in=variant.equivalent_variants,
                                                 tag_id=settings.TAG_REQUIRES_CLASSIFICATION):
        if variant_tag.sample_id:
            in_case = variant_tag.sample_id in sample_ids
        else:
            in_case = bool(variant_tag.analysis) and \
                bool(sample_ids.intersection(s.pk for s in variant_tag.analysis.get_samples()))
        if in_case and variant_tag.can_write(user):
            variant_tags.append(variant_tag)
    return _retire_variant_tags(variant_tags, classification, user)

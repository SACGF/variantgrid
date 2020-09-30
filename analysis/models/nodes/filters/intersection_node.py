from typing import Optional

from django.conf import settings
from django.db import models
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q
from django.db.models.signals import post_delete
from django.dispatch.dispatcher import receiver
import logging
import os
import subprocess

from analysis.models.nodes.analysis_node import AnalysisNode
from annotation.vcf_files.variants_to_vcf import write_qs_to_vcf_file_sort_alphabetically
from genes.hgvs import get_hgvs_variant
from snpdb.models import GenomicIntervalsCollection, GenomicInterval, Sample, \
    VCFBedIntersection, Cohort, VariantCollection
from snpdb.models.models_variant import Variant


class IntersectionNode(AnalysisNode):
    SELECTED_INTERVALS = 0
    CUSTOM_INTERVAL = 1
    HGVS = 2
    BACKEND_ENRICHMENT_KIT = 3

    genomic_intervals_collection = models.ForeignKey(GenomicIntervalsCollection, null=True, blank=True, on_delete=SET_NULL)
    genomic_interval = models.OneToOneField(GenomicInterval, null=True, on_delete=SET_NULL)
    # HGVS linking done in both form and save() method below
    hgvs_name = models.TextField(null=True, blank=True)
    hgvs_variant = models.ForeignKey(Variant, null=True, blank=True, on_delete=SET_NULL)
    left = models.IntegerField(default=0)
    right = models.IntegerField(default=0)
    accordion_panel = models.IntegerField(default=0)

    def valid_custom_genomic_interval(self):
        return self.accordion_panel == self.CUSTOM_INTERVAL and self.genomic_interval

    def valid_selected_genomic_intervals_collection(self):
        return self.accordion_panel == self.SELECTED_INTERVALS and self.genomic_intervals_collection

    def valid_hgvs(self):
        return self.accordion_panel == self.HGVS and self.hgvs_name

    def valid_backend_enrichment_kit(self):
        pbi, _ = self.get_vcf_bed_intersection_and_enrichment_kit()
        return self.accordion_panel == self.BACKEND_ENRICHMENT_KIT and pbi is not None

    def modifies_parents(self):
        return any([self.valid_custom_genomic_interval(),
                    self.valid_selected_genomic_intervals_collection(),
                    self.valid_hgvs(),
                    self.valid_backend_enrichment_kit()])

    @property
    def use_cache(self):
        return super().use_cache or self.valid_selected_genomic_intervals_collection()

    def _get_node_q(self) -> Optional[Q]:
        q = None
        if self.accordion_panel == self.CUSTOM_INTERVAL:
            q_chrom = Variant.get_chrom_q(self.genomic_interval.chrom)
            q = q_chrom & Q(locus__position__gte=self.genomic_interval.start,
                            locus__position__lte=self.genomic_interval.end)
        elif self.accordion_panel == self.HGVS:
            # if hgvs_variant doesn't exist - then it's not in the system - so will always be nothing
            if self.hgvs_variant:
                q = Q(pk=self.hgvs_variant_id)
            else:
                q = self.q_none()
        else:
            variant_collection = None
            if self.accordion_panel == self.SELECTED_INTERVALS:
                raise ValueError("Should never be here - NodeCache should have been generated!")
            if self.accordion_panel == self.BACKEND_ENRICHMENT_KIT:
                pbi, _ = self.get_vcf_bed_intersection_and_enrichment_kit()
                variant_collection = pbi.variant_collection

            if variant_collection:
                q = Q(variantcollectionrecord__variant_collection=variant_collection)
        return q

    def get_vcf_bed_intersection_and_enrichment_kit(self):
        input_sample_ids = self.get_sample_ids()
        num_samples = len(input_sample_ids)

        pbi = None
        enrichment_kit = None
        if num_samples == 1:
            sample_id = input_sample_ids[0]
            sample = Sample.objects.get(pk=sample_id)
            pbi, enrichment_kit = VCFBedIntersection.get_with_enrichment_kit_for_sample(sample)
        elif num_samples > 1:
            all_samples = set(input_sample_ids)
            try:
                extra_cohort_filter_kwargs = {'cohortcount__collection__isnull': False}
                containing_cohort = Cohort.get_cohort_containing_all_samples(all_samples, extra_cohort_filter_kwargs=extra_cohort_filter_kwargs)
                vcf = containing_cohort.vcf
                backend_vcf = vcf.uploadedvcf.backendvcf
                enrichment_kits = list(backend_vcf.sample_sheet.get_sample_enrichment_kits())
                if len(enrichment_kits) == 1:
                    enrichment_kit = enrichment_kits[0]
                    pbi = VCFBedIntersection.get_for_vcf_and_enrichment_kit(vcf, enrichment_kit)
            except:
                pass

        return pbi, enrichment_kit

    def handle_ancestor_input_samples_changed(self):
        pbi, _ = self.get_vcf_bed_intersection_and_enrichment_kit()
        if pbi:
            logging.info("Setting to backend enrichment_kit")
            self.accordion_panel = self.BACKEND_ENRICHMENT_KIT

    def _get_method_summary(self):
        method_summary = 'No filtering applied.'
        if self.modifies_parents():
            if self.accordion_panel == self.CUSTOM_INTERVAL:
                method_summary = f"Filtering to custom interval {self.genomic_interval}"
            elif self.accordion_panel == self.SELECTED_INTERVALS:
                method_summary = f"Filtering to selected interval {self.genomic_intervals_collection.name}"
            elif self.accordion_panel == self.BACKEND_ENRICHMENT_KIT:
                _, enrichment_kit = self.get_vcf_bed_intersection_and_enrichment_kit()
                method_summary = f"Filtering to enrichment_kit {enrichment_kit}"
        return method_summary

    def get_node_name(self):
        name = ''
        if self.modifies_parents():
            if self.accordion_panel == self.CUSTOM_INTERVAL:
                name = str(self.genomic_interval)
            elif self.accordion_panel == self.SELECTED_INTERVALS:
                name = self.genomic_intervals_collection.name
            elif self.accordion_panel == self.BACKEND_ENRICHMENT_KIT:
                _, enrichment_kit = self.get_vcf_bed_intersection_and_enrichment_kit()
                name = f"Enrichment Kit: {enrichment_kit}"
        return name

    def save(self, **kwargs):
        # HGVS name is validated in IntersectionNodeForm, and linked to a variant if one is found
        # But it's possible a variant isn't there at form save time, but will appear later
        # Thus if hgvs_name is set, but hgvs_variant isn't - recheck in save (eg re-connected to diff source node)
        if self.hgvs_name and not self.hgvs_variant:
            self.hgvs_variant = get_hgvs_variant(self.hgvs_name, self.analysis.genome_build)
        return super().save(**kwargs)

    def save_clone(self):
        orig_genomic_interval = self.genomic_interval

        # genomic_interval is a 1-to-1 field, so don't want to copy it in super().save_clone()
        if self.genomic_interval:
            self.genomic_interval = self.genomic_interval.clone()

        copy = super().save_clone()
        self.genomic_interval = orig_genomic_interval
        return copy

    def write_cache(self, variant_collection: VariantCollection):
        if self.valid_selected_genomic_intervals_collection():
            bed_file = self.genomic_intervals_collection.processed_file
            if bed_file is None or not os.path.exists(bed_file):
                msg = f"BED file: {bed_file} does not exist"
                raise ValueError(msg)

            # Open a pipe to intersectBed, which also uploads into the VariantCollection
            args = [settings.INTERSECT_BED_SCRIPT, bed_file, str(variant_collection.pk)]
            intercept_bed_pipe = subprocess.Popen(args, stdin=subprocess.PIPE)  # , stdout=subprocess.PIPE)
            parent_node = self.get_single_parent()
            parent_queryset = parent_node.get_queryset()
            write_qs_to_vcf_file_sort_alphabetically(parent_queryset, intercept_bed_pipe.stdin)
            intercept_bed_pipe.communicate()
        else:
            super().write_cache(variant_collection)

    @staticmethod
    def get_node_class_label():
        return "Intervals intersection"


@receiver(post_delete, sender=IntersectionNode)
def post_delete_intersection_node(sender, instance, *args, **kwargs):
    if instance.genomic_interval is not None:
        instance.genomic_interval.delete()

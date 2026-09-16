import logging
import operator
import os
import subprocess
from functools import cached_property, reduce
from io import TextIOWrapper
from typing import Optional

from auditlog.registry import auditlog
from django.conf import settings
from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.db.models import CASCADE
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q

from analysis.models.nodes.analysis_node import AnalysisNode, NodeAuditLogMixin
from analysis.models.nodes.node_display import NodeIcon
from analysis.variant_text import get_variant_text_summary, parse_region, resolve_variant_text
from snpdb.models import (
    Cohort,
    GenomicIntervalsCollection,
    Sample,
    VariantCollection,
    VCFBedIntersection,
)
from snpdb.models.models_genome import Contig
from snpdb.models.models_variant import Variant
from snpdb.variants_to_vcf import write_qs_to_vcf_file_sort_alphabetically


class IntersectionNode(AnalysisNode):
    # accordion_panel is the accordion index in intersectionnode_editor.html
    SELECTED_INTERVALS = 0
    VARIANTS = 1
    BACKEND_ENRICHMENT_KIT = 2
    CONTIG = 3

    # A pk__in list stops being the cheap option somewhere around here - above it we let
    # NodeCache write the result out to a VariantCollection instead
    VARIANTS_CACHE_MIN = 1000

    genomic_intervals_collection = models.ForeignKey(GenomicIntervalsCollection, null=True, blank=True, on_delete=SET_NULL)
    # Free text entries (one per line) resolved to variants/regions in the form and save() below
    variant_text = models.TextField(null=True, blank=True)
    variant_ids = ArrayField(models.IntegerField(), default=list, blank=True)
    variant_regions = ArrayField(models.TextField(), default=list, blank=True)
    variant_text_unresolved = ArrayField(models.TextField(), default=list, blank=True)
    left = models.IntegerField(default=0)
    right = models.IntegerField(default=0)
    accordion_panel = models.IntegerField(default=0)

    @cached_property
    def contig_ids(self) -> list[int]:
        if self.pk is None:
            return []
        qs = self.intersectionnodecontig_set.filter(
            contig__genomebuildcontig__genome_build=self.analysis.genome_build
        ).order_by("contig__genomebuildcontig__order")
        return list(qs.values_list("contig_id", flat=True))

    def valid_selected_genomic_intervals_collection(self):
        return self.accordion_panel == self.SELECTED_INTERVALS and self.genomic_intervals_collection

    def valid_variants(self):
        # Text that resolved to nothing still filters (to nothing) - that's the honest answer
        return self.accordion_panel == self.VARIANTS and bool(self.variant_text)

    def valid_backend_enrichment_kit(self):
        pbi, _ = self.get_vcf_bed_intersection_and_enrichment_kit()
        return self.accordion_panel == self.BACKEND_ENRICHMENT_KIT and pbi is not None

    def valid_contig(self):
        return self.accordion_panel == self.CONTIG and bool(self.contig_ids)

    def modifies_parents(self):
        return any([self.valid_selected_genomic_intervals_collection(),
                    self.valid_variants(),
                    self.valid_backend_enrichment_kit(),
                    self.valid_contig()])

    @property
    def use_cache(self):
        if self.valid_variants() and len(self.variant_ids) > self.VARIANTS_CACHE_MIN:
            return True
        return super().use_cache or self.valid_selected_genomic_intervals_collection()

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if self.genomic_intervals_collection:
            errors.extend(self._get_genome_build_errors("genomic_intervals_collection",
                                                        self.genomic_intervals_collection.genome_build))
        return errors

    def _get_variants_q(self) -> Q:
        q_list = []
        if self.variant_ids:
            q_list.append(Q(pk__in=self.variant_ids))
        for region in self.variant_regions:
            chrom, start, end = parse_region(region)
            q_list.append(Variant.get_chrom_q(chrom) & Q(locus__position__gte=start, locus__position__lte=end))
        if not q_list:
            # Entries that aren't in the system will always match nothing
            return self.q_none()
        return reduce(operator.or_, q_list)

    def _get_node_q(self) -> Optional[Q]:
        q = None
        if self.accordion_panel == self.VARIANTS:
            q = self._get_variants_q()
        elif self.accordion_panel == self.CONTIG:
            if self.contig_ids:
                q = reduce(operator.or_, [Q(locus__contig_id=c) for c in self.contig_ids])
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
            except Exception:
                pass

        return pbi, enrichment_kit

    def handle_ancestor_input_samples_changed(self):
        AUTO_SWITCH_TO_PANEL_KIT = False
        if AUTO_SWITCH_TO_PANEL_KIT:
            pbi, _ = self.get_vcf_bed_intersection_and_enrichment_kit()
            if pbi:
                logging.info("Setting to backend enrichment_kit")
                self.accordion_panel = self.BACKEND_ENRICHMENT_KIT

    def _get_method_summary(self):
        method_summary = 'No filtering applied.'
        if self.modifies_parents():
            if self.accordion_panel == self.VARIANTS:
                method_summary = f"Filtering to variant entries: {self.variant_text}"
            elif self.accordion_panel == self.SELECTED_INTERVALS:
                method_summary = f"Filtering to selected interval {self.genomic_intervals_collection.name}"
            elif self.accordion_panel == self.BACKEND_ENRICHMENT_KIT:
                _, enrichment_kit = self.get_vcf_bed_intersection_and_enrichment_kit()
                method_summary = f"Filtering to enrichment_kit {enrichment_kit}"
            elif self.accordion_panel == self.CONTIG:
                method_summary = f"Filtering to contigs {self._get_contig_names()}"
        return method_summary

    def _get_contig_names(self) -> str:
        contig_names = list(Contig.objects.filter(pk__in=self.contig_ids,
                                                  genomebuildcontig__genome_build=self.analysis.genome_build)
                            .order_by("genomebuildcontig__order")
                            .values_list("name", flat=True))
        return ", ".join(contig_names)

    def get_node_name(self):
        name = ''
        if self.modifies_parents():
            if self.accordion_panel == self.VARIANTS:
                name = self._get_variant_text_name()
            elif self.accordion_panel == self.SELECTED_INTERVALS:
                name = self.genomic_intervals_collection.name
            elif self.accordion_panel == self.BACKEND_ENRICHMENT_KIT:
                _, enrichment_kit = self.get_vcf_bed_intersection_and_enrichment_kit()
                name = f"Enrichment Kit: {enrichment_kit}"
            elif self.accordion_panel == self.CONTIG:
                name = self._get_contig_names()
        return name

    def _get_variant_text_name(self) -> str:
        entries = self.get_variant_text_entries()
        if len(entries) == 1:
            return entries[0]
        return f"{len(entries)} variant entries"

    def get_variant_text_entries(self) -> list[str]:
        return [line for line in (self.variant_text or "").splitlines() if line.strip()]

    def _get_variant_text_summary(self) -> Optional[str]:
        return get_variant_text_summary(len(self.get_variant_text_entries()), self.variant_ids,
                                        self.variant_regions, self.variant_text_unresolved)

    def get_live_data_notes(self) -> list[str]:
        notes = super().get_live_data_notes()
        if self.valid_variants():
            notes.append(self._get_variant_text_summary())
        return notes

    def get_warnings(self) -> list[str]:
        warnings = super().get_warnings()
        if self.valid_variants() and self.variant_text_unresolved:
            unresolved = ", ".join(self.variant_text_unresolved)
            warnings.append(f"Not found in this database: {unresolved} - "
                            f"use manual variant entry if you want them created.")
        return warnings

    @staticmethod
    def get_help_text() -> str:
        return "Filter to BED files, chromosomes, regions or a list of variants"

    def save(self, *args, **kwargs):
        # Entries are resolved in IntersectionNodeForm, but a clone (eg an analysis template run against
        # a different genome build) carries variant PKs that don't belong to this analysis - re-resolve those
        if self.variant_text and self._variant_text_needs_resolving():
            self._resolve_variant_text()
        return super().save(*args, **kwargs)

    def _variant_text_needs_resolving(self) -> bool:
        if not (self.variant_ids or self.variant_regions or self.variant_text_unresolved):
            return True  # Never resolved - eg a template imported from JSON
        return not self._variant_ids_match_genome_build()

    def _variant_ids_match_genome_build(self) -> bool:
        if not self.variant_ids:
            return True
        genome_build = self.analysis.genome_build
        num_in_build = Variant.objects.filter(Variant.get_contigs_q(genome_build),
                                              pk__in=self.variant_ids).count()
        return num_in_build == len(self.variant_ids)

    def _resolve_variant_text(self):
        resolution = resolve_variant_text(self.analysis.user, self.analysis.genome_build, self.variant_text)
        self.variant_ids = resolution.variant_ids
        self.variant_regions = resolution.regions
        self.variant_text_unresolved = resolution.unresolved + resolution.errors

    def save_clone(self):
        contig_ids = self.contig_ids  # Save before clone
        copy = super().save_clone()
        for contig_id in contig_ids:
            copy.intersectionnodecontig_set.create(contig_id=contig_id)
        return copy

    def write_cache(self, variant_collection: VariantCollection):
        if self.valid_selected_genomic_intervals_collection():
            bed_file = self.genomic_intervals_collection.processed_file
            if bed_file is None or not os.path.exists(bed_file):
                msg = f"BED file: {bed_file} does not exist"
                raise ValueError(msg)

            # Open a pipe to intersectBed, which also uploads into the VariantCollection
            args = [settings.INTERSECT_BED_SCRIPT, bed_file, str(variant_collection.pk)]
            with subprocess.Popen(args, stdin=subprocess.PIPE) as intercept_bed_pipe:
                parent_node = self.get_single_parent()
                parent_queryset = parent_node.get_queryset()
                # VCF is written as text; wrap the binary pipe so encoding happens at the handle
                text_stdin = TextIOWrapper(intercept_bed_pipe.stdin, encoding="utf-8", write_through=True)
                write_qs_to_vcf_file_sort_alphabetically(parent_queryset, text_stdin)
                text_stdin.flush()
                text_stdin.detach()
                intercept_bed_pipe.communicate()
        else:
            super().write_cache(variant_collection)

    @staticmethod
    def get_node_class_label():
        return "Intervals and variants"

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(symbol="node-icon-intervals")

    @classmethod
    def get_node_class_label_short(cls) -> str:
        return "Intervals"


class IntersectionNodeContig(NodeAuditLogMixin, models.Model):
    """ Stores multi-select contig values """
    intersection_node = models.ForeignKey(IntersectionNode, on_delete=CASCADE)
    contig = models.ForeignKey(Contig, on_delete=CASCADE)

    class Meta:
        unique_together = ("intersection_node", "contig")

    def _get_node(self):
        return self.intersection_node

    def __str__(self):
        return f"IntersectionNodeContig {self.intersection_node_id}: {self.contig}"


auditlog.register(IntersectionNode)
auditlog.register(IntersectionNodeContig)

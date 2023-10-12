from typing import List

from django.contrib.auth.models import User
from django.db import models
from django.db.models import CASCADE, SET_NULL, PROTECT, Q
from django_extensions.db.models import TimeStampedModel

from analysis.models.enums import TagLocation
from analysis.models.models_analysis import Analysis
from analysis.models.nodes.analysis_node import AnalysisNode
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsAutoInitialSaveMixin
from snpdb.models import Variant, GenomeBuild, Tag, VariantAllele, Allele


class VariantTagsImport(TimeStampedModel):
    user = models.ForeignKey(User, on_delete=CASCADE)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)

    def __str__(self):
        num_tags = self.importedvarianttag_set.count()
        return f"{num_tags} tags ({self.genome_build}) by {self.user} on {self.created}"


class ImportedVariantTag(models.Model):
    """ Contains all data as not everything will match etc """
    variant_tags_import = models.ForeignKey(VariantTagsImport, on_delete=CASCADE)
    variant_string = models.TextField()
    genome_build_string = models.TextField()
    gene_symbol_string = models.TextField(null=True)
    tag_string = models.TextField()
    variant_id = models.IntegerField()
    analysis_id = models.IntegerField(null=True)
    node_id = models.IntegerField(null=True)
    analysis_name = models.TextField(null=True)
    username = models.TextField()
    created = models.DateTimeField()  # Time on original server


class VariantTag(GuardianPermissionsAutoInitialSaveMixin, TimeStampedModel):
    """ A tag in an analysis. Has create create/delete signal handlers:
        @see analysis.signals.signal_handlers._update_analysis_on_variant_tag_change """
    variant = models.ForeignKey(Variant, on_delete=PROTECT)
    allele = models.ForeignKey(Allele, null=True, on_delete=PROTECT)  # moved to own field as an optimisation
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    tag = models.ForeignKey(Tag, on_delete=CASCADE)
    analysis = models.ForeignKey(Analysis, null=True, on_delete=SET_NULL)
    location = models.CharField(max_length=1, choices=TagLocation.choices, default=TagLocation.ANALYSIS)
    imported_variant_tag = models.ForeignKey(ImportedVariantTag, null=True, on_delete=CASCADE)
    # Most recent node where it was added
    node = models.ForeignKey(AnalysisNode, null=True, on_delete=SET_NULL)  # Keep even if node deleted
    user = models.ForeignKey(User, on_delete=CASCADE)

    def can_view(self, user) -> bool:
        """ Delegate to Analysis if set """
        if self.analysis:
            return self.analysis.can_view(user)
        return super().can_view(user)

    def can_write(self, user) -> bool:
        """ Delegate to Analysis if set """
        if self.analysis:
            return self.analysis.can_write(user)
        return super().can_write(user)

    @property
    def canonical_c_hgvs(self):
        return self.variant.get_canonical_c_hgvs(self.genome_build)

    @property
    def gene_symbol(self):
        gs = None
        if cta := self.variant.get_canonical_transcript_annotation(self.genome_build):
            if tv := cta.transcript_version:
                gs = tv.gene_version.gene_symbol
        return gs

    @staticmethod
    def get_for_build(genome_build: GenomeBuild, tags_qs=None, variant_qs=None):
        """ Returns tags visible within a build
            tags_qs - set to filter - default (None) = all tags
            variant_qs - set to filter - default (None) = all variants """
        if tags_qs is None:
            tags_qs = VariantTag.objects.all()

        va_kwargs = {
            "genome_build": genome_build,
            "allele__in": tags_qs.values_list("allele")
        }
        if variant_qs is not None:
            va_kwargs["variant__in"] = variant_qs

        va_qs = VariantAllele.objects.filter(**va_kwargs)
        # Do distinct as we used to get 2x results with MT (which share contigs across builds)
        tags_qs = VariantTag.objects.filter(allele__in=va_qs.values_list("allele", flat=True))
        return tags_qs.distinct()

    @staticmethod
    def variants_for_build_q(genome_build, tags_qs, tag_ids: List[str]) -> Q:
        if tags_qs is None:
            tags_qs = VariantTag.objects.all()
        if tag_ids:
            tags_qs = tags_qs.filter(tag__in=tag_ids)

        return Q(variantallele__genome_build=genome_build,
                 variantallele__allele__in=tags_qs.values_list("allele"))

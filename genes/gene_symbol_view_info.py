import operator
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime, timedelta
from functools import cached_property, reduce
from typing import Optional, Union

from django.conf import settings
from django.db.models import QuerySet
from django.db.models.query_utils import Q

from analysis.models import VariantTag
from annotation.models import Citation
from annotation.models.models import (
    AnnotationVersion,
    DBNSFPGeneAnnotation,
    DBNSFPGeneAnnotationVersion,
)
from classification.models import ClassificationModification
from classification.models.classification_utils import classification_gene_symbol_filter
from genes.models import (
    HGNC,
    Gene,
    GeneCoverage,
    GeneCoverageCanonicalTranscript,
    GeneList,
    GeneSymbol,
    GeneVersion,
    GnomADGeneConstraint,
    PanelAppServer,
)
from genes.models_enums import AnnotationConsortium
from library.utils import defaultdict_to_dict
from ontology.models import OntologyService, OntologySnake, OntologyTerm
from snpdb.models import Sample, VariantGridColumn, VariantZygosityCountCollection
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_user_settings import UserSettings
from snpdb.variant_queries import get_has_classifications_q, get_variant_queryset_for_gene_symbol


def _get_omim_and_hpo_for_gene_symbol(gene_symbol: GeneSymbol) -> list[tuple[OntologyTerm, list[OntologyTerm]]]:
    omim_and_hpo_for_gene = []
    try:
        # max_depth = 0 for direct links only
        for omim in OntologySnake.terms_for_gene_symbol(gene_symbol, OntologyService.OMIM, max_depth=0).leafs():
            hpo_list = OntologySnake.snake_from(omim, OntologyService.HPO, max_depth=0).leafs()
            omim_and_hpo_for_gene.append((omim, hpo_list))
    except ValueError:  # in case we don't have this gene symbol available
        pass

    return omim_and_hpo_for_gene


@dataclass(frozen=True)
class HasVariants:
    has_tagged_variants: bool
    has_observed_variants: bool
    has_classified_variants: bool

    @property
    def has_variants(self) -> bool:
        return self.has_tagged_variants or self.has_observed_variants or self.has_classified_variants

    def __bool__(self):
        return self.has_variants


class GeneSymbolViewInfo:

    def __init__(self, gene_symbol: GeneSymbol, desired_genome_build: Optional[GenomeBuild], user,
                 tool_tips=None):
        self.gene_symbol = gene_symbol
        self.desired_genome_build = desired_genome_build
        self.tool_tips = tool_tips
        if user is not None:
            self.user = user
            if self.tool_tips is None:
                user_settings = UserSettings.get_for_user(user)
                self.tool_tips = user_settings.tool_tips

    @cached_property
    def omim_and_hpo_for_gene(self) -> list[tuple[OntologyTerm, list[OntologyTerm]]]:
        return _get_omim_and_hpo_for_gene_symbol(self.gene_symbol)

    @cached_property
    def hgnc(self) -> Optional[HGNC]:
        return self.gene_symbol.hgnc_set.order_by("status").first()

    @cached_property
    def citations_ids(self) -> list[str]:
        return sorted(set(Citation.objects.filter(genesymbolcitation__gene_symbol=self.gene_symbol).values_list('id', flat=True)))

    @cached_property
    def dbnsfp_gene_annotation(self) -> Optional[DBNSFPGeneAnnotation]:
        dga = None
        if dbnsfp_gene_version := DBNSFPGeneAnnotationVersion.latest():
            dga = self.gene_symbol.dbnsfpgeneannotation_set.filter(version=dbnsfp_gene_version).first()
        return dga

    @cached_property
    def gene_constraint(self) -> Optional[GnomADGeneConstraint]:
        return GnomADGeneConstraint.objects.filter(gene_symbol=self.gene_symbol).order_by("-mane_select").first()

    @cached_property
    def gene_summary(self) -> Optional[str]:
        refseq_gene: Gene
        if refseq_gene := Gene.objects.filter(annotation_consortium=AnnotationConsortium.REFSEQ,
                                              geneversion__gene_symbol=self.gene_symbol).first():
            return refseq_gene.summary
        return None

    @cached_property
    def gene_version(self) -> Optional[GeneVersion]:
        gene_versions = self.gene_symbol.geneversion_set.all()
        if gene_versions.exists():
            # This page is shown using the users default genome build
            # However - it's possible the gene doesn't exist for a particular genome build.
            # If gene has a version for a build in user settings, use that and show a warning message
            gene_version = gene_versions.filter(genome_build=self.desired_genome_build).first()
            if not gene_version:
                gene_version = gene_versions.first()  # Try another build
            return gene_version
        return None

    @cached_property
    def genome_build(self) -> GenomeBuild:
        if gene_version := self.gene_version:
            return gene_version.genome_build
        return self.desired_genome_build

    def warnings(self) -> list[str]:
        warnings = []
        if self.gene_version:
            # This page is shown using the users default genome build
            # However - it's possible the gene doesn't exist for a particular genome build.
            # If gene has a version for a build in user settings, use that and show a warning message
            if self.genome_build != self.desired_genome_build:
                warnings.append(f"This symbol is not associated with any genes in build {self.desired_genome_build}, viewing in build {self.genome_build}")
        else:
            warnings.append("There are no genes linked against this symbol")
        return warnings

    @cached_property
    def has_variants(self) -> HasVariants:

        has_tagged_variants = False
        has_observed_variants = False
        has_classified_variants = False

        if self.gene_version:
            annotation_version = AnnotationVersion.latest(self.genome_build)
            gene_variant_qs = get_variant_queryset_for_gene_symbol(self.gene_symbol, annotation_version,
                                                                   traverse_aliases=True)
            gene_variant_qs, vzcc = VariantZygosityCountCollection.annotate_global_germline_counts(gene_variant_qs)
            has_observed_variants = gene_variant_qs.filter(**{f"{vzcc.non_ref_call_alias}__gt": 0}).exists()

            has_tagged_variants = VariantTag.get_for_build(self.genome_build, variant_qs=gene_variant_qs).exists()

            # has classifications isn't 100% in sync with the classification table: this code looks at VariantAlleles
            # wheras the classification table will filter on gene symbol and transcript evidence keys
            q = get_has_classifications_q(self.genome_build)
            has_classified_variants = gene_variant_qs.filter(q).exists()

        return HasVariants(
            has_tagged_variants=has_tagged_variants,
            has_observed_variants=has_observed_variants,
            has_classified_variants=has_classified_variants)

    @property
    def has_classified_variants(self):
        return self.has_variants.has_classified_variants

    @cached_property
    def consortium_genes_and_aliases(self) -> dict[str, set[str]]:
        consortium_genes_and_aliases = defaultdict(lambda: defaultdict(set))
        gene: Gene
        for gene in self.gene_symbol.genes:
            aliases = consortium_genes_and_aliases[gene.get_annotation_consortium_display()][gene.identifier]
            aliases.update(gene.get_symbols().exclude(symbol=self.gene_symbol))
        return defaultdict_to_dict(consortium_genes_and_aliases)

    @cached_property
    def gene_external_urls(self) -> dict[str, str]:
        gene_external_urls: dict[str, str] = {}
        for gene in self.gene_symbol.genes:
            gene_external_urls[gene.identifier] = gene.get_external_url()
        return gene_external_urls

    @cached_property
    def annotation_description(self):
        descriptions = {}
        if self.tool_tips:
            descriptions = VariantGridColumn.get_column_descriptions()
            descriptions["gnomad_gene_constraint"] = """
            constraint score shown in gnomAD is the ratio of the observed / expected (oe) number of loss-of-function
            variants in that gene. The expected counts are based on a mutational model that takes sequence context,
            coverage and methylation into account. Low oe values are indicative of strong intolerance. Range is 90%
            confidence interval. <a href='http://gnomad-sg.org/help/constraint'>Details at gnomAD</a> """
            descriptions["essential_gene"] = f"""
                <p><b>CRISPR:</b> {descriptions['essential_gene_crispr']}</p>
                <p><b>CRISPR2:</b>{descriptions['essential_gene_crispr2']}</p>
                <p><b>Gene Trap:</b>{descriptions['essential_gene_gene_trap']}</p>
            """
        return descriptions

    @cached_property
    def has_gene_coverage(self) -> bool:
        if not settings.VIEW_GENE_SYMBOL_SHOW_GENE_COVERAGE:
            return False
        has_gene_coverage = GeneCoverage.get_for_symbol(self.genome_build, self.gene_symbol).exists()
        if has_gene_coverage:
            return True
        has_canonical_gene_coverage = GeneCoverageCanonicalTranscript.get_for_symbol(self.genome_build, self.gene_symbol).exists()
        return has_canonical_gene_coverage

    @property
    def has_samples_in_other_builds(self) -> bool:
        return Sample.objects.exclude(vcf__genome_build=self.genome_build).exists()

    @cached_property
    def gene_in_gene_lists(self) -> bool:
        gene_lists_qs = GeneList.filter_for_user(self.user)
        gene_in_gene_lists = GeneList.visible_gene_lists_containing_gene_symbol(gene_lists_qs, self.gene_symbol).exists()
        return gene_in_gene_lists

    @cached_property
    def classifications(self) -> QuerySet[ClassificationModification]:
        # Note this is loaded in Ajax
        classifications_qs = ClassificationModification.objects.none()
        if filters := classification_gene_symbol_filter(self.gene_symbol):
            classifications = ClassificationModification.objects.filter(filters).filter(
                is_last_published=True).exclude(classification__withdrawn=True)
            classifications_qs = ClassificationModification.filter_for_user(user=self.user, queryset=classifications)

        classifications_qs = classifications_qs.select_related('classification', 'classification__lab', 'classification__allele', 'classification__allele__clingen_allele')
        return classifications_qs

    @cached_property
    def unmatched_classifications(self) -> QuerySet[ClassificationModification]:
        "Only return classifications with no grouping"
        evidence_q_list = []
        for symbol in self.gene_symbol.alias_meta.alias_symbol_strs:
            evidence_q_list.append(Q(published_evidence__gene_symbol__value__iexact=symbol))
        classifications_qs = ClassificationModification.objects.filter(
            classification__allele_info__allele__isnull=True,
            classification__withdrawn=False,
            is_last_published=True,
            classification__created__lte=datetime.now() - timedelta(minutes=1),
        ).filter(reduce(operator.or_, evidence_q_list))
        classifications_qs = ClassificationModification.filter_for_user(user=self.user, queryset=classifications_qs)
        classifications_qs = classifications_qs.select_related('classification', 'classification__lab')
        return sorted(classifications_qs[0:100], key=lambda c: c.curated_date_check, reverse=True)

    @cached_property
    def unmatched_classifications_title(self):
        if count := len(self.unmatched_classifications):
            return f"{count} Unmatched Classification{'s' if count > 1 else ''} for {self.gene_symbol}"

    def panel_app_servers(self) -> Union[QuerySet, Iterable[PanelAppServer]]:
        return PanelAppServer.objects.order_by("pk")

    def show_classifications_hotspot_graph(self) -> bool:
        return settings.VIEW_GENE_HOTSPOT_GRAPH_CLASSIFICATIONS and self.has_variants.has_classified_variants

    def show_hotspot_graph(self) -> bool:
        return settings.VIEW_GENE_HOTSPOT_GRAPH and self.has_variants.has_observed_variants

import operator
from collections import defaultdict
from datetime import date
from functools import reduce
from typing import Optional, Set, List, Dict

from django.db import models
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.query_utils import Q

from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.cohort_mixin import AncestorSampleMixin
from annotation.models import VariantTranscriptAnnotation, OntologyTerm
from genes.models import GeneSymbol
from ontology.models import GeneDiseaseClassification, OntologyTermRelation
from patients.models_enums import Zygosity
from snpdb.models import Contig, Sample


class MOINode(AncestorSampleMixin, AnalysisNode):
    PANEL_CUSTOM = 0
    PANEL_PATIENT = 1

    # Sample isn't mandatory, but if you supply it, you can use the zygosity
    # Probably want to be able to swap the panel out if sample changes as per cohort node and the zyg filter
    sample = models.ForeignKey(Sample, null=True, blank=True, on_delete=SET_NULL)
    require_zygosity = models.BooleanField(default=True)
    min_classification = models.CharField(max_length=1, choices=GeneDiseaseClassification.choices,
                                          default=GeneDiseaseClassification.MODERATE)
    min_date = models.DateField(null=True, blank=True)
    max_date = models.DateField(null=True, blank=True)
    accordion_panel = models.IntegerField(default=0)

    @property
    def use_patient(self):
        return self.accordion_panel == self.PANEL_PATIENT

    def modifies_parents(self):
        return bool(self.get_gene_disease_relations())

    def get_gene_disease_relations(self) -> List[OntologyTermRelation]:
        """ Filtered by node settings """
        gene_disease_relations = []
        ontology_terms = self._get_all_ontology_term_ids()
        gene_disease_qs = self.analysis.annotation_version.ontology_version.get_gene_disease_relations_qs()
        gene_disease_qs = gene_disease_qs.filter(source_term__in=ontology_terms)

        valid_classifications = set(GeneDiseaseClassification.get_above_min(self.min_classification))
        mode_of_inheritance = set(self.moinodemodeofinheritance_set.all().values_list("mode_of_inheritance", flat=True))
        submitters = set(self.moinodesubmitter_set.all().values_list("submitter", flat=True))

        for otr in gene_disease_qs:
            # Only include the sources that match, so we can use that to build docs / grids
            filtered_sources = []
            for source in otr.extra["sources"]:
                if source["gencc_classification"] not in valid_classifications:
                    continue
                if mode_of_inheritance:
                    if source["mode_of_inheritance"] not in mode_of_inheritance:
                        continue
                if submitters:
                    if source["submitter"] not in submitters:
                        continue
                if self.min_date or self.max_date:
                    evaluated_date = date.fromisoformat(source["evaluated_date"])
                    if self.min_date and self.min_date > evaluated_date:
                        continue
                    if self.max_date and self.max_date < evaluated_date:
                        continue
                filtered_sources.append(source)

            if filtered_sources:
                otr.extra["sources"] = filtered_sources
                # Override the save method, so it can't be saved by accident
                otr.save = None
                gene_disease_relations.append(otr)

        return gene_disease_relations

    def get_gene_symbols_qs(self):
        hgnc_names = set()
        for otr in self.get_gene_disease_relations():
            hgnc_names.add(otr.dest_term.name)
        gene_symbols_qs = GeneSymbol.objects.filter(symbol__in=hgnc_names)
        return gene_symbols_qs

    def _get_all_ontology_term_ids(self):
        """ Returns all terms (not filtered by any settings) """
        ontology_term_ids = []
        if self.accordion_panel == self.PANEL_PATIENT:
            if self.sample:
                if patient := self.sample.patient:
                    ontology_term_ids = patient.get_ontology_term_ids()
        else:
            ontology_term_ids = self.moinodeontologyterm_set.values_list("ontology_term", flat=True)
        return ontology_term_ids

    def _get_gene_qs(self):
        gene_symbols_qs = self.get_gene_symbols_qs()
        return self.analysis.gene_annotation_release.genes_for_symbols(gene_symbols_qs)

    def _get_zygosities(self, moi: str) -> Set[str]:
        het_or_hom = {Zygosity.HET, Zygosity.HOM_ALT}
        recessive_zygosities = {Zygosity.HOM_ALT}
        if self.require_zygosity:
            dominant_zygosities = {Zygosity.HET}
        else:
            dominant_zygosities = {Zygosity.HET, Zygosity.HOM_ALT}

        MOI_ZYGOSITY_Q = {
            'Autosomal dominant': dominant_zygosities,
            'Autosomal dominant inheritance with maternal imprinting HP:0012275': dominant_zygosities,
            'Autosomal dominant inheritance with paternal imprinting': dominant_zygosities,
            'Autosomal recessive': dominant_zygosities,
            # Digenic = A type of multifactorial inheritance governed by the simultaneous action of two gene loci.
            'Digenic inheritance HP:0010984': het_or_hom,
            'Mitochondrial': het_or_hom,
            # Semidominant - A/A homozygote has a mutant phenotype
            #                A/a heterozygote has a less severe phenotype
            #                a/a homozygote (REF) is wild-type.
            'Semidominant': het_or_hom,
            'Somatic mosaicism': None,  # Could be called as anything, even Ref
            'Unknown': het_or_hom,
            'X-linked': het_or_hom,
            'X-linked dominant': dominant_zygosities,
            'X-linked recessive': recessive_zygosities,
            'Y-linked inheritance': None,  # Anything on Y will be HOM anyway, don't want to exclude due to
        }
        return MOI_ZYGOSITY_Q[moi]

    def _get_zygosity_q(self, moi: str) -> Optional[Q]:
        q = None
        if zygosities := self._get_zygosities(moi):
            _alias, field = self.sample.get_cohort_genotype_alias_and_field("zygosity")
            q = Q(**{f"{field}__in": zygosities})
        return q

    def _get_genes_q_from_hgnc(self, hgnc_names: Set[str]):
        gene_symbols_qs = GeneSymbol.objects.filter(symbol__in=hgnc_names)
        gene_qs = self.analysis.gene_annotation_release.genes_for_symbols(gene_symbols_qs)
        variant_annotation_version = self.analysis.annotation_version.variant_annotation_version
        return VariantTranscriptAnnotation.get_overlapping_genes_q(variant_annotation_version, gene_qs)

    def _get_moi_genes(self):
        moi_genes = defaultdict(set)
        for otr in self.get_gene_disease_relations():
            for source in otr.extra["sources"]:
                moi_genes[source["mode_of_inheritance"]].add(otr.dest_term.name)
        return moi_genes

    def _get_node_arg_q_dict(self) -> Dict[Optional[str], Dict[str, Q]]:
        arg_q_dict = {}
        if self.sample:
            moi_genes = self._get_moi_genes()
            or_filters = []
            for moi, hgnc_names in moi_genes.items():
                q_genes = self._get_genes_q_from_hgnc(hgnc_names)
                if q_zygosity := self._get_zygosity_q(moi):
                    or_filters.append(q_zygosity & q_genes)
                else:
                    or_filters.append(q_genes)  # Any zygosity
            q = reduce(operator.or_, or_filters)
            arg_q_dict[self.sample.zygosity_alias] = {str(q): q}
        else:
            gene_qs = self._get_gene_qs()
            variant_annotation_version = self.analysis.annotation_version.variant_annotation_version
            q_genes = VariantTranscriptAnnotation.get_overlapping_genes_q(variant_annotation_version, gene_qs)
            arg_q_dict[None] = {str(q_genes): q_genes}
        return arg_q_dict

    def _get_node_contigs(self) -> Optional[Set[Contig]]:
        contig_qs = Contig.objects.filter(transcriptversion__genome_build=self.analysis.genome_build,
                                          transcriptversion__gene_version__gene__in=self._get_gene_qs())
        return set(contig_qs.distinct())

    def _get_method_summary(self):
        if self.modifies_parents():
            method_list = []
            if self.sample:
                moi_genes = self._get_moi_genes()
                li_list = []
                for moi, hgnc_names in moi_genes.items():
                    if zygosities := self._get_zygosities(moi):
                        zygosity = ", ".join(sorted([Zygosity.display(z) for z in zygosities]))
                    else:
                        zygosity = "Any Zygosity"
                    li_list.append(f"<li>{moi} - {zygosity} - {', '.join(sorted(hgnc_names))}</li>")

                method_list.append(f"<ul>{''.join(li_list)}</ul>")
                method_summary = "".join([f"<p>{s}</p>" for s in method_list])
            else:
                gene_symbols = self.get_gene_symbols_qs()
                method_summary = f"{', '.join(sorted(gene_symbols))}"
        else:
            method_summary = 'No filters applied.'

        return method_summary

    def get_node_name(self):
        name = ''
        if ontology_term_ids := self._get_all_ontology_term_ids():
            terms = sorted([ot.name for ot in OntologyTerm.objects.filter(pk__in=ontology_term_ids)])
            terms_text = ", ".join(terms)
            if self.accordion_panel == self.PANEL_PATIENT:
                name = f"{terms_text or 'No terms'} (from patient)"
            else:
                name = terms_text
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Filter to curated Gene/Disease relationships"

    def save_clone(self):
        ontology_term_list = list(self.moinodeontologyterm_set.all())
        moi_list = list(self.moinodemodeofinheritance_set.all())
        submitter_list = list(self.moinodesubmitter_set.all())
        copy = super().save_clone()

        for moi_ot in ontology_term_list:
            copy.moinodeontologyterm_set.create(ontology_term=moi_ot.ontology_term)

        for moi_moi in moi_list:
            copy.moinodemodeofinheritance_set.create(mode_of_inheritance=moi_moi.mode_of_inheritance)

        for moi_submitter in submitter_list:
            copy.moinodesubmitter_set.create(submitter=moi_submitter.moi_submitter)

        return copy

    @staticmethod
    def get_node_class_label():
        return "Mode of Inheritance"


class MOINodeOntologyTerm(models.Model):
    node = models.ForeignKey(MOINode, on_delete=CASCADE)
    ontology_term = models.ForeignKey(OntologyTerm, on_delete=CASCADE)

    class Meta:
        unique_together = ("node", "ontology_term")


class MOINodeModeOfInheritance(models.Model):
    """ If you have none of these, use any """
    node = models.ForeignKey(MOINode, on_delete=CASCADE)
    mode_of_inheritance = models.TextField()

    class Meta:
        unique_together = ("node", "mode_of_inheritance")


class MOINodeSubmitter(models.Model):
    """ If you have none of these, use any """
    node = models.ForeignKey(MOINode, on_delete=CASCADE)
    submitter = models.TextField()

    class Meta:
        unique_together = ("node", "submitter")

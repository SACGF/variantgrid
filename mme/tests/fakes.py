""" Lightweight fakes for building MME patient profiles without the heavy
    Classification/variant-matching machinery, plus a real-DB builder for the tests that
    need eligibility (a derived queryset filter) to actually resolve. """
from dataclasses import dataclass, field
from typing import Optional

from classification.enums import SubmissionSource
from classification.enums.classification_enums import ClinicalSignificance, ShareLevel, SpecialEKeys
from classification.models.classification import (
    Classification,
    ClassificationModification,
    ConditionResolved,
)
from genes.models import Gene, GeneAnnotationImport, GeneSymbol, GeneVersion
from ontology.models import OntologyService, OntologyTerm
from snpdb.models import GenomeBuild


def make_classification(lab, user, gene_symbol: str = "BRCA1",
                        share_level: str = ShareLevel.PUBLIC.value,
                        clinical_significance: str = ClinicalSignificance.VUS,
                        withdrawn: bool = False,
                        is_last_published: bool = True) -> Classification:
    """ A real Classification with a deterministic "latest published" modification, so
        mme_eligible_classifications() (a queryset filter on that modification) resolves.
        Defaults produce an ELIGIBLE record; vary one argument per test to break one layer. """
    vc = Classification.create(
        user=user, lab=lab, data={SpecialEKeys.GENE_SYMBOL: {'value': gene_symbol}},
        save=True, source=SubmissionSource.API)
    if withdrawn:
        vc.withdrawn = True
        vc.save()
    ClassificationModification.objects.filter(classification=vc).update(is_last_published=False)
    ClassificationModification.objects.create(
        classification=vc, user=user, source="TEST", delta={},
        share_level=share_level, clinical_significance=clinical_significance,
        is_last_published=is_last_published, published=True)
    return vc


def make_gene_version(gene_id: str, symbol: str, annotation_consortium: str,
                      version: int = 1) -> GeneVersion:
    """ One (gene, symbol) row of the gene annotation MME resolves identifiers through.
        Call twice with the same gene_id and different symbols/versions to model a rename. """
    genome_build = GenomeBuild.grch38()
    gene, _ = Gene.objects.get_or_create(identifier=gene_id,
                                         annotation_consortium=annotation_consortium)
    import_source, _ = GeneAnnotationImport.objects.get_or_create(
        url="fake", genome_build=genome_build, annotation_consortium=annotation_consortium)
    gene_symbol, _ = GeneSymbol.objects.get_or_create(symbol=symbol)
    gene_version, _ = GeneVersion.objects.get_or_create(
        gene=gene, gene_symbol=gene_symbol, version=version,
        genome_build=genome_build, import_source=import_source)
    return gene_version


def make_term(term_id: str, service: OntologyService, index: int, name: str) -> OntologyTerm:
    """ Unsaved OntologyTerm with full control over routing-relevant fields. """
    return OntologyTerm(id=term_id, ontology_service=service, index=index, name=name)


@dataclass
class FakeCoordinate:
    chrom: str
    position: int
    ref: str
    alt: str
    is_symbolic: bool = False


@dataclass
class FakeVariant:
    coordinate: FakeCoordinate


@dataclass
class FakeGenomeBuild:
    name: str


@dataclass
class FakeOrganization:
    name: str = ""


@dataclass
class FakeLab:
    name: str = ""
    contact_name: str = ""
    contact_email: str = ""
    url: str = ""
    organization: Optional[FakeOrganization] = None


@dataclass
class FakeResolvedVariantInfo:
    gene_symbol: Optional[GeneSymbol] = None


@dataclass
class FakeAlleleInfo:
    """ Stands in for ImportedAlleleInfo, whose only role here is carrying the resolved
        (curated-transcript) gene symbol MME prefers over the evidence key. """
    grch37: Optional[FakeResolvedVariantInfo] = None
    grch38: Optional[FakeResolvedVariantInfo] = None


@dataclass
class FakeClassification:
    pk: int = 1
    terms: list = field(default_factory=list)
    gene_symbol: Optional[str] = None
    variant: Optional[FakeVariant] = None
    genome_build_name: str = "GRCh38"
    sample = None
    has_build: bool = True
    lab: Optional[FakeLab] = None
    lab_id: Optional[int] = None
    allele_info: Optional[FakeAlleleInfo] = None

    @property
    def condition_resolution_obj(self) -> Optional[ConditionResolved]:
        if self.terms:
            return ConditionResolved(terms=list(self.terms))
        return None

    def get(self, key):
        if key == SpecialEKeys.GENE_SYMBOL:
            return self.gene_symbol
        return None

    def get_genome_build(self) -> FakeGenomeBuild:
        if not self.has_build:
            raise ValueError("Classification does not have a value for genome build")
        return FakeGenomeBuild(self.genome_build_name)

    def get_variant_for_build(self, genome_build) -> Optional[FakeVariant]:
        return self.variant


@dataclass
class FakeSubmission:
    classification: FakeClassification
    external_patient_id: str = "vg:1"
    node_id: str = "testnode"
    classification_id: int = 1

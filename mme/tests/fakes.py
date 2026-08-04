""" Lightweight fakes for building MME patient profiles without the heavy
    Classification/variant-matching machinery, plus a real-DB builder for the tests that
    need eligibility (a derived queryset filter) to actually resolve. """
from dataclasses import dataclass, field
from typing import Optional

from classification.enums import SubmissionSource
from classification.enums.classification_enums import SpecialEKeys, ShareLevel, ClinicalSignificance
from classification.models.classification import (
    ConditionResolved, Classification, ClassificationModification,
)
from ontology.models import OntologyTerm, OntologyService


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

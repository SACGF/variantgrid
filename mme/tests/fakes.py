""" Lightweight fakes for building MME patient profiles without the heavy
    Classification/variant-matching machinery. """
from dataclasses import dataclass, field
from typing import Optional

from classification.enums.classification_enums import SpecialEKeys
from classification.models.classification import ConditionResolved
from ontology.models import OntologyTerm, OntologyService


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

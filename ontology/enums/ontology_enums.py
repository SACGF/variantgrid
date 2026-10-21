"""
Enums for ontology terms: OntologyService (the ontology prefixes - MONDO, OMIM, HPO, HGNC stored locally, the rest
linked externally - with their id lengths, URLs, display order and prefix aliases), OntologyTermStatus,
and the gene-disease strength scales PanelAppClassification and GeneDiseaseClassification (GenCC).
Stored codes: add values, never change one.
"""
import re
from typing import Optional, Union

from django.db import models
from django.db.models import TextChoices

from library.utils import Constant


class OntologyService(models.TextChoices):
    MONDO = "MONDO", "MONDO"
    OMIM = "OMIM", "OMIM"
    HPO = "HP", "HP"
    HGNC = "HGNC", "HGNC"

    #Not stored locally
    DOID = "DOID", "DOID"
    ORPHANET = "ORPHA", "ORPHA"
    MEDGEN = "MedGen", "MedGen"
    MeSH = "MeSH", "MeSH"

    EXPECTED_LENGTHS: dict[str, int] = Constant({
        MONDO[0]: 7,
        OMIM[0]: 6,
        HPO[0]: 7,
        HGNC[0]: 1,  # HGNC ids aren't typically 0 padded, because they're not monsters
        DOID[0]: None,  # variable length with padded 0s
        ORPHANET[0]: 1,  # ORPHANET ids aren't typically 0 padded
        MEDGEN[0]: None,
        MeSH[0]: None
    })

    IMPORTANCE: dict[str, int] = Constant({
        MONDO[0]: 2,
        OMIM[0]: 3,
        HPO[0]: 4,  # put HPO relationships last as they occasionally spam OMIM
        DOID[0]: 5,
        ORPHANET[0]: 6,
        HGNC[0]: 1,  # show gene relationships first
        MEDGEN[0]: 7,
        MeSH[0]: 8
    })

    URLS: dict[str, str] = Constant({
        MONDO[0]: "https://monarchinitiative.org/disease/MONDO:${1}",
        OMIM[0]: "http://www.omim.org/entry/${1}",
        HPO[0]: "https://hpo.jax.org/app/browse/term/HP:${1}",
        HGNC[0]: "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:${1}",
        DOID[0]: "https://www.ebi.ac.uk/ols/ontologies/doid/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FDOID_${1}",
        ORPHANET[0]: "https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=EN&Expert=${1}",
        MEDGEN[0]: "https://www.ncbi.nlm.nih.gov/medgen/${1}",
        MeSH[0]: "https://meshb.nlm.nih.gov/record/ui?ui=${1}"
    })

    LOCAL_ONTOLOGY_PREFIXES: set[str] = Constant({
        MONDO[0],
        OMIM[0],
        HPO[0],
        HGNC[0]
    })

    CONDITION_ONTOLOGIES: set[str] = Constant({
        MONDO[0],
        OMIM[0],
        HPO[0],
        DOID[0],
        ORPHANET[0],
        MEDGEN[0],
        MeSH[0]
    })

    # Alternate external spellings that differ from our canonical value (e.g. Monarch returns
    # "Orphanet" for "ORPHA"). These can't be derived from the enum, so they're listed here.
    # Prefixes that only differ by case are handled automatically by iterating the enum.
    PREFIX_ALIASES: dict[str, str] = Constant({
        "ORPHANET": ORPHANET[0],
        "MIM": OMIM[0],
        "HPO": HPO[0],
    })

    @classmethod
    def resolve_prefix(cls, prefix: str) -> Optional['OntologyService']:
        """ Map a (possibly aliased or differently-cased) prefix to a supported OntologyService,
            or None if it isn't one we support (e.g. MPATH from the external Monarch search). """
        upper = prefix.strip().upper()
        # Canonical values matched case-insensitively - derived from the enum so a newly added
        # OntologyService is recognised without touching PREFIX_ALIASES.
        by_value = {service.value.upper(): service for service in cls}
        if service := by_value.get(upper):
            return service
        if canonical := cls.PREFIX_ALIASES.get(upper):
            return cls(canonical)
        return None

    @classmethod
    def is_supported_id(cls, term_id: str) -> bool:
        """ True when term_id is a well-formed <prefix>:<postfix> id whose prefix is a supported
            OntologyService. Unsupported ontologies (e.g. MPATH from the external Monarch search)
            and malformed ids both return False. """
        parts = re.split("[:|_]", term_id or "")
        if len(parts) != 2:
            return False
        return cls.resolve_prefix(parts[0]) is not None

    @staticmethod
    def index_to_id(ontology_service: 'OntologyService', index: Union[int|str]):
        num_part = str(index)
        if expected_length := OntologyService.EXPECTED_LENGTHS.get(ontology_service):
            num_part = str(index).rjust(expected_length, '0')
        return f"{ontology_service}:{num_part}"


class PanelAppClassification(models.TextChoices):
    RED = "1", "Expert Review Red"
    AMBER = "2", "Expert Review Amber"
    GREEN = "3", "Expert Review Green"

    @property
    def is_strong_enough(self) -> bool:
        return self == PanelAppClassification.GREEN

    def __lt__(self, other) -> bool:
        return int(self.value) < int(other.value)

    @staticmethod
    def get_by_label_pac(label: str) -> 'PanelAppClassification':
        for pac in PanelAppClassification:
            if pac.label == label:
                return pac
        raise ValueError(f"No PanelAppClassification for {label}")

    @staticmethod
    def get_above_min(min_classification: str) -> set[str]:
        classifications = set()
        for e in reversed(PanelAppClassification):
            classifications.add(e.label)
            if e.value == min_classification:
                break
        return classifications


class GeneDiseaseClassification(models.TextChoices):
    # @see https://thegencc.org/faq.html#validity-termsdelphi-survey - where sort order comes from
    # Not using e.g. "GENCC:100009" (= Supportive) as that has different sort order than page above
    REFUTED = "1", "Refuted Evidence"
    NO_KNOWN = "2", "No Known Disease Relationship"
    ANIMAL = "3", "Animal Model Only"
    DISPUTED = "4", "Disputed Evidence"
    LIMITED = "5", "Limited"
    SUPPORTIVE = "6", "Supportive"
    MODERATE = "7", "Moderate"
    STRONG = "8", "Strong"
    DEFINITIVE = "9", "Definitive"

    @property
    def is_strong_enough(self) -> bool:
        return self in {GeneDiseaseClassification.STRONG, GeneDiseaseClassification.DEFINITIVE}

    @staticmethod
    def get_by_label(label: str) -> 'GeneDiseaseClassification':
        for gdc in GeneDiseaseClassification:
            if gdc.label == label:
                return gdc
        raise ValueError(f"No GeneDiseaseClassification for {label}")

    @staticmethod
    def get_above_min(min_classification: str) -> set[str]:
        classifications = set()
        for e in reversed(GeneDiseaseClassification):
            classifications.add(e.label)
            if e.value == min_classification:
                break
        return classifications


class OntologyTermStatus(TextChoices):
    CONDITION = 'C'  # Also phenotypes are mixed up in there right now
    DEPRECATED = 'D'
    NON_CONDITION = 'N'
    STUB = 'S'

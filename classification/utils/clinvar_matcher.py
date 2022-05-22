import csv
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Set, Iterator

from lazy import lazy
from pyhgvs import InvalidHGVSName

from annotation.models import ClinVar
from classification.enums import SpecialEKeys
from classification.models import Classification, ClinVarExport, ClinVarAllele
from genes.hgvs import CHGVS
from library.guardian_utils import admin_bot
from ontology.models import OntologyTerm, OntologySnake
from snpdb.models import GenomeBuild, Variant, Allele, ClinVarKey
from variantopedia.search import search_hgvs, SearchResult, ClassifyVariant
import re

C_HGVS_AND_P_DOT = re.compile(r"(?P<c_hgvs>,*?) ?(?P<p_hgvs>/(p[.].*))")


class ClinVarLegacyMatchType(str, Enum):
    CLINVAR_VARIANT_ID = "clinvar_variant_id"
    VARIANT_PREFERRED_VARIANT = "variant_preferred_variant"
    VARIANT_PREFERRED_IMPORTED_C_HGVS = "variant_preferred_imported_c_hgvs"


@dataclass
class ClinVarLegacyMatch:
    clinvar_export: ClinVarExport
    c_hgvs_matches: bool
    condition_matches: bool
    clin_sig_matches: bool


@dataclass
class ClinVarLegacyMatches:
    allele: Allele
    match_types: Set[ClinVarLegacyMatchType]
    classifications: List[Classification]
    clinvar_export_matches: List[ClinVarLegacyMatch]


CLINICAL_SIGNIFICANCE_MAP = {
    "Benign": "B",
    "Likely benign": "LB",
    "Uncertain significance": "VUS",
    "Likely pathogenic": "LP",
    "Pathogenic": "P"
}


class ClinVarLegacyColumn(str, Enum):
    VariationID = "VariationID"
    AlleleID = "AlleleID"
    Your_record_id = "Your_record_id"
    SCV = "SCV"
    RCV = "RCV"
    Your_variant_description_HGVS = "Your_variant_description_HGVS"
    Your_variant_description_chromosome_coordinates = "Your_variant_description_chromosome_coordinates"
    Preferred_variant_name = "Preferred_variant_name"
    Your_condition_name = "Your_condition_name"
    Your_condition_identifier = "Your_condition_identifier"
    ClinVar_condition_name = "ClinVar_condition_name"
    Clinical_significance = "Clinical_significance"
    Date_last_evaluated = "Date_last_evaluated"
    Assertion_criteria = "Assertion_criteria"
    Submitted_gene = "Submitted_gene"


@dataclass
class ClinVarLegacyRow:
    clinvar_key: ClinVarKey
    data: Dict[str, str]

    def get_column(self, field_name: str):
        return self.data[field_name]

    @property
    def labs(self):
        return list(self.clinvar_key.lab_set.all())

    @property
    def scv(self):
        return self.get_column(ClinVarLegacyColumn.SCV)

    @property
    def variant_clinvar_id(self):
        return int(self.get_column(ClinVarLegacyColumn.VariationID))

    @property
    def variant_description(self):
        return self.get_column(ClinVarLegacyColumn.Your_variant_description_HGVS)

    @property
    def variant_preferred(self):
        return self.get_column(ClinVarLegacyColumn.Preferred_variant_name)

    @property
    def condition_identifier(self):
        return self.get_column(ClinVarLegacyColumn.Your_condition_identifier)

    @property
    def clinical_significance(self):
        return self.get_column(ClinVarLegacyColumn.Clinical_significance)

    @property
    def submitted_gene(self):
        return self.get_column(ClinVarLegacyColumn.Submitted_gene)

    @staticmethod
    def load_file(file, clinvar_key: ClinVarKey) -> Iterator['ClinVarLegacyRow']:
        csv_f = csv.reader(file, delimiter='\t')
        header_indexes: Dict[str, int] = dict()
        for row in csv_f:
            if not header_indexes:
                if len(row) > 1 and row[1] == 'VariationID':
                    # work out weird header row halfway down the file (and one column over)
                    # don't process the header row itself as data, now skip to the next row
                    header_indexes = {label: index for index, label in enumerate(row)}
            else:
                row_dict: Dict[str, str] = dict()
                for header_key, index in header_indexes.items():
                    row_dict[header_key] = row[index]
                yield ClinVarLegacyRow(
                    clinvar_key=clinvar_key,
                    data=row_dict
                )
        if not header_indexes:
            raise ValueError("Did not find the header row")

    @property
    def c_hgvs_preferred_str(self) -> str:
        if match := C_HGVS_AND_P_DOT.match(self.variant_preferred):
            return match.group('c_hgvs')
        else:
            return self.variant_preferred

    @property
    def c_hgvs_with_gene_symbol(self) -> Optional[CHGVS]:
        if c_hgvs_str := self.c_hgvs_preferred_str:
            c_hgvs = CHGVS(c_hgvs_str)
            c_hgvs.gene_symbol = self.submitted_gene
            return c_hgvs
        return None

    @property
    def clinical_significance_code(self) -> str:
        return CLINICAL_SIGNIFICANCE_MAP.get(self.clinical_significance)

    @lazy
    def ontology_terms(self) -> List[OntologyTerm]:
        terms: List[OntologyTerm] = list()
        if individual_terms := self.condition_identifier.split(';'):
            for individual_term in individual_terms:
                if individual_term.startswith("Human Phenotype Ontology:"):
                    individual_term = individual_term.split(":", maxsplit=1)[1]
                terms.append(OntologyTerm.get_or_stub(individual_term))
        return terms

    def find_variant_grid_allele(self) -> List[ClinVarLegacyMatches]:
        allele_to_match_types: Dict[Allele, Set[ClinVarLegacyMatchType]] = defaultdict(set)

        if clinvar_annotation := ClinVar.objects.filter(clinvar_variation_id=self.variant_clinvar_id).order_by('-version').first():
            if allele := clinvar_annotation.variant.allele:
                allele_to_match_types[allele].add(ClinVarLegacyMatchType.CLINVAR_VARIANT_ID)

        try:
            if results := search_hgvs(self.c_hgvs_preferred_str, user=admin_bot(), genome_build=GenomeBuild.grch38(), variant_qs=Variant.objects.all()):
                for result in results:
                    variant: Optional[Variant] = None
                    if isinstance(result, SearchResult):
                        result = result.record
                    if isinstance(result, Variant):
                        variant = result
                    elif isinstance(result, ClassifyVariant):
                        variant = result.variant
                    if variant:
                        if allele := variant.allele:
                            allele_to_match_types[allele].add(ClinVarLegacyMatchType.VARIANT_PREFERRED_VARIANT)
        except InvalidHGVSName:
            pass
        except NotImplementedError:
            pass

        if c_hgvs := self.c_hgvs_with_gene_symbol:
            # allow for some transcript version increases
            for attempt_increase in range(0, 3):
                if attempt_increase and c_hgvs.transcript_parts.version:
                    c_hgvs = c_hgvs.with_transcript_version(c_hgvs.transcript_parts.version + attempt_increase)
                if c_hgvs:
                    if allele_ids := set(Classification.objects.filter(evidence__c_hgvs__value=str(c_hgvs), lab__in=self.labs).values_list('allele', flat=True)):
                        alleles = Allele.objects.filter(pk__in=allele_ids)
                        for allele in alleles:
                            allele_to_match_types[allele].add(ClinVarLegacyMatchType.VARIANT_PREFERRED_IMPORTED_C_HGVS)

        all_matches: [ClinVarLegacyMatches] = list()
        for allele, match_types in allele_to_match_types.items():
            classifications = list(Classification.objects.filter(lab__in=self.labs, allele=allele))
            clinvar_export_matches: List[ClinVarLegacyMatch] = list()
            if clinvar_allele := ClinVarAllele.objects.filter(allele=allele, clinvar_key=self.clinvar_key):
                if clinvar_exports := ClinVarExport.objects.filter(clinvar_allele=clinvar_allele):
                    for clinvar_export in clinvar_exports:

                        condition_matches = False
                        clin_sig_matches = False
                        c_hgvs_matches = False

                        if classification_based_on := clinvar_export.classification_based_on:
                            clinical_significance = classification_based_on.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
                            clin_sig_matches = clinical_significance == self.clinical_significance_code
                            c_hgvs_matches = CHGVS(classification_based_on.get(SpecialEKeys.C_HGVS)) == self.c_hgvs_with_gene_symbol

                        umbrella = clinvar_export.condition_resolved
                        for term in self.ontology_terms:
                            if OntologySnake.snake_from(term, umbrella, max_depth=3):
                                condition_matches = True

                        clinvar_export_match = ClinVarLegacyMatch(
                            clinvar_export=clinvar_export,
                            c_hgvs_matches=c_hgvs_matches,
                            condition_matches=condition_matches,
                            clin_sig_matches=clin_sig_matches
                        )
                        clinvar_export_matches.append(clinvar_export_match)
            ccwm = ClinVarLegacyMatches(
                allele=allele,
                match_types=match_types,
                classifications=classifications,
                clinvar_export_matches=clinvar_export_matches
            )
            all_matches.append(ccwm)
        return all_matches

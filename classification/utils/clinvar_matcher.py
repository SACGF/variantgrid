"""

This file is to do with matching records that were manually submitting to ClinVar with records in VariantGrid.

The typical scenario is a lab would have submitted records in Excel to ClinVar, and later signed up with Shariant
and started submitting to ClinVar via Shariant. So now we want to avoid submitting those records as new records.

When you see the word "Legacy" below, that's an in the record was submitting to ClinVar using the legacy Excel format.

"""

import csv
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Set, Iterator

from lazy import lazy
from pyhgvs import InvalidHGVSName

from annotation.models import ClinVar
from classification.enums import SpecialEKeys
from classification.models import Classification, ClinVarExport, ClinVarAllele, EvidenceKeyMap, \
    ConditionResolved
from genes.hgvs import CHGVS
from library.guardian_utils import admin_bot
from ontology.models import OntologyTerm, OntologySnake, OntologyTermRelation
from snpdb.models import GenomeBuild, Variant, Allele, ClinVarKey
from variantopedia.search import search_hgvs, SearchResult, ClassifyVariant

C_HGVS_AND_P_DOT = re.compile(r"^(?P<c_hgvs>.+?)( \((?P<p_hgvs>p[.].+)\))?$")


class ClinVarLegacyAlleleMatchType(str, Enum):
    CLINVAR_VARIANT_ID = "ClinVar VariantID matches"
    VARIANT_PREFERRED_VARIANT = "ClinVar c.HGVS matches"
    VARIANT_PREFERRED_IMPORTED_C_HGVS = "imported c.HGVS matches"

    def __str__(self):
        return self.value

    def __repr__(self):
        return self.value


class ClinVarLegacyExportMatchType(str, Enum):
    CONDITION_MATCHES = "condition matches"
    CLINICAL_SIGNIFICANCE_MATCHES = "clin sig matches"
    SCV_MATCHES = "SCV matches"

    def __str__(self):
        return self.value

    def __repr__(self):
        return self.value


@dataclass
class ClinVarLegacyMatch:
    """
    This represents the match between a legacy ClinVar (as in records that were submitted to ClinVar prior to Shariant)
    Specifically between a LegacyMatch (which will reference many of these) with a ClinVarExport in Shariant
    This will allow us to match up the two
    """
    clinvar_export: ClinVarExport
    match_types: Set[ClinVarLegacyExportMatchType]

    @property
    def clinical_significance(self) -> str:
        if based_on := self.clinvar_export.classification_based_on:
            return EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(based_on.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))
        return ''

    @lazy
    def condition(self) -> Optional[ConditionResolved]:
        if cm := self.clinvar_export.classification_based_on:
            if cons := cm.classification.condition_resolution_obj:
                return cons
        return self.clinvar_export.condition_resolved

@dataclass
class ClinVarLegacyMatches:
    allele: Allele
    match_types: Set[ClinVarLegacyAlleleMatchType]
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
    """
    Legacy refers to this is the export from ClinVar for a lab prior to using Shariant
    We want to match up the reocrds with ones within
    """
    clinvar_key: ClinVarKey
    data: Dict[str, str]

    @lazy
    def clinvar_export(self) -> Optional[ClinVarExport]:
        return ClinVarExport.objects.filter(scv=self.scv_no_version).first()

    @staticmethod
    def from_data_str(clinvar_key: ClinVarKey, data_str: str):
        return ClinVarLegacyRow(clinvar_key=clinvar_key, data=json.loads(data_str) if data_str else {})

    @property
    def data_str(self) -> str:
        return json.dumps(self.data)

    def get_column(self, field_name: str):
        return self.data[field_name]

    @property
    def labs(self):
        return list(self.clinvar_key.lab_set.all())

    @property
    def scv(self):
        return self.get_column(ClinVarLegacyColumn.SCV)

    @property
    def scv_no_version(self):
        return self.scv.split('.')[0]

    @property
    def scv_version(self):
        return self.scv.split('.')[1]

    @property
    def variant_clinvar_id(self):
        return int(self.get_column(ClinVarLegacyColumn.VariationID))

    @property
    def variant_description_str(self) -> Optional[str]:
        return self.get_column(ClinVarLegacyColumn.Your_variant_description_HGVS)

    @property
    def variant_description(self) -> Optional[CHGVS]:
        if cd := self.variant_description_str:
            c_hgvs = CHGVS(cd)
            c_hgvs.gene_symbol = self.submitted_gene
            return c_hgvs

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
        header_indexes: Dict[str, int] = {}
        for row in csv_f:
            if not header_indexes:
                if len(row) > 1 and row[1] == 'VariationID':
                    # work out weird header row halfway down the file (and one column over)
                    # don't process the header row itself as data, now skip to the next row
                    header_indexes = {label: index for index, label in enumerate(row)}
            else:
                row_dict: Dict[str, str] = {}
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
        if not self.variant_preferred:
            return self.variant_preferred
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
    def c_hgvs_differs(self) -> bool:
        return self.c_hgvs_with_gene_symbol != self.variant_description

    @property
    def clinical_significance_code(self) -> str:
        return CLINICAL_SIGNIFICANCE_MAP.get(self.clinical_significance)

    @lazy
    def ontology_terms(self) -> List[OntologyTerm]:
        terms: List[OntologyTerm] = []
        if condition_identifier := self.condition_identifier.strip():
            if individual_terms := condition_identifier.split(';'):
                for individual_term in individual_terms:
                    # Fun with ClinVar prefixes
                    if individual_term.startswith("Human Phenotype Ontology:"):
                        individual_term = individual_term.split(":", maxsplit=1)[1]
                    elif individual_term.startswith("MONDO:MONDO"):
                        individual_term = individual_term.split(":", maxsplit=1)[1]

                    try:
                        terms.append(OntologyTerm.get_or_stub(individual_term))
                    except ValueError:
                        # might be a bunch of ontologies that we don't support
                        pass
        return terms

    def find_variant_grid_allele(self) -> List[ClinVarLegacyMatches]:

        allele_to_match_types: Dict[Allele, Set[ClinVarLegacyAlleleMatchType]] = defaultdict(set)

        if clinvar_annotation := ClinVar.objects.filter(clinvar_variation_id=self.variant_clinvar_id).order_by('-version').first():
            if allele := clinvar_annotation.variant.allele:
                allele_to_match_types[allele].add(ClinVarLegacyAlleleMatchType.CLINVAR_VARIANT_ID)

        try:
            if c_hgvs_preferred_str := self.c_hgvs_preferred_str:
                if results := search_hgvs(c_hgvs_preferred_str, user=admin_bot(), genome_build=GenomeBuild.grch38(), variant_qs=Variant.objects.all()):
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
                                allele_to_match_types[allele].add(ClinVarLegacyAlleleMatchType.VARIANT_PREFERRED_VARIANT)
        except InvalidHGVSName:
            pass
        except NotImplementedError:
            pass
        except Variant.MultipleObjectsReturned:
            # feel this really needs to be handled
            pass

        if c_hgvs := self.c_hgvs_with_gene_symbol:
            # allow for some transcript version increases
            if c_hgvs.transcript_parts.version:
                test_c_hgvs: CHGVS
                c_hgvs_strs: List[str] = []
                for attempt_increase in range(-3, 3):
                    if c_hgvs.transcript_parts.version + attempt_increase >= 1:
                        if test_c_hgvs := c_hgvs.with_transcript_version(c_hgvs.transcript_parts.version + attempt_increase):
                            c_hgvs_strs.append(str(test_c_hgvs))

                if c_hgvs_strs:
                    if allele_ids := set(Classification.objects.filter(evidence__c_hgvs__value__in=c_hgvs_strs, lab__in=self.labs).values_list('allele', flat=True)):
                        alleles = Allele.objects.filter(pk__in=allele_ids)
                        for allele in alleles:
                            allele_to_match_types[allele].add(ClinVarLegacyAlleleMatchType.VARIANT_PREFERRED_IMPORTED_C_HGVS)

        all_matches: [ClinVarLegacyMatches] = []
        for allele, match_types in allele_to_match_types.items():
            classifications = list(Classification.objects.filter(lab__in=self.labs, allele=allele, withdrawn=False))
            clinvar_export_matches: List[ClinVarLegacyMatch] = []
            if clinvar_allele := ClinVarAllele.objects.filter(allele=allele, clinvar_key=self.clinvar_key).first():
                if clinvar_exports := ClinVarExport.objects.filter(clinvar_allele=clinvar_allele):
                    for clinvar_export in clinvar_exports:

                        export_match_types: Set[ClinVarLegacyExportMatchType] = set()

                        if classification_based_on := clinvar_export.classification_based_on:
                            clinical_significance = classification_based_on.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
                            if clinical_significance == self.clinical_significance_code:
                                export_match_types.add(ClinVarLegacyExportMatchType.CLINICAL_SIGNIFICANCE_MATCHES)
                            #if CHGVS(classification_based_on.get(SpecialEKeys.C_HGVS)) == self.c_hgvs_with_gene_symbol:
                            #    export_match_types.add(ClinVarLegacyExportMatchType.C_HGVS_MATCHES)

                        umbrella = clinvar_export.condition_resolved.terms
                        for term in self.ontology_terms:
                            if mondo_term := OntologyTermRelation.as_mondo(term):
                                term = mondo_term
                            for umbrella_term in umbrella:
                                if umbrella_term.ontology_service == term.ontology_service:
                                    if term == umbrella_term or OntologySnake.check_if_ancestor(descendant=term, ancestor=umbrella_term, max_levels=2):
                                        export_match_types.add(ClinVarLegacyExportMatchType.CONDITION_MATCHES)
                                        break

                        if clinvar_export.scv == self.scv_no_version:
                            export_match_types.add(ClinVarLegacyExportMatchType.SCV_MATCHES)

                        clinvar_export_match = ClinVarLegacyMatch(
                            clinvar_export=clinvar_export,
                            match_types=export_match_types
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

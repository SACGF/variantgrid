import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from enum import StrEnum
from functools import cached_property
from typing import Iterator, Any, Optional, TypeVar, Generic
from django.utils.timezone import now
from annotation.models import ClinVarVersion, AnnotationVersion, ClinVar
from classification.enums import SpecialEKeys
from classification.models import ClinVarExport, ConditionResolved
from classification.models.clinvar_legacy import ClinVarLegacy, ClinVarLegacyAlleleMatchType
from classification.utils.clinvar_matcher import ClinVarLegacyRow
from genes.hgvs import CHGVS
from library.guardian_utils import admin_bot
from library.utils import first
from ontology.models import OntologySnake, OntologyTermRelation
from snpdb.models import ClinVarKey, GenomeBuild, Allele, Variant
from snpdb.search import SearchInput
from snpdb.signals.variant_search import search_hgvs


class ClinVarLegacyColumn(StrEnum):
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


T = TypeVar('T')


@dataclass(frozen=True)
class ClinVarLegacyMatchProperty(Generic[T]):
    value: T
    matches: bool


@dataclass(frozen=True)
class ClinVarLegacyMatch:
    clinvar_legacy: ClinVarLegacy
    clinvar_export: ClinVarExport

    @cached_property
    def scv_match(self) -> ClinVarLegacyMatchProperty[str]:
        scv = self.clinvar_export.scv
        matches = self.clinvar_export.scv == self.clinvar_legacy.scv
        return ClinVarLegacyMatchProperty(value=scv, matches=matches)

    @cached_property
    def c_hgvs_match(self) -> ClinVarLegacyMatchProperty[Optional[CHGVS]]:
        c_hgvs: Optional[CHGVS] = None
        matches = False
        if classification_based_on := self.clinvar_export.classification_based_on:
            c_hgvs = classification_based_on.classification.c_parts
            if legacy_c_hgvs := self.clinvar_legacy.your_hgvs_obj:
                if c_hgvs.transcript_parts.identifier == legacy_c_hgvs.transcript_parts.identifier and c_hgvs.c_dot == legacy_c_hgvs.c_dot:
                    matches = True
        return ClinVarLegacyMatchProperty(value=c_hgvs, matches=matches)

    @cached_property
    def condition_match(self) -> ClinVarLegacyMatchProperty[Optional[ConditionResolved]]:
        export_condition: Optional[ConditionResolved] = None
        matches = False
        if classification_based_on := self.clinvar_export.classification_based_on:
            export_condition = classification_based_on.classification.condition_resolution_obj
            if legacy_condition := self.clinvar_legacy.condition_obj:
                mondo_condition = OntologyTermRelation.as_mondo(export_condition.terms[0]) or export_condition.terms[0]
                mondo_legacy_condition = OntologyTermRelation.as_mondo(legacy_condition.terms[0]) or legacy_condition
                matches = mondo_condition == mondo_legacy_condition
        return ClinVarLegacyMatchProperty(value=export_condition, matches=matches)

    @cached_property
    def classification_match(self) -> ClinVarLegacyMatchProperty[Optional[str]]:
        classification: Optional[str] = None
        matches = False
        if classification_based_on := self.clinvar_export.classification_based_on:
            classification = classification_based_on.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            legacy_classification = self.clinvar_legacy.classification_code
            if classification == legacy_classification or classification.startswith("VUS") and legacy_classification.startswith("VUS"):
                matches = True
        return ClinVarLegacyMatchProperty(value=classification, matches=matches)

    @cached_property
    def allele_match(self) -> ClinVarLegacyMatchProperty[Allele]:
        export_allele = self.clinvar_export.clinvar_allele.allele
        matches = export_allele == self.clinvar_legacy.allele
        return ClinVarLegacyMatchProperty(
            value=export_allele,
            matches=matches
        )

    @cached_property
    def clinvar_key_match(self) -> ClinVarLegacyMatchProperty[ClinVarKey]:
        clinvar_key = self.clinvar_export.clinvar_allele.clinvar_key
        return ClinVarLegacyMatchProperty(
            value=clinvar_key,
            matches=clinvar_key == self.clinvar_legacy.clinvar_key
        )



class ClinVarLegacyService:

    def __init__(self, clinvar_key: ClinVarKey):
        self.clinvar_key = clinvar_key

    @cached_property
    def existing_legacy_records(self) -> dict[str, ClinVarLegacy]:
        legacy_dict: dict[str, ClinVarLegacy] = {}
        for existing in ClinVarLegacy.objects.filter(clinvar_key=self.clinvar_key).iterator():
            legacy_dict[existing.scv] = existing
        return legacy_dict

    def _process_row(self, row: dict[str, Any]):
        scv = row.get(ClinVarLegacyColumn.SCV)
        # don't worry about SCV version numbers, they don't affect matching
        if "." in scv:
            scv = scv.split(".")[0]

        update_me: ClinVarLegacy
        if existing := self.existing_legacy_records.get(row.get('SCV')):
            update_me = existing
        else:
            update_me = ClinVarLegacy(scv=scv, clinvar_key=self.clinvar_key)
            self.existing_legacy_records[scv] = update_me

        # TODO - really should check if thse values change
        # and if they do, mark the record as dirty
        # and if they never change, move this logic to only when we get a new row
        update_me.clinvar_variation_id = row.get(ClinVarLegacyColumn.VariationID)
        update_me.clinvar_allele_id = row.get(ClinVarLegacyColumn.AlleleID)
        update_me.your_variant_description_hgvs = row.get(ClinVarLegacyColumn.Your_variant_description_HGVS)
        update_me.your_variant_description_chromosome_coordinates = row.get(ClinVarLegacyColumn.Your_variant_description_chromosome_coordinates)
        update_me.preferred_variant_name = row.get(ClinVarLegacyColumn.Preferred_variant_name)
        update_me.your_condition_name = row.get(ClinVarLegacyColumn.Your_condition_name)
        update_me.your_condition_identifier = row.get(ClinVarLegacyColumn.Your_condition_identifier)
        update_me.clinical_significance = row.get(ClinVarLegacyColumn.Clinical_significance)

        date_last_evaluated = None
        if raw_date := row.get(ClinVarLegacyColumn.Date_last_evaluated):
            date_last_evaluated = datetime.strptime(raw_date, "%b %d, %Y")

        update_me.date_last_evaluated = date_last_evaluated
        update_me.assertion_criteria = row.get(ClinVarLegacyColumn.Assertion_criteria)
        update_me.submitted_gene = row.get(ClinVarLegacyColumn.Submitted_gene)

    def load_file(self, file, delimiter='\t') -> Iterator['ClinVarLegacyRow']:
        csv_f = csv.reader(file, delimiter=delimiter)
        header_indexes: dict[str, int] = {}
        # Clinvar cummalative summary file
        # it has a list of rows describing the columns (all these descriptors appear in the first column)
        # then finally the actual columns (the first column is then left blank - only used for the descriptors above)
        for row in csv_f:
            if not header_indexes:
                # test to see if this row is the actual header
                if len(row) > 1 and row[1] == 'VariationID':
                    # work out weird header row halfway down the file (and one column over)
                    # don't process the header row itself as data, now skip to the next row
                    header_indexes = {label: index for index, label in enumerate(row)}
            else:
                # we now have a header
                row_dict: dict[str, str] = {}
                for header_key, index in header_indexes.items():
                    row_dict[header_key] = row[index]

                self._process_row(row_dict)

        if not header_indexes:
            raise ValueError("Did not find the header row")

        ClinVarLegacy.objects.bulk_create(
            objs=self.existing_legacy_records.values(),
            update_conflicts=True,
            unique_fields=['scv'],
            update_fields=[
                'clinvar_variation_id',
                'clinvar_allele_id',
                'your_variant_description_hgvs',
                'your_variant_description_chromosome_coordinates',
                'preferred_variant_name',
                'your_condition_name',
                'your_condition_identifier',
                'clinical_significance',
                'date_last_evaluated',
                'assertion_criteria',
                'submitted_gene'
            ]
        )

    @cached_property
    def clinvar_version(self) -> ClinVarVersion:
        # is it better just to select latest ClinVarVersion
        return AnnotationVersion.latest(GenomeBuild.grch38()).clinvar_version

    @staticmethod
    def find_matches(clinvar_legacy: ClinVarLegacy) -> list[ClinVarExport]:
        query = ClinVarExport.objects.filter(scv=clinvar_legacy.scv)
        if allele := clinvar_legacy.allele:
            query = query.union(ClinVarExport.objects.filter(
                clinvar_allele__clinvar_key=clinvar_legacy.clinvar_key,
                clinvar_allele__allele=allele
            ))
        return [ClinVarLegacyMatch(clinvar_legacy=clinvar_legacy, clinvar_export=clinvar_export) for clinvar_export in query]

    def match_allele(self, clinvar_legacy: ClinVarLegacy) -> None:
        clinvar_legacy.allele = None
        clinvar_legacy.allele_match_type = ClinVarLegacyAlleleMatchType.MATCHED_NOT_FOUND
        clinvar_legacy.allele_match_date = now()

        if clinvar_annotation := ClinVar.objects.filter(
                clinvar_variation_id=clinvar_legacy.clinvar_variation_id,
                version=self.clinvar_version
        ).first():
            if variant := clinvar_annotation.variant:
                if allele := variant.allele:
                    clinvar_legacy.allele = allele
                    clinvar_legacy.allele_match_type = ClinVarLegacyAlleleMatchType.MATCHED_ON_CLINVAR_ALLELE

        if not clinvar_legacy.allele:
            # use search, this is v slow
            if c_hgvs_preferred_str := clinvar_legacy.hgvs_obj.full_c_hgvs:
                search_input = SearchInput(user=admin_bot(), search_string=c_hgvs_preferred_str,
                                           genome_build_preferred=GenomeBuild.grch38())
                try:
                    response = search_hgvs(sender=None, search_input=search_input)
                    # note, the method signature is correct, the annotation on search_hgvs makes it take a search_input not a search_input_instance
                    matched_alleles = set()
                    for result_entry in response.results:
                        result = result_entry.preview.obj
                        allele: Optional[Allele] = None
                        if isinstance(result, Variant):
                            allele = result.allele
                        elif isinstance(result, Allele):
                            allele = result
                        if allele:
                            matched_alleles.add(allele)

                    if len(matched_alleles) == 1:
                        clinvar_legacy.allele = first(matched_alleles)
                        clinvar_legacy.allele_match_type = ClinVarLegacyAlleleMatchType.MATCHED_ON_HGVS
                except ValueError:
                    pass
        clinvar_legacy.save()
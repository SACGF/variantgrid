import csv
from dataclasses import dataclass
from datetime import datetime
from enum import StrEnum
from functools import cached_property
from typing import Iterator, Any, Optional, TypeVar, Generic

from celery import Task
from django.db.models import QuerySet
from django.utils.timezone import now
from annotation.models import ClinVarVersion, AnnotationVersion, ClinVar
from classification.enums import SpecialEKeys
from classification.models import ClinVarExport, ConditionResolved
from classification.models.clinvar_legacy import ClinVarLegacy, ClinVarLegacyAlleleMatchType
from genes.hgvs import CHGVS
from library.guardian_utils import admin_bot
from library.utils import first
from ontology.models import OntologyTermRelation
from snpdb.models import ClinVarKey, GenomeBuild, Allele, Variant, ClinVarExportTypeBucket
from snpdb.search import SearchInput
from snpdb.signals.variant_search import search_hgvs
from variantgrid.celery import app


class ClinVarLegacyColumn(StrEnum):
    """
    Lists the columns as provided by ClinVar for parsing
    Note that with the introduction of Somatic/Clinical Impact records, the column
    Clinical_significance has been renamed to Germline_classification with no type indicator or column for Tier/Somatic value etc
    """
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
    Germline_classification = "Germline_classification"
    Date_last_evaluated = "Date_last_evaluated"
    Assertion_criteria = "Assertion_criteria"
    Submitted_gene = "Submitted_gene"


T = TypeVar('T')


@dataclass(frozen=True)
class ClinVarLegacyMatchProperty(Generic[T]):
    """
    When matching between a ClinVarLegacy and ClinVarExport record, keep the value of the ClinVarExport
    and provide a boolean indicating if the value matched the Legacy (sometimes a little more complex than just == )
    """
    value: T
    matches: bool


@dataclass(frozen=True)
class ClinVarLegacyMatch:
    """
    Represents a possible match between ClinVarLegacy and ClinVarExport.
    Might be matched based on allele & export type or SCV
    """

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


class ClinVarLegacyImporter:

    def __init__(self, clinvar_key: ClinVarKey):
        self.clinvar_key = clinvar_key

    @cached_property
    def existing_legacy_records(self) -> dict[str, ClinVarLegacy]:
        """
        Get all ClinVarlegacy records for the clinvar_key, in a dict keyed by SCV
        """
        legacy_dict: dict[str, ClinVarLegacy] = {}
        for existing in ClinVarLegacy.objects.filter(clinvar_key=self.clinvar_key).iterator():
            legacy_dict[existing.scv] = existing
        return legacy_dict

    @staticmethod
    def clean_condition_identifier(identifier: str) -> str:
        """
        ClinVar has conditions a big screwed up (also need to check for HPO as I believe that's expanded for some reason)
        :param identifier: The condition identifier as provided by ClinVar
        :return: The standard condition identifier as can be parsed by VariantGrid
        """
        if identifier:
            if identifier.startswith("MONDO:MONDO:"):
                identifier = identifier.replace("MONDO:MONDO:", "MONDO:")
            elif identifier.isnumeric():
                identifier = f"OMIM:{identifier}"
        return identifier

    def _process_row(self, row: dict[str, Any]) -> ClinVarLegacy:
        """
        Pass in a row of the ClinVarLegacy import, and update an existing ClinVarLegacy record (matching on SCV) or create
        a new one.
        Updates the dictionary returned by existing_legacy_records, so need to bulk_update the contents of that dict after
        """

        scv = row.get(ClinVarLegacyColumn.SCV)
        # remove version number from SCV, so we can match on the number (version number just increments with each
        # submission and doesn't alter variant matching)
        if "." in scv:
            scv = scv.split(".")[0]

        update_me: ClinVarLegacy
        if existing := self.existing_legacy_records.get(row.get('SCV')):
            update_me = existing
        else:
            update_me = ClinVarLegacy(scv=scv, clinvar_key=self.clinvar_key, clinvar_bucket=ClinVarExportTypeBucket.GERMLINE)
            self.existing_legacy_records[scv] = update_me

        # and if they do, mark the record as dirty
        # and if they never change, move this logic to only when we get a new row
        update_me.clinvar_variation_id = int(row.get(ClinVarLegacyColumn.VariationID))
        update_me.clinvar_allele_id = int(row.get(ClinVarLegacyColumn.AlleleID))
        update_me.your_variant_description_hgvs = row.get(ClinVarLegacyColumn.Your_variant_description_HGVS)
        update_me.your_variant_description_chromosome_coordinates = row.get(ClinVarLegacyColumn.Your_variant_description_chromosome_coordinates)
        update_me.preferred_variant_name = row.get(ClinVarLegacyColumn.Preferred_variant_name)
        update_me.your_condition_name = row.get(ClinVarLegacyColumn.Your_condition_name)
        update_me.your_condition_identifier = ClinVarLegacyImporter.clean_condition_identifier(row.get(ClinVarLegacyColumn.Your_condition_identifier))
        update_me.germline_classification = row.get(ClinVarLegacyColumn.Germline_classification) or row.get(ClinVarLegacyColumn.Clinical_significance)

        date_last_evaluated = None
        if raw_date := row.get(ClinVarLegacyColumn.Date_last_evaluated):
            date_last_evaluated = datetime.strptime(raw_date, "%b %d, %Y")

        update_me.date_last_evaluated = date_last_evaluated
        update_me.assertion_criteria = row.get(ClinVarLegacyColumn.Assertion_criteria)
        update_me.submitted_gene = row.get(ClinVarLegacyColumn.Submitted_gene)
        return update_me

    def load_file(self, file, delimiter='\t') -> None:
        """
        ClinVar cumulative summary file
        it has a list of rows describing the columns (all these descriptors appear in the first column).
        then finally the actual columns (the first column is then left blank - only used for the descriptors above)
        :param file: File object that can be put into a CSV reader
        :param delimiter: Delimiter for the file, as the file seems to come in tsv and csv
        :return:
        """
        csv_f = csv.reader(file, delimiter=delimiter)
        header_indexes: dict[str, int] = {}

        updated_rows: list[ClinVarLegacy] = []

        for row in csv_f:
            if not header_indexes:
                # test to see if this row is the actual header
                if len(row) > 1 and row[1] == 'VariationID':
                    # work out the unusual header row halfway down the file (and one column over)
                    # don't process the header row itself as data, now skip to the next row
                    header_indexes = {label: index for index, label in enumerate(row)}
            else:
                # we now have a header
                row_dict: dict[str, str] = {}
                for header_key, index in header_indexes.items():
                    row_dict[header_key] = row[index]

                updated_rows.append(self._process_row(row_dict))

        if not header_indexes:
            raise ValueError("Did not find the header row")

        # note we don't delete rows that weren't present in the file, htat has to be done manually for now

        ClinVarLegacy.objects.bulk_create(
            objs=updated_rows,
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
                'germline_classification',
                'date_last_evaluated',
                'assertion_criteria',
                'submitted_gene'
            ]
        )


@staticmethod
def clinvar_legacy_find_matches(clinvar_legacy: ClinVarLegacy) -> list[ClinVarLegacyMatch]:
    """
    Returns ClinVarLegacyMatch objects for the given ClinVarLegacy
    Matches on SCV unioned with Allele and ClinVar Export type (assuming the clinvar_legacy has matched to an allele)
    """
    query = ClinVarExport.objects.filter(scv=clinvar_legacy.scv)
    if allele := clinvar_legacy.allele:
        allele_matches_qs = ClinVarExport.objects.filter(
            clinvar_allele__clinvar_key=clinvar_legacy.clinvar_key,
            clinvar_allele__allele=allele,
            clinvar_allele__clinvar_export_bucket=clinvar_legacy.clinvar_bucket
        )
        query = query.union(allele_matches_qs)
    return [ClinVarLegacyMatch(clinvar_legacy=clinvar_legacy, clinvar_export=clinvar_export) for clinvar_export in
            query]


def clinvar_legacy_match_alleles_for_clinvar_key(clinvar_key: str):
    clinvar_key = ClinVarKey.objects.get(pk=clinvar_key)
    return clinvar_legacy_match_alleles_for_query(clinvar_key.clinvarlegacy_set.all())


def clinvar_legacy_match_alleles_for_query(qs: QuerySet[ClinVarLegacy]):
    """
    Matches all ClinVarLegacy records to alleles (as best we can).
    Is a long running so best done as a task
    :param qs: String key of a ClinVarKey
    """

    variation_id_to_legacy: dict[int, ClinVarLegacy] = {}
    for legacy in qs.iterator():
        variation_id_to_legacy[legacy.clinvar_variation_id] = legacy

    match_date = now()

    latest_clinvar = AnnotationVersion.latest(GenomeBuild.grch38()).clinvar_version
    for cv in ClinVar.objects.filter(
        clinvar_variation_id__in=variation_id_to_legacy.keys(),
        version=latest_clinvar
    ).select_related('variant'):
        legacy = variation_id_to_legacy[cv.clinvar_variation_id]
        legacy.allele = cv.variant.allele
        legacy.allele_match_type = ClinVarLegacyAlleleMatchType.MATCHED_ON_CLINVAR_ALLELE
        legacy.allele_match_date = match_date

    for cv in variation_id_to_legacy.values():
        if cv.allele_match_date != match_date:
            found_allele_in_search = False
            if c_hgvs_preferred_str := cv.hgvs_obj.full_c_hgvs:
                search_input = SearchInput(user=admin_bot(), search_string=c_hgvs_preferred_str,
                                           genome_build_preferred=GenomeBuild.grch38())
                try:
                    response = search_hgvs(sender=None, search_input=search_input)
                    # note, the method signature is correct, the annotation on search_hgvs makes it take a search_input not a search_input_instance
                    matched_alleles: set[Allele] = set()
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
                        cv.allele = first(matched_alleles)
                        cv.allele_match_type = ClinVarLegacyAlleleMatchType.MATCHED_ON_HGVS
                        cv.allele_match_date = match_date
                        found_allele_in_search = True
                except:
                    pass

            if not found_allele_in_search:
                cv.allele = None
                cv.allele_match_type = ClinVarLegacyAlleleMatchType.MATCHED_NOT_FOUND
                cv.allele_match_date = match_date

    ClinVarLegacy.objects.bulk_update(variation_id_to_legacy.values(), fields=["allele", "allele_match_type", "allele_match_date"])


class ClinVarLegacyAlleleMatchTask(Task):

    def run(self, clinvar_key_id: str):
        clinvar_legacy_match_alleles_for_clinvar_key(clinvar_key_id)


ClinVarLegacyAlleleMatchTask = app.register_task(ClinVarLegacyAlleleMatchTask())
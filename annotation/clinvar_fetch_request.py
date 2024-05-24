import logging
import time
from dataclasses import dataclass, field
from datetime import timedelta
from typing import Type, Optional, Collection
from urllib.error import HTTPError

from django.db import transaction
from django.utils import timezone

from annotation.clinvar_xml_parser import CLINVAR_RECORD_CACHE_DAYS, ClinVarXmlParser
from annotation.clinvar_xml_parser_via_vcv import ClinVarXmlParserViaVCV
from annotation.models import ClinVarRecordCollection, ClinVarVersion, ClinVar
from snpdb.models import GenomeBuild, Allele


@dataclass
class ClinVarFetchRequest:
    """
    Is used to retrieve individual ClinVar records from the ClinVar web service for a given clinvar_variation_id
    """

    clinvar_variation_id: int
    """
    The API currently only allows us to ask for all records for that clinvar_variation_id and then we can cache
    only the ones we want. if we've already cached records (but with a lower min_stars, we can just re-use that)
    """

    max_cache_age: timedelta = field(default=timedelta(days=CLINVAR_RECORD_CACHE_DAYS))
    """
    How old until the cache is considered stale, provide seconds=0 if you want to force a refresh
    """

    parser: Type[ClinVarXmlParser] = ClinVarXmlParserViaVCV

    clinvar_versions: Optional[Collection[ClinVarVersion]] = None

    def fetch(self) -> ClinVarRecordCollection:
        fetch_date = timezone.now()

        clinvar_versions = self.clinvar_versions
        if not clinvar_versions:
            clinvar_versions = []
            for genome_build in GenomeBuild.builds_with_annotation():
                if clinvar_version := ClinVarVersion.objects.filter(genome_build=genome_build).order_by('-created').first():
                    clinvar_versions.append(clinvar_version)

        with transaction.atomic():
            clinvar_record_collection, created = ClinVarRecordCollection.objects.get_or_create(clinvar_variation_id=self.clinvar_variation_id)
            # the select_for_update() stops two simultaneous requests for updating the same clinvar_record_collection
            clinvar_record_collection = ClinVarRecordCollection.objects.select_for_update().get(pk=clinvar_record_collection.pk)
            fetch_from_clinvar = True
            if not created:
                if \
                        (clinvar_record_collection.parser_version == self.parser.PARSER_VERSION) and \
                        clinvar_record_collection.last_loaded and \
                        ((fetch_date - clinvar_record_collection.last_loaded) <= self.max_cache_age):
                    # if all the above is true, then our cache is fine
                    fetch_from_clinvar = False

            allele_id = clinvar_record_collection.allele_id
            if not allele_id:
                clinvar_record: ClinVar
                if clinvar_record := ClinVar.objects.filter(clinvar_variation_id=self.clinvar_variation_id).first():
                    if variant := clinvar_record.variant:
                        if allele := variant.allele:
                            allele_id = allele.pk
                clinvar_record_collection.allele_id = allele_id

            if not allele_id:
                raise ValueError(f"Couldn't determine Allele for clinvar_variation_id {self.clinvar_variation_id}")

            if fetch_from_clinvar:
                # so while Entrez does automatically retry on 500s, ClinVar has been providing 400s (Bad Request) when
                # the request is fine
                attempt_count = 2
                while True:
                    # loop is broken out of if it works, or raise if it fails after attempt_count
                    try:
                        response = self.parser.load_from_clinvar_id(clinvar_variation_id=self.clinvar_variation_id)

                        # update our cache
                        clinvar_record_collection.last_loaded = fetch_date
                        clinvar_record_collection.urls = response.urls
                        clinvar_record_collection.parser_version = self.parser.PARSER_VERSION

                        clinvar_record_collection.update_with_records_and_save(response.all_records)
                        break

                    except HTTPError as http_ex:
                        if http_ex.code == 400:
                            attempt_count -= 1
                            if attempt_count > 0:
                                logging.warning("400 from Entrez when fetching ClinVarRecord, waiting then trying again")
                                time.sleep(10)
                                continue
                        # out of attempts or not 400
                        raise

        return clinvar_record_collection

import re
from typing import List, Optional

from django.core.exceptions import ObjectDoesNotExist
from django.db import models
from django_extensions.db.models import TimeStampedModel
from requests import request

from snpdb.models import GenomeBuild

DBSNP_PATTERN = re.compile(r"^rs(\d+)$")


class DbSNP(TimeStampedModel):
    """ Wrapper around DbSNP API service, @see https://api.ncbi.nlm.nih.gov/variation/v0/ """
    # id = "rs" stripped off, eg rs121964969 => 121964969
    id = models.IntegerField(primary_key=True)
    api_response = models.JSONField(null=False)  # returned

    @staticmethod
    def get_id_from_code(dbsnp_rs_id: str) -> int:
        m = DBSNP_PATTERN.match(dbsnp_rs_id)
        if not m:
            raise ValueError(f"{dbsnp_rs_id} didn't match dbSNP rs pattern")
        return int(m.group(1))

    @staticmethod
    def get(dbsnp_rs_id: str) -> 'DbSNP':
        dbsnp_id = DbSNP.get_id_from_code(dbsnp_rs_id)
        try:
            dbsnp = DbSNP.objects.get(pk=dbsnp_id)
        except DbSNP.DoesNotExist:
            url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{dbsnp_id}"
            r = request(url=url, method='get')
            api_response = r.json()
            dbsnp = DbSNP.objects.create(pk=dbsnp_id, api_response=api_response)
        return dbsnp

    @staticmethod
    def get_for_variant(variant: 'Variant', variant_annotation_version) -> 'DbSNP':
        dbsnp = None
        try:
            variant_annotation = variant.variantannotation_set.get(version=variant_annotation_version)
            if variant_annotation.dbsnp_rs_id:
                dbsnp = DbSNP.get(variant_annotation.dbsnp_rs_id)
        except ObjectDoesNotExist:
            pass
        return dbsnp

    def get_alleles_for_genome_build(self, genome_build: GenomeBuild) -> List:
        placements_with_allele = self.api_response["primary_snapshot_data"]["placements_with_allele"]
        build_with_patch = genome_build.get_build_with_patch()

        for pwa in placements_with_allele:
            pa = pwa["placement_annot"]
            if pa["seq_type"] == "refseq_chromosome":
                for sit in pa["seq_id_traits_by_assembly"]:
                    if sit.get("assembly_name") == build_with_patch:
                        return pwa["alleles"]
        return []

    def get_g_hgvs(self, genome_build: GenomeBuild, alt=None) -> Optional[str]:
        g_hgvs = None
        alleles = self.get_alleles_for_genome_build(genome_build)
        if alleles:
            num_alleles = len(alleles)
            if num_alleles > 1:
                if alt is None:
                    msg = f"{self} had {num_alleles} alleles for {genome_build} - need alt to distinguish"
                    raise ValueError(msg)
                alleles_with_alt = []
                for allele in alleles:
                    spdi = allele["allele"]["spdi"]
                    if alt == spdi["inserted_sequence"]:
                        alleles_with_alt.append(allele["hgvs"])

                num_alleles = len(alleles_with_alt)
                if num_alleles:
                    if num_alleles == 1:
                        g_hgvs = alleles_with_alt[0]
                    else:
                        msg = f"{self} had {num_alleles} alleles for {genome_build} w/alt={alt}"
                        raise ValueError(msg)
                else:
                    raise ValueError(f"{self}/{genome_build} has no alleles w/alt={alt}")
        return g_hgvs

    def __str__(self):
        return f"rs{self.pk}"

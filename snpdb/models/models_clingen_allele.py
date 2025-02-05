import logging
import re
from functools import cached_property, lru_cache
from typing import Optional

from django.conf import settings
from django.db import models
from django_extensions.db.models import TimeStampedModel

from snpdb.models.models_enums import SequenceRole
from snpdb.models.models_genome import GenomeBuild, Contig


class ClinGenAllele(TimeStampedModel):
    """ Canonical Allele - @see http://reg.clinicalgenome.org """
    # id = "CA" stripped off, e.g. CA7169043 => 7169043
    id = models.BigIntegerField(primary_key=True)
    api_response = models.JSONField(null=False)  # returned

    CLINGEN_ALLELE_SERVER_ERROR_TYPE = "ServerError"  # Set fake API response of errorType to indicate server error
    CLINGEN_ALLELE_URL_PATTERN = re.compile(r"http.*/allele/CA(\d+)$")
    CLINGEN_ALLELE_CODE_PATTERN = re.compile(r"^CA(\d+)$", re.IGNORECASE)
    CLINGEN_ALLELE_MAX_ALLELE_SIZE = 10_000

    class ClinGenAlleleRegistryException(ValueError):
        """ Base exception """

    class ClinGenBuildNotInResponseError(ClinGenAlleleRegistryException):
        pass

    class ClinGenNonChromosomeLiftoverError(ClinGenAlleleRegistryException):
        pass

    class ClinGenMissingAlleleID(ClinGenAlleleRegistryException):
        """ Coordinate is not yet assigned ID and stored on server """

    class ClinGenHGVSReferenceBaseUnavailableError(ClinGenAlleleRegistryException):
        """ HGVS doesn't have reference """

    @staticmethod
    @lru_cache
    def refseq_valid_and_invalid_contigs(genome_build) -> tuple[set,set]:
        invalid_contigs = Contig.objects.none()
        valid_contigs = genome_build.contigs
        if settings.LIFTOVER_TO_CHROMOSOMES_ONLY:
            valid_contigs = valid_contigs.filter(role=SequenceRole.ASSEMBLED_MOLECULE)
            invalid_contigs = genome_build.contigs.exclude(role=SequenceRole.ASSEMBLED_MOLECULE)

        refseq_contigs = set(valid_contigs.values_list("refseq_accession", flat=True))
        refseq_invalid_contigs = set(invalid_contigs.values_list("refseq_accession", flat=True))

        return refseq_contigs, refseq_invalid_contigs

    @staticmethod
    def _strip_transcript_version(transcript_id):
        """ strip dot version ie NM_198798.2 => NM_198798 """
        return transcript_id.split(".")[0]

    @staticmethod
    def get_id_from_response(api_response):
        # Looks like "@id": "http://reg.test.genome.network/allele/CA7169043",
        url_id = api_response["@id"]
        if m := ClinGenAllele.CLINGEN_ALLELE_URL_PATTERN.match(url_id):
            return int(m.group(1))
        # logging.error("Bad @id in ClinGen response:\n%s", api_response)
        raise ClinGenAllele.ClinGenMissingAlleleID(f"Couldn't retrieve ClinGen AlleleID from @id '{url_id}'")

    @staticmethod
    def get_id_from_code(code):
        """ returns -1 on error """
        if m := ClinGenAllele.CLINGEN_ALLELE_CODE_PATTERN.match(code):
            return int(m.group(1))
        return -1

    @staticmethod
    def looks_like_id(code):
        return ClinGenAllele.get_id_from_code(code) >= 0

    @cached_property
    def transcript_alleles_by_transcript_accession(self) -> dict[str, dict]:
        ta_by_tv = {}
        if transcript_alleles := self.api_response.get("transcriptAlleles"):
            for ta in transcript_alleles:
                for t_hgvs in ta["hgvs"]:
                    transcript_accession = t_hgvs.split(":")[0]
                    ta_by_tv[transcript_accession] = ta
        return ta_by_tv

    @cached_property
    def transcript_alleles_by_transcript(self) -> dict[str, dict]:
        """ Version stripped off """
        ta_by_t = {}
        for transcript_accession, ta in self.transcript_alleles_by_transcript_accession.items():
            transcript_id = ClinGenAllele._strip_transcript_version(transcript_accession)
            ta_by_t[transcript_id] = ta
        return ta_by_t

    def _get_transcript_allele(self, transcript_accession, match_version=True) -> Optional[dict]:
        if match_version:
            ta = self.transcript_alleles_by_transcript_accession.get(transcript_accession)
        else:
            transcript_id = ClinGenAllele._strip_transcript_version(transcript_accession)
            ta = self.transcript_alleles_by_transcript.get(transcript_id)
        return ta

    def get_p_hgvs(self, transcript_accession, match_version=True):
        if ta := self._get_transcript_allele(transcript_accession, match_version):
            if protein_effect := ta.get("proteinEffect"):
                if p_hgvs := protein_effect.get("hgvs"):
                    return p_hgvs
        return None

    def get_c_hgvs_variant(self, hgvs_converter, transcript_accession):
        """ c.HGVS has reference bases on it """
        from genes.models import TranscriptVersionSequenceInfo

        hgvs_variant = None
        raw_hgvs_string, t_data = self._get_raw_hgvs_and_data(transcript_accession)
        if raw_hgvs_string:  # Has for this transcript version
            hgvs_variant = hgvs_converter.create_hgvs_variant(raw_hgvs_string)
            # Sometimes ClinGen return "n." on NM transcripts - reported as a bug 22/9/21
            if hgvs_variant.kind == "n":
                if transcript_accession.startswith("NM_") or "proteinEffect" in t_data:
                    hgvs_variant.kind = 'c'

            if not hgvs_variant.gene:  # Ref/Ens HGVSs have transcript no gene, LRG is set as gene
                hgvs_variant.gene = t_data.get("geneSymbol")

            if hgvs_variant.mutation_type in {"dup", "del", "delins"}:
                # We want to add reference bases onto HGVS but ClinGen reference sequence is wrong (see issue #493)
                coord = t_data["coordinates"][0]
                if "startIntronOffset" in coord:
                    # We have to accept these - as we get so many but we can't add on the bases
                    msg = "A coding DNA reference sequence does not contain intron or 5' and 3' gene "\
                          "flanking sequences and can therefore not be used as a reference to describe " \
                          "variants in these regions"
                    if settings.CLINGEN_ALLELE_REGISTRY_REQUIRE_REF_ALLELE:
                        raise ClinGenAllele.ClinGenHGVSReferenceBaseUnavailableError(msg)
                    else:
                        logging.warning(msg)
                elif transcript_accession.startswith("LRG_"):
                    logging.warning("Don't have sequence for LRGs, relying on ClinGenAlleleRegistry reference "
                                    "bases which may be wrong")
                else:
                    from genes.models import NoTranscript

                    try:
                        if tvsi := TranscriptVersionSequenceInfo.get(transcript_accession):
                            if hgvs_variant.mutation_type == "dup":
                                ref_end = coord["end"]
                                ref_start = ref_end - len(coord["allele"])
                            else:
                                ref_start = coord["start"]
                                ref_end = coord["end"]
                            hgvs_variant.ref_allele = tvsi.sequence[ref_start:ref_end]
                        else:
                            msg = f"Could not retrieve reference allele for '{transcript_accession}'"
                            if settings.CLINGEN_ALLELE_REGISTRY_REQUIRE_REF_ALLELE:
                                raise ClinGenAllele.ClinGenHGVSReferenceBaseUnavailableError(msg)
                            else:
                                logging.warning(msg)
                    except NoTranscript as e_no_transcript:
                        if settings.CLINGEN_ALLELE_REGISTRY_REQUIRE_REF_ALLELE:
                            raise ClinGenAllele.ClinGenHGVSReferenceBaseUnavailableError() from e_no_transcript
                        else:
                            logging.warning(e_no_transcript)

        return hgvs_variant

    def _get_raw_hgvs_and_data(self, transcript_accession, match_version=True) -> tuple[Optional[str],
                                                                                        Optional[dict]]:
        if ta := self._get_transcript_allele(transcript_accession, match_version):
            transcript_id = ClinGenAllele._strip_transcript_version(transcript_accession)
            for t_hgvs in ta["hgvs"]:
                hgvs_transcript_accession = t_hgvs.split(":")[0]
                if match_version:
                    if transcript_accession == hgvs_transcript_accession:
                        return t_hgvs, ta
                else:
                    hgvs_transcript_id = ClinGenAllele._strip_transcript_version(hgvs_transcript_accession)
                    if transcript_id == hgvs_transcript_id:
                        return t_hgvs, ta
        return None, None

    def get_g_hgvs(self, genome_build: GenomeBuild):
        refseq_contigs, refseq_invalid_contigs = self.refseq_valid_and_invalid_contigs(genome_build)
        genomic_alleles = self.api_response["genomicAlleles"]

        invalid_contigs = set()
        unknown_contigs = set()
        # Work on the contig level as AlleleRegistry only has MT for "GRCh38" (even though contigs are the same)
        for ga in genomic_alleles:
            for h in ga["hgvs"]:
                contig = h.split(":", 1)[0]
                if contig in refseq_contigs:
                    return h

                # Only report bad ones for this genome build
                if ga.get("referenceGenome") == genome_build.name:
                    if contig in refseq_invalid_contigs:
                        invalid_contigs.add(contig)
                    else:
                        unknown_contigs.add(contig)

        if invalid_contigs:
            ic_msg = "settings.LIFTOVER_TO_CHROMOSOMES_ONLY=True disabled liftover to non-chrom contigs"
            msg = f"{self}/{genome_build} - {ic_msg}: {','.join(invalid_contigs)}"
            raise Contig.ContigNotInBuildError(msg)

        if unknown_contigs:
            msg = f"{self}/{genome_build} - unknown contigs: {','.join(unknown_contigs)}"
            raise Contig.ContigNotInBuildError(msg)

        msg = f"{self}/{genome_build} not in ClinGenAllele genomicAlleles response"
        raise ClinGenAllele.ClinGenBuildNotInResponseError(msg)

    def get_variant_coordinate(self, genome_build: GenomeBuild) -> 'VariantCoordinate':
        from genes.hgvs import get_hgvs_variant_coordinate
        g_hgvs = self.get_g_hgvs(genome_build)
        return get_hgvs_variant_coordinate(g_hgvs, genome_build)

    def get_variant_string(self, genome_build: GenomeBuild, abbreviate=False):
        from snpdb.models import Variant
        return Variant.format_tuple(*self.get_variant_coordinate(genome_build), abbreviate=abbreviate)

    @property
    def human_url(self):
        params = (settings.CLINGEN_ALLELE_REGISTRY_DOMAIN, self)
        return "%s/redmine/projects/registry/genboree_registry/by_caid?caid=%s" % params

    @staticmethod
    def format_clingen_allele(pk):
        if pk:
            return f"CA{pk:06}"
        return "Unregistered Allele"

    def __str__(self):
        """ ClinGen has minimum length of 6 (0 padded) but can go longer """
        return self.format_clingen_allele(self.pk)

import logging
import re
from typing import Dict, Optional, Tuple

from django.conf import settings
from django.db import models
from django_extensions.db.models import TimeStampedModel
from lazy import lazy
from pyhgvs import HGVSName

from snpdb.models.models_enums import SequenceRole
from snpdb.models.models_genome import GenomeBuild, Contig


class ClinGenAllele(TimeStampedModel):
    """ Canonical Allele - @see http://reg.clinicalgenome.org """
    # id = "CA" stripped off, eg CA7169043 => 7169043
    id = models.BigIntegerField(primary_key=True)
    api_response = models.JSONField(null=False)  # returned

    CLINGEN_ALLELE_SERVER_ERROR_TYPE = "ServerError"  # Set fake API response of errorType to indicate server error
    CLINGEN_ALLELE_URL_PATTERN = re.compile(r"http.*/allele/CA(\d+)$")
    CLINGEN_ALLELE_CODE_PATTERN = re.compile(r"^CA(\d+)")
    CLINGEN_ALLELE_MAX_REPRESENTATION_SIZE = 10000

    class ClinGenBuildNotInResponseError(ValueError):
        pass

    class ClinGenNonChromosomeLiftoverError(ValueError):
        pass

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
        msg = f"Couldn't retrieve ClinGen AlleleID from @id '{url_id}'"
        raise ValueError(msg)

    @staticmethod
    def get_id_from_code(code):
        """ returns -1 on error """
        if m := ClinGenAllele.CLINGEN_ALLELE_CODE_PATTERN.match(code):
            return int(m.group(1))
        return -1

    @staticmethod
    def looks_like_id(code):
        return ClinGenAllele.get_id_from_code(code) >= 0

    @lazy
    def transcript_alleles_by_transcript_accession(self) -> Dict[str, Dict]:
        ta_by_tv = {}
        if transcript_alleles := self.api_response.get("transcriptAlleles"):
            for ta in transcript_alleles:
                for t_hgvs in ta["hgvs"]:
                    transcript_accession = t_hgvs.split(":")[0]
                    ta_by_tv[transcript_accession] = ta
        return ta_by_tv

    @lazy
    def transcript_alleles_by_transcript(self) -> Dict[str, Dict]:
        """ Version stripped off """
        ta_by_t = {}
        for transcript_accession, ta in self.transcript_alleles_by_transcript_accession.items():
            transcript_id = ClinGenAllele._strip_transcript_version(transcript_accession)
            ta_by_t[transcript_id] = ta
        return ta_by_t

    def _get_transcript_allele(self, transcript_accession, match_version=True) -> Optional[Dict]:
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

    def get_c_hgvs(self, transcript_accession):
        """ c.HGVS has reference bases on it """
        from genes.models import TranscriptVersionSequenceInfo

        hgvs_string = None
        raw_hgvs_string, t_data = self._get_raw_c_hgvs_and_data(transcript_accession)
        if raw_hgvs_string:  # Has for this transcript version
            hgvs_name = HGVSName(raw_hgvs_string)
            if not hgvs_name.gene:  # Ref/Ens HGVSs have transcript no gene, LRG is set as gene
                hgvs_name.gene = t_data.get("geneSymbol")
            if hgvs_name.mutation_type in {"dup", "del", "delins"}:
                # We want to add reference bases onto HGVS but ClinGen reference sequence is wrong (see issue #493)
                coord = t_data["coordinates"][0]
                if "startIntronOffset" in coord:
                    # We have to accept these - as we get so many but we can't add on the bases
                    logging.warning("A coding DNA reference sequence does not contain intron or 5' and 3' gene "
                                    "flanking sequences and can therefore not be used as a reference to describe "
                                    "variants in these regions")
                else:
                    if hgvs_name.mutation_type == "dup":
                        ref_end = coord["end"]
                        ref_start = ref_end - len(coord["allele"])
                    else:
                        ref_start = coord["start"]
                        ref_end = coord["end"]

                    tvsi = TranscriptVersionSequenceInfo.get(transcript_accession)
                    hgvs_name.ref_allele = tvsi.sequence[ref_start:ref_end]
            hgvs_string = hgvs_name.format()
        return hgvs_string

    def _get_raw_c_hgvs_and_data(self, transcript_accession, match_version=True) -> Tuple[Optional[str],
                                                                                          Optional[Dict]]:
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

    @staticmethod
    def filtered_genomic_alleles(genomic_alleles, genome_build: GenomeBuild):
        return [ga for ga in genomic_alleles if ga.get("referenceGenome") == genome_build.name]

    def get_g_hgvs(self, genome_build: GenomeBuild):
        invalid_contigs = Contig.objects.none()
        valid_contigs = genome_build.contigs
        if settings.LIFTOVER_TO_CHROMOSOMES_ONLY:
            valid_contigs = valid_contigs.filter(role=SequenceRole.ASSEMBLED_MOLECULE)
            invalid_contigs = genome_build.contigs.exclude(role=SequenceRole.ASSEMBLED_MOLECULE)

        refseq_contigs = set(valid_contigs.values_list("refseq_accession", flat=True))
        refseq_invalid_contigs = set(invalid_contigs.values_list("refseq_accession", flat=True))
        genomic_alleles = self.api_response["genomicAlleles"]

        invalid_contigs = set()
        unknown_contigs = set()
        for ga in self.filtered_genomic_alleles(genomic_alleles, genome_build):
            for h in ga["hgvs"]:
                contig = h.split(":", 1)[0]
                if contig in refseq_contigs:
                    return h
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

    def get_variant_tuple(self, genome_build: GenomeBuild) -> 'VariantCoordinate':
        from genes.hgvs import get_hgvs_variant_tuple
        g_hgvs = self.get_g_hgvs(genome_build)
        return get_hgvs_variant_tuple(g_hgvs, genome_build)

    def get_variant_string(self, genome_build: GenomeBuild, abbreviate=False):
        from snpdb.models import Variant
        return Variant.format_tuple(*self.get_variant_tuple(genome_build), abbreviate=abbreviate)

    @property
    def human_url(self):
        params = (settings.CLINGEN_ALLELE_REGISTRY_DOMAIN, self)
        return "%s/redmine/projects/registry/genboree_registry/by_caid?caid=%s" % params

    @staticmethod
    def format_clingen_allele(pk):
        return f"CA{pk:06}"

    def __str__(self):
        """ ClinGen has minimum length of 6 (0 padded) but can go longer """
        return self.format_clingen_allele(self.pk)

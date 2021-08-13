from django.conf import settings
from django.db import models
from django_extensions.db.models import TimeStampedModel
import re

from snpdb.models.models_enums import SequenceRole
from snpdb.models.models_genome import GenomeBuild, Contig


class ClinGenAllele(TimeStampedModel):
    """ Canonical Allele - @see http://reg.clinicalgenome.org """
    # id = "CA" stripped off, eg CA7169043 => 7169043
    id = models.BigIntegerField(primary_key=True)
    api_response = models.JSONField(null=False)  # returned

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

    def get_p_hgvs(self, transcript_id, match_version=True):
        if not match_version:
            transcript_id = ClinGenAllele._strip_transcript_version(transcript_id)

        transcript_alleles = self.api_response.get("transcriptAlleles")
        if transcript_alleles:
            for ta in transcript_alleles:
                for t_hgvs in ta["hgvs"]:
                    hgvs_transcript = t_hgvs.split(":")[0]
                    if not match_version:
                        hgvs_transcript = ClinGenAllele._strip_transcript_version(hgvs_transcript)
                    if transcript_id == hgvs_transcript:
                        protein_effect = ta.get("proteinEffect")
                        if protein_effect:
                            p_hgvs = protein_effect.get("hgvs")
                            if p_hgvs:
                                return p_hgvs
        return None

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

    def __str__(self):
        return f"CA{self.pk}"

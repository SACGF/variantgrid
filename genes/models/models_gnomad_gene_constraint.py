from typing import Optional

from django.db import models
from django.db.models.deletion import CASCADE

from genes.models.models_gene import GeneSymbol, Transcript, TranscriptVersion


class GnomADGeneConstraint(models.Model):
    """ @see https://gnomad.broadinstitute.org/downloads#v4-constraint """
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, null=True, on_delete=CASCADE)
    transcript_version = models.ForeignKey(TranscriptVersion, null=True, on_delete=CASCADE)
    cached_web_resource = models.ForeignKey('annotation.CachedWebResource', on_delete=CASCADE)

    # Boolean indicator for MANE Select transcript
    mane_select = models.BooleanField(null=True, blank=True)

    # Number of observed high and low confidence loss-of-function (pLoF) variants
    lof_hc_lc_obs = models.IntegerField(null=True, blank=True)

    # Number of expected high and low confidence pLoF variants
    lof_hc_lc_exp = models.IntegerField(null=True, blank=True)

    # Number of possible high and low confidence pLoF variants
    lof_hc_lc_possible = models.IntegerField(null=True, blank=True)

    # Observed over expected ratio for high and low confidence pLoF variants
    lof_hc_lc_oe = models.FloatField(null=True, blank=True)

    # Mutation rate for high and low confidence pLoF variants
    lof_hc_lc_mu = models.FloatField(null=True, blank=True)

    # Probability of loss-of-function intolerance
    lof_hc_lc_pli = models.FloatField(null=True, blank=True)

    # Probability for recessive genes
    lof_hc_lc_prec = models.FloatField(null=True, blank=True)

    # Probability for unconstrained genes
    lof_hc_lc_pnull = models.FloatField(null=True, blank=True)

    # Number of observed high confidence pLoF variants
    lof_obs = models.IntegerField(null=True, blank=True)

    # Number of expected high confidence pLoF variants
    lof_exp = models.FloatField(null=True, blank=True)

    # Number of possible high confidence pLoF variants
    lof_possible = models.IntegerField(null=True, blank=True)

    # Observed over expected ratio for high confidence pLoF variants
    lof_oe = models.FloatField(null=True)

    # Mutation rate for high confidence pLoF variants
    lof_mu = models.FloatField(null=True, blank=True)

    # Probability of loss-of-function intolerance (high confidence)
    lof_pli = models.FloatField(null=True, blank=True)

    # Probability for recessive genes (high confidence)
    lof_prec = models.FloatField(null=True, blank=True)

    # Probability for unconstrained genes (high confidence)
    lof_pnull = models.FloatField(null=True, blank=True)

    # Lower bound of 90% CI for o/e ratio (high confidence pLoF variants)
    lof_oe_ci_lower = models.FloatField(null=True)

    # Upper bound of 90% CI for o/e ratio (high confidence pLoF variants)
    lof_oe_ci_upper = models.FloatField(null=True)

    # Raw Z-score for pLoF variants
    lof_z_raw = models.FloatField(null=True, blank=True)

    # Z-score for pLoF variants
    lof_z_score = models.FloatField(null=True, blank=True)

    # Number of observed missense variants
    mis_obs = models.IntegerField(null=True, blank=True)

    # Number of expected missense variants
    mis_exp = models.FloatField(null=True, blank=True)

    # Number of possible missense variants
    mis_possible = models.IntegerField(null=True, blank=True)

    # Observed over expected ratio for missense variants
    mis_oe = models.FloatField(null=True, blank=True)

    # Mutation rate for missense variants
    mis_mu = models.FloatField(null=True, blank=True)

    # Lower bound of 90% CI for o/e ratio (missense variants)
    mis_oe_ci_lower = models.FloatField(null=True, blank=True)

    # Upper bound of 90% CI for o/e ratio (missense variants)
    mis_oe_ci_upper = models.FloatField(null=True, blank=True)

    # Raw Z-score for missense variants
    mis_z_raw = models.FloatField(null=True, blank=True)

    # Z-score for missense variants
    mis_z_score = models.FloatField(null=True, blank=True)

    # Number of observed "probably damaging" missense variants (PolyPhen-2)
    mis_pphen_obs = models.IntegerField(null=True, blank=True)

    # Number of expected "probably damaging" missense variants (PolyPhen-2)
    mis_pphen_exp = models.IntegerField(null=True, blank=True)

    # Number of possible "probably damaging" missense variants (PolyPhen-2)
    mis_pphen_possible = models.IntegerField(null=True, blank=True)

    # O/E ratio for "probably damaging" missense variants (PolyPhen-2)
    mis_pphen_oe = models.FloatField(null=True, blank=True)

    # Number of observed synonymous variants
    syn_obs = models.IntegerField(null=True, blank=True)

    # Number of expected synonymous variants
    syn_exp = models.FloatField(null=True, blank=True)

    # Number of possible synonymous variants
    syn_possible = models.IntegerField(null=True, blank=True)

    # Observed over expected ratio for synonymous variants
    syn_oe = models.FloatField(null=True, blank=True)

    # Mutation rate for synonymous variants
    syn_mu = models.FloatField(null=True, blank=True)

    # Lower bound of 90% CI for o/e ratio (synonymous variants)
    syn_oe_ci_lower = models.FloatField(null=True, blank=True)

    # Upper bound of 90% CI for o/e ratio (synonymous variants)
    syn_oe_ci_upper = models.FloatField(null=True, blank=True)

    # Raw Z-score for synonymous variants
    syn_z_raw = models.FloatField(null=True, blank=True)

    # Z-score for synonymous variants
    syn_z_score = models.FloatField(null=True, blank=True)

    # Reason transcript does not have constraint metrics
    constraint_flags = models.JSONField(default=list)

    GNOMAD_BASE_URL = "https://gnomad.broadinstitute.org"

    @property
    def transcript_url(self):
        if tv := self.transcript_version:
            transcript_id = tv.transcript_id
        else:
            transcript_id = self.transcript_id
        return f"{self.GNOMAD_BASE_URL}/transcript/{transcript_id}"

    @property
    def gene_url(self):
        return f"{self.GNOMAD_BASE_URL}/gene/{self.gene_id}"

    @property
    def gene_symbol_url(self):
        return f"{self.GNOMAD_BASE_URL}/gene/{self.gene_symbol_id}"

    def _get_oe_summary(self, prefix: str) -> str:
        oe = getattr(self, f"{prefix}_oe")
        summary = "n/a"
        if oe is not None:
            oe_ci_lower = getattr(self, f"{prefix}_oe_ci_lower")
            oe_ci_upper = getattr(self, f"{prefix}_oe_ci_upper")
            summary = f"{oe:.2f} ({oe_ci_lower:.2f} - {oe_ci_upper:.2f})"
        return summary

    @property
    def mis_oe_summary(self) -> str:
        return self._get_oe_summary("mis")

    @property
    def syn_oe_summary(self) -> str:
        return self._get_oe_summary("syn")

    @property
    def lof_oe_summary(self) -> str:
        return self._get_oe_summary("lof")

    @property
    def transcript_accession(self) -> str:
        if self.transcript_version:
            return str(self.transcript_version)
        return str(self.transcript)

    @staticmethod
    def get_for_transcript_version(transcript_version: TranscriptVersion) -> Optional['GnomADGeneConstraint']:
        """ GnomADGeneConstraint uses Ensembl gene/transcripts - so load the most specific
            possible (transcript, gene, then symbol) """
        return GnomADGeneConstraint.get_for_transcript_version_with_method_and_url(transcript_version)[0]

    @staticmethod
    def get_for_transcript_version_with_method_and_url(transcript_version: TranscriptVersion) \
            -> tuple[Optional['GnomADGeneConstraint'], Optional[str], Optional[str]]:
        """ GnomADGeneConstraint uses Ensembl gene/transcripts - so load the most specific
            possible (transcript version, transcript then symbol) """
        ggc = None
        method = None
        url = None
        qs = GnomADGeneConstraint.objects.all()
        # May be able to get via transcript version or just transcript
        if ggc := qs.filter(transcript_version=transcript_version).first():
            url = ggc.transcript_url
            method = "transcript version"
        elif ggc := qs.filter(transcript=transcript_version.transcript).first():
            url = ggc.transcript_url
            method = "transcript"

        if ggc is None:
            if ggc := qs.filter(gene_symbol=transcript_version.gene_version.gene_symbol).first():
                url = ggc.gene_symbol_url
                method = "gene symbol"
        return ggc, method, url

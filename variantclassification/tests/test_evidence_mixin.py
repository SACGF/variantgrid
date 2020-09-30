from variantclassification.models import EvidenceMixin
from variantclassification.models.evidence_mixin import VCStore


class BasicEvidence(EvidenceMixin):

    def __init__(self, evidence: VCStore):
        self.evidence = evidence

    @property
    def _evidence(self) -> VCStore:
        return self.evidence

# doesn't work without Transcripts loaded now
# class EvidenceMixinTest(TestCase):
#
#     @override_settings(VARIANT_ANNOTATION_TRANSCRIPT_PREFERENCES=['refseq_transcript_accession'])
#     def test_get_transcript(self):
#         # if transcript version is in c.hgvs use it
#         be = BasicEvidence({
#             SpecialEKeys.C_HGVS: "NM_020975.5(RET):c.867+48A>G",
#             SpecialEKeys.REFSEQ_TRANSCRIPT_ID: "NM_020975",
#             SpecialEKeys.GENOME_BUILD: "GRCh37"
#         })
#         self.assertEqual(be.transcript, "NM_020975.5")
#
#         # if transcript version is in c.hgvs but transcript doesn't match
#         # value in transcript field, use the raw transcript value
#         be = BasicEvidence({
#             SpecialEKeys.C_HGVS: "NM_020975.5(RET):c.867+48A>G",
#             SpecialEKeys.REFSEQ_TRANSCRIPT_ID: "NM_033333",
#             SpecialEKeys.GENOME_BUILD: "GRCh37"
#         })
#         self.assertEqual(be.transcript, "NM_033333")
#
#         # if there is no transcript field, use the contents of c.hgvs
#         be = BasicEvidence({
#             SpecialEKeys.C_HGVS: "NM_020975.5(RET):c.867+48A>G",
#             SpecialEKeys.GENOME_BUILD: "GRCh37"
#         })
#         self.assertEqual(be.transcript, "NM_020975.5")

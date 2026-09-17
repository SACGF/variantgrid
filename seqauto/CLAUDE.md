# seqauto — agent notes
Owns: SequencingRun, SampleSheet/SequencingSample, the sequencing file records (UnalignedReads, BamFile, SingleSampleVCF,
JointCalledVCF), QC models, EnrichmentKit and gold coverage, and the seqauto REST API.
API:
- A client-visible API change needs a name in `variantgrid/views_rest.py:API_FEATURES` (see
  claude/guides/operations.md#authentication-surface).

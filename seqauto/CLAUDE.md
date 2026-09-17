# seqauto — agent notes
Owns: SequencingRun, SampleSheet/SequencingSample, the sequencing file records (UnalignedReads, BamFile, SingleSampleVCF,
JointCalledVCF), QC models, EnrichmentKit and gold coverage, and the seqauto REST API.
API:
- `seqauto/views_rest.py:API_FEATURES` is the contract clients read from */seqauto/api/v1/capabilities* to decide which calls
  a server accepts (VG3 and VG4 run side by side). Add a name in the same change as a client-visible feature, and keep
  names once added. `upload_file_types` is derived from the import task factories, so it needs no upkeep.

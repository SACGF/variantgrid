# Sparse test reference fastas

Stand-ins for the real reference genome fastas, used by the unit test suite on
GitHub Actions (`variantgrid/settings/env/github_actions.py` points each build's
`reference_fasta` here). Each file has the same contig names, order and lengths as
the real NCBI fasta it mirrors, but the sequence is all `N` except:

- the regions the test suite actually fetches (plus 10kb padding each side), which
  carry the real reference sequence, and
- contigs under 100kb that the suite touches (eg mitochondria), included in full.

That keeps them ~12MB each instead of ~3GB while behaving identically for every
sequence read the tests perform.

## Symptom that these are stale

A test that passes on a dev machine (which uses the real fastas in
`/data/annotation/fasta/`) but fails on CI with a sequence-related assertion - eg
reference bases coming back as `N` - means the test fetches a region these files
don't carry yet.

## Regenerating

On a machine with the real reference fastas:

```bash
python3 manage.py test --keepdb --testrunner=variantgrid.test_runner.FastaRecordingRunner
python3 scripts/generate_sparse_test_fastas.py /tmp/vg_fasta_regions.jsonl \
    variantgrid/data/reference/sparse_test_fastas
```

The recording runner logs every fasta region the suite reads (override the output
path with `VG_FASTA_REGIONS_FILE`); the generator rebuilds the sparse fastas from
those regions and verifies each one round-trips byte-identically against the real
fasta before it finishes.

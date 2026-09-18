# NGS-pipelines: one vg_api_full.py for VG3 and VG4, runnable off the TAU box

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-17; revised by Claude Opus 5 (claude-opus-5), 2026-09-17
Status: in progress - both PRs implemented on NGS-pipelines branch `sapath443_vg_capabilities`, awaiting test against VG3/VG4

Pipeline half of [sapath#443](https://github.com/SACGF/variantgrid_sapath/issues/443), for the
*GMP-TAU/NGS-pipelines* repo (checked out at `../NGS-pipelines`, script *scripts/vg_api_full.py*,
config *misc/vg_api_config.yaml*, tests *tests/test_vg_api_full.py*). It goes to TAU as two pull
requests: the first is independent of any server or client change, the second uses the capability probe.

## What it builds on (done)

**Server.** `variantgrid/views_rest.py:CapabilitiesView` at `GET /api/v1/capabilities` returns
`version`, `git_hash`, `features` (`variantgrid/views_rest.py:API_FEATURES`) and `upload_file_types`
(derived from the import task factories). The contract and how a VG3 server answers are in
`claude/guides/operations.md` (authentication surface).

**Client.** variantgrid_api 1.5.0 is on PyPI (repo `../variantgrid_api`, CHANGELOG `[1.5.0]`,
SACGF/variantgrid_api#20 and #21). The parts this script relies on:

```python
@dataclass(frozen=True)
class ServerCapabilities:          # variantgrid_api.data_models
    version: str                   # "vg4.0-12-gc174556", or "legacy"
    git_hash: Optional[str] = None
    features: FrozenSet[str] = frozenset()
    upload_file_types: FrozenSet[str] = frozenset()

class UnsupportedFeaturePolicy(Enum):   # variantgrid_api.api_client
    SKIP = "skip"     # log a warning, return None
    ERROR = "error"   # raise UnsupportedFeatureError (the default)
```

- `VariantGridAPI(server, token, unsupported_feature_policy=...)`. `vg_api.capabilities` is fetched
  once on first use; a 404 or a redirect (VG3 sending `/api/` to its login page) gives
  `ServerCapabilities.LEGACY`. Ungated calls make no probe.
- `supports(feature)`, `accepts_upload(file_type)` for branches a skip can't express.
- Gated: `create_patient` / `create_specimen` / `create_extraction` (`patients`),
  `create_specimen_measure(s)` (`specimen_measures`), `link_sequencing_sample_extraction`
  (`link_extraction`), `poll_upload_status` / `wait_for_annotation` / `download_annotated` /
  `annotate_vcf` (`upload_status`).
- `upload_file(filename, metadata=None, file_type=None)`: `file_type` isn't sent, it gates the upload
  on `accepts_upload(file_type)` (SKIP returns `None`). Non-empty `metadata` needs `upload_metadata`;
  under SKIP the file is still uploaded, without the metadata.
- `MockVariantGridAPI(capabilities=None, unsupported_feature_policy=ERROR)`: default capabilities are
  a current server with every gated feature and the TSO 500 upload types;
  `capabilities=ServerCapabilities.LEGACY` behaves as VG3. A call skipped under SKIP is not recorded
  in `mock.calls`.
- *examples/example_tso500.py* in the client repo is the reference flow for a TSO 500 pair.

## PR 1: run the script anywhere

Today the script can only run on the TAU box: the config path is a literal
`/tau/ngs_pipelines/shared_repo/misc/vg_api_config.yaml`, and `--server` is `choices=['VG_prod', 'VG_test']`,
a key into that config. Neither obstacle is about the pipeline logic.

- `--server` accepts either a key from the config's `servers` map or a URL. A value containing `://`
  is used as-is and then `--api-token` (or the `VARIANTGRID_API_TOKEN` environment variable, which
  the `vg_api` CLI already honours) is required; a key still resolves through the config as now.
  *tso500_run_wrapper.py* and *reanalysis_and_rerun.py* keep passing `VG_prod` and see no change.
- `--config` (default the current literal path, overridable by `VG_API_CONFIG`) replaces the
  module-level `open()` so importing the module in tests needs no `/tau`. The tests already inject
  `MockVariantGridAPI` into `run_api()`; with this change they also stop depending on the file
  existing.
- The `/tau/data/...` literals for the samplesheet directory and the WGS report base move into the
  yaml under a `paths:` key, read the same way `variant_callers` and `aligners` are. Values stay
  exactly what they are now.

With that, a local VG on the `vg3_sapath_prod` branch (own database, since its migrations differ from
master) is a valid target: build a fake run directory the way `_make_haem_run` in the tests does and
run:

```bash
scripts/vg_api_full.py --server http://localhost:8000 --api-token ... --batchid TSO_26_001_260101_NB551_0001 --results_dir /tmp/fake_runs
```

## PR 2: capability-gated TSO 500 upload

The shared TAU environment's variantgrid_api pin moves to `>=1.5.0`.

Construct the client with `unsupported_feature_policy=UnsupportedFeaturePolicy.SKIP` in `__main__`
so every new call below is a logged no-op on a server that lacks it.

**Files per TSO 500 pair.** Beside the per-arm small-variant VCF the script finds today
(`vcfsuffix = {'.hard-filtered.vcf', '_SpliceVariants.vcf'}` under `Results/<pairID>/<sample_id>/`),
add to the `'TSO' in experiment` branch:

| File | Where | Upload |
|---|---|---|
| `<pairID>_CombinedVariantOutput.tsv` | `Results/<pairID>/` | `upload_file(path, file_type="dragen_tso500_combined_variant_output", metadata={"genome_build": "GRCh37"})` |
| `<sample_id>_SpliceVariants.vcf` | RNA arm | only when `not vg_api.accepts_upload("dragen_tso500_combined_variant_output")` |
| *AllFusions.csv*, CNV VCFs | as today | unchanged |

The splice line is the one real branch: VG4 takes the splice calls from the CVO and would double
them if the VCF also arrived, VG3 has no CVO importer and needs the VCF. So the RNA arm's
`vcfsuffix` choice reads the capability. A new `upload_combined_variant_output` entry in `API_STEPS`
does the CVO upload for every pair in `samp_to_pairIDs`, after `upload_single_sample_vcf_file`.
The RNA arm only becomes a `SequencingFile` with `--include_rna`, which *tso500_run_wrapper.py* doesn't
pass, so a default run sends no splice VCF to either server; the CVO upload and both extraction links
don't depend on it.

**Patient chain.** For each pair, one `link_sequencing_sample_extraction(SequencingSampleLookup(...), extraction_reference)`
per arm, with the extraction reference being the container suffix the CVO names (the DNA and RNA
sample IDs' trailing `C` / `B` containers). On VG4 the server parks the claim until the extraction
exists (202 with `match_status` `Pending`), on VG3 the call is skipped and returns `None`. Name, DOB, sex, tissue and dates keep arriving from Mocha
through the patient / specimen / extraction API, so this script posts none of them; TMB, MSI and
the other specimen measures come out of the CVO on the server, so it posts no `specimen_measure`
either.

**Log line.** After constructing the client, print `vg_api.capabilities.version` and whether the CVO
is accepted, so a run's log says which shape it took.

## Tests

In *tests/test_vg_api_full.py*, a `_make_tso_run(tmp_path, batchid)` helper alongside
`_make_haem_run` that lays out a two-arm pair (`Results/<pairID>/<dna>/`, `<rna>/`, the CVO, the
splice VCF, `Logs_Intermediates/...` BAM and FASTQ paths), then:

- with the default mock (VG4 shape): `upload_file` is called with the CVO and never with the splice
  VCF; `link_sequencing_sample_extraction` is called once per arm;
- with `MockVariantGridAPI(capabilities=ServerCapabilities.LEGACY)`: the splice VCF is uploaded, the
  CVO and link calls are absent from `mock.calls`;
- PR 1: `--server http://x` with no token exits with a usage error, `--server VG_test` resolves
  through the config.

## Rollout

1. PR 1 merges and the TAU environment installs nothing new.
2. The VG4 test server (`VG_test` in their config) is deployed from a master that has
   `/api/v1/capabilities`, and TAU's environment installs variantgrid_api 1.5.0.
3. PR 2 merges. Against VG3 prod the runs are byte-for-byte what they are today plus one capabilities
   probe per run (a 404 or login redirect). When prod moves to VG4 the same script starts sending the CVO and
   stops sending the splice VCF with no further change on TAU's side.

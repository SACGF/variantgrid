# variantgrid_api 1.5: capability probe and feature-gated calls

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-17
Status: draft

Client half of [sapath#443](https://github.com/SACGF/variantgrid_sapath/issues/443), for the
*SACGF/variantgrid_api* repo (checked out at `../variantgrid_api`). It consumes the endpoint in
`claude/plans/sapath443_api_capabilities_server_plan.md` and is consumed by
`claude/plans/sapath443_ngs_pipelines_plan.md`. It also carries the patient / specimen / extraction
methods asked for in [variantgrid_api#20](https://github.com/SACGF/variantgrid_api/issues/20),
since those are exactly the calls that need gating.

## Data

A frozen dataclass in *src/variantgrid_api/data_models.py*, the parsed capabilities response:

```python
@dataclass(frozen=True)
class ServerCapabilities:
    version: str                       # "vg4.0-12-gc174556", or "legacy" when the endpoint 404s
    git_hash: Optional[str]
    features: frozenset[str]
    upload_file_types: frozenset[str]
```

`ServerCapabilities.LEGACY` is the instance for a 404: `version="legacy"`, both sets empty.

Plus the dataclasses from variantgrid_api#20 §1 (`Patient`, `Specimen`, `Extraction`,
`SpecimenMeasure`, `ExternalPK`, `ExternalReference`), as specified there.

An enum beside `EmptyInputPolicy` in *src/variantgrid_api/api_client.py*:

```python
class UnsupportedFeaturePolicy(Enum):
    SKIP = "skip"     # log a warning, return None
    ERROR = "error"   # raise UnsupportedFeatureError
```

## Client behaviour

`VariantGridAPI.__init__` gains `unsupported_feature_policy=UnsupportedFeaturePolicy.ERROR`. The
default is ERROR because a silent skip against a server the caller expected to be VG4 hides a
misconfiguration; TAU passes SKIP explicitly for the migration window (pipeline plan).

- `capabilities` is a lazily fetched, cached property: one `GET seqauto/api/v1/capabilities` on
  first use, `ServerCapabilities.LEGACY` on 404. Any other non-2xx raises as every other call does.
  Lazy so that constructing the client stays free of network I/O (tests, dry runs).
- `supports(feature: str) -> bool` and `accepts_upload(file_type: str) -> bool` read the two sets.
  These are public because the pipeline has one genuine branch (CVO versus splice VCF) that a
  no-op cannot express.
- `_require(feature)` applies the policy. Every new method calls it first:

| Method | Requires |
|---|---|
| `create_patient`, `create_specimen`, `create_extraction` | `patients` |
| `create_specimen_measure`, `create_specimen_measures` | `specimen_measures` |
| `link_sequencing_sample_extraction` | `link_extraction` |
| `poll_upload_status`, `wait_for_annotation`, `download_annotated`, `annotate_vcf` | `upload_status` |
| `upload_file(..., metadata={...})` with a non-empty `metadata` | `upload_metadata` |

Existing methods (`create_sequencing_run`, `create_sample_sheet`, `upload_file` with no metadata and
so on) stay ungated: they work on both servers.

`upload_file` gains `metadata: Optional[dict] = None`, sent as extra query params (variantgrid_api#20
§3), and `file_type: Optional[str] = None`. When `file_type` is given the call is gated on
`accepts_upload(file_type)`; with SKIP that means a CVO posted to VG3 is logged and skipped rather
than stored as a broken gene-list import. Existing callers pass neither and see no change.

Under SKIP a gated method returns `None`. The mock returns the same so pipeline tests see one shape.

## Mock

`MockVariantGridAPI(capabilities: Optional[ServerCapabilities] = None)`. Default is a VG4-shaped
instance with every feature and upload type in the tables above, so existing tests keep passing.
`MockVariantGridAPI(capabilities=ServerCapabilities.LEGACY)` is the VG3 shape. `supports()` /
`accepts_upload()` mirror the real client; the new methods record and return canned dicts like the
existing ones. `unsupported_feature_policy` is honoured so a pipeline test can assert that a VG3
run skipped `create_patient` and posted no CVO.

## Tests

*tests/test_api_client_capabilities.py* with `responses`, in the style of *tests/test_api_client.py*:

- 200 parses into `ServerCapabilities` and is fetched once across several `supports()` calls;
- 404 gives `LEGACY`, and a gated method raises `UnsupportedFeatureError` under ERROR, returns
  `None` and logs under SKIP;
- `upload_file(file_type="dragen_tso500_combined_variant_output")` posts on VG4 shape and skips on
  LEGACY under SKIP;
- an ungated method makes no capabilities request at all.

*tests/test_mock_variantgrid_api.py*: the LEGACY mock skips and the default mock records.

## Release

- *pyproject.toml* version `1.5.0`; CHANGELOG entry under `## [1.5.0]` naming both this plan's
  issue and variantgrid_api#20.
- README gains a short "Talking to more than one VariantGrid version" section: the SKIP policy, and
  `api.supports("...")` for the one branch a caller has to write.
- Requires-python stays at 3.8; `frozenset[str]` in annotations needs `from __future__ import annotations`
  at the top of the module.

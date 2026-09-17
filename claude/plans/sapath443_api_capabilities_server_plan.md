# API capabilities endpoint, so one client can talk to VG3 and VG4

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-17
Status: landed 6d736fe3e, c17455640 (VG3 port on local branch vg3_sapath_prod_capabilities)

Server half of [sapath#443](https://github.com/SACGF/variantgrid_sapath/issues/443): TAU's pipeline
(*GMP-TAU/NGS-pipelines*, *scripts/vg_api_full.py*) posts to VG3 prod today and will post to VG4 once
it is deployed. During the changeover the same client code runs against both, and the two servers
accept different things. This plan gives the client a way to ask which. The client and pipeline halves
are `claude/plans/sapath443_api_capabilities_client_plan.md` and
`claude/plans/sapath443_ngs_pipelines_plan.md`.

## Why a feature list and not a version number

VG3 is the `vg3_sapath_prod` branch. It forked from master in 2021 and has had features cherry-picked
onto it since (the seqauto `bulk_create` endpoints are there), so "3 or 4" says little about what a
given deployment accepts. What a client needs to know is which endpoints and which uploaded file
types exist, and the answer to that keeps changing on both sides. So the endpoint reports features by
name, and a version string only as information for humans and logs.

What differs today, from `seqauto/urls.py`, `patients/urls.py`, `upload/urls.py` on master against
the same files on `vg3_sapath_prod`:

| Master has, VG3 lacks | Client consequence |
|---|---|
| *patients/api/v1/patient*, *specimen*, *extraction*, *specimen_measure* (+ *bulk_create*) | skip on VG3 |
| *seqauto/api/v1/sequencing_sample/link_extraction* | skip on VG3 |
| *upload/api/v1/upload_status*, *download* | skip on VG3 |
| *seqauto/api/v1/single_sample_vcf* (VG3 only has the *vcf_file* alias) | alias works on both |
| *file_upload* processes *_CombinedVariantOutput.tsv* (`UploadedFileTypes.DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT`) | on VG3 the `.tsv` is claimed by the gene-list importer and fails; the client just sends the TSO500 variants vcf |
| *file_upload* accepts upload metadata query params (`genome_build`, `extraction`) | VG3 ignores unknown query params |
| */api/schema* (drf-spectacular) | lists endpoints but cannot say which file types `file_upload` will process |

The last row is why the OpenAPI schema is not enough on its own: the CVO goes through the same URL on both servers, and only the processing differs.

## The endpoint

`GET seqauto/api/v1/capabilities`, a read-only `rest_framework.views.APIView` in `seqauto/views_rest.py`,
registered in `seqauto/urls.py` next to `SequencingSampleExtractionLinkView`. Token auth as for every
other API call.

It lives under *seqauto/api/* rather than */api/* because VG3's `PUBLIC_PATHS` lists `^/seqauto/api/.*`
but not `^/api/.*`: on VG3 a request to `/api/v1/capabilities` is redirected to the login page (a 200
of HTML after `requests` follows the redirect), whereas `/seqauto/api/v1/capabilities` is a clean 404.
The client treats 404 as "legacy server, no capabilities".

Response:

```json
{
  "version": "vg4-gc17455640",
  "git_hash": "2130cffe0",
  "features": [
    "patients",
    "specimen_measures",
    "link_extraction",
    "upload_status",
    "joint_called_vcf_cross_run",
    "upload_metadata"
  ],
  "upload_file_types": [
    "vcf",
    "gene_coverage",
    "dragen_tso500_all_fusions",
    "dragen_tso500_combined_variant_output",
    "gene_level_cnv_vcf"
  ]
}
```

- `version` is the `VARIANTGRID_VERSION` setting (`git describe` against the latest `vg<major>.*` tag),
  `git_hash` from `library.git.Git` as `variantgrid/views.py:version` already does.
- `features` is a module-level tuple in `seqauto/views_rest.py`. Each name is a fact about this
  codebase, appended when the feature lands and never removed while the endpoint exists. The names
  above are the ones the client plan consumes; the list is the contract, so a new client-visible
  feature adds a name here in the same change.
- `upload_file_types` is derived, not hand-kept: the `UploadedFileTypes` value's name in lower case
  for every non-abstract `ImportTaskFactory` from
  `upload/import_task_factories/import_task_factory.py:get_import_task_factories`, so a deployment
  whose `IMPORT_TASK_FACTORY_IMPORTS` setting leaves a factory out reports honestly. Only the
  factories a client can upload to are worth listing; the internal ones (`ANALYSIS`, `LIFTOVER`,
  `MANUAL_VARIANT_ENTRY`, `WIKI_*`, `VARIANT_TAGS`, `CLINVAR`) are excluded by a small denylist
  beside the tuple.

Decorate with `@extend_schema` like the neighbouring views so `/api/docs` shows it.

## VG3 branch

Cherry-pick a reduced view onto `vg3_sapath_prod`: same URL and shape, `version` from the `vg3.*` tag (`vg3.0-622-g2eaf26c28`),
`features` empty, `upload_file_types` from the same derivation (which on VG3 yields `vcf`,
`gene_coverage`, `gene_list`, `bed` and so on). Then the client's 404 fallback is only for
deployments older than that cherry-pick, and a VG3 answer is positive evidence rather than an
absence. This is optional for the client to work; it is worth doing because it makes the TAU
integration test against VG3 assert something.

## Tests

A new *seqauto/tests/test_capabilities.py*, an `APITestCase` in the style of
`seqauto/tests/test_extraction_link.py`:

- authenticated GET returns 200 with the three keys, `patients` in `features` and
  `dragen_tso500_combined_variant_output` in `upload_file_types`;
- anonymous GET is 401 (the path is public to Django so DRF answers, and `IsAuthenticated` refuses).

## Docs

- A new `seqauto/CLAUDE.md` (the app has none yet) with one API note: the feature list is the client
  contract, add a name when a client-visible feature lands.
- `claude/guides/operations.md#authentication-surface`: the endpoint is under *seqauto/api/* for the
  VG3 public-path reason above, so nobody later "tidies" it to */api/v1/*.

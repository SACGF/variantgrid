# NGS-pipelines: one vg_api_full.py for VG3 and VG4, runnable off the TAU box

Written by Claude Fable 5.1 (claude-fable-5-1), 2026-09-17; revised by Claude Opus 5 (claude-opus-5), 2026-09-18
Status: in progress - implemented on NGS-pipelines branch `sapath443_vg_capabilities` (davmlaw fork), awaiting test against VG3/VG4 and a PR to GMP-TAU

Pipeline half of [sapath#443](https://github.com/SACGF/variantgrid_sapath/issues/443), for the
*GMP-TAU/NGS-pipelines* repo (checked out at `../NGS-pipelines`).

## Done

- **Server:** `variantgrid/views_rest.py:CapabilitiesView` (`GET /api/v1/capabilities`), contract in
  `claude/guides/operations.md` (authentication surface).
- **Client:** variantgrid_api 1.6.0 on PyPI; *examples/example_tso500.py* in `../variantgrid_api` is the reference TSO 500 flow.
- **Pipeline:** commit `d68ebc4` on `sapath443_vg_capabilities`. *scripts/vg_api_full.py* takes `--server` as a
  config key or URL (URL needs `--api-token` / `VARIANTGRID_API_TOKEN`), `--config` / `VG_API_CONFIG`, and reads
  the `/tau/data` paths from the yaml's `paths:`. Runs with `UnsupportedFeaturePolicy.SKIP`: uploads the TSO 500
  CombinedVariantOutput and links each arm's extraction when the server supports it, and falls back to the splice
  VCF on VG3. The config gains a `VG4_test` server; the TAU env pin moves to `variantgrid-api>=1.6.0`.
  Tests: `_make_tso_run` and the `test_run_api_tso_*` / `test_resolve_server_*` tests.

## Remaining

1. **Test locally against VG3:** a local VG on the `vg3_sapath_prod` branch (own database, its migrations differ
   from master). Build a fake run directory the way `_make_tso_run` does and run
   ```bash
   scripts/vg_api_full.py --server http://localhost:8000 --api-token ... --batchid TSO_26_001_260101_NB551_0001 --results_dir /tmp/fake_runs
   ```
   Expect the log's capabilities line to say `legacy`, and the uploads to match what the script sends today.
2. **Test against VG4:** the same run against `VG4_test`, deployed from a master that has `/api/v1/capabilities`.
   Expect the CVO upload, no splice VCF, and one extraction link per arm (202 `Pending` until Mocha creates the extraction).
3. **PR to GMP-TAU** from the fork branch. The PR description carries the rollout: TAU's environment installs
   variantgrid_api 1.6.0 with the merge. Against VG3 prod, runs are unchanged apart from one capabilities probe;
   when prod moves to VG4, the same script switches to the CVO with no further change on TAU's side.

Once the PR is merged, delete this plan.
